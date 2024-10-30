import csv
import base64
import zlib
import numpy as np
import h5py
from lxml import etree
import os
import re

def sanitize_filename(filename):
    """
    Removes or replaces characters that are invalid in Windows filenames.
    
    Parameters:
        filename (str): The original filename.
    
    Returns:
        str: A sanitized filename safe for use in Windows.
    """
    # Replace invalid characters with underscore
    return re.sub(r'[<>:"/\\|?*]', '_', filename)

def parse_mzml(file_path, output_dir):
    """
    Parses an mzML file, extracts metadata and numerical data, and stores them in CSV and HDF5 files.
    
    Parameters:
        file_path (str): Path to the mzML file.
        output_dir (str): Directory to store the output CSV and HDF5 files.
        
    Returns:
        tuple: Paths to the generated CSV and HDF5 files.
    """
    # Define namespaces
    namespaces = {
        'mzml': 'http://psi.hupo.org/ms/mzml',
        'cv': 'http://psi.hupo.org/ms/mzml/cv'
    }
    
    # Extract base filename without extension and sanitize it
    base_filename = os.path.splitext(os.path.basename(file_path))[0]
    base_filename = sanitize_filename(base_filename)
    
    # Define output file paths
    metadata_csv = os.path.join(output_dir, f"{base_filename}_metadata.csv")
    data_h5 = os.path.join(output_dir, f"{base_filename}_data.h5")
    
    try:
        # Open CSV for writing metadata
        with open(metadata_csv, mode='w', newline='', encoding='utf-8') as csvfile:
            fieldnames = [
                'Index Number',
                'Cycle Number',
                'Experiment Number',
                'MS Level',
                'Total Ion Current',
                'Scan Start Time (min)',
                'Selected Ion m/z'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Open HDF5 file for storing numerical data
            with h5py.File(data_h5, 'w') as h5file:
                # Initialize iterparse context for 'spectrum' elements
                context = etree.iterparse(file_path, events=('end',), tag='{http://psi.hupo.org/ms/mzml}spectrum')
                
                for event, elem in context:
                    # Initialize a dictionary to hold metadata
                    spectrum_data = {
                        'Index Number': None,
                        'Cycle Number': None,
                        'Experiment Number': None,
                        'MS Level': None,
                        'Total Ion Current': None,
                        'Scan Start Time (min)': None,
                        'Selected Ion m/z': None
                    }
                    
                    # Extract 'index' attribute
                    spectrum_index = elem.get('index')
                    spectrum_data['Index Number'] = int(spectrum_index) if spectrum_index else None
                    
                    # Extract and parse 'id' attribute
                    spectrum_id = elem.get('id')  # e.g., "sample=1 period=1 cycle=2 experiment=2"
                    if spectrum_id:
                        id_parts = spectrum_id.split()
                        for part in id_parts:
                            key, value = part.split('=')
                            if key == 'cycle':
                                spectrum_data['Cycle Number'] = int(value)
                            elif key == 'experiment':
                                spectrum_data['Experiment Number'] = int(value)
                    
                    # Initialize placeholders for numerical data
                    mz_array = None
                    intensity_array = None
                    
                    # Iterate through child elements to extract metadata and data
                    for child in elem:
                        tag = etree.QName(child).localname
                        
                        if tag == 'cvParam':
                            accession = child.get('accession')
                            name = child.get('name')
                            value = child.get('value', '')
                            
                            if accession == 'MS:1000511':  # MS level
                                try:
                                    spectrum_data['MS Level'] = int(value)
                                except ValueError:
                                    spectrum_data['MS Level'] = None
                            elif accession == 'MS:1000285':  # Total Ion Current
                                try:
                                    spectrum_data['Total Ion Current'] = float(value)
                                except ValueError:
                                    spectrum_data['Total Ion Current'] = None
                        
                        elif tag == 'scanList':
                            # Assuming count="1", process the first <scan>
                            scan = child.find('{http://psi.hupo.org/ms/mzml}scan', namespaces)
                            if scan is not None:
                                for scan_child in scan:
                                    scan_tag = etree.QName(scan_child).localname
                                    
                                    if scan_tag == 'cvParam':
                                        scan_accession = scan_child.get('accession')
                                        scan_name = scan_child.get('name')
                                        scan_value = scan_child.get('value', '')
                                        
                                        if scan_accession == 'MS:1000016':  # Scan start time
                                            try:
                                                spectrum_data['Scan Start Time (min)'] = float(scan_value)
                                            except ValueError:
                                                spectrum_data['Scan Start Time (min)'] = None
                                        
                        elif tag == 'precursorList':
                            # Process the first <precursor>
                            precursor = child.find('{http://psi.hupo.org/ms/mzml}precursor', namespaces)
                            if precursor is not None:
                                # Process <selectedIonList>
                                selected_ion_list = precursor.find('{http://psi.hupo.org/ms/mzml}selectedIonList', namespaces)
                                if selected_ion_list is not None:
                                    selected_ion = selected_ion_list.find('{http://psi.hupo.org/ms/mzml}selectedIon', namespaces)
                                    if selected_ion is not None:
                                        for si_child in selected_ion:
                                            si_tag = etree.QName(si_child).localname
                                            
                                            if si_tag == 'cvParam':
                                                si_accession = si_child.get('accession')
                                                si_name = si_child.get('name')
                                                si_value = si_child.get('value', '')
                                                
                                                if si_accession == 'MS:1000744':  # Selected ion m/z
                                                    try:
                                                        spectrum_data['Selected Ion m/z'] = float(si_value)
                                                    except ValueError:
                                                        spectrum_data['Selected Ion m/z'] = None
                                        
                        elif tag == 'binaryDataArrayList':
                            # Process all <binaryDataArray> elements
                            for bda in child.findall('{http://psi.hupo.org/ms/mzml}binaryDataArray', namespaces):
                                # Initialize variables
                                data_type = None
                                compression = False
                                array_name = None
                                binary_data = None
                                
                                # Extract cvParams and binary data
                                for bda_child in bda:
                                    bda_tag = etree.QName(bda_child).localname
                                    
                                    if bda_tag == 'cvParam':
                                        bda_accession = bda_child.get('accession')
                                        bda_name = bda_child.get('name')
                                        bda_value = bda_child.get('value', '')
                                        
                                        if bda_accession == 'MS:1000523':  # 64-bit float
                                            data_type = np.float64
                                        elif bda_accession == 'MS:1000522':  # 64-bit integer
                                            data_type = np.int64
                                        elif bda_accession == 'MS:1000574':  # zlib compression
                                            compression = True
                                        elif 'm/z array' in bda_name.lower():
                                            array_name = 'm/z'
                                        elif 'intensity array' in bda_name.lower():
                                            array_name = 'intensity'
                                        elif 'time array' in bda_name.lower():
                                            array_name = 'time'
                                        # Add more conditions if needed
                                    
                                    elif bda_tag == 'binary':
                                        binary_data = bda_child.text
                                
                                # Decode Base64
                                if binary_data:
                                    try:
                                        decoded_data = base64.b64decode(binary_data)
                                    except base64.binascii.Error as e:
                                        print(f"Base64 decoding error in spectrum {spectrum_data['Index Number']}: {e}")
                                        continue
                                    
                                    # Decompress if needed
                                    if compression:
                                        try:
                                            decompressed_data = zlib.decompress(decoded_data)
                                        except zlib.error as e:
                                            print(f"Zlib decompression error in spectrum {spectrum_data['Index Number']}: {e}")
                                            decompressed_data = decoded_data
                                    else:
                                        decompressed_data = decoded_data
                                    
                                    # Convert to numerical array
                                    try:
                                        numerical_array = np.frombuffer(decompressed_data, dtype=data_type)
                                    except ValueError as e:
                                        print(f"Error converting binary data to numpy array in spectrum {spectrum_data['Index Number']}: {e}")
                                        numerical_array = None
                                    
                                    # Assign to the appropriate array
                                    if array_name == 'm/z':
                                        mz_array = numerical_array
                                    elif array_name == 'intensity':
                                        intensity_array = numerical_array
                                    # Add handling for other arrays if needed

                    # After processing all child elements, store data
                    # Write metadata to CSV
                    writer.writerow(spectrum_data)
                    
                    # Store numerical data in HDF5
                    if mz_array is not None and intensity_array is not None:
                        try:
                            # Create a group for the spectrum index
                            grp = h5file.create_group(str(spectrum_data['Index Number']))
                            grp.create_dataset('m/z', data=mz_array)
                            grp.create_dataset('intensity', data=intensity_array)
                        except Exception as e:
                            print(f"Failed to create HDF5 datasets for spectrum {spectrum_data['Index Number']}: {e}")
                    
                    # Clear the processed element to free memory
                    elem.clear()
                    while elem.getprevious() is not None:
                        del elem.getparent()[0]
        
    except OSError as e:
        print(f"Failed to create or write to output files: {e}")
        raise

    print(f"Parsing complete for '{file_path}'.")
    return metadata_csv, data_h5

if __name__ == "__main__":
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Parse mzML files and extract metadata and numerical data.")
    parser.add_argument('mzml_files', metavar='mzml', type=str, nargs='+',
                        help='Path(s) to mzML file(s) to be parsed.')
    parser.add_argument('--output_dir', type=str, default='parsed_data',
                        help='Directory to store the output CSV and HDF5 files.')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse each mzML file
    for mzml_file in args.mzml_files:
        try:
            parse_mzml(mzml_file, args.output_dir)
        except Exception as e:
            print(f"Error parsing file '{mzml_file}': {e}")
