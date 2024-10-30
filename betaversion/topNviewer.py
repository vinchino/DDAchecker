import csv
import base64
import zlib
import numpy as np
import h5py
from lxml import etree
import os
import re
import argparse

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
                            if '=' in part:
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
                                if binary_data and data_type:
                                    try:
                                        decoded_data = base64.b64decode(binary_data)
                                        if compression:
                                            decompressed_data = zlib.decompress(decoded_data)
                                        else:
                                            decompressed_data = decoded_data
                                        numerical_array = np.frombuffer(decompressed_data, dtype=data_type)
                                        
                                        # Assign to the appropriate array
                                        if array_name == 'm/z':
                                            mz_array = numerical_array
                                        elif array_name == 'intensity':
                                            intensity_array = numerical_array
                                        # Add handling for other arrays if needed
                                    except (base64.binascii.Error, zlib.error, ValueError) as e:
                                        print(f"Error decoding binary data in spectrum {spectrum_data['Index Number']}: {e}")
                    
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
                            
                            # Store Cycle Number and Experiment Number as attributes
                            grp.attrs['Cycle Number'] = spectrum_data['Cycle Number']
                            grp.attrs['Experiment Number'] = spectrum_data['Experiment Number']
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

def create_cycle_experiment_mapping(metadata_csv):
    """
    Creates a mapping from cycle number to experiment 1 Index Numbers.
    
    Parameters:
        metadata_csv (str): Path to the metadata CSV file.
    
    Returns:
        dict: Mapping from cycle number to list of experiment 1 Index Numbers.
    """
    cycle_to_exp1_indices = {}
    
    with open(metadata_csv, mode='r', newline='', encoding='utf-8') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            cycle_number = int(row['Cycle Number']) if row['Cycle Number'] else None
            experiment_number = int(row['Experiment Number']) if row['Experiment Number'] else None
            index_number = int(row['Index Number']) if row['Index Number'] else None
            
            if cycle_number is not None and experiment_number == 1 and index_number is not None:
                if cycle_number not in cycle_to_exp1_indices:
                    cycle_to_exp1_indices[cycle_number] = []
                cycle_to_exp1_indices[cycle_number].append(index_number)
    
    return cycle_to_exp1_indices

def annotate_csv(metadata_csv, data_h5, output_dir):
    """
    Annotates the metadata CSV with Matched Peak Rank and Matched Peak Intensity based on HDF5 data.
    
    Parameters:
        metadata_csv (str): Path to the metadata CSV file.
        data_h5 (str): Path to the HDF5 data file.
        output_dir (str): Directory to store the annotated CSV file.
    
    Returns:
        str: Path to the annotated CSV file.
    """
    import sys

    # Define output annotated CSV path
    base_filename = os.path.splitext(os.path.basename(metadata_csv))[0]
    annotated_csv = os.path.join(output_dir, f"{base_filename}_annotated.csv")
    
    try:
        # Create a mapping from cycle number to experiment 1 Index Numbers
        cycle_to_exp1_indices = create_cycle_experiment_mapping(metadata_csv)
        
        # Open the metadata CSV for reading
        with open(metadata_csv, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            rows = list(reader)  # Read all rows into memory
        
        # Open the HDF5 file for reading
        with h5py.File(data_h5, 'r') as h5file:
            # Create a mapping from Index Number to (m/z, intensity)
            index_to_data = {}
            for index in h5file.keys():
                try:
                    mz = h5file[index]['m/z'][:]
                    intensity = h5file[index]['intensity'][:]
                    index_to_data[int(index)] = (mz, intensity)
                except KeyError as e:
                    print(f"Missing dataset in HDF5 for spectrum {index}: {e}", file=sys.stderr)
        
        # Open the annotated CSV for writing
        with open(annotated_csv, mode='w', newline='', encoding='utf-8') as outfile:
            fieldnames = reader.fieldnames + ['Matched Peak Rank', 'Matched Peak Intensity']
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for row in rows:
                ms_level = int(row['MS Level']) if row['MS Level'] else None
                experiment_number = int(row['Experiment Number']) if row['Experiment Number'] else None
                cycle_number = int(row['Cycle Number']) if row['Cycle Number'] else None
                selected_mz = float(row['Selected Ion m/z']) if row['Selected Ion m/z'] else None
                index_number = int(row['Index Number']) if row['Index Number'] else None
                
                # Initialize annotation fields
                row['Matched Peak Rank'] = None
                row['Matched Peak Intensity'] = None
                
                # Proceed only if MS Level == 2 and Selected Ion m/z is available
                if ms_level == 2 and selected_mz is not None and cycle_number is not None:
                    # Get experiment 1 Index Numbers for the same cycle number
                    exp1_indices = cycle_to_exp1_indices.get(cycle_number, [])
                    
                    best_rank = None
                    best_intensity = None
                    best_mz_diff = None
                    
                    # Iterate over all experiment 1 spectra in the same cycle
                    for exp1_index in exp1_indices:
                        mz_intensity = index_to_data.get(exp1_index)
                        if mz_intensity is None:
                            continue
                        exp1_mz, exp1_intensity = mz_intensity
                        
                        # Calculate tolerance in m/z
                        tolerance = selected_mz * 50e-6  # 50 ppm
                        
                        # Find matches within 50 ppm
                        mz_diff = np.abs(exp1_mz - selected_mz)
                        matches = mz_diff <= tolerance
                        
                        if np.any(matches):
                            # If multiple matches, select the one with the smallest m/z difference
                            matched_indices = np.where(matches)[0]
                            closest_idx = matched_indices[np.argmin(mz_diff[matches])]
                            
                            matched_intensity = exp1_intensity[closest_idx]
                            matched_mz_diff = mz_diff[closest_idx]
                            
                            # Rank the intensities in descending order
                            sorted_indices = np.argsort(-exp1_intensity)
                            
                            # Find the rank (1-based)
                            rank = np.where(sorted_indices == closest_idx)[0][0] + 1
                            
                            # Update best match if this is better
                            if (best_rank is None) or (matched_mz_diff < best_mz_diff):
                                best_rank = rank
                                best_intensity = matched_intensity
                                best_mz_diff = matched_mz_diff
                                
                    # Annotate the spectrum data if a match was found
                    if best_rank is not None:
                        row['Matched Peak Rank'] = best_rank
                        row['Matched Peak Intensity'] = best_intensity
                
                # Write the updated row to the annotated CSV
                writer.writerow(row)
    
    except Exception as e:
        print(f"Error during annotation: {e}")
        raise
    
    print(f"Annotated CSV created at '{annotated_csv}'.")
    return annotated_csv

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Parse mzML files, extract metadata and numerical data, and annotate matched peaks.")
    parser.add_argument('mzml_files', metavar='mzml', type=str, nargs='+',
                        help='Path(s) to mzML file(s) to be parsed.')
    parser.add_argument('--output_dir', type=str, default='parsed_data',
                        help='Directory to store the output CSV and HDF5 files.')
    parser.add_argument('--annotate', action='store_true',
                        help='Flag to perform annotation after parsing.')
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Process each mzML file
    for mzml_file in args.mzml_files:
        try:
            print(f"Parsing mzML file '{mzml_file}'...")
            metadata_csv, data_h5 = parse_mzml(mzml_file, args.output_dir)
            print(f"Metadata CSV: '{metadata_csv}', Data HDF5: '{data_h5}'")
            
            if args.annotate:
                print(f"Starting annotation for '{mzml_file}'...")
                annotated_csv = annotate_csv(metadata_csv, data_h5, args.output_dir)
                print(f"Annotated CSV: '{annotated_csv}'")
        
        except Exception as e:
            print(f"Error processing file '{mzml_file}': {e}")

if __name__ == "__main__":
    main()
