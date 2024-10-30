# gui.py

import sys
import os
import subprocess
import numpy as np
import h5py
import csv
import re  # For parsing strings
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QPushButton, QFileDialog, QWidget, QAbstractItemView, QListWidget, QListWidgetItem,
    QLabel, QLineEdit, QHBoxLayout
)
from PyQt5.QtCore import Qt, QEvent
import pyqtgraph as pg  # For plotting
from PyQt5 import QtWidgets, QtGui

class FileListWidget(QListWidget):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setAcceptDrops(True)
        self.setSelectionMode(QAbstractItemView.SingleSelection)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()
    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()
    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
            file_paths = []
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    file_paths.append(str(url.toLocalFile()))
            self.main_window.load_mzml_files(file_paths)
        else:
            event.ignore()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("mzML File Processor")
        self.setGeometry(100, 100, 1200, 600)
        
        # Create splitter for the three panes
        splitter = QSplitter(Qt.Horizontal)
        
        # Left pane: List of mzML filenames and controls
        self.left_pane_widget = QWidget()
        left_layout = QVBoxLayout()
        
        # Instructions label
        instructions = QLabel("Drag and drop mzML files here or click 'Add Files' to select.")
        left_layout.addWidget(instructions)
        
        # File list widget
        self.left_pane = FileListWidget(main_window=self, parent=self)
        left_layout.addWidget(self.left_pane)
        
        # Add Files button
        add_files_button = QPushButton("Add Files")
        add_files_button.clicked.connect(self.select_files)
        left_layout.addWidget(add_files_button)
        
        # Processing parameters
        parameters_layout = QHBoxLayout()
        
        ppm_label = QLabel("PPM Tolerance:")
        self.ppm_input = QLineEdit("50")
        self.ppm_input.setFixedWidth(50)
        
        topx_label = QLabel("Top X Peaks:")
        self.topx_input = QLineEdit("30")
        self.topx_input.setFixedWidth(50)
        
        offset_label = QLabel("Cycle Offset Range:")
        self.offset_input = QLineEdit("-8 to +8")
        self.offset_input.setFixedWidth(70)
        
        parameters_layout.addWidget(ppm_label)
        parameters_layout.addWidget(self.ppm_input)
        parameters_layout.addWidget(topx_label)
        parameters_layout.addWidget(self.topx_input)
        parameters_layout.addWidget(offset_label)
        parameters_layout.addWidget(self.offset_input)
        
        left_layout.addLayout(parameters_layout)
        
        # Process Data button at the bottom of the left pane
        self.process_button = QPushButton("Process Data")
        left_layout.addWidget(self.process_button)
        self.left_pane_widget.setLayout(left_layout)
        
        # Middle pane: Metadata table
        self.middle_pane = QTableWidget()
        self.middle_pane.setEditTriggers(QAbstractItemView.NoEditTriggers)
        
        # Right pane: Plotting area
        self.right_pane = QWidget()
        right_layout = QVBoxLayout()
        # Plot button
        self.plot_button = QPushButton("No MS2 for Top Peaks")
        right_layout.addWidget(self.plot_button)
        # Plot widget for the scatter plot
        self.plot_widget = pg.PlotWidget()
        right_layout.addWidget(self.plot_widget)
        # Spectrum plot widget for mass spectra
        self.spectrum_plot_widget = pg.PlotWidget()
        right_layout.addWidget(self.spectrum_plot_widget)
        self.right_pane.setLayout(right_layout)
        
        # Add panes to splitter
        splitter.addWidget(self.left_pane_widget)
        splitter.addWidget(self.middle_pane)
        splitter.addWidget(self.right_pane)
        
        self.setCentralWidget(splitter)
        
        # Connect signals
        self.process_button.clicked.connect(self.process_data)
        self.plot_button.clicked.connect(self.plot_no_ms2_for_top_peaks)
        self.left_pane.itemSelectionChanged.connect(self.on_file_selected)
        
        # Variables to store file paths and data
        self.mzml_file_paths = []
        self.metadata_dict = {}  # Dictionary to hold metadata for each file
        self.headers = []



    
    def select_files(self):
        options = QFileDialog.Options()
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select mzML File(s)", "", "mzML Files (*.mzML);;All Files (*)", options=options
        )
        if files:
            self.load_mzml_files(files)
    
    def load_mzml_files(self, file_paths):
        for file_path in file_paths:
            if file_path not in self.mzml_file_paths:
                self.mzml_file_paths.append(file_path)
                filename = os.path.basename(file_path)
                item = QListWidgetItem(filename)
                item.setToolTip(file_path)
                self.left_pane.addItem(item)
    
    def on_file_selected(self):
        selected_items = self.left_pane.selectedItems()
        if selected_items:
            item = selected_items[0]
            file_path = item.toolTip()
            if file_path in self.metadata_dict:
                self.display_metadata_for_file(file_path)
            else:
                # Parse the file if not already parsed
                self.parse_and_display_metadata_for_file(file_path)
        else:
            # Clear the table if no file is selected
            self.middle_pane.clear()
    
    def parse_and_display_metadata_for_file(self, file_path):
        output_dir = os.path.join(os.getcwd(), 'parsed')  # Save files in 'parsed' folder
        os.makedirs(output_dir, exist_ok=True)
        
        # Call the parse_mzml.py script
        try:
            subprocess.run(['python', 'parse_mzml.py', file_path, '--output_dir', output_dir], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running parse_mzml.py on {file_path}: {e}")
            return
        
        # Determine the output file names
        base_filename = os.path.splitext(os.path.basename(file_path))[0]
        sanitized_filename = self.sanitize_filename(base_filename)
        metadata_csv = os.path.join(output_dir, f"{sanitized_filename}_metadata.csv")
        data_h5 = os.path.join(output_dir, f"{sanitized_filename}_data.h5")
        
        # Load metadata from CSV
        metadata = self.load_metadata_from_csv(metadata_csv)
        if metadata:
            self.metadata_dict[file_path] = {
                'metadata': metadata,
                'data_h5': data_h5
            }
            self.display_metadata_for_file(file_path)
    
    def sanitize_filename(self, filename):
        """
        Removes or replaces characters that are invalid in Windows filenames.
        """
        import re
        return re.sub(r'[<>:"/\\|?*]', '_', filename)
    
    def load_metadata_from_csv(self, csv_path):
        metadata = []
        with open(csv_path, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            if not self.headers:
                self.headers = reader.fieldnames
            for row in reader:
                # Convert appropriate fields to numeric types
                for key in row:
                    if key in ['Index Number', 'Cycle Number', 'Experiment Number', 'MS Level']:
                        row[key] = int(row[key]) if row[key] else None
                    elif key in ['Total Ion Current', 'Scan Start Time (min)', 'Selected Ion m/z']:
                        row[key] = float(row[key]) if row[key] else None
                metadata.append(row)
        return metadata
    
    def display_metadata_for_file(self, file_path):
        data = self.metadata_dict.get(file_path)
        if data:
            metadata = data['metadata']
            self.middle_pane.setRowCount(len(metadata))
            self.middle_pane.setColumnCount(len(self.headers))
            self.middle_pane.setHorizontalHeaderLabels(self.headers)
            
            for row_idx, data_row in enumerate(metadata):
                for col_idx, key in enumerate(self.headers):
                    value = data_row.get(key, '')
                    item = QTableWidgetItem(str(value) if value is not None else '')
                    self.middle_pane.setItem(row_idx, col_idx, item)
        else:
            self.middle_pane.clear()
    
    def process_data(self):
        selected_items = self.left_pane.selectedItems()
        if not selected_items:
            print("No file selected.")
            return
        item = selected_items[0]
        file_path = item.toolTip()
        
        # Check if the metadata for this file has been parsed
        if file_path not in self.metadata_dict:
            print("Metadata not parsed for this file. Parsing now...")
            self.parse_and_display_metadata_for_file(file_path)
        
        # Get processing parameters
        try:
            ppm_tolerance = float(self.ppm_input.text())
        except ValueError:
            ppm_tolerance = 50.0
            self.ppm_input.setText("50")
        
        try:
            top_x_peaks = int(self.topx_input.text())
        except ValueError:
            top_x_peaks = 30
            self.topx_input.setText("30")
        
        offset_text = self.offset_input.text()
        offset_parts = offset_text.replace('to', '').split()
        if len(offset_parts) == 2:
            try:
                offset_start = int(offset_parts[0])
                offset_end = int(offset_parts[1])
            except ValueError:
                offset_start, offset_end = -8, 8
                self.offset_input.setText("-8 to +8")
        else:
            offset_start, offset_end = -8, 8
            self.offset_input.setText("-8 to +8")
        
        # Get the metadata and h5 file for the selected file
        data_dict = self.metadata_dict[file_path]
        metadata = data_dict['metadata']
        data_h5_path = data_dict['data_h5']
        
        # Open the HDF5 file
        with h5py.File(data_h5_path, 'r') as h5file:
            # Step 1: For each cycle number, get mz-intensity data from experiment=1
            cycle_data = {}
            for data in metadata:
                if data['Experiment Number'] == 1:
                    cycle = data['Cycle Number']
                    index = str(data['Index Number'])
                    if index in h5file:
                        mz_array = h5file[index]['m/z'][:]
                        intensity_array = h5file[index]['intensity'][:]
                        cycle_data[cycle] = {
                            'mz': mz_array,
                            'intensity': intensity_array
                        }
            
            # Step 2: Rank the top X peaks
            top_peaks = {}
            for cycle, data in cycle_data.items():
                intensity = data['intensity']
                mz = data['mz']
                if len(intensity) == 0:
                    continue
                indices = intensity.argsort()[-top_x_peaks:][::-1]
                top_peaks[cycle] = {
                    'mz': mz[indices],
                    'intensity': intensity[indices],
                    'rank': list(range(1, len(indices)+1))
                }
            
            # Step 2 continued: Collect selected ion m/z from experiments != 1
            selected_ions = {}
            for data in metadata:
                if data['Experiment Number'] != 1:
                    cycle = data['Cycle Number']
                    selected_mz = data.get('Selected Ion m/z')
                    if selected_mz:
                        if cycle not in selected_ions:
                            selected_ions[cycle] = []
                        selected_ions[cycle].append(selected_mz)
            
            # Step 3 and 4: Perform peak matching
            ms1top_matched_MS2_col = []
            no_match_ms1top_col = []
            for data in metadata:
                if data['Experiment Number'] == 1:
                    cycle = data['Cycle Number']
                    top_peak_data = top_peaks.get(cycle)
                    if not top_peak_data:
                        ms1top_matched_MS2_col.append('')
                        no_match_ms1top_col.append('')
                        continue
                    mz_values = top_peak_data['mz']
                    ranks = top_peak_data['rank']
                    matched = []
                    unmatched = []
                    for mz_value, rank in zip(mz_values, ranks):
                        match_found = False
                        for offset in range(offset_start, offset_end + 1):
                            offset_cycle = cycle + offset
                            ions = selected_ions.get(offset_cycle, [])
                            for ion_mz in ions:
                                ppm_diff = abs(mz_value - ion_mz) / ion_mz * 1e6
                                if ppm_diff <= ppm_tolerance:
                                    match_str = f"{mz_value:.2f}mz_cycle{offset:+}_rank{rank}"
                                    matched.append(match_str)
                                    match_found = True
                                    break
                            if match_found:
                                break
                        if not match_found:
                            unmatched_str = f"{mz_value:.2f}mz_rank{rank}"
                            unmatched.append(unmatched_str)
                    ms1top_matched_MS2_col.append('; '.join(matched))
                    no_match_ms1top_col.append('; '.join(unmatched))
                else:
                    ms1top_matched_MS2_col.append('')
                    no_match_ms1top_col.append('')
            
            # Add new columns to the metadata
            for idx, data_row in enumerate(metadata):
                data_row['ms1top_matched_MS2'] = ms1top_matched_MS2_col[idx]
                data_row['no_match_ms1top'] = no_match_ms1top_col[idx]
            
            # Update headers if new columns are added
            if 'ms1top_matched_MS2' not in self.headers:
                self.headers.extend(['ms1top_matched_MS2', 'no_match_ms1top'])
            
            # Refresh the display
            self.display_metadata_for_file(file_path)

            # Store cycle_data and rt_cycle_map for later use
            self.metadata_dict[file_path]['cycle_data'] = cycle_data
            self.metadata_dict[file_path]['rt_cycle_map'] = {}
            for data in metadata:
                if data['Experiment Number'] == 1:
                    rt = data['Scan Start Time (min)']
                    cycle = data['Cycle Number']
                    self.metadata_dict[file_path]['rt_cycle_map'][rt] = cycle            

    
    def plot_no_ms2_for_top_peaks(self):
        selected_items = self.left_pane.selectedItems()
        if not selected_items:
            print("No file selected.")
            return
        item = selected_items[0]
        file_path = item.toolTip()
        self.current_file_path = file_path

        # Get the metadata for the selected file
        data = self.metadata_dict.get(file_path)
        if not data:
            print("Metadata not available for this file.")
            return
        metadata = data['metadata']
        
        # Check if 'no_match_ms1top' column is available
        if 'no_match_ms1top' not in metadata[0]:
            print("Please process data before plotting.")
            return
        
        # Extract retention times and mz values
        spots = []
        
        for row in metadata:
            rt = row.get('Scan Start Time (min)')
            no_match_str = row.get('no_match_ms1top')
            if rt is None or not no_match_str:
                continue
            # Parse the no_match_ms1top string
            entries = no_match_str.split('; ')
            for entry in entries:
                if not entry:
                    continue
                # Entry format: "303.23mz_rank6"
                match = re.match(r'([0-9.]+)mz_rank([0-9]+)', entry)
                if match:
                    mz = float(match.group(1))
                    rank = int(match.group(2))
                    spots.append({'pos': (rt, mz), 'data': {'rank': rank, 'mz': mz, 'rt': rt}})
        
        if not spots:
            print("No data to plot.")
            return
        
        # Clear previous plot
        self.plot_widget.clear()

        # Set the background color
        self.plot_widget.setBackground('black')  # Change to your preferred color        
        
        # Create scatter plot with hover tooltips
        self.scatter = pg.ScatterPlotItem(
            size=7,
            pen=pg.mkPen(None),
            brush=pg.mkBrush(255, 255, 255, 180),
            hoverable=True,
            hoverSymbol='o',
            hoverSize=15,
            hoverBrush=pg.mkBrush('orange'),
        )
        self.scatter.setData(spots)
        self.scatter.hoverEvent = self.scatterHoverEvent
        self.plot_widget.addItem(self.scatter)
        
        self.plot_widget.setLabel('bottom', 'Retention Time (min)')
        self.plot_widget.setLabel('left', 'm/z')
        self.plot_widget.setTitle('No MS2 for Top Peaks')
        self.plot_widget.showGrid(x=True, y=True)

        # Connect the click event handler
        self.scatter.sigClicked.connect(self.on_scatter_clicked)    

    def scatterHoverEvent(self, ev):
        if ev.isExit():
            QtWidgets.QToolTip.hideText()
            return
        pos = ev.pos()
        points = self.scatter.pointsAt(pos)
        if len(points) > 0:
            point = points[0]
            data = point.data()
            rank = data['rank']
            mz = data['mz']
            rt = data['rt']
            screen_pos = ev.screenPos()
            QtWidgets.QToolTip.showText(screen_pos.toPoint(),
                                        f"Rank: {rank}\nRT: {rt:.2f} min\nm/z: {mz:.4f}",
                                        widget=self.plot_widget)
        else:
            QtWidgets.QToolTip.hideText()

    def on_scatter_clicked(self, scatter_plot_item, points):
        if not points:
            return
        point = points[0]
        data = point.data()
        rank = data['rank']
        mz_clicked = data['mz']
        rt = data['rt']
        file_path = self.current_file_path
        # Retrieve the rt_cycle_map
        rt_cycle_map = self.metadata_dict[file_path]['rt_cycle_map']
        # Find the closest rt to the clicked point's rt
        rts = np.array(list(rt_cycle_map.keys()))
        idx = (np.abs(rts - rt)).argmin()
        rt_closest = rts[idx]
        cycle_number = rt_cycle_map[rt_closest]
        # Retrieve the cycle_data
        cycle_data = self.metadata_dict[file_path]['cycle_data']
        spectrum_data = cycle_data.get(cycle_number)
        if spectrum_data is None:
            print("Spectrum data not found for cycle number:", cycle_number)
            return
        mz_array = spectrum_data['mz']
        intensity_array = spectrum_data['intensity']
        # Plot the mass spectrum
        self.plot_spectrum(mz_array, intensity_array, mz_clicked, rt)  # Pass rt as an argument

    def plot_spectrum(self, mz_array, intensity_array, mz_clicked, rt):  # Add rt parameter
        # Clear previous plot
        self.spectrum_plot_widget.clear()
        # Create arrays for x and y values to plot vertical lines
        x_values = np.repeat(mz_array, 3)
        y_values = np.zeros(len(x_values))
        y_values[1::3] = intensity_array
        y_values[2::3] = np.nan  # To break the lines between peaks
        # Plot the vertical lines (mass spectrum)
        self.spectrum_plot_widget.plot(x_values, y_values, pen='w')
        # Highlight the clicked m/z value
        idx = (np.abs(mz_array - mz_clicked)).argmin()
        mz_highlight = mz_array[idx]
        intensity_highlight = intensity_array[idx]
        # Plot the highlighted line in red
        self.spectrum_plot_widget.plot([mz_highlight, mz_highlight], [0, intensity_highlight], pen=pg.mkPen('r', width=2))
        self.spectrum_plot_widget.setLabel('bottom', 'm/z')
        self.spectrum_plot_widget.setLabel('left', 'Intensity')
        self.spectrum_plot_widget.setTitle('Mass Spectrum at RT {:.2f} min'.format(rt))  # Use rt here
        


    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
