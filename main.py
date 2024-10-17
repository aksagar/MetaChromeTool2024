"""
Author: Md Abdul Kader Sagar
Email: sagarm2@nih.gov
Institute: National Cancer Institute/NIH

This code is designed for analyzing metaphase chromosomes using the Napari platform.
It facilitates the visualization and segmentation of chromosome images, enabling users 
to efficiently assess chromosome structures and perform quantitative analysis.
The code integrates tools for detecting centromeres and measuring CENP-A levels 
within metaphase chromosome regions, enhancing the accuracy of chromosome analysis.
"""


import os
import napari
import numpy as np
import pandas as pd
from napari.utils.notifications import show_info
from magicgui import magicgui
from qtpy.QtWidgets import QFileDialog, QVBoxLayout, QWidget, QLineEdit, QLabel, QHBoxLayout, QCheckBox, QPushButton
from superqt import QLabeledSlider
from skimage.draw import line
from image_processor import ImageProcessor
from batch_processor import BatchProcessor

# Initialize viewer and processor
viewer = napari.Viewer()
processor = ImageProcessor()
batch_processor = None
current_folder_path = ""  # Global variable to store the current folder path
segment_done = False

# Global list to store images
images = [None, None, None]  # To hold DAPI, DNA-FISH, and CENPC images

# Flags to ensure each process runs once
detect_dna_fish_done = False
detect_cenpc_done = False

# Create a QWidget to hold the sliders and the detect buttons
class ControlWidgetDNAFISH(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)  # Range is 0 to 100 to work with integer values
        self.slider.setValue(40)  # Default value
        self.slider.setSingleStep(1)  # Step size
        self.slider.valueChanged.connect(self.reset_dna_fish_flag)  # Connect slider change to reset function
        self.layout.addWidget(self.slider)

        self.detect_dna_fish_spots_button = detect_dna_fish_spots.native
        self.layout.addWidget(self.detect_dna_fish_spots_button)

        self.setLayout(self.layout)

    def reset_dna_fish_flag(self):
        global detect_dna_fish_done
        detect_dna_fish_done = False
        show_info("Threshold slider changed, reset spot detection flag for DNA-FISH channel")

class ControlWidgetCENPC(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)  # Range is 0 to 100 to work with integer values
        self.slider.setValue(40)  # Default value
        self.slider.setSingleStep(1)  # Step size
        self.slider.valueChanged.connect(self.reset_cenpc_flag)  # Connect slider change to reset function
        self.layout.addWidget(self.slider)

        self.detect_cenpc_spots_button = detect_cenpc_spots.native
        self.layout.addWidget(self.detect_cenpc_spots_button)

        self.setLayout(self.layout)

    def reset_cenpc_flag(self):
        global detect_cenpc_done
        detect_cenpc_done = False
        show_info("Threshold slider changed, reset spot detection flag for CENPC channel")

# Create a QWidget to hold the channel identifier text boxes
class ChannelIdentifiers(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()

        self.dapi_label = QLabel("DAPI Channel Identifier:")
        self.dapi_text = QLineEdit()
        self.layout.addWidget(self.dapi_label)
        self.layout.addWidget(self.dapi_text)

        self.dna_fish_label = QLabel("DNA-FISH Channel Identifier:")
        self.dna_fish_text = QLineEdit()
        self.layout.addWidget(self.dna_fish_label)
        self.layout.addWidget(self.dna_fish_text)

        self.cenpc_label = QLabel("CENPC Channel Identifier:")
        self.cenpc_text = QLineEdit()
        self.layout.addWidget(self.cenpc_label)
        self.layout.addWidget(self.cenpc_text)

        self.setLayout(self.layout)

# Add a checkbox beside the segment DAPI button
@magicgui(call_button="Segment DAPI Image")
def segment_image():
    global segment_done, images
    if segment_done:
        show_info("Segmentation has already been done.")
        return

    if images[0] is not None:
        try:
            masks = processor.segment_image(images[0])
            viewer.add_labels(masks, name="Cellpose Segmented DAPI")
            show_info("Segmented DAPI image using Cellpose")
            segment_done = True
        except Exception as e:
            show_info(f"Error segmenting image: {e}")

def merge_chromosomes():
    try:
        shapes_layer = viewer.layers['Shapes']
        for coords in shapes_layer.data:
            updated_nuclei = processor.merge_nuclei_with_line(coords)
        processor.nuclei_merged = updated_nuclei.copy()
        viewer.add_labels(processor.nuclei_merged, name="Merged Chromosomes")
        show_info("Merged chromosomes based on drawn lines.")
    except Exception as e:
        show_info(f"Error merging chromosomes: {e}")

def remove_chromosome():
    try:
        shapes_layer = viewer.layers['Shapes']
        for coords in shapes_layer.data:
            updated_nuclei = processor.remove_nuclei_with_line(coords)
        processor.nuclei_merged = updated_nuclei.copy()
        viewer.add_labels(processor.nuclei_merged, name="Updated Chromosomes")
        show_info("Removed chromosomes based on drawn lines.")
    except Exception as e:
        show_info(f"Error removing chromosomes: {e}")


class SegmentDAPIWidget(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()  # Changed to QVBoxLayout to stack elements vertically

        # Segment button and checkbox
        self.segment_layout = QHBoxLayout()
        self.segment_button = segment_image.native
        self.checkbox = QCheckBox("Skip Segmentation")
        self.segment_layout.addWidget(self.segment_button)
        self.segment_layout.addWidget(self.checkbox)
        self.layout.addLayout(self.segment_layout)

        # Buttons for merging and removing chromosomes
        self.buttons_layout = QHBoxLayout()
        self.merge_button = QPushButton("Merge Chromosomes")
        self.merge_button.clicked.connect(merge_chromosomes)
        self.buttons_layout.addWidget(self.merge_button)

        self.remove_button = QPushButton("Remove Chromosome")
        self.remove_button.clicked.connect(remove_chromosome)
        self.buttons_layout.addWidget(self.remove_button)
        
        self.layout.addLayout(self.buttons_layout)

        self.setLayout(self.layout)

    def is_checked(self):
        return self.checkbox.isChecked()


segment_dapi_widget = SegmentDAPIWidget()

@magicgui(call_button="Load Images")
def load_images():
    global segment_done, detect_dna_fish_done, detect_cenpc_done, current_folder_path, images
    folder_path = QFileDialog.getExistingDirectory(caption='Select Image Folder')
    if folder_path:
        current_folder_path = folder_path
        try:
            dapi_id = channel_identifiers.dapi_text.text()
            dna_fish_id = channel_identifiers.dna_fish_text.text()
            cenpc_id = channel_identifiers.cenpc_text.text()
            images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
            
            # Store existing shapes layer data
            shapes_data = []
            for layer in viewer.layers:
                if isinstance(layer, napari.layers.Shapes):
                    shapes_data.append(layer.data)
            
            viewer.layers.clear()
            
            # Restore shapes layer
            shapes_layer = viewer.add_shapes(name='Shapes', edge_color='yellow', edge_width=2)
            if shapes_data:
                shapes_layer.data = shapes_data[0]
            
            if segment_dapi_widget.is_checked():
                images.insert(0, None)  # Add None for DAPI to images
                viewer.add_image(images[1], name='DNA-FISH')
                viewer.add_image(images[2], name='CENPC')
            else:
                viewer.add_image(images[0], name='DAPI')
                viewer.add_image(images[1], name='DNA-FISH')
                viewer.add_image(images[2], name='CENPC')
            segment_done = False
            detect_dna_fish_done = False
            detect_cenpc_done = False
            show_info(f"Loaded {len(images)} images from: {folder_path}")
        except Exception as e:
            show_info(f"Error loading images: {e}")

@magicgui(call_button="Detect DNA-FISH Spots")
def detect_dna_fish_spots():
    global detect_dna_fish_done, images
    if detect_dna_fish_done:
        show_info("Spot detection for DNA-FISH has already been done.")
        return

    threshold = control_widget_dna_fish.slider.value() / 100  # Convert slider value to 0-1 range
    if images[1] is not None:
        try:
            if segment_dapi_widget.is_checked():
                spots, labels = processor.detect_spots_no_segmentation(images[1], threshold)
            else:
                spots = processor.detect_spots(images[1], 'DNA-FISH', threshold)
            viewer.add_labels(spots, name="Spots in DNA-FISH")
            show_info(f"Detected and labeled spots in DNA-FISH image with threshold {threshold}")
            detect_dna_fish_done = True
        except Exception as e:
            show_info(f"Error detecting spots: {e}")

@magicgui(call_button="Detect CENPC Spots")
def detect_cenpc_spots():
    global detect_cenpc_done, images
    if detect_cenpc_done:
        show_info("Spot detection for CENPC has already been done.")
        return

    threshold = control_widget_cenpc.slider.value() / 100  # Convert slider value to 0-1 range
    if images[2] is not None:
        try:
            if segment_dapi_widget.is_checked():
                spots, labels = processor.detect_spots_no_segmentation(images[2], threshold)
            else:
                spots = processor.detect_spots(images[2], 'CENPC', threshold)
            viewer.add_labels(2 * spots, name="Spots in CENPC")
            show_info(f"Detected and labeled spots in CENPC image with threshold {threshold}")
            detect_cenpc_done = True
        except Exception as e:
            show_info(f"Error detecting spots: {e}")

@magicgui(call_button="Find Common")
def find_common():
    try:
        threshold_dna_fish = control_widget_dna_fish.slider.value() / 100  # Convert slider value to 0-1 range
        threshold_cenpc = control_widget_cenpc.slider.value() / 100  # Convert slider value to 0-1 range

        if segment_dapi_widget.is_checked():
            show_info("Skipping find common due to checkbox selection.")
            return

        common_nuclei = processor.find_common(threshold_dna_fish, threshold_cenpc)
        if common_nuclei is None:
            show_info("No common labels found.")
            return
        viewer.add_labels(common_nuclei, name="Matched Nuclei")
        show_info("Found common labels and updated the view.")
    except Exception as e:
        show_info(f"Error finding common labels: {e}")

@magicgui(call_button="Get Intensity at CENPC Location")

def get_intensity_at_cenpc_location():
    try:
        if segment_dapi_widget.is_checked():
            show_info("Calculating intensity at all DNA-FISH locations without segmentation.")
            df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
        else:
            df_with_cenpc_inten = processor.gen_intensity_from_df(processor.img_cenpc, processor.df_centroid_dna_fish)
        print(df_with_cenpc_inten)
        # Save the DataFrame in the current folder using the folder's name
        folder_name = os.path.basename(current_folder_path)
        save_path = os.path.join(current_folder_path, f"{folder_name}_intensity.csv")
        df_with_cenpc_inten.to_csv(save_path, index=False)
        show_info(f"DataFrame saved to {save_path}")
    except Exception as e:
        show_info(f"Error getting intensity: {e}")

@magicgui(call_button="Run All")
def run_all():
    global images
    try:
        threshold_dna_fish = control_widget_dna_fish.slider.value() / 100  # Convert slider value to 0-1 range
        threshold_cenpc = control_widget_cenpc.slider.value() / 100  # Convert slider value to 0-1 range

        if all(img is not None for img in images if img is not None):
            if segment_dapi_widget.is_checked():
                spots_dna_fish, labels_dna_fish = processor.detect_spots_no_segmentation(images[1], threshold_dna_fish)
                #spots_cenpc, labels_cenpc = processor.detect_spots_no_segmentation(images[2], threshold_cenpc)
                
                df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
            else:
                masks = processor.segment_image(images[0])
                viewer.add_labels(masks, name="Cellpose Segmented")
                spots_dna_fish = processor.detect_spots(images[1], 'DNA-FISH', threshold_dna_fish)
                spots_cenpc = processor.detect_spots(images[2], 'CENPC', threshold_cenpc)

                common_nuclei = processor.find_common(threshold_dna_fish, threshold_cenpc)
                if common_nuclei is None:
                    show_info("No common labels found.")
                    return
                viewer.add_labels(common_nuclei, name="Matched Nuclei")

                df_with_cenpc_inten = processor.gen_intensity_from_df(processor.img_cenpc, processor.df_centroid_dna_fish)

            viewer.add_labels(spots_dna_fish, name="Spots in DNA-FISH")
            viewer.add_labels(2 * spots_cenpc, name="Spots in CENPC")

            # Save the DataFrame in the current folder using the folder's name
            folder_name = os.path.basename(current_folder_path)
            save_path = os.path.join(current_folder_path, f"{folder_name}_intensity.csv")
            df_with_cenpc_inten.to_csv(save_path, index=False)
            show_info(f"DataFrame saved to {save_path}")

            print(df_with_cenpc_inten)
            show_info("Run all processing completed.")
        else:
            show_info("Ensure that all images are loaded.")
    except Exception as e:
        show_info(f"Error during run all processing: {e}")

@magicgui(call_button="Batch Processing")
def batch_processing():
    root_folder = QFileDialog.getExistingDirectory(caption='Select Root Folder')
    if root_folder:
        dapi_id = channel_identifiers.dapi_text.text()
        dna_fish_id = channel_identifiers.dna_fish_text.text()
        cenpc_id = channel_identifiers.cenpc_text.text()
        batch_processor.batch_processing(root_folder, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())

# Add the widgets to the viewer
control_widget_dna_fish = ControlWidgetDNAFISH()
control_widget_cenpc = ControlWidgetCENPC()
channel_identifiers = ChannelIdentifiers()
batch_processor = BatchProcessor(processor, control_widget_dna_fish, control_widget_cenpc)

viewer.window.add_dock_widget(channel_identifiers, area='right', name='Channel Identifiers')
viewer.window.add_dock_widget(load_images, area='right', name='')
viewer.window.add_dock_widget(segment_dapi_widget, area='right', name='Segment DAPI Control')
viewer.window.add_dock_widget(control_widget_dna_fish, area='right', name='Detect DNA-FISH Spot Control')
viewer.window.add_dock_widget(control_widget_cenpc, area='right', name='Detect CENPC Spot Control')
viewer.window.add_dock_widget(find_common, area='right', name='')
viewer.window.add_dock_widget(get_intensity_at_cenpc_location, area='right', name='')

# Create Run All button and set its color to green
run_all_button = run_all.native
run_all_button.setStyleSheet("background-color: green; color: white;")
viewer.window.add_dock_widget(run_all_button, area='right', name='Run All')

batch_processing_button = batch_processing.native
batch_processing_button.setStyleSheet("background-color: blue; color: white;")
viewer.window.add_dock_widget(batch_processing_button, area='right', name='Batch Processing')

# Create a shapes layer for drawing lines
shapes_layer = viewer.add_shapes(name='Shapes', edge_color='yellow', edge_width=2)

# Add a callback to handle line drawing for merging/removing nuclei
@shapes_layer.mouse_drag_callbacks.append
def shapes_layer_callback(viewer, event):
    if event.type == 'mouse_press' and 'Shift' in event.modifiers:
        yield
        while event.type == 'mouse_move':
            yield
        # Get the line coordinates from the shapes layer
        line_coords = np.array(shapes_layer.data[-1])
        if 'merge' in event.modifiers:
            processor.merge_nuclei_with_line(line_coords)
        elif 'remove' in event.modifiers:
            processor.remove_nuclei_with_line(line_coords)

# Start the napari event loop
napari.run()
