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
from qtpy.QtWidgets import QFileDialog, QVBoxLayout, QWidget, QLineEdit, QLabel, QHBoxLayout, QCheckBox, QPushButton, QListWidget, QListWidgetItem
from superqt import QLabeledSlider
from skimage.draw import line
from .image_processor import ImageProcessor
from .batch_processor import BatchProcessor
from .segmentation_postprocessing import SegmentationPostprocessing
from datetime import datetime

from PIL import Image
import matplotlib.pyplot as plt

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


# Apply to control widget buttons
# Apply to control widget buttons
class ControlWidgetDNAFISH(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        
        # Slider
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)
        self.slider.setValue(40)
        self.slider.setSingleStep(1)
        self.slider.valueChanged.connect(self.reset_dna_fish_flag)
        self.layout.addWidget(self.slider)

        # Detect button
        self.detect_dna_fish_spots_button = detect_dna_fish_spots.native
        self.detect_dna_fish_spots_button.setStyleSheet(BUTTON_STYLE)
        self.layout.addWidget(self.detect_dna_fish_spots_button)

        # Create horizontal layout for delete and save buttons
        self.button_layout = QHBoxLayout()
        
        # Delete button
        self.delete_dna_fish_spots_button = delete_dna_fish_spots.native
        self.delete_dna_fish_spots_button.setStyleSheet("""
            QPushButton {
                background-color: #ff6b6b;
                color: white;
                border-radius: 5px;
                padding: 5px;
                min-height: 25px;
            }
            QPushButton:hover {
                background-color: #ff8787;
            }
            QPushButton:pressed {
                background-color: #fa5252;
            }
        """)
        self.button_layout.addWidget(self.delete_dna_fish_spots_button)

        # Save button
        self.save_dna_fish_spots_button = save_dna_fish_spots.native
        self.save_dna_fish_spots_button.setStyleSheet("""
            QPushButton {
                background-color: #40c057;
                color: white;
                border-radius: 5px;
                padding: 5px;
                min-height: 25px;
            }
            QPushButton:hover {
                background-color: #51cf66;
            }
            QPushButton:pressed {
                background-color: #37b24d;
            }
        """)
        self.button_layout.addWidget(self.save_dna_fish_spots_button)

        # Add button layout to main layout
        self.layout.addLayout(self.button_layout)
        self.setLayout(self.layout)

    def reset_dna_fish_flag(self):
        global detect_dna_fish_done
        detect_dna_fish_done = False
        show_info("Threshold slider changed, reset spot detection flag for Channel 1")


class ControlWidgetDNAFISH1(QWidget):
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
        show_info("Threshold slider changed, reset spot detection flag for Channel 1")


class ControlWidgetCENPC(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        
        # Slider
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)
        self.slider.setValue(40)
        self.slider.setSingleStep(1)
        self.slider.valueChanged.connect(self.reset_cenpc_flag)
        self.layout.addWidget(self.slider)

        # Detect button
        self.detect_cenpc_spots_button = detect_cenpc_spots.native
        self.detect_cenpc_spots_button.setStyleSheet(BUTTON_STYLE)
        self.layout.addWidget(self.detect_cenpc_spots_button)

        # Create horizontal layout for delete and save buttons
        self.button_layout = QHBoxLayout()
        
        # Delete button
        self.delete_cenpc_spots_button = delete_cenpc_spots.native
        self.delete_cenpc_spots_button.setStyleSheet("""
            QPushButton {
                background-color: #ff6b6b;
                color: white;
                border-radius: 5px;
                padding: 5px;
                min-height: 25px;
            }
            QPushButton:hover {
                background-color: #ff8787;
            }
            QPushButton:pressed {
                background-color: #fa5252;
            }
        """)
        self.button_layout.addWidget(self.delete_cenpc_spots_button)

        # Save button
        self.save_cenpc_spots_button = save_cenpc_spots.native
        self.save_cenpc_spots_button.setStyleSheet("""
            QPushButton {
                background-color: #40c057;
                color: white;
                border-radius: 5px;
                padding: 5px;
                min-height: 25px;
            }
            QPushButton:hover {
                background-color: #51cf66;
            }
            QPushButton:pressed {
                background-color: #37b24d;
            }
        """)
        self.button_layout.addWidget(self.save_cenpc_spots_button)

        # Add button layout to main layout
        self.layout.addLayout(self.button_layout)
        self.setLayout(self.layout)

    def reset_cenpc_flag(self):
        global detect_cenpc_done
        detect_cenpc_done = False
        show_info("Threshold slider changed, reset spot detection flag for Channel 2")

class ControlWidgetCENPC1(QWidget):
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

from qtpy.QtCore import Qt  # Add this import at the top with other imports

class ChromosomeCountWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QHBoxLayout()
        self.label = QLabel("Chromosome Count: --")  # Start with -- instead of 0
        self.label.setStyleSheet("""
            QLabel {
                background-color: rgba(0, 0, 0, 0.7);
                color: white;
                padding: 5px;
                border-radius: 3px;
                min-width: 150px;  # Set minimum width
                max-width: 150px;  # Set maximum width
                min-height: 25px;  # Set minimum height
                max-height: 25px;  # Set maximum height
            }
        """)
        self.label.setAlignment(Qt.AlignCenter)  # Center the text
        self.layout.addWidget(self.label)
        self.setLayout(self.layout)

    def update_count(self, count):
        self.label.setText(f"Chromosome Count: {count}")
# Add these new widget classes after ChromosomeCountWidget
class SpotCountWidget(QWidget):
    def __init__(self, spot_type):
        super().__init__()
        self.layout = QHBoxLayout()
        self.spot_type = spot_type
        self.label = QLabel(f"{spot_type} Spot Count: --")  # Start with -- instead of 0
        self.label.setStyleSheet("""
            QLabel {
                background-color: rgba(0, 0, 0, 0.7);
                color: white;
                padding: 5px;
                border-radius: 3px;
                min-width: 150px;  # Set minimum width
                max-width: 150px;  # Set maximum width
                min-height: 25px;  # Set minimum height
                max-height: 25px;  # Set maximum height
            }
        """)
        self.label.setAlignment(Qt.AlignCenter)  # Center the text
        self.layout.addWidget(self.label)
        self.setLayout(self.layout)

    def update_count(self, count):
        self.label.setText(f"{self.spot_type} Spot Count: {count}")
# Create a QWidget to hold the channel identifier text boxes

import json
import os


def save_channel_settings(identifiers):
    """Save channel identifiers to settings file."""
    settings = {
        'segmentation_channel': identifiers.dapi_text.text(),
        'channel1': identifiers.dna_fish_text.text(),
        'channel2': identifiers.cenpc_text.text()
    }
    
    # Save to user's home directory or application directory
    settings_dir = os.path.expanduser('~/.napari_chromosome')
    os.makedirs(settings_dir, exist_ok=True)
    settings_file = os.path.join(settings_dir, 'channel_settings.json')
    
    with open(settings_file, 'w') as f:
        json.dump(settings, f)

def load_channel_settings():
    """Load channel identifiers from settings file."""
    settings_file = os.path.join(os.path.expanduser('~/.napari_chromosome'), 'channel_settings.json')
    
    # Default values
    settings = {
        'segmentation_channel': '435',
        'channel1': '525',
        'channel2': '679'
    }
    
    if os.path.exists(settings_file):
        try:
            with open(settings_file, 'r') as f:
                loaded_settings = json.load(f)
                settings.update(loaded_settings)
        except Exception as e:
            show_info(f"Error loading settings: {str(e)}")
    
    return settings

# Modify the ChannelIdentifiers class to use settings
class ChannelIdentifiers(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.layout.setSpacing(1)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # Load saved settings
        settings = load_channel_settings()
        
        # Segmentation Channel
        seg_layout = QHBoxLayout()
        seg_layout.setContentsMargins(0, 0, 0, 0)
        self.dapi_text = QLineEdit()
        self.dapi_text.setText(settings['segmentation_channel'])
        self.dapi_text.setFixedWidth(50)
        self.dapi_label = QLabel("Segmentation Channel Identifier")
        seg_layout.addWidget(self.dapi_text)
        seg_layout.addWidget(self.dapi_label)
        seg_layout.addStretch()
        self.layout.addLayout(seg_layout)

        # Channel 1
        ch1_layout = QHBoxLayout()
        ch1_layout.setContentsMargins(0, 0, 0, 0)
        self.dna_fish_text = QLineEdit()
        self.dna_fish_text.setText(settings['channel1'])
        self.dna_fish_text.setFixedWidth(50)
        self.dna_fish_label = QLabel("Channel 1 Identifier")
        ch1_layout.addWidget(self.dna_fish_text)
        ch1_layout.addWidget(self.dna_fish_label)
        ch1_layout.addStretch()
        self.layout.addLayout(ch1_layout)

        # Channel 2
        ch2_layout = QHBoxLayout()
        ch2_layout.setContentsMargins(0, 0, 0, 0)
        self.cenpc_text = QLineEdit()
        self.cenpc_text.setText(settings['channel2'])
        self.cenpc_text.setFixedWidth(50)
        self.cenpc_label = QLabel("Channel 2 Identifier")
        ch2_layout.addWidget(self.cenpc_text)
        ch2_layout.addWidget(self.cenpc_label)
        ch2_layout.addStretch()
        self.layout.addLayout(ch2_layout)

        self.setLayout(self.layout)



# Add a checkbox beside the segment DAPI button

def segment_image_BU1():
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

@magicgui(call_button="Segment (DAPI) Image")
def segment_image():
    global segment_done, images


    if images[0] is not None:
        try:
            # Pass current folder path to save results
            masks = processor.segment_image(images[0], save_dir=current_folder_path)
            viewer.add_labels(masks, name="Cellpose Segmented (DAPI)")
            
            # Update chromosome count
            unique_chromosomes = len(np.unique(masks)) - 1
            chromosome_counter.update_count(unique_chromosomes)
            
            show_info("Segmented DAPI image using Cellpose")
            segment_done = True
        except Exception as e:
            show_info(f"Error segmenting image: {e}")




from skimage import draw, morphology
import scipy.ndimage as ndi
import numpy as np


BUTTON_STYLE = """
    QPushButton {
        background-color: #4a4a4a;
        color: white;
        border-radius: 5px;
        padding: 5px;
        min-height: 25px;
    }
    QPushButton:hover {
        background-color: #5a5a5a;
    }
    QPushButton:pressed {
        background-color: #3a3a3a;
    }
"""
folder_list_widget = QListWidget()
folder_list_widget.setMinimumHeight(200)  # Make it taller
#folder_list_widget.setMinimumWidth(200)   # Make it wider
folder_list_widget.setStyleSheet("""
    QListWidget {
        background-color: #2d2d2d;
        color: white;
        border: 1px solid #3d3d3d;
        border-radius: 5px;
    }
    QListWidget::item {
        padding: 5px;
    }
    QListWidget::item:selected {
        background-color: #4a4a4a;
    }
    QListWidget::item:hover {
        background-color: #3a3a3a;
    }
""")


# Apply to SegmentDAPIWidget buttons
class SegmentDAPIWidget(QWidget):
    def __init__(self, postprocessing=None):
        super().__init__()
        self.postprocessing = postprocessing
        self.layout = QVBoxLayout()

        # Segment button and checkbox
        self.segment_layout = QHBoxLayout()
        self.segment_button = segment_image.native
        self.segment_button.setStyleSheet(BUTTON_STYLE)
        self.checkbox = QCheckBox("Skip Segmentation")
        self.segment_layout.addWidget(self.segment_button)
        self.segment_layout.addWidget(self.checkbox)
        self.layout.addLayout(self.segment_layout)

        # Buttons for merging, removing, and splitting chromosomes
        self.buttons_layout = QHBoxLayout()
        
        self.merge_button = QPushButton("Merge")
        self.merge_button.setStyleSheet(BUTTON_STYLE)
        self.merge_button.clicked.connect(self.postprocessing.merge_chromosomes)
        self.buttons_layout.addWidget(self.merge_button)

        self.remove_button = QPushButton("Remove")
        self.remove_button.setStyleSheet(BUTTON_STYLE)
        self.remove_button.clicked.connect(self.postprocessing.remove_chromosome)
        self.buttons_layout.addWidget(self.remove_button)

        self.split_button = QPushButton("Split")
        self.split_button.setStyleSheet(BUTTON_STYLE)
        self.split_button.clicked.connect(self.postprocessing.split_chromosome)
        self.buttons_layout.addWidget(self.split_button)

        self.save_button = QPushButton("Save")
        self.save_button.setStyleSheet(BUTTON_STYLE)
        self.save_button.clicked.connect(self.postprocessing.save_segmentation)  # Changed this line
        self.buttons_layout.addWidget(self.save_button)
        
        self.layout.addLayout(self.buttons_layout)
        self.setLayout(self.layout)

    def is_checked(self):
        return self.checkbox.isChecked()





@magicgui(call_button="Load Images")
def load_images():
    save_channel_settings(channel_identifiers)

    selected_folder = QFileDialog.getExistingDirectory(caption='Select Folder')
    if selected_folder:
        try:
            # Clear previous items
            folder_list_widget.clear()
            
            # If selected folder contains images, go up one level
            if any(f.lower().endswith(('.tif', '.tiff')) for f in os.listdir(selected_folder)):
                root_folder = os.path.dirname(selected_folder)
            else:
                root_folder = selected_folder
            
            # Store full paths but display only folder names
            for folder_name in os.listdir(root_folder):
                # Skip the "Exported Segmentation" folder
                if folder_name == "Exported Segmentation":
                    continue
                    
                folder_path = os.path.join(root_folder, folder_name)
                if os.path.isdir(folder_path):
                    item = QListWidgetItem(os.path.basename(folder_path))
                    # Store full path as item data
                    item.setData(Qt.UserRole, folder_path)
                    folder_list_widget.addItem(item)
            
            # Select the originally chosen folder if it's in the list
            if root_folder != selected_folder:
                for i in range(folder_list_widget.count()):
                    item = folder_list_widget.item(i)
                    if item.data(Qt.UserRole) == selected_folder:
                        folder_list_widget.setCurrentItem(item)
                        break
            
            show_info(f"Found {folder_list_widget.count()} folders")
        except Exception as e:
            show_info(f"Error loading folders: {str(e)}")



def load_images_BU2():
    global segment_done, detect_dna_fish_done, detect_cenpc_done, current_folder_path, images
    folder_path = QFileDialog.getExistingDirectory(caption='Select Image Folder')
    if folder_path:
        current_folder_path = folder_path
        postprocessing.set_current_folder(current_folder_path) 
        try:
            dapi_id = channel_identifiers.dapi_text.text()
            dna_fish_id = channel_identifiers.dna_fish_text.text()
            cenpc_id = channel_identifiers.cenpc_text.text()
            
            images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
            
            viewer.layers.clear()
            
            if segment_dapi_widget.is_checked():
                images.insert(0, None)
                viewer.add_image(images[1], name='DNA-FISH')
                viewer.add_image(images[2], name='CENPC')
            else:
                viewer.add_image(images[0], name='DAPI')
                viewer.add_image(images[1], name='DNA-FISH')
                viewer.add_image(images[2], name='CENPC')
            
            # Check for and load intermediate results
            intermediate_path = os.path.join(folder_path, "intermediate_results")
            if os.path.exists(intermediate_path):
                # Load segmentation
                seg_file = os.path.join(intermediate_path, "segmentation.npy")
                if os.path.exists(seg_file):
                    processor.nuclei = np.load(seg_file)
                    viewer.add_labels(processor.nuclei, name="Cellpose Segmented")
                    segment_done = True
                    unique_chromosomes = len(np.unique(processor.nuclei)) - 1
                    chromosome_counter.update_count(unique_chromosomes)
                    show_info("Loaded existing segmentation")
                
                # Load DNA-FISH spots
                dna_fish_file = os.path.join(intermediate_path, "dna_fish_spots.npy")
                dna_fish_centroids_file = os.path.join(intermediate_path, "dna_fish_centroids.npy")
                if os.path.exists(dna_fish_file) and os.path.exists(dna_fish_centroids_file):
                    processor.labels_dna_fish = np.load(dna_fish_file)
                    processor.dna_fish_centroids = np.load(dna_fish_centroids_file)
                    if processor.dna_fish_centroids is not None and len(processor.dna_fish_centroids) > 0:
                        squares = [
                            [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                            for x, y in processor.dna_fish_centroids
                        ]
                        viewer.add_shapes(
                            squares,
                            shape_type='polygon',
                            edge_color="yellow",
                            face_color=[1, 1, 0, 0.2],
                            edge_width=2,
                            name="Channel 1 Spots",
                            opacity=0.8

                        )
                        dna_fish_counter.update_count(len(processor.dna_fish_centroids))
                    detect_dna_fish_done = True
                
                # Load CENPC spots
                cenpc_file = os.path.join(intermediate_path, "cenpc_spots.npy")
                cenpc_centroids_file = os.path.join(intermediate_path, "cenpc_centroids.npy")
                if os.path.exists(cenpc_file) and os.path.exists(cenpc_centroids_file):
                    processor.labels_cenpc = np.load(cenpc_file)
                    processor.cenpc_centroids = np.load(cenpc_centroids_file)
                    if processor.cenpc_centroids is not None and len(processor.cenpc_centroids) > 0:
                        squares = [
                            [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                            for x, y in processor.cenpc_centroids
                        ]
                        viewer.add_shapes(
                            squares,
                            shape_type='polygon',
                            edge_color="skyblue",
                            face_color=[0, 0.5, 1, 0.2],
                            edge_width=2,
                            name="Channel 2 Spots",
                            opacity=0.8
                        )
                        cenpc_counter.update_count(len(processor.cenpc_centroids))
                    detect_cenpc_done = True

            else:
                segment_done = False
                detect_dna_fish_done = False
                detect_cenpc_done = False
            
            show_info(f"Loaded images from: {os.path.basename(folder_path)}")
        except Exception as e:
            show_info(f"Error loading images: {e}")


def save_segmentation_for_folder(folder_path):
    """Helper function to save segmentation for a single folder"""
    try:
        # Check for intermediate results folder
        intermediate_path = os.path.join(folder_path, "intermediate_results")
        seg_file = os.path.join(intermediate_path, "segmentation.npy")
        
        if os.path.exists(seg_file):
            # Load the segmentation
            segmentation = np.load(seg_file)
            
            # Convert to uint16 for saving as PNG
            segmentation_uint16 = ((segmentation > 0) * segmentation).astype(np.uint16)
            
            # Create "Exported Segmentation" folder in root directory
            root_dir = os.path.dirname(folder_path)
            export_dir = os.path.join(root_dir, "Exported Segmentation")
            os.makedirs(export_dir, exist_ok=True)
            
            # Save as PNG in export directory
            folder_name = os.path.basename(folder_path)
            png_path = os.path.join(export_dir, f"{folder_name}_segmentation.png")
            Image.fromarray(segmentation_uint16).save(png_path)
            
    except Exception as e:
        show_info(f"Error processing folder {os.path.basename(folder_path)}: {str(e)}")

@magicgui(call_button="Save Segmentation as PNG")
def save_segmentation():
    try:
        # Get the checkbox state
        save_all = save_all_checkbox.isChecked()
        
        if save_all:
            # Process all folders in the list
            for i in range(folder_list_widget.count()):
                folder_path = folder_list_widget.item(i).data(Qt.UserRole)
                save_segmentation_for_folder(folder_path)
            show_info("Saved segmentation PNGs for all folders in 'Exported Segmentation' directory")
        else:
            # Process only current folder
            if current_folder_path:
                save_segmentation_for_folder(current_folder_path)
                show_info("Saved segmentation PNG in 'Exported Segmentation' directory")
            else:
                show_info("No folder selected")
    except Exception as e:
        show_info(f"Error saving segmentation: {str(e)}")

# Create a horizontal layout for the save segmentation button and checkbox
class SaveSegmentationWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # Save button
        self.save_button = save_segmentation.native
        self.save_button.setStyleSheet(BUTTON_STYLE)
        self.layout.addWidget(self.save_button)
        
        # Checkbox for "All"
        global save_all_checkbox  # Make it accessible to the save_segmentation function
        save_all_checkbox = QCheckBox("All")
        save_all_checkbox.setStyleSheet("""
            QCheckBox {
                color: white;
                spacing: 5px;
            }
            QCheckBox::indicator {
                width: 13px;
                height: 13px;
            }
            QCheckBox::indicator:unchecked {
                background-color: #2d2d2d;
                border: 1px solid #4a4a4a;
            }
            QCheckBox::indicator:checked {
                background-color: #40c057;
                border: 1px solid #2d2d2d;
            }
        """)
        self.layout.addWidget(save_all_checkbox)
        
        self.setLayout(self.layout)



def load_images0_BU():
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


from napari.utils.colormaps import DirectLabelColormap

@magicgui(call_button="Detect Channel 1 spots")
def detect_dna_fish_spots():
    global detect_dna_fish_done, images, current_folder_path
    if detect_dna_fish_done:
        show_info("Spot detection for channel 1 has already been done.")
        return

    threshold = control_widget_dna_fish.slider.value() / 100
    if images[1] is not None:
        try:
            # Get centroids based on segmentation mode
            centroids = processor.detect_spots_cent(images[1], 'Channel 1', threshold, save_dir=current_folder_path)
            
            # Visualize spots if centroids were found
            if centroids is not None and len(centroids) > 0:
                # Update spot count
                spot_count = len(centroids)
                dna_fish_counter.update_count(spot_count)

                # Create squares around centroids
                squares = [
                    [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                    for x, y in centroids
                ]
                
                # Remove existing DNA-FISH layers
                for layer in list(viewer.layers):
                    if layer.name in ["Channel 1 Spots", "Centroids in Channel 1"]:  
                        viewer.layers.remove(layer)

                # Add new visualization layers
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="yellow",
                    face_color=[1, 1, 0, 0.2],
                    edge_width=2,
                    name="Channel 1 Spots",
                    opacity=0.8

                )
                detect_dna_fish_done = True
                show_info(f"Detected {spot_count} spots in channel 1 image with threshold {threshold}")
            else:
                show_info("No spots detected in channel 1 image")

                detect_dna_fish_done = False
                dna_fish_counter.update_count(0)
                

        except Exception as e:
            show_info(f"Error detecting spots: {str(e)}")
            detect_dna_fish_done = False
            dna_fish_counter.update_count(0)
    else:
        show_info("Channel 1 image not loaded")



@magicgui(call_button="Delete channel 1 Spots")
def delete_dna_fish_spots():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.dna_fish_centroids is None:
            show_info("No channel 1 spots detected yet")
            return
            

        if len(shapes_layer.data) == 0:
            show_info("Please draw a line or select points to remove spots.")
            return
            
        # Delete spots
        processor.delete_dna_fish_spots_with_line(viewer)
        
        # Update spot counter
        if processor.dna_fish_centroids is not None:
            dna_fish_counter.update_count(len(processor.dna_fish_centroids))
            
        # Clear the shapes layer after deletion
        shapes_layer.data = []
        
    except Exception as e:
        show_info(f"Error during spot deletion: {str(e)}")


@magicgui(call_button="Save channel 1 Spots")
def save_dna_fish_spots():
    if processor.dna_fish_centroids is None:
        show_info("No Channel 1 spots to save")
        return
        

    if current_folder_path:
        intermediate_path = os.path.join(current_folder_path, "intermediate_results")
        os.makedirs(intermediate_path, exist_ok=True)
        np.save(os.path.join(intermediate_path, "dna_fish_centroids.npy"), processor.dna_fish_centroids) # check this
        show_info("Channel 1 spots saved successfully")
    else:
        show_info("No folder selected")





@magicgui(call_button="Delete Channel 2 Spots")
def delete_cenpc_spots():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.cenpc_centroids is None:
            show_info("No Channel 2 spots detected yet")
            return
            
        if len(shapes_layer.data) == 0:
            show_info("Please draw a line or select points to remove spots.")
            return
            
        # Delete spots
        processor.delete_cenpc_spots_with_line(viewer)
        
        # Update spot counter
        if processor.cenpc_centroids is not None:
            cenpc_counter.update_count(len(processor.cenpc_centroids))
            
        # Clear the shapes layer after deletion
        shapes_layer.data = []
        
    except Exception as e:
        show_info(f"Error during spot deletion: {str(e)}")

@magicgui(call_button="Save Channel 2 Spots")
def save_cenpc_spots():
    if processor.cenpc_centroids is None:
        show_info("No Channel 2 spots to save")
        return
        
    if current_folder_path:
        intermediate_path = os.path.join(current_folder_path, "intermediate_results")
        os.makedirs(intermediate_path, exist_ok=True)
        np.save(os.path.join(intermediate_path, "cenpc_centroids.npy"), processor.cenpc_centroids)
        show_info("Channel 2 spots saved successfully")
    else:
        show_info("No folder selected")




def delete_dna_fish_spots_with_line(self, line_coords):
    """Delete DNA-FISH spots that intersect with the drawn line."""
    if self.dna_fish_centroids is None or len(self.dna_fish_centroids) == 0:
        return

    # Convert line coordinates to pixel coordinates
    start_point = line_coords[0]
    end_point = line_coords[1]
    
    # Create a mask of the line
    img_shape = self.img_dna_fish.shape if self.img_dna_fish is not None else (1024, 1024)
    line_mask = np.zeros(img_shape, dtype=bool)
    rr, cc = line(int(start_point[0]), int(start_point[1]), 
                  int(end_point[0]), int(end_point[1]))
    valid_points = (rr >= 0) & (rr < img_shape[0]) & (cc >= 0) & (cc < img_shape[1])
    rr, cc = rr[valid_points], cc[valid_points]
    line_mask[rr, cc] = True
    
    # Buffer the line to make it easier to select spots
    from scipy.ndimage import binary_dilation
    line_mask = binary_dilation(line_mask, iterations=3)
    
    # Find spots that don't intersect with the line
    kept_spots = []
    square_size = 5  # Half size of the square around each spot
    
    for spot in self.dna_fish_centroids:
        spot_y, spot_x = int(spot[0]), int(spot[1])
        
        # Check if any part of the square around the spot intersects with the line
        y_min = max(0, spot_y - square_size)
        y_max = min(img_shape[0], spot_y + square_size + 1)
        x_min = max(0, spot_x - square_size)
        x_max = min(img_shape[1], spot_x + square_size + 1)
        
        square_region = line_mask[y_min:y_max, x_min:x_max]
        if not np.any(square_region):  # If no intersection with the line
            kept_spots.append(spot)
    
    # Update centroids
    self.dna_fish_centroids = np.array(kept_spots) if kept_spots else np.array([])
    
    # Update the viewer
    if len(viewer.layers) > 0:
        # Remove existing DNA-FISH spots layer
        for layer in viewer.layers:
            if 'Channel 1 Spots' in layer.name:
                viewer.layers.remove(layer)
        

        # Add updated squares for remaining spots
        if len(self.dna_fish_centroids) > 0:
            squares = [
                [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                for x, y in self.dna_fish_centroids
            ]
            viewer.add_shapes(
                squares,
                shape_type='polygon',
                edge_color="yellow",
                face_color=[1, 1, 0, 0.2],
                edge_width=2,
                name="Channel 1 Spots",
                opacity=0.8

            )


@magicgui(call_button="Detect channel 2 Spots")
def detect_cenpc_spots():
    global detect_cenpc_done, images
    if detect_cenpc_done:
        show_info("Spot detection for channel 2 has already been done.")
        return


    threshold = control_widget_cenpc.slider.value() / 100
    if images[2] is not None:
        try:
            centroids = None
            # Get centroids based on segmentation mode
            if segment_dapi_widget.is_checked():
                processor.detect_spots_cent(images[2], 'Channel 2', threshold, save_dir=current_folder_path)
                centroids = processor.cenpc_centroids
            else:

                processor.detect_spots_cent(images[2], 'Channel 2', threshold, save_dir=current_folder_path)
                centroids = processor.cenpc_centroids


            # Visualize spots if centroids were found
            if centroids is not None and len(centroids) > 0:
                # Update spot count
                spot_count = len(centroids)
                cenpc_counter.update_count(spot_count)

                # Create squares around centroids
                squares = [
                    [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                    for x, y in centroids
                ]
                
                # Remove existing CENPC layers
                for layer in list(viewer.layers):
                    if layer.name in ["Channel 2 Spots", "Centroids in Channel 2"]:
                        viewer.layers.remove(layer)

                # Add new visualization layer
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="skyblue",
                    face_color=[0, 0.5, 1, 0.2],
                    edge_width=2,
                    name="Channel 2 Spots",
                    opacity=0.8

                )

            show_info(f"Detected and labeled spots in channel 2 image with threshold {threshold}")
            detect_cenpc_done = True
        except Exception as e:
            show_info(f"Error detecting spots: {e}")



@magicgui(call_button="Find Common")
def find_common():
    try:
        if segment_dapi_widget.is_checked():
            show_info("Skipping find common due to checkbox selection.")
            return

        common_nuclei = processor.find_common()
        if common_nuclei is None:
            show_info("No common labels found.")
            return
            
        # Save common nuclei to intermediate_results directory
        if current_folder_path:
            intermediate_path = os.path.join(current_folder_path, "intermediate_results")
            os.makedirs(intermediate_path, exist_ok=True)
            common_file = os.path.join(intermediate_path, "common_nuclei.npy")
            np.save(common_file, common_nuclei)
            
        viewer.add_labels(common_nuclei, name="Matched Chromosome")
        show_info("Found common labels and updated the view.")
        return common_nuclei  # Add this return statement
    except Exception as e:
        show_info(f"Error finding common labels: {e}")
        return None



@magicgui(call_button="Get Intensity at Spots Locations")
def get_intensity_at_cenpc_location():
    try:
        # Create intermediate_results directory if it doesn't exist
        if not current_folder_path:
            show_info("No folder selected")
            return
            
        intermediate_path = os.path.join(current_folder_path, "intermediate_results")
        os.makedirs(intermediate_path, exist_ok=True)

        # Basic validation checks with spot counts
        if processor.dna_fish_centroids is None or len(processor.dna_fish_centroids) == 0:
            show_info("Please detect Channel 1 spots first")
            return
        else:
            print(f"Number of Channel 1 spots: {len(processor.dna_fish_centroids)}")
            
        if processor.cenpc_centroids is None or len(processor.cenpc_centroids) == 0:
            show_info("Please detect Channel 2 spots first")
            return
        else:
            print(f"Number of Channel 2 spots: {len(processor.cenpc_centroids)}")
            
        if images[1] is None or images[2] is None:
            show_info("Both channel images are required")
            return

        # Set the images in processor
        processor.img_dna_fish = images[1]
        processor.img_cenpc = images[2]

        # Initialize DataFrames
        df_common_dna_fish = None
        df_common_cenpc = None

        # Process based on segmentation mode
        if segment_dapi_widget.is_checked():
            show_info("Calculating intensities at both channel locations without segmentation.")
            df_common_dna_fish = pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X'])
            df_common_cenpc = pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X'])
            
            print(f"Channel 1 spots DataFrame shape: {df_common_dna_fish.shape}")
            print(f"Channel 2 spots DataFrame shape: {df_common_cenpc.shape}")
            
        else:
            show_info("Calculating intensities at both channel locations in segmented regions.")
            if processor.nuclei is None:
                show_info("Please segment the image first")
                return
                
            common_nuclei = processor.find_common()
            if common_nuclei is None:
                show_info("No common regions found for intensity calculation")
                return

            # Get spots only in common regions
            df_common_dna_fish = processor.get_spots_in_common_regions(
                pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X']), 
                common_nuclei
            )
            df_common_cenpc = processor.get_spots_in_common_regions(
                pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X']), 
                common_nuclei
            )
            
            print(f"Channel 1 spots in common regions: {len(df_common_dna_fish) if df_common_dna_fish is not None else 0}")
            print(f"Channel 2 spots in common regions: {len(df_common_cenpc) if df_common_cenpc is not None else 0}")
            
            if df_common_dna_fish is None or df_common_cenpc is None:
                show_info("No spots found in common regions")
                return

        # Calculate intensities for both channels independently
        df_ch2_at_ch1 = processor.measure_intensity_at_spots(
            intensity_image=processor.img_cenpc,  # Measure Channel 2 intensity
            spots_df=df_common_dna_fish,         # at Channel 1 spots
            channel_name='Channel1'
        )

        df_ch1_at_ch2 = processor.measure_intensity_at_spots(
            intensity_image=processor.img_dna_fish,  # Measure Channel 1 intensity
            spots_df=df_common_cenpc,               # at Channel 2 spots
            channel_name='Channel2'
        )

        folder_name = os.path.basename(current_folder_path)
        saved_files = []

        # Save and process Channel 2 intensity at Channel 1 spots
        if df_ch2_at_ch1 is not None and not df_ch2_at_ch1.empty:
            df_ch2_at_ch1['Skip_Segmentation'] = segment_dapi_widget.is_checked()
            df_ch2_at_ch1['Folder'] = folder_name
            
            # Remove any rows with NaN values
            df_ch2_at_ch1 = df_ch2_at_ch1.dropna()
            
            if len(df_ch2_at_ch1) > 0:
                save_path = os.path.join(intermediate_path, f"{folder_name}_ch2_intensity_at_ch1_spots.csv")
                df_ch2_at_ch1.to_csv(save_path, index=False)
                saved_files.append(save_path)
                
                print("\nChannel 2 intensity at Channel 1 spots:")
                print(f"Mean intensity: {df_ch2_at_ch1['Intensity'].mean():.2f}")
                print(f"Median intensity: {df_ch2_at_ch1['Intensity'].median():.2f}")
                print(f"Standard deviation: {df_ch2_at_ch1['Intensity'].std():.2f}")

        # Save and process Channel 1 intensity at Channel 2 spots
        if df_ch1_at_ch2 is not None and not df_ch1_at_ch2.empty:
            df_ch1_at_ch2['Skip_Segmentation'] = segment_dapi_widget.is_checked()
            df_ch1_at_ch2['Folder'] = folder_name
            
            # Remove any rows with NaN values
            df_ch1_at_ch2 = df_ch1_at_ch2.dropna()
            
            if len(df_ch1_at_ch2) > 0:
                save_path = os.path.join(intermediate_path, f"{folder_name}_ch1_intensity_at_ch2_spots.csv")
                df_ch1_at_ch2.to_csv(save_path, index=False)
                saved_files.append(save_path)
                
                print("\nChannel 1 intensity at Channel 2 spots:")
                print(f"Mean intensity: {df_ch1_at_ch2['Intensity'].mean():.2f}")
                print(f"Median intensity: {df_ch1_at_ch2['Intensity'].median():.2f}")
                print(f"Standard deviation: {df_ch1_at_ch2['Intensity'].std():.2f}")

        if saved_files:
            show_info("Intensity measurements saved to:\n" + "\n".join(saved_files))
        else:
            show_info("No intensity measurements were saved")
            
    except Exception as e:
        show_info(f"Error measuring intensities: {str(e)}")
        print(f"Detailed error: {str(e)}")  # More detailed error in console


def get_spots_in_common_regions(self, spots_df, common_nuclei):
    """
    Filter spots DataFrame to only include spots in common regions.
    
    Args:
        spots_df: DataFrame with spot coordinates
        common_nuclei: Label image of common regions
        
    Returns:
        DataFrame with only spots in common regions
    """
    if spots_df is None or common_nuclei is None:
        return None
        
    # Filter spots to only those in common regions
    common_spots = []
    for _, row in spots_df.iterrows():
        y, x = int(row['Y']), int(row['X'])
        if 0 <= y < common_nuclei.shape[0] and 0 <= x < common_nuclei.shape[1]:
            if common_nuclei[y, x] > 0:  # spot is in a common region
                common_spots.append(row)
    
    return pd.DataFrame(common_spots)


@magicgui(call_button="Run All")
def run_all():
    global images, current_folder_path
    try:
        if not current_folder_path:
            show_info("No folder selected")
            return
            
        threshold_dna_fish = control_widget_dna_fish.slider.value() / 100
        threshold_cenpc = control_widget_cenpc.slider.value() / 100

        if all(img is not None for img in images if img is not None):
            # Create intermediate_results directory
            intermediate_path = os.path.join(current_folder_path, "intermediate_results")
            os.makedirs(intermediate_path, exist_ok=True)

            # Remove existing spot layers before processing
            for layer in list(viewer.layers):
                if layer.name in ["Channel 1 Spots", "Channel 2 Spots", "Cellpose Segmented", "Matched Chromosome"]:
                    viewer.layers.remove(layer)

            # Case 1: No Segmentation (checkbox checked)
            if segment_dapi_widget.is_checked():
                # Detect spots without segmentation
                processor.detect_spots_cent(images[1], 'Channel 1', threshold_dna_fish, save_dir=current_folder_path)
                processor.detect_spots_cent(images[2], 'Channel 2', threshold_cenpc, save_dir=current_folder_path)
                
                df_common_dna_fish = pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X'])
                df_common_cenpc = pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X'])
                
            # Case 2: With Segmentation (checkbox unchecked)
            else:
                try:
                    # Segment DAPI
                    masks = processor.segment_image(images[0], save_dir=current_folder_path)
                    if masks is None:
                        show_info("Segmentation failed")
                        return
                        
                    viewer.add_labels(masks, name="Cellpose Segmented")
                    chromosome_counter.update_count(len(np.unique(masks)) - 1)
                except Exception as e:
                    show_info(f"Error during segmentation: {str(e)}")
                    return
                
                # Detect spots
                processor.detect_spots_cent(images[1], 'Channel 1', threshold_dna_fish, save_dir=current_folder_path)
                processor.detect_spots_cent(images[2], 'Channel 2', threshold_cenpc, save_dir=current_folder_path)

                # Find common regions and get spots in those regions
                common_nuclei = processor.find_common()
                if common_nuclei is None:
                    show_info("No common regions found")
                    return

                viewer.add_labels(common_nuclei, name="Matched Chromosome")

                df_common_dna_fish = processor.get_spots_in_common_regions(
                    pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X']), 
                    common_nuclei
                )
                df_common_cenpc = processor.get_spots_in_common_regions(
                    pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X']), 
                    common_nuclei
                )

            # Update spot visualization and counters
            if processor.dna_fish_centroids is not None:
                dna_fish_counter.update_count(len(processor.dna_fish_centroids))
                squares = [
                    [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                    for x, y in processor.dna_fish_centroids
                ]
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="yellow",
                    face_color=[1, 1, 0, 0.2],
                    edge_width=2,
                    name="Channel 1 Spots",
                    opacity=0.8
                )

            if processor.cenpc_centroids is not None:
                cenpc_counter.update_count(len(processor.cenpc_centroids))
                squares = [
                    [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                    for x, y in processor.cenpc_centroids
                ]
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="skyblue",
                    face_color=[0, 0.5, 1, 0.2],
                    edge_width=2,
                    name="Channel 2 Spots",
                    opacity=0.8
                )

            # Calculate intensities for both channels independently
            df_ch2_at_ch1 = processor.measure_intensity_at_spots(
                intensity_image=processor.img_cenpc,  # Measure Channel 2 intensity
                spots_df=df_common_dna_fish,         # at Channel 1 spots
                channel_name='Channel1'
            )

            df_ch1_at_ch2 = processor.measure_intensity_at_spots(
                intensity_image=processor.img_dna_fish,  # Measure Channel 1 intensity
                spots_df=df_common_cenpc,               # at Channel 2 spots
                channel_name='Channel2'
            )

            folder_name = os.path.basename(current_folder_path)
            saved_files = []

            # Save and process Channel 2 intensity at Channel 1 spots
            if df_ch2_at_ch1 is not None and not df_ch2_at_ch1.empty:
                df_ch2_at_ch1['Skip_Segmentation'] = segment_dapi_widget.is_checked()
                df_ch2_at_ch1['Folder'] = folder_name
                df_ch2_at_ch1['Channel1_Threshold'] = threshold_dna_fish
                df_ch2_at_ch1['Channel2_Threshold'] = threshold_cenpc
                
                # Remove any rows with NaN values
                df_ch2_at_ch1 = df_ch2_at_ch1.dropna()
                
                if len(df_ch2_at_ch1) > 0:
                    save_path = os.path.join(intermediate_path, f"{folder_name}_ch2_intensity_at_ch1_spots.csv")
                    df_ch2_at_ch1.to_csv(save_path, index=False)
                    saved_files.append(save_path)
                    
                    print("\nChannel 2 intensity at Channel 1 spots:")
                    print(f"Mean intensity: {df_ch2_at_ch1['Intensity'].mean():.2f}")
                    print(f"Median intensity: {df_ch2_at_ch1['Intensity'].median():.2f}")
                    print(f"Standard deviation: {df_ch2_at_ch1['Intensity'].std():.2f}")

            # Save and process Channel 1 intensity at Channel 2 spots
            if df_ch1_at_ch2 is not None and not df_ch1_at_ch2.empty:
                df_ch1_at_ch2['Skip_Segmentation'] = segment_dapi_widget.is_checked()
                df_ch1_at_ch2['Folder'] = folder_name
                df_ch1_at_ch2['Channel1_Threshold'] = threshold_dna_fish
                df_ch1_at_ch2['Channel2_Threshold'] = threshold_cenpc
                
                # Remove any rows with NaN values
                df_ch1_at_ch2 = df_ch1_at_ch2.dropna()
                
                if len(df_ch1_at_ch2) > 0:
                    save_path = os.path.join(intermediate_path, f"{folder_name}_ch1_intensity_at_ch2_spots.csv")
                    df_ch1_at_ch2.to_csv(save_path, index=False)
                    saved_files.append(save_path)
                    
                    print("\nChannel 1 intensity at Channel 2 spots:")
                    print(f"Mean intensity: {df_ch1_at_ch2['Intensity'].mean():.2f}")
                    print(f"Median intensity: {df_ch1_at_ch2['Intensity'].median():.2f}")
                    print(f"Standard deviation: {df_ch1_at_ch2['Intensity'].std():.2f}")

            if saved_files:
                show_info("Intensity measurements saved to:\n" + "\n".join(saved_files))
            else:
                show_info("No intensity measurements were saved")

            show_info("Run all processing completed")
        else:
            show_info("Ensure that all images are loaded")
    except Exception as e:
        show_info(f"Error during run all processing: {str(e)}")

@magicgui(call_button="Batch Load")
def batch_load():
    root_folder = QFileDialog.getExistingDirectory(caption='Select Root Folder for Batch Loading')
    if root_folder:
        try:
            # Clear previous items
            folder_list_widget.clear()
            
            # Store full paths but display only folder names
            for folder_name in os.listdir(root_folder):
                # Skip the "Exported Segmentation" folder
                if folder_name == "Exported Segmentation":
                    continue
                    
                folder_path = os.path.join(root_folder, folder_name)
                if os.path.isdir(folder_path):
                    item = QListWidgetItem(os.path.basename(folder_path))
                    # Store full path as item data
                    item.setData(Qt.UserRole, folder_path)
                    folder_list_widget.addItem(item)
            
            show_info(f"Found {folder_list_widget.count()} folders")
        except Exception as e:
            show_info(f"Error loading folders: {str(e)}")

def on_folder_selected():
    try:
        current_item = folder_list_widget.currentItem()
        if current_item:
            # Get the full path from item data
            selected_folder = current_item.data(Qt.UserRole)
            if selected_folder:
                postprocessing.set_current_folder(selected_folder)

                dapi_id = channel_identifiers.dapi_text.text()
                dna_fish_id = channel_identifiers.dna_fish_text.text()
                cenpc_id = channel_identifiers.cenpc_text.text()
                
                global segment_done, detect_dna_fish_done, detect_cenpc_done, current_folder_path, images
                current_folder_path = selected_folder
                
                images = processor.load_images(current_folder_path, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
                
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
                    viewer.add_image(images[1], name='Channel 1')
                    viewer.add_image(images[2], name='Channel 2')
                else:
                    viewer.add_image(images[0], name='DAPI')
                    viewer.add_image(images[1], name='Channel 1')
                    viewer.add_image(images[2], name='Channel 2')

                
                # Check for and load intermediate results
                intermediate_path = os.path.join(current_folder_path, "intermediate_results")
                if os.path.exists(intermediate_path):
                    # Load segmentation if it exists
                    seg_file = os.path.join(intermediate_path, "segmentation.npy")
                    if os.path.exists(seg_file):
                        processor.nuclei = np.load(seg_file)
                        viewer.add_labels(processor.nuclei, name="Cellpose Segmented")
                        segment_done = True
                        
                        # Update chromosome count
                        unique_chromosomes = len(np.unique(processor.nuclei)) - 1
                        chromosome_counter.update_count(unique_chromosomes)
                        show_info("Loaded existing segmentation")
                        
                        # Add loading of common nuclei here
                        common_file = os.path.join(intermediate_path, "common_nuclei.npy")
                        if os.path.exists(common_file):
                            try:
                                common_nuclei = np.load(common_file)
                                viewer.add_labels(common_nuclei, name="Matched Chromosome")
                                show_info("Loaded matched chromosomes")
                            except Exception as e:
                                show_info(f"Error loading matched chromosomes: {str(e)}")
                    else:
                        segment_done = False
                        show_info("No existing segmentation found")
                    
                    # Load DNA-FISH spots if they exist
                    dna_fish_file = os.path.join(intermediate_path, "dna_fish_spots.npy")
                    dna_fish_centroids_file = os.path.join(intermediate_path, "dna_fish_centroids.npy")
                    if os.path.exists(dna_fish_file) and os.path.exists(dna_fish_centroids_file):
                        processor.labels_dna_fish = np.load(dna_fish_file)
                        processor.dna_fish_centroids = np.load(dna_fish_centroids_file)
                        
                        # Visualize DNA-FISH spots
                        if processor.dna_fish_centroids is not None and len(processor.dna_fish_centroids) > 0:
                            squares = [
                                [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                                for x, y in processor.dna_fish_centroids
                            ]
                            viewer.add_shapes(
                                squares,
                                shape_type='polygon',
                                edge_color="yellow",
                                face_color=[1, 1, 0, 0.2],
                                edge_width=2,
                                name="Channel 1 Spots",
                                opacity=0.8
                            )
                            dna_fish_counter.update_count(len(processor.dna_fish_centroids))

                            detect_dna_fish_done = True
                            show_info("Loaded existing Channel 1 spots")
                    
                    # Load CENPC spots if they exist
                    cenpc_file = os.path.join(intermediate_path, "cenpc_spots.npy")
                    cenpc_centroids_file = os.path.join(intermediate_path, "cenpc_centroids.npy")
                    if os.path.exists(cenpc_file) and os.path.exists(cenpc_centroids_file):
                        processor.labels_cenpc = np.load(cenpc_file)
                        processor.cenpc_centroids = np.load(cenpc_centroids_file)
                        
                        # Visualize CENPC spots
                        if processor.cenpc_centroids is not None and len(processor.cenpc_centroids) > 0:
                            squares = [
                                [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                                for x, y in processor.cenpc_centroids
                            ]
                            viewer.add_shapes(
                                squares,
                                shape_type='polygon',
                                edge_color="skyblue",
                                face_color=[0, 0.5, 1, 0.2],
                                edge_width=2,
                                name="Channel 2 Spots", 
                                opacity=0.8
                            )
                            cenpc_counter.update_count(len(processor.cenpc_centroids))

                            detect_cenpc_done = True
                            show_info("Loaded existing Channel 2 spots")
                else:
                    segment_done = False
                    detect_dna_fish_done = False
                    detect_cenpc_done = False
                
                show_info(f"Loaded images from: {os.path.basename(current_folder_path)}")
                
    except Exception as e:
        show_info(f"Error loading selected folder: {str(e)}")



# Connect the selection signal to our handler
folder_list_widget.itemClicked.connect(on_folder_selected)

# Create a container widget for the list
folder_list_container = QWidget()
layout = QVBoxLayout()
layout.addWidget(folder_list_widget)
folder_list_container.setLayout(layout)


def process_folder_data(folder_path, folder_name, intermediate_path, dapi_id, dna_fish_id, cenpc_id, processor, segment_dapi_widget):
    """Process a single folder's data and return intensity results and spot counts"""
    try:
        # Load saved data from intermediate_results directory
        seg_file = os.path.join(intermediate_path, "segmentation.npy")
        dna_fish_file = os.path.join(intermediate_path, "dna_fish_spots.npy")
        dna_fish_centroids_file = os.path.join(intermediate_path, "dna_fish_centroids.npy")
        cenpc_file = os.path.join(intermediate_path, "cenpc_spots.npy")
        cenpc_centroids_file = os.path.join(intermediate_path, "cenpc_centroids.npy")
        common_file = os.path.join(intermediate_path, "common_nuclei.npy")

        # Check if required files exist
        if not all(os.path.exists(f) for f in [dna_fish_file, dna_fish_centroids_file, cenpc_file, cenpc_centroids_file]):
            raise FileNotFoundError("Missing required spot detection files")

        # Load the saved data
        processor.labels_dna_fish = np.load(dna_fish_file)
        processor.dna_fish_centroids = np.load(dna_fish_centroids_file)
        processor.labels_cenpc = np.load(cenpc_file)
        processor.cenpc_centroids = np.load(cenpc_centroids_file)

        # Load the original images for intensity calculation
        skip_segmentation = segment_dapi_widget.is_checked()
        if skip_segmentation:
            images = processor.load_images(folder_path, None, dna_fish_id, cenpc_id, True)
            if images is None:
                raise ValueError("Couldn't load images")
            processor.img_cenpc = images[1]  # CENPC is second image in skip mode
            processor.img_dna_fish = images[0]  # DNA-FISH is first image in skip mode
            df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
        else:
            images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, False)
            if images is None:
                raise ValueError("Couldn't load images")
            processor.img_cenpc = images[2]  # CENPC is third image in regular mode
            processor.img_dna_fish = images[1]  # DNA-FISH is second image in regular mode
            
            # Load segmentation and common nuclei
            if not os.path.exists(seg_file):
                raise FileNotFoundError("Missing segmentation file")
            if not os.path.exists(common_file):
                raise FileNotFoundError("Missing common nuclei file")
                
            processor.nuclei = np.load(seg_file)
            processor.common_nuclei = np.load(common_file)
            
            # Calculate intensities using saved data
            df_with_cenpc_inten = processor.gen_intensity_from_df(
                processor.img_cenpc,
                processor.df_centroid_dna_fish
            )

        if df_with_cenpc_inten is None or df_with_cenpc_inten.empty:
            raise ValueError(f"No intensity measurements found for folder: {folder_name}")

        # Save intensity results
        intensity_save_path = os.path.join(folder_path, f"{folder_name}_intensity.csv")
        df_with_cenpc_inten.to_csv(intensity_save_path, index=False)

        return {
            'df': df_with_cenpc_inten,
            'channel_1_count': len(processor.dna_fish_centroids) if processor.dna_fish_centroids is not None else 0,
            'channel_2_count': len(processor.cenpc_centroids) if processor.cenpc_centroids is not None else 0,
            'mean_intensity': df_with_cenpc_inten['Channel2_Intensity'].mean(),

            'common_regions': True if not skip_segmentation else False

        }

    except Exception as e:
        raise Exception(f"Error processing folder {folder_name}: {str(e)}")

def process_folder_data_BU(folder_path, folder_name, intermediate_path, dapi_id, dna_fish_id, cenpc_id, processor, segment_dapi_widget):
    """Process a single folder's data and return intensity results and spot counts"""
    
    # Load images based on segmentation mode
    skip_segmentation = segment_dapi_widget.is_checked()
    if skip_segmentation:
        images = processor.load_images(folder_path, None, dna_fish_id, cenpc_id, True)
        if images is None:
            raise ValueError("Couldn't load images")
        processor.img_cenpc = images[1]  # CENPC is second image in skip mode
        processor.img_dna_fish = images[0]  # DNA-FISH is first image in skip mode
    else:
        images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, False)
        if images is None:
            raise ValueError("Couldn't load images")
        processor.img_cenpc = images[2]  # CENPC is third image in regular mode
        processor.img_dna_fish = images[1]  # DNA-FISH is second image in regular mode

    # Load spot detection data
    dna_fish_file = os.path.join(intermediate_path, "dna_fish_spots.npy")
    cenpc_file = os.path.join(intermediate_path, "cenpc_spots.npy")
    
    if not (os.path.exists(dna_fish_file) and os.path.exists(cenpc_file)):
        raise FileNotFoundError("Missing spot detection files")
        
    dna_fish_spots = np.load(dna_fish_file)
    cenpc_spots = np.load(cenpc_file)
    
    # Calculate intensities based on segmentation mode
    if skip_segmentation:
        processor.dna_fish_centroids = np.argwhere(dna_fish_spots > 0)
        processor.cenpc_centroids = np.argwhere(cenpc_spots > 0)
        df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
    else:
        df_with_cenpc_inten = processor.gen_intensity_from_df(
            processor.img_cenpc,
            processor.df_centroid_dna_fish
        )
    
    if df_with_cenpc_inten is None or df_with_cenpc_inten.empty:
        raise ValueError("No intensity measurements found")
        
    # Save intensity results
    intensity_save_path = os.path.join(intermediate_path, f"{folder_name}_intensity.csv")
    df_with_cenpc_inten.to_csv(intensity_save_path, index=False)
    
    return {
        'df': df_with_cenpc_inten,
        'dna_fish_count': len(df_with_cenpc_inten),
        'cenpc_count': len(np.argwhere(cenpc_spots > 0)),
        'mean_intensity': df_with_cenpc_inten['CENPC_Intensity'].mean()
    }


@magicgui(
    call_button="Batch Processing",
    use_current_settings={'widget_type': 'CheckBox', 'text': 'Use Current UI Settings'}
)
def batch_processing(use_current_settings: bool):
    try:
        summary_data = []
        all_ch1_intensities = []
        all_ch2_intensities = []
        
        folders_to_process = [
            folder_list_widget.item(i).data(Qt.UserRole)
            for i in range(folder_list_widget.count())
        ]
        if not folders_to_process:
            show_info("No folders in the list to process")
            return
            
        if use_current_settings:
            # Using current UI settings
            threshold_dna_fish = control_widget_dna_fish.slider.value() / 100
            threshold_cenpc = control_widget_cenpc.slider.value() / 100
            skip_segmentation = segment_dapi_widget.is_checked()
            dapi_id = channel_identifiers.dapi_text.text()
            dna_fish_id = channel_identifiers.dna_fish_text.text()
            cenpc_id = channel_identifiers.cenpc_text.text()
            show_info(f"Using current UI settings:\n"
                     f"Channel 1 Threshold: {threshold_dna_fish}\n"
                     f"Channel 2 Threshold: {threshold_cenpc}\n"
                     f"Skip Segmentation: {skip_segmentation}\n"
                     f"DAPI ID: {dapi_id}\n"
                     f"Channel 1 ID: {dna_fish_id}\n"
                     f"Channel 2 ID: {cenpc_id}")

            for folder_path in folders_to_process:
                folder_name = os.path.basename(folder_path)
                try:
                    print(f"\nProcessing folder: {folder_name}")
                    global current_folder_path, images
                    current_folder_path = folder_path
                    
                    # Create intermediate_results directory
                    intermediate_path = os.path.join(folder_path, "intermediate_results")
                    os.makedirs(intermediate_path, exist_ok=True)
                    
                    # Load images based on segmentation mode
                    images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, skip_segmentation)
                    if images is None or any(img is None for img in images if not skip_segmentation or images.index(img) > 0):
                        raise ValueError("Failed to load one or more images")
                    
                    chromosome_count = 0
                    matched_count = 0
                    
                    # Set the images in processor
                    processor.img_dna_fish = images[1]
                    processor.img_cenpc = images[2]

                    # Case 1: No Segmentation (checkbox checked)
                    if skip_segmentation:
                        # Detect spots without segmentation
                        processor.detect_spots_cent(images[1], 'Channel 1', threshold_dna_fish, save_dir=current_folder_path)
                        processor.detect_spots_cent(images[2], 'Channel 2', threshold_cenpc, save_dir=current_folder_path)
                        
                        df_common_dna_fish = pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X'])
                        df_common_cenpc = pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X'])
                        
                    # Case 2: With Segmentation (checkbox unchecked)
                    else:
                        # Segment DAPI
                        masks = processor.segment_image(images[0], save_dir=current_folder_path)
                        if masks is None:
                            raise ValueError("Segmentation failed")
                        chromosome_count = len(np.unique(masks)) - 1
                        
                        # Detect spots
                        processor.detect_spots_cent(images[1], 'Channel 1', threshold_dna_fish, save_dir=current_folder_path)
                        processor.detect_spots_cent(images[2], 'Channel 2', threshold_cenpc, save_dir=current_folder_path)
                        
                        # Find common regions
                        common_nuclei = processor.find_common()
                        if common_nuclei is None:
                            raise ValueError("No common regions found")
                        matched_count = len(np.unique(common_nuclei)) - 1
                        
                        # Get spots in common regions
                        df_common_dna_fish = processor.get_spots_in_common_regions(
                            pd.DataFrame(processor.dna_fish_centroids, columns=['Y', 'X']), 
                            common_nuclei
                        )
                        df_common_cenpc = processor.get_spots_in_common_regions(
                            pd.DataFrame(processor.cenpc_centroids, columns=['Y', 'X']), 
                            common_nuclei
                        )

                    # Calculate intensities for both channels independently
                    df_ch2_at_ch1 = processor.measure_intensity_at_spots(
                        intensity_image=processor.img_cenpc,  # Measure Channel 2 intensity
                        spots_df=df_common_dna_fish,         # at Channel 1 spots
                        channel_name='Channel1'
                    )

                    df_ch1_at_ch2 = processor.measure_intensity_at_spots(
                        intensity_image=processor.img_dna_fish,  # Measure Channel 1 intensity
                        spots_df=df_common_cenpc,               # at Channel 2 spots
                        channel_name='Channel2'
                    )

                    if df_ch2_at_ch1 is not None and df_ch1_at_ch2 is not None:
                        # Add metadata to Channel 2 at Channel 1 spots
                        df_ch2_at_ch1['Skip_Segmentation'] = skip_segmentation
                        df_ch2_at_ch1['Folder'] = folder_name
                        df_ch2_at_ch1['Channel1_Threshold'] = threshold_dna_fish
                        df_ch2_at_ch1['Channel2_Threshold'] = threshold_cenpc
                        df_ch2_at_ch1['Chromosome_Count'] = chromosome_count
                        df_ch2_at_ch1['Matched_Nuclei_Count'] = matched_count

                        # Add metadata to Channel 1 at Channel 2 spots
                        df_ch1_at_ch2['Skip_Segmentation'] = skip_segmentation
                        df_ch1_at_ch2['Folder'] = folder_name
                        df_ch1_at_ch2['Channel1_Threshold'] = threshold_dna_fish
                        df_ch1_at_ch2['Channel2_Threshold'] = threshold_cenpc
                        df_ch1_at_ch2['Chromosome_Count'] = chromosome_count
                        df_ch1_at_ch2['Matched_Nuclei_Count'] = matched_count

                        # Save intensity measurements
                        ch2_at_ch1_path = os.path.join(intermediate_path, f"{folder_name}_ch2_intensity_at_ch1_spots.csv")
                        ch1_at_ch2_path = os.path.join(intermediate_path, f"{folder_name}_ch1_intensity_at_ch2_spots.csv")
                        
                        df_ch2_at_ch1.to_csv(ch2_at_ch1_path, index=False)
                        df_ch1_at_ch2.to_csv(ch1_at_ch2_path, index=False)

                        # Add to summary data
                        summary_row = {
                            'Folder': folder_name,
                            'Channel1_Spots': len(df_ch2_at_ch1),
                            'Channel2_Spots': len(df_ch1_at_ch2),
                            'Mean_Ch2_at_Ch1': df_ch2_at_ch1['Intensity'].mean(),
                            'Mean_Ch1_at_Ch2': df_ch1_at_ch2['Intensity'].mean(),
                            'Chromosome_Count': chromosome_count,
                            'Matched_Nuclei_Count': matched_count,
                            'Skip_Segmentation': skip_segmentation,
                            'Channel1_Threshold': threshold_dna_fish,
                            'Channel2_Threshold': threshold_cenpc
                        }
                        summary_data.append(summary_row)
                        all_ch1_intensities.append(df_ch1_at_ch2)
                        all_ch2_intensities.append(df_ch2_at_ch1)
                        
                        print(f"\nResults for folder: {folder_name}")
                        print(f"Channel 1 spots: {len(df_ch2_at_ch1)}")
                        print(f"Channel 2 spots: {len(df_ch1_at_ch2)}")
                        print(f"Mean Ch2 at Ch1 intensity: {summary_row['Mean_Ch2_at_Ch1']:.2f}")
                        print(f"Mean Ch1 at Ch2 intensity: {summary_row['Mean_Ch1_at_Ch2']:.2f}")
                        if not skip_segmentation:
                            print(f"Chromosome count: {chromosome_count}")
                            print(f"Matched nuclei count: {matched_count}")

                except Exception as e:
                    print(f"Error processing folder {folder_name}: {str(e)}")
                    continue

        else: 
        
            # Process existing saved data
            for folder_path in folders_to_process:
                folder_name = os.path.basename(folder_path)
                try:
                    print(f"\nAnalyzing saved data for folder: {folder_name}")
                    
                    # Check for required files in intermediate_results
                    intermediate_path = os.path.join(folder_path, "intermediate_results")
                    ch2_at_ch1_file = os.path.join(intermediate_path, f"{folder_name}_ch2_intensity_at_ch1_spots.csv")
                    ch1_at_ch2_file = os.path.join(intermediate_path, f"{folder_name}_ch1_intensity_at_ch2_spots.csv")
                    
                    if not (os.path.exists(ch2_at_ch1_file) and os.path.exists(ch1_at_ch2_file)):
                        print(f"Skipping {folder_name}: Missing intensity CSV files")
                        continue
                        
                    # Load intensity data
                    df_ch2_at_ch1 = pd.read_csv(ch2_at_ch1_file)
                    df_ch1_at_ch2 = pd.read_csv(ch1_at_ch2_file)
                    
                    # Load spot counts from intermediate results if they exist
                    dna_fish_centroids_file = os.path.join(intermediate_path, "dna_fish_centroids.npy")
                    cenpc_centroids_file = os.path.join(intermediate_path, "cenpc_centroids.npy")
                    
                    channel_1_count = 0
                    channel_2_count = 0
                    
                    if os.path.exists(dna_fish_centroids_file):
                        channel_1_spots = np.load(dna_fish_centroids_file)
                        channel_1_count = len(channel_1_spots)
                        
                    if os.path.exists(cenpc_centroids_file):
                        channel_2_spots = np.load(cenpc_centroids_file)
                        channel_2_count = len(channel_2_spots)
                    
                    # Check for segmentation data
                    seg_file = os.path.join(intermediate_path, "segmentation.npy")
                    common_file = os.path.join(intermediate_path, "common_nuclei.npy")
                    
                    chromosome_count = 0
                    matched_count = 0
                    
                    if os.path.exists(seg_file):
                        segmentation = np.load(seg_file)
                        chromosome_count = len(np.unique(segmentation)) - 1
                        
                    if os.path.exists(common_file):
                        common_nuclei = np.load(common_file)
                        matched_count = len(np.unique(common_nuclei)) - 1
                    
                    # Create summary row
                    summary_row = {
                        'Folder': folder_name,
                        'Channel1_Spots': channel_1_count,
                        'Channel2_Spots': channel_2_count,
                        'Mean_Ch2_at_Ch1': df_ch2_at_ch1['Intensity'].mean(),
                        'Mean_Ch1_at_Ch2': df_ch1_at_ch2['Intensity'].mean(),
                        'Chromosome_Count': chromosome_count,
                        'Matched_Nuclei_Count': matched_count,
                        'Skip_Segmentation': df_ch2_at_ch1['Skip_Segmentation'].iloc[0] if 'Skip_Segmentation' in df_ch2_at_ch1.columns else None,
                        'Channel1_Threshold': df_ch2_at_ch1['Channel1_Threshold'].iloc[0] if 'Channel1_Threshold' in df_ch2_at_ch1.columns else None,
                        'Channel2_Threshold': df_ch2_at_ch1['Channel2_Threshold'].iloc[0] if 'Channel2_Threshold' in df_ch2_at_ch1.columns else None
                    }
                    summary_data.append(summary_row)
                    
                    # Add to all intensities lists
                    all_ch1_intensities.append(df_ch1_at_ch2)
                    all_ch2_intensities.append(df_ch2_at_ch1)
                    
                    print(f"\nResults for folder: {folder_name}")
                    print(f"Channel 1 spots: {channel_1_count}")
                    print(f"Channel 2 spots: {channel_2_count}")
                    print(f"Mean Ch2 at Ch1 intensity: {summary_row['Mean_Ch2_at_Ch1']:.2f}")
                    print(f"Mean Ch1 at Ch2 intensity: {summary_row['Mean_Ch1_at_Ch2']:.2f}")
                    if chromosome_count > 0:
                        print(f"Chromosome count: {chromosome_count}")
                        print(f"Matched nuclei count: {matched_count}")
                    
                except Exception as e:
                    print(f"Error processing saved data for {folder_name}: {str(e)}")
                    continue



        if summary_data:
            try:
                # Create summary DataFrame
                df_summary = pd.DataFrame(summary_data)
                
                # Save summary to root directory
                root_dir = os.path.dirname(folders_to_process[0])
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                summary_path = os.path.join(root_dir, f"batch_summary_{timestamp}.csv")
                df_summary.to_csv(summary_path, index=False)
                
                # Save combined intensity data for each channel
                if all_ch1_intensities and all_ch2_intensities:
                    df_all_ch1 = pd.concat(all_ch1_intensities, ignore_index=True)
                    df_all_ch2 = pd.concat(all_ch2_intensities, ignore_index=True)
                    
                    ch1_path = os.path.join(root_dir, f"all_ch1_intensities_{timestamp}.csv")
                    ch2_path = os.path.join(root_dir, f"all_ch2_intensities_{timestamp}.csv")
                    
                    df_all_ch1.to_csv(ch1_path, index=False)
                    df_all_ch2.to_csv(ch2_path, index=False)
                
                show_info(f"Batch processing completed. Summary saved to:\n{summary_path}")

            except Exception as e:
                    show_info(f"Error saving summary data: {str(e)}")   
        else:
                show_info("No summary data generated")


            
    except Exception as e:
        show_info(f"Error during batch processing: {str(e)}")



chromosome_counter = ChromosomeCountWidget()
dna_fish_counter = SpotCountWidget("Channel 1")
cenpc_counter = SpotCountWidget("Channel 2")



postprocessing = SegmentationPostprocessing(viewer, processor, chromosome_counter)

segment_dapi_widget = SegmentDAPIWidget(postprocessing=postprocessing)








# Add the widgets to the viewer
control_widget_dna_fish = ControlWidgetDNAFISH()
control_widget_cenpc = ControlWidgetCENPC()
channel_identifiers = ChannelIdentifiers()
batch_processor = BatchProcessor(processor, control_widget_dna_fish, control_widget_cenpc)

# Add the widgets to the viewer
viewer.window.add_dock_widget(channel_identifiers, area='right', name='Channel Identifiers')
viewer.window.add_dock_widget(load_images, area='right', name='')
viewer.window.add_dock_widget(segment_dapi_widget, area='right', name='Segment DAPI Control')
viewer.window.add_dock_widget(control_widget_dna_fish, area='right', name='Detect Channel 1 Spot Control')
viewer.window.add_dock_widget(control_widget_cenpc, area='right', name='Detect Channel 2 Spot Control')
viewer.window.add_dock_widget(find_common, area='right', name='')
viewer.window.add_dock_widget(get_intensity_at_cenpc_location, area='right', name='')



toggle_button = QPushButton("Toggle All Layers")
toggle_button.setStyleSheet(BUTTON_STYLE)

def toggle_all_layers():
    current_state = next(iter(viewer.layers)).visible if viewer.layers else True
    for layer in viewer.layers:
        layer.visible = not current_state

toggle_button.clicked.connect(toggle_all_layers)
viewer.window.add_dock_widget(toggle_button, area='left', name='Toggle Layers')


# Create Run All button and set its color to green
run_all_button = run_all.native
run_all_button.setStyleSheet("background-color: green; color: white;")
viewer.window.add_dock_widget(run_all_button, area='right', name='Run All')

# Create Save Segmentation widget with button and checkbox
save_seg_widget = SaveSegmentationWidget()
viewer.window.add_dock_widget(save_seg_widget, area='right', name='Save Segmentation')


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
        if shapes_layer.data:
            line_coords = np.array(shapes_layer.data[-1])
            if 'merge' in event.modifiers:
                processor.merge_nuclei_with_line(line_coords)
            elif 'remove' in event.modifiers:
                processor.remove_nuclei_with_line(line_coords)
            elif 'split' in event.modifiers:
                processor.split_chromosome_with_line(line_coords)
        else:
            print("No line found for merging, removing, or splitting nuclei.")

# Start the napari event loop

viewer.window.add_dock_widget(
    chromosome_counter, 
    area='bottom',
    name='Chromosome Counter'
)

viewer.window.add_dock_widget(
    dna_fish_counter,
    area='bottom',
    name='DNA-FISH Counter'
)

viewer.window.add_dock_widget(
    cenpc_counter,
    area='bottom',
    name='CENPC Counter'
)




batch_load_button = batch_load.native
batch_load_button.setStyleSheet(BUTTON_STYLE)
viewer.window.add_dock_widget(batch_load_button, area='left', name='Batch Load')
viewer.window.add_dock_widget(folder_list_container, area='left', name='Folder List')

# Apply to other magicgui buttons
load_images.native.setStyleSheet(BUTTON_STYLE)
find_common.native.setStyleSheet(BUTTON_STYLE)
get_intensity_at_cenpc_location.native.setStyleSheet(BUTTON_STYLE)
run_all.native.setStyleSheet(BUTTON_STYLE)
batch_processing.native.setStyleSheet(BUTTON_STYLE)

def main():
    """Main entry point for the metachrome application."""
    napari.run()

if __name__ == "__main__":
    main()