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
from image_processor import ImageProcessor
from batch_processor import BatchProcessor
from segmentation_postprocessing import SegmentationPostprocessing

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
        show_info("Threshold slider changed, reset spot detection flag for DNA-FISH channel")


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
        show_info("Threshold slider changed, reset spot detection flag for DNA-FISH channel")
class ControlWidgetCENPC(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)
        self.slider.setValue(40)
        self.slider.setSingleStep(1)
        self.slider.valueChanged.connect(self.reset_cenpc_flag)
        self.layout.addWidget(self.slider)

        self.detect_cenpc_spots_button = detect_cenpc_spots.native
        self.detect_cenpc_spots_button.setStyleSheet(BUTTON_STYLE)
        self.layout.addWidget(self.detect_cenpc_spots_button)

        self.setLayout(self.layout)

    def reset_cenpc_flag(self):
        global detect_cenpc_done
        detect_cenpc_done = False
        show_info("Threshold slider changed, reset spot detection flag for CENPC channel")

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
                            name="DNA-FISH Spots",
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
                            name="CENPC Spots",
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
        show_info("Spot detection for DNA-FISH has already been done.")
        return

    threshold = control_widget_dna_fish.slider.value() / 100
    if images[1] is not None:
        try:
            # Get centroids based on segmentation mode
            centroids = processor.detect_spots_cent(images[1], 'DNA-FISH', threshold, save_dir=current_folder_path)
            
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
                    if layer.name in ["DNA-FISH Spots", "Centroids in DNA-FISH"]:
                        viewer.layers.remove(layer)

                # Add new visualization layers
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="yellow",
                    face_color=[1, 1, 0, 0.2],
                    edge_width=2,
                    name="DNA-FISH Spots",
                    opacity=0.8
                )
                detect_dna_fish_done = True
                show_info(f"Detected {spot_count} spots in DNA-FISH image with threshold {threshold}")
            else:
                show_info("No spots detected in DNA-FISH image")
                detect_dna_fish_done = False
                dna_fish_counter.update_count(0)
                
        except Exception as e:
            show_info(f"Error detecting spots: {str(e)}")
            detect_dna_fish_done = False
            dna_fish_counter.update_count(0)
    else:
        show_info("DNA-FISH image not loaded")


@magicgui(call_button="Delete channel 1 Spots")
def delete_dna_fish_spots():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.dna_fish_centroids is None:
            show_info("No DNA-FISH spots detected yet")
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
        show_info("No DNA-FISH spots to save")
        return
        
    if current_folder_path:
        intermediate_path = os.path.join(current_folder_path, "intermediate_results")
        os.makedirs(intermediate_path, exist_ok=True)
        np.save(os.path.join(intermediate_path, "dna_fish_centroids.npy"), processor.dna_fish_centroids)
        show_info("DNA-FISH spots saved successfully")
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
            if 'DNA-FISH Spots' in layer.name:
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
                name="DNA-FISH Spots",
                opacity=0.8
            )


@magicgui(call_button="Detect channel 2 Spots")
def detect_cenpc_spots():
    global detect_cenpc_done, images
    if detect_cenpc_done:
        show_info("Spot detection for CENPC has already been done.")
        return

    threshold = control_widget_cenpc.slider.value() / 100
    if images[2] is not None:
        try:
            centroids = None
            # Get centroids based on segmentation mode
            if segment_dapi_widget.is_checked():
                processor.detect_spots_cent(images[2], 'CENPC', threshold, save_dir=current_folder_path)
                centroids = processor.cenpc_centroids
            else:
                processor.detect_spots_cent(images[2], 'CENPC', threshold, save_dir=current_folder_path)
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
                    if layer.name in ["CENPC Spots", "Centroids in CENPC"]:
                        viewer.layers.remove(layer)

                # Add new visualization layer
                viewer.add_shapes(
                    squares,
                    shape_type='polygon',
                    edge_color="skyblue",
                    face_color=[0, 0.5, 1, 0.2],
                    edge_width=2,
                    name="CENPC Spots",
                    opacity=0.8
                )

            show_info(f"Detected and labeled spots in CENPC image with threshold {threshold}")
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
    except Exception as e:
        show_info(f"Error finding common labels: {e}")

@magicgui(call_button="Get Intensity at channel 2 Location")
def get_intensity_at_cenpc_location():
    try:
        if segment_dapi_widget.is_checked():
            # Case 1: No segmentation (checkbox is checked)
            show_info("Calculating intensity at all DNA-FISH locations without segmentation.")
            df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
        else:
            # Case 2: With segmentation (checkbox is unchecked)
            show_info("Calculating intensity at DNA-FISH locations in segmented regions.")
            if processor.nuclei is None:
                show_info("Please segment the image first")
                return
                
            if processor.df_centroid_dna_fish is None:
                show_info("Please detect DNA-FISH spots first")
                return
                
            if processor.img_cenpc is None:
                show_info("CENPC image not found")
                return
                
            df_with_cenpc_inten = processor.gen_intensity_from_df(
                processor.img_cenpc, 
                processor.df_centroid_dna_fish
            )

        if df_with_cenpc_inten is not None and not df_with_cenpc_inten.empty:
            # Save results
            folder_name = os.path.basename(current_folder_path)
            save_path = os.path.join(current_folder_path, f"{folder_name}_intensity.csv")
            df_with_cenpc_inten.to_csv(save_path, index=False)
            show_info(f"Intensity data saved to: {save_path}")
            print("\nIntensity measurements:")
            print(df_with_cenpc_inten)
        else:
            show_info("No intensity measurements found")
            
    except Exception as e:
        show_info(f"Error calculating intensity: {str(e)}")

        
@magicgui(call_button="Run All")
def run_all():
    global images
    try:
        threshold_dna_fish = control_widget_dna_fish.slider.value() / 100
        threshold_cenpc = control_widget_cenpc.slider.value() / 100

        if all(img is not None for img in images if img is not None):
            # Case 1: No Segmentation (checkbox checked)
            if segment_dapi_widget.is_checked():
                # Detect spots without segmentation
                processor.detect_spots_cent(images[0], 'DNA-FISH', threshold_dna_fish, save_dir=current_folder_path)
                processor.detect_spots_cent(images[1], 'CENPC', threshold_cenpc, save_dir=current_folder_path)
                
                # Update spot counters
                if processor.dna_fish_centroids is not None:
                    dna_fish_counter.update_count(len(processor.dna_fish_centroids))
                if processor.cenpc_centroids is not None:
                    cenpc_counter.update_count(len(processor.cenpc_centroids))
                
                # Calculate intensities
                df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()

            # Case 2: With Segmentation (checkbox unchecked)
            else:
                # Segment DAPI
                masks = processor.segment_image(images[0], save_dir=current_folder_path)
                viewer.add_labels(masks, name="Cellpose Segmented")
                chromosome_counter.update_count(len(np.unique(masks)) - 1)
                
                # Detect spots
                processor.detect_spots_cent(images[1], 'DNA-FISH', threshold_dna_fish, save_dir=current_folder_path)
                processor.detect_spots_cent(images[2], 'CENPC', threshold_cenpc, save_dir=current_folder_path)

                # Update spot counters
                if processor.dna_fish_centroids is not None:
                    dna_fish_counter.update_count(len(processor.dna_fish_centroids))
                if processor.cenpc_centroids is not None:
                    cenpc_counter.update_count(len(processor.cenpc_centroids))

                # Find common regions (saving is handled within find_common)
                find_common()

                # Calculate intensities
                df_with_cenpc_inten = processor.gen_intensity_from_df(
                    processor.img_cenpc,
                    processor.df_centroid_dna_fish
                )

            # Save intensity results
            if df_with_cenpc_inten is not None and not df_with_cenpc_inten.empty:
                folder_name = os.path.basename(current_folder_path)
                save_path = os.path.join(current_folder_path, f"{folder_name}_intensity.csv")
                df_with_cenpc_inten.to_csv(save_path, index=False)
                show_info(f"Intensity data saved to: {save_path}")
                print("\nIntensity measurements:")
                print(df_with_cenpc_inten)
            else:
                show_info("No intensity measurements found")

            show_info("Run all processing completed")
        else:
            show_info("Ensure that all images are loaded")
    except Exception as e:
        show_info(f"Error during run all processing: {e}")


@magicgui(call_button="Batch Load")
def batch_load():
    root_folder = QFileDialog.getExistingDirectory(caption='Select Root Folder for Batch Loading')
    if root_folder:
        try:
            # Clear previous items
            folder_list_widget.clear()
            
            # Store full paths but display only folder names
            for folder_name in os.listdir(root_folder):
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
                    viewer.add_image(images[1], name='DNA-FISH')
                    viewer.add_image(images[2], name='CENPC')
                else:
                    viewer.add_image(images[0], name='DAPI')
                    viewer.add_image(images[1], name='DNA-FISH')
                    viewer.add_image(images[2], name='CENPC')
                
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
                                name="DNA-FISH Spots",
                                opacity=0.8
                            )
                            dna_fish_counter.update_count(len(processor.dna_fish_centroids))
                            detect_dna_fish_done = True
                            show_info("Loaded existing DNA-FISH spots")
                    
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
                                name="CENPC Spots",
                                opacity=0.8
                            )
                            cenpc_counter.update_count(len(processor.cenpc_centroids))
                            detect_cenpc_done = True
                            show_info("Loaded existing CENPC spots")
                else:
                    segment_done = False
                    detect_dna_fish_done = False
                    detect_cenpc_done = False
                
                show_info(f"Loaded images from: {os.path.basename(current_folder_path)}")
                
    except Exception as e:
        show_info(f"Error loading selected folder: {str(e)}")

def on_folder_selected_BU2():
    try:
        current_item = folder_list_widget.currentItem()
        if current_item:
            # Get the full path from item data
            selected_folder = current_item.data(Qt.UserRole)
            if selected_folder:
                # Rest of your existing loading code...
                dapi_id = channel_identifiers.dapi_text.text()
                dna_fish_id = channel_identifiers.dna_fish_text.text()
                cenpc_id = channel_identifiers.cenpc_text.text()
                
                global segment_done, detect_dna_fish_done, detect_cenpc_done, current_folder_path, images
                current_folder_path = selected_folder
                
                images = processor.load_images(selected_folder, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
                
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
                show_info(f"Loaded images from: {os.path.basename(selected_folder)}")
                
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
            'dna_fish_count': len(processor.dna_fish_centroids) if processor.dna_fish_centroids is not None else 0,
            'cenpc_count': len(processor.cenpc_centroids) if processor.cenpc_centroids is not None else 0,
            'mean_intensity': df_with_cenpc_inten['CENPC_Intensity'].mean(),
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
        all_intensities = []
        
        folders_to_process = [
            folder_list_widget.item(i).data(Qt.UserRole)
            for i in range(folder_list_widget.count())
        ]
        if not folders_to_process:
            show_info("No folders in the list to process")
            return
            
        if use_current_settings:
            threshold_dna_fish = control_widget_dna_fish.slider.value() / 100
            threshold_cenpc = control_widget_cenpc.slider.value() / 100
            
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
                    dapi_id = channel_identifiers.dapi_text.text()
                    dna_fish_id = channel_identifiers.dna_fish_text.text()
                    cenpc_id = channel_identifiers.cenpc_text.text()
                    images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
                    
                    chromosome_count = 0
                    matched_count = 0
                    df_with_cenpc_inten = None
                    
                    if segment_dapi_widget.is_checked():  # Skip Segmentation case
                        # Detect spots without segmentation
                        processor.detect_spots_cent(images[0], 'DNA-FISH', threshold_dna_fish, save_dir=current_folder_path)
                        processor.detect_spots_cent(images[1], 'CENPC', threshold_cenpc, save_dir=current_folder_path)
                        
                        # Calculate intensities
                        df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()
                        
                    else:  # With Segmentation case
                        # Segment DAPI
                        masks = processor.segment_image(images[0], save_dir=current_folder_path)
                        chromosome_count = len(np.unique(masks)) - 1
                        
                        # Detect spots
                        processor.detect_spots_cent(images[1], 'DNA-FISH', threshold_dna_fish, save_dir=current_folder_path)
                        processor.detect_spots_cent(images[2], 'CENPC', threshold_cenpc, save_dir=current_folder_path)
                        
                        # Find common regions
                        find_common()
                        if processor.common_nuclei is not None:
                            matched_count = len(np.unique(processor.common_nuclei)) - 1
                            df_with_cenpc_inten = processor.gen_intensity_from_df(
                                processor.img_cenpc,
                                processor.df_centroid_dna_fish
                            )
                    
                    # Add to summaries
                    summary_row = {
                        'Folder': folder_name,
                        'DNA-FISH Spots': len(processor.dna_fish_centroids) if processor.dna_fish_centroids is not None else 0,
                        'CENPC Spots': len(processor.cenpc_centroids) if processor.cenpc_centroids is not None else 0,
                        'Mean CENPC Intensity': df_with_cenpc_inten['CENPC_Intensity'].mean() if df_with_cenpc_inten is not None else None,
                        'Chromosome Count': chromosome_count,
                        'Matched Nuclei Count': matched_count,
                        'Skip Segmentation': segment_dapi_widget.is_checked(),
                        'Using UI Settings': True
                    }
                    summary_data.append(summary_row)
                    
                    # Save intensity results if available
                    if df_with_cenpc_inten is not None and not df_with_cenpc_inten.empty:
                        df_with_cenpc_inten['Folder'] = folder_name
                        all_intensities.append(df_with_cenpc_inten)
                        save_path = os.path.join(folder_path, f"{folder_name}_intensity.csv")
                        df_with_cenpc_inten.to_csv(save_path, index=False)
                    
                except Exception as e:
                    print(f"Error processing {folder_name}: {str(e)}")
                    summary_row = {
                        'Folder': folder_name,
                        'DNA-FISH Spots': 0,
                        'CENPC Spots': 0,
                        'Mean CENPC Intensity': None,
                        'Chromosome Count': 0,
                        'Matched Nuclei Count': 0,
                        'Skip Segmentation': segment_dapi_widget.is_checked(),
                        'Using UI Settings': True,
                        'Error': str(e)
                    }
                    summary_data.append(summary_row)
                    continue
        
        # Save summary files
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_path = os.path.join(os.path.dirname(folders_to_process[0]), "batch_processing_summary.csv")
            summary_df.to_csv(summary_path, index=False)
            
            if all_intensities:
                combined_df = pd.concat(all_intensities, ignore_index=True)
                combined_path = os.path.join(os.path.dirname(folders_to_process[0]), "all_intensities_summary.csv")
                combined_df.to_csv(combined_path, index=False)
                
                show_info(f"Processing complete. Summaries saved to:\n{summary_path}\n{combined_path}")
            else:
                show_info("No intensity data found to combine")
        else:
            show_info("No folders were successfully processed")
            
    except Exception as e:
        show_info(f"Error during batch processing: {str(e)}")               

chromosome_counter = ChromosomeCountWidget()
dna_fish_counter = SpotCountWidget("DNA-FISH")
cenpc_counter = SpotCountWidget("CENPC")



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
viewer.window.add_dock_widget(control_widget_dna_fish, area='right', name='Detect DNA-FISH Spot Control')
viewer.window.add_dock_widget(control_widget_cenpc, area='right', name='Detect CENPC Spot Control')
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

napari.run()