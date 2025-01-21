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
class ControlWidgetDNAFISH(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.slider = QLabeledSlider(orientation='horizontal')
        self.slider.setRange(0, 100)
        self.slider.setValue(40)
        self.slider.setSingleStep(1)
        self.slider.valueChanged.connect(self.reset_dna_fish_flag)
        self.layout.addWidget(self.slider)

        self.detect_dna_fish_spots_button = detect_dna_fish_spots.native
        self.detect_dna_fish_spots_button.setStyleSheet(BUTTON_STYLE)
        self.layout.addWidget(self.detect_dna_fish_spots_button)

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
class ChannelIdentifiers(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()  # Initialize the layout first
        
        self.dapi_label = QLabel("DAPI Channel Identifier:")
        self.dapi_text = QLineEdit()
        self.dapi_text.setText("435")  # Default wavelength for DAPI
        self.layout.addWidget(self.dapi_label)
        self.layout.addWidget(self.dapi_text)

        self.dna_fish_label = QLabel("DNA-FISH Channel Identifier:")
        self.dna_fish_text = QLineEdit()
        self.dna_fish_text.setText("525")  # Default wavelength for DNA-FISH
        self.layout.addWidget(self.dna_fish_label)
        self.layout.addWidget(self.dna_fish_text)

        self.cenpc_label = QLabel("CENPC Channel Identifier:")
        self.cenpc_text = QLineEdit()
        self.cenpc_text.setText("679")  # Default wavelength for CENPC
        self.layout.addWidget(self.cenpc_label)
        self.layout.addWidget(self.cenpc_text)

        self.setLayout(self.layout)  # Set the layout for the widget

# Add a checkbox beside the segment DAPI button
@magicgui(call_button="Segment DAPI Image")
def segment_image1():
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

@magicgui(call_button="Segment DAPI Image")
def segment_image():
    global segment_done, images
    if segment_done:
        show_info("Segmentation has already been done.")
        return

    if images[0] is not None:
        try:
            # Pass current folder path to save results
            masks = processor.segment_image(images[0], save_dir=current_folder_path)
            viewer.add_labels(masks, name="Cellpose Segmented DAPI")
            
            # Update chromosome count
            unique_chromosomes = len(np.unique(masks)) - 1
            chromosome_counter.update_count(unique_chromosomes)
            
            show_info("Segmented DAPI image using Cellpose")
            segment_done = True
        except Exception as e:
            show_info(f"Error segmenting image: {e}")


def segment_image2():
    global segment_done, images
    if segment_done:
        show_info("Segmentation has already been done.")
        return

    if images[0] is not None:
        try:
            masks = processor.segment_image(images[0])
            viewer.add_labels(masks, name="Cellpose Segmented DAPI")
            
            # Update chromosome count
            unique_chromosomes = len(np.unique(masks)) - 1  # Subtract 1 to exclude background
            chromosome_counter.update_count(unique_chromosomes)
            
            show_info("Segmented DAPI image using Cellpose")
            segment_done = True
        except Exception as e:
            show_info(f"Error segmenting image: {e}")
def remove_chromosome():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.nuclei is None:
            show_info("No segmented chromosomes found. Please segment chromosomes first.")
            return
            
        current_labels = processor.nuclei.copy()
        
        # Process each shape
        for shape_coords in shapes_layer.data:
            # Create point mask
            mask = np.zeros_like(current_labels, dtype=bool)
            y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
            mask[y, x] = True
            
            # Get label at clicked point
            label_value = current_labels[y, x]
            
            if label_value > 0:
                # Remove the chromosome
                current_labels[current_labels == label_value] = 0
                print(f"Removed chromosome with label: {label_value}")
        
        # Update processor and display
        processor.nuclei = current_labels.copy()
        processor.nuclei_split = current_labels.copy()
        
        # Update or create the labels layer
        if 'Segmentation updated' in viewer.layers:
            viewer.layers['Segmentation updated'].data = processor.nuclei
        else:
            viewer.add_labels(processor.nuclei, name='Segmentation updated')
        
        # Update chromosome count
        unique_chromosomes = len(np.unique(processor.nuclei)) - 1
        chromosome_counter.update_count(unique_chromosomes)
        
        # Clear the shapes layer
        shapes_layer.data = []
        
        show_info(f"Removed chromosome. New total: {unique_chromosomes}")
        
    except Exception as e:
        show_info(f"Error removing chromosome: {str(e)}")

def merge_chromosomes():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.nuclei is None:
            show_info("No segmented chromosomes found. Please segment chromosomes first.")
            return
            
        current_labels = processor.nuclei.copy()
        
        if len(shapes_layer.data) == 0:
            show_info("Please draw a line or select points to merge chromosomes.")
            return
            
        # Get all unique labels along the line or at points
        labels_to_merge = set()
        for shape_coords in shapes_layer.data:
            # For each shape (line or point)
            if len(shape_coords) == 1:  # Single point
                y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                label_value = current_labels[y, x]
                if label_value > 0:
                    labels_to_merge.add(label_value)
            else:  # Line or multiple points
                # Create line mask
                mask = np.zeros_like(current_labels, dtype=bool)
                for i in range(len(shape_coords) - 1):
                    start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                    end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                    
                    rr, cc = draw.line(start_y, start_x, end_y, end_x)
                    valid_coords = (rr < current_labels.shape[0]) & (cc < current_labels.shape[1])
                    rr, cc = rr[valid_coords], cc[valid_coords]
                    mask[rr, cc] = True
                
                # Get unique labels along the line
                line_labels = set(current_labels[mask])
                line_labels.discard(0)  # Remove background label
                labels_to_merge.update(line_labels)
        
        if len(labels_to_merge) < 2:
            show_info("Please select different chromosomes to merge.")
            return
            
        # Merge to the smallest label value
        target_label = min(labels_to_merge)
        for label in labels_to_merge:
            current_labels[current_labels == label] = target_label
            
        print(f"Merged chromosomes with labels: {labels_to_merge} to label {target_label}")
        
        # Update processor and display
        processor.nuclei = current_labels.copy()
        processor.nuclei_split = current_labels.copy()
        
        # Update or create the labels layer
        if 'Segmentation updated' in viewer.layers:
            viewer.layers['Segmentation updated'].data = processor.nuclei
        else:
            viewer.add_labels(processor.nuclei, name='Segmentation updated')
        
        # Update chromosome count
        unique_chromosomes = len(np.unique(processor.nuclei)) - 1
        chromosome_counter.update_count(unique_chromosomes)
        
        # Clear the shapes layer
        shapes_layer.data = []
        
        show_info(f"Merged chromosomes. New total: {unique_chromosomes}")
        
    except Exception as e:
        show_info(f"Error merging chromosomes: {str(e)}")
        
from skimage import draw, morphology
import scipy.ndimage as ndi
import numpy as np

def split_chromosome():
    try:
        shapes_layer = viewer.layers['Shapes']
        
        if processor.nuclei is None:
            show_info("No segmented chromosomes found. Please segment chromosomes first.")
            return
            
        current_labels = processor.nuclei.copy()
        
        # Process each shape drawn
        for shape_coords in shapes_layer.data:
            print(f"Shape has {len(shape_coords)} points")
            
            # Create mask for the entire path
            mask = np.zeros_like(current_labels, dtype=bool)
            
            if len(shape_coords) == 2:  # Straight line (2 points)
                start_y, start_x = int(shape_coords[0][0]), int(shape_coords[0][1])
                end_y, end_x = int(shape_coords[1][0]), int(shape_coords[1][1])
                
                rr, cc = draw.line(start_y, start_x, end_y, end_x)
                valid_coords = (rr < current_labels.shape[0]) & (cc < current_labels.shape[1])
                rr, cc = rr[valid_coords], cc[valid_coords]
                mask[rr, cc] = True
                print("Processing straight line")
                
            else:  # Polygon or complex path (more than 2 points)
                print("Processing polygon/complex path")
                for i in range(len(shape_coords) - 1):
                    start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                    end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                    
                    rr, cc = draw.line(start_y, start_x, end_y, end_x)
                    valid_coords = (rr < current_labels.shape[0]) & (cc < current_labels.shape[1])
                    rr, cc = rr[valid_coords], cc[valid_coords]
                    mask[rr, cc] = True
            
            # Use cross-shaped structuring element for minimal separation
            cross = np.array([[0,1,0],
                            [1,1,1],
                            [0,1,0]], dtype=bool)
            mask = morphology.dilation(mask, cross)
            
            # Get all labels along the path
            path_labels = set(current_labels[mask])
            path_labels.discard(0)
            print(f"Labels crossed by the path: {path_labels}")
            
            if len(path_labels) > 0:
                label_value = list(path_labels)[0]
                print(f"Attempting to split chromosome with label: {label_value}")
                
                chromosome_mask = current_labels == label_value
                temp_mask = chromosome_mask.copy()
                temp_mask[mask] = False
                
                split_regions, num_regions = ndi.label(temp_mask, structure=np.ones((3,3)))
                
                if num_regions >= 2:
                    current_labels[chromosome_mask] = 0
                    max_label = current_labels.max()
                    
                    for i in range(1, num_regions + 1):
                        region_mask = split_regions == i
                        new_label = max_label + i
                        current_labels[region_mask] = new_label
                        print(f"Created new region with label: {new_label}")
                    
                    print(f"Successfully split chromosome {label_value} into {num_regions} parts")
                else:
                    print("Path did not split the chromosome into multiple parts")
            else:
                print("Path does not cross any labeled chromosomes")
        
        # Update processor and display
        processor.nuclei = current_labels.copy()
        processor.nuclei_split = current_labels.copy()
        
        # Update or create the labels layer with new name
        if 'Segmentation updated' in viewer.layers:
            viewer.layers['Segmentation updated'].data = processor.nuclei
        else:
            viewer.add_labels(processor.nuclei, name='Segmentation updated')
        
        unique_chromosomes = len(np.unique(processor.nuclei)) - 1
        chromosome_counter.update_count(unique_chromosomes)
        
        shapes_layer.data = []
        
        show_info(f"Split chromosomes based on drawn path. New total: {unique_chromosomes}")
        
    except Exception as e:
        show_info(f"Error splitting chromosomes: {str(e)}")
        raise e

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
        
        self.merge_button = QPushButton("Merge Chromosomes")
        self.merge_button.setStyleSheet(BUTTON_STYLE)
        self.merge_button.clicked.connect(self.postprocessing.merge_chromosomes)
        self.buttons_layout.addWidget(self.merge_button)

        self.remove_button = QPushButton("Remove Chromosome")
        self.remove_button.setStyleSheet(BUTTON_STYLE)
        self.remove_button.clicked.connect(self.postprocessing.remove_chromosome)
        self.buttons_layout.addWidget(self.remove_button)

        self.split_button = QPushButton("Split Chromosome")
        self.split_button.setStyleSheet(BUTTON_STYLE)
        self.split_button.clicked.connect(self.postprocessing.split_chromosome)
        self.buttons_layout.addWidget(self.split_button)
        
        self.layout.addLayout(self.buttons_layout)
        self.setLayout(self.layout)


    def is_checked(self):
        return self.checkbox.isChecked()




class SegmentDAPIWidget1(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()  # Changed to QVBoxLayout to stack elements vertically
        self.postprocessing = SegmentationPostprocessing(self.viewer, self.processor, self.chromosome_counter)

        # Segment button and checkbox
        self.segment_layout = QHBoxLayout()
        self.segment_button = segment_image.native
        self.checkbox = QCheckBox("Skip Segmentation")
        self.segment_layout.addWidget(self.segment_button)
        self.segment_layout.addWidget(self.checkbox)
        self.layout.addLayout(self.segment_layout)

        # Buttons for merging, removing, and splitting chromosomes
        self.buttons_layout = QHBoxLayout()
        self.merge_button = QPushButton("Merge Chromosomes")
        self.merge_button.clicked.connect(merge_chromosomes)
        self.buttons_layout.addWidget(self.merge_button)

        self.remove_button = QPushButton("Remove Chromosome")
        self.remove_button.clicked.connect(remove_chromosome)
        self.buttons_layout.addWidget(self.remove_button)

        # Split Chromosome button
        self.split_button = QPushButton("Split Chromosome")
        self.split_button.clicked.connect(split_chromosome)  # Connect to the external function
        self.buttons_layout.addWidget(self.split_button)
        
        self.layout.addLayout(self.buttons_layout)
        self.setLayout(self.layout)

    def is_checked(self):
        return self.checkbox.isChecked()




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

@magicgui(call_button="Detect DNA-FISH Spots")
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


@magicgui(call_button="Detect CENPC Spots")
def detect_cenpc_spots():
    global detect_cenpc_done, images, current_folder_path
    if detect_cenpc_done:
        show_info("Spot detection for CENPC has already been done.")
        return

    threshold = control_widget_cenpc.slider.value() / 100
    if images[2] is not None:
        try:
            # Get centroids based on segmentation mode
            centroids = processor.detect_spots_cent(images[2], 'CENPC', threshold, save_dir=current_folder_path)
            
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
                detect_cenpc_done = True
                show_info(f"Detected {spot_count} spots in CENPC image with threshold {threshold}")
            else:
                show_info("No spots detected in CENPC image")
                detect_cenpc_done = False
                cenpc_counter.update_count(0)
                
        except Exception as e:
            show_info(f"Error detecting spots: {str(e)}")
            detect_cenpc_done = False
            cenpc_counter.update_count(0)
    else:
        show_info("CENPC image not loaded")



@magicgui(call_button="Detect CENPC Spots")
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




def detect_dna_fish_spots_BU2():
    global detect_dna_fish_done, images
    if detect_dna_fish_done:
        show_info("Spot detection for DNA-FISH has already been done.")
        return

    threshold = control_widget_dna_fish.slider.value() / 100
    if images[1] is not None:
        try:
            # Get centroids based on segmentation mode
            centroids = None
            if segment_dapi_widget.is_checked():
                spots, labels = processor.detect_spots_no_segmentation(images[1], threshold, channel='DNA-FISH')
                centroids = processor.dna_fish_centroids
            else:
                processor.detect_spots_cent(images[1], 'DNA-FISH', threshold)
                centroids = processor.dna_fish_centroids

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

            show_info(f"Detected and labeled spots in DNA-FISH image with threshold {threshold}")
            detect_dna_fish_done = True
        except Exception as e:
            show_info(f"Error detecting spots: {e}")
#@magicgui(call_button="Detect DNA-FISH Spots")
def detect_dna_fish_spots_BU():
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

def detect_cenpc_spots_BU2():
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
                spots, labels = processor.detect_spots_no_segmentation(images[2], threshold, channel='CENPC')
                centroids = processor.cenpc_centroids
            else:
                processor.detect_spots_cent(images[2], 'CENPC', threshold)
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

                print(f"CENPC spots detected: {spot_count}")

            show_info(f"Detected and labeled spots in CENPC image with threshold {threshold}")
            detect_cenpc_done = True
        except Exception as e:
            show_info(f"Error detecting spots: {e}")
            print(f"Exception details: {str(e)}")


def detect_cenpc_spots_BU():
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
                processor.detect_spots_cent(images[1], 'DNA-FISH', threshold_dna_fish, save_dir=current_folder_path)
                processor.detect_spots_cent(images[2], 'CENPC', threshold_cenpc, save_dir=current_folder_path)
                
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

                # Process results
                common_nuclei = processor.find_common(threshold_dna_fish, threshold_cenpc)
                if common_nuclei is not None:
                    viewer.add_labels(common_nuclei, name="Matched Nuclei")
                df_with_cenpc_inten = processor.gen_intensity_from_df(processor.img_cenpc, processor.df_centroid_dna_fish)

            # Save results
            folder_name = os.path.basename(current_folder_path)
            save_path = os.path.join(current_folder_path, f"{folder_name}_intensity.csv")
            df_with_cenpc_inten.to_csv(save_path, index=False)

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
    # ... existing loading code ...
    
    # Check for and load intermediate results
    intermediate_path = os.path.join(selected_folder, "intermediate_results")
    if os.path.exists(intermediate_path):
        # Load segmentation (existing code)
        
        # Load DNA-FISH spots
        dna_fish_file = os.path.join(intermediate_path, "dna_fish_spots.npy")
        dna_fish_centroids_file = os.path.join(intermediate_path, "dna_fish_centroids.npy")
        if os.path.exists(dna_fish_file) and os.path.exists(dna_fish_centroids_file):
            processor.labels_dna_fish = np.load(dna_fish_file)
            processor.dna_fish_centroids = np.load(dna_fish_centroids_file)
            detect_dna_fish_done = True
            
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
        
        # Load CENPC spots
        cenpc_file = os.path.join(intermediate_path, "cenpc_spots.npy")
        cenpc_centroids_file = os.path.join(intermediate_path, "cenpc_centroids.npy")
        if os.path.exists(cenpc_file) and os.path.exists(cenpc_centroids_file):
            processor.labels_cenpc = np.load(cenpc_file)
            processor.cenpc_centroids = np.load(cenpc_centroids_file)
            detect_cenpc_done = True
            
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

def on_folder_selected():
    try:
        current_item = folder_list_widget.currentItem()
        if current_item:
            # Get the full path from item data
            selected_folder = current_item.data(Qt.UserRole)
            if selected_folder:
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



@magicgui(call_button="Batch Processing")
def batch_processing():
    root_folder = QFileDialog.getExistingDirectory(caption='Select Root Folder')
    if root_folder:
        try:
            # Get channel identifiers and thresholds
            dapi_id = channel_identifiers.dapi_text.text()
            dna_fish_id = channel_identifiers.dna_fish_text.text()
            cenpc_id = channel_identifiers.cenpc_text.text()
            threshold_dna_fish = control_widget_dna_fish.slider.value() / 100
            threshold_cenpc = control_widget_cenpc.slider.value() / 100

            # Process each subfolder
            for folder_name in os.listdir(root_folder):
                folder_path = os.path.join(root_folder, folder_name)
                if not os.path.isdir(folder_path):
                    continue

                try:
                    print(f"Processing folder: {folder_path}")
                    # Load images
                    images = processor.load_images(folder_path, dapi_id, dna_fish_id, cenpc_id, segment_dapi_widget.is_checked())
                    if images is None:
                        print(f"Skipping {folder_path} - missing required images")
                        continue

                    # Create intermediate results directory
                    intermediate_path = os.path.join(folder_path, "intermediate_results")
                    os.makedirs(intermediate_path, exist_ok=True)

                    # Case 1: No Segmentation (checkbox checked)
                    if segment_dapi_widget.is_checked():
                        # Detect spots without segmentation
                        processor.detect_spots_cent(images[1], 'DNA-FISH', threshold_dna_fish, save_dir=folder_path)
                        processor.detect_spots_cent(images[2], 'CENPC', threshold_cenpc, save_dir=folder_path)
                        
                        # Calculate intensities
                        df_with_cenpc_inten = processor.calculate_intensity_all_dna_fish()

                    # Case 2: With Segmentation (checkbox unchecked)
                    else:
                        # Segment DAPI
                        masks = processor.segment_image(images[0], save_dir=folder_path)
                        
                        # Detect spots
                        processor.detect_spots_cent(images[1], 'DNA-FISH', threshold_dna_fish, save_dir=folder_path)
                        processor.detect_spots_cent(images[2], 'CENPC', threshold_cenpc, save_dir=folder_path)
                        
                        # Process results
                        common_nuclei = processor.find_common(threshold_dna_fish, threshold_cenpc)
                        df_with_cenpc_inten = processor.gen_intensity_from_df(processor.img_cenpc, processor.df_centroid_dna_fish)

                    # Save results
                    save_path = os.path.join(folder_path, f"{folder_name}_intensity.csv")
                    df_with_cenpc_inten.to_csv(save_path, index=False)
                    print(f"Processed {folder_path}")
                    print(f"Found {len(df_with_cenpc_inten)} spots")

                except Exception as e:
                    print(f"Error processing {folder_path}: {str(e)}")
                    continue

            show_info("Batch processing completed")
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
