
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



from skimage import draw, morphology
import scipy.ndimage as ndi
import numpy as np
import os
from scipy import ndimage  # Add this import



from napari.utils.notifications import show_info


class SegmentationPostprocessing:
    def __init__(self, viewer, processor, chromosome_counter):
        self.viewer = viewer
        self.processor = processor
        self.chromosome_counter = chromosome_counter
        self.segmentation_modified = False
        self.current_folder_path = None  # Add this line


    def set_current_folder(self, folder_path):
        """Set the current folder path"""
        self.current_folder_path = folder_path

    def _save_updated_segmentation(self, updated_labels):
        """Helper function to save updated segmentation"""
        try:
            if self.current_folder_path:
                intermediate_path = os.path.join(self.current_folder_path, "intermediate_results")
                os.makedirs(intermediate_path, exist_ok=True)
                np.save(os.path.join(intermediate_path, "segmentation.npy"), updated_labels)
                print("Saved updated segmentation")
            else:
                print("No folder path available to save segmentation")
        except Exception as e:
            print(f"Error saving segmentation: {str(e)}")


    def save_segmentation(self):
        try:
            if self.processor.nuclei is None:
                self.show_info("No segmentation to save. Please segment chromosomes first.")
                return
                
            if not self.current_folder_path:
                self.show_info("No folder selected. Please load images first.")
                return
                
            # Create intermediate_results directory if it doesn't exist
            intermediate_path = os.path.join(self.current_folder_path, "intermediate_results")
            os.makedirs(intermediate_path, exist_ok=True)
            
            # Save segmentation
            seg_file = os.path.join(intermediate_path, "segmentation.npy")
            np.save(seg_file, self.processor.nuclei)
            
            self.show_info(f"Segmentation saved to {seg_file}")
            
        except Exception as e:
            self.show_info(f"Error saving segmentation: {str(e)}")

    def remove_chromosome(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
            if len(shapes_layer.data) == 0:
                show_info("Please draw a line or select points to remove chromosomes.")
                return
                
            # Get all unique labels along the line or at points
            labels_to_remove = set()
            for shape_coords in shapes_layer.data:
                # For each shape (line or point)
                if len(shape_coords) == 1:  # Single point
                    y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                    label_value = current_labels[y, x]
                    if label_value > 0:
                        labels_to_remove.add(label_value)
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
                    labels_to_remove.update(line_labels)
            
            if len(labels_to_remove) == 0:
                show_info("No chromosomes selected for removal.")
                return
                
            # Remove selected chromosomes
            for label in labels_to_remove:
                current_labels[current_labels == label] = 0
                
            print(f"Removed chromosomes with labels: {labels_to_remove}")
            
            # Update processor and display
            self.processor.nuclei = current_labels.copy()
            self.processor.nuclei_split = current_labels.copy()
            
            # Update or create the labels layer
            if 'Segmentation updated' in self.viewer.layers:
                self.viewer.layers['Segmentation updated'].data = self.processor.nuclei
            else:
                self.viewer.add_labels(self.processor.nuclei, name='Segmentation updated')
            
            # Update chromosome count
            unique_chromosomes = len(np.unique(self.processor.nuclei)) - 1
            self.chromosome_counter.update_count(unique_chromosomes)
            
            # Set modified flag
            self.segmentation_modified = True
            
            # Clear the shapes layer
            shapes_layer.data = []
            
            show_info(f"Removed chromosomes. New total: {unique_chromosomes}")
            
        except Exception as e:
            show_info(f"Error removing chromosomes: {str(e)}")

    def split_chromosome(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
            if len(shapes_layer.data) == 0:
                show_info("Please draw a line to split chromosome.")
                return
                
            # Create line mask for the drawn line
            shape_coords = shapes_layer.data[0]  # Use first line
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
            
            if len(line_labels) != 1:
                show_info("Please draw a line through a single chromosome.")
                return
                
            label_to_split = line_labels.pop()
            
            # Create a mask of the selected chromosome
            chrom_mask = current_labels == label_to_split
            
            # Dilate the line mask slightly to ensure separation
            split_line = morphology.binary_dilation(mask, morphology.disk(1))
            
            # Remove the line from the chromosome
            chrom_mask[split_line] = False
            
            # Label the separated regions
            labeled_regions, num_regions = ndi.label(chrom_mask)
            
            if num_regions < 2:
                show_info("Split failed - please draw the line completely across the chromosome.")
                return
                
            # Update the labels
            new_label = np.max(current_labels) + 1
            region_mask = labeled_regions == 2  # Use second region
            current_labels[region_mask] = new_label
            
            # Update processor and display
            self._update_display(current_labels)
            
            print(f"Split chromosome {label_to_split} into {label_to_split} and {new_label}")
            show_info(f"Split chromosome. New total: {len(np.unique(current_labels)) - 1}")
            
            # Set modified flag
            self.segmentation_modified = True
            
        except Exception as e:
            show_info(f"Error splitting chromosome: {str(e)}")

    def merge_chromosomes(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
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
            self.processor.nuclei = current_labels.copy()
            self.processor.nuclei_split = current_labels.copy()
            
            # Update or create the labels layer
            if 'Segmentation updated' in self.viewer.layers:
                self.viewer.layers['Segmentation updated'].data = self.processor.nuclei
            else:
                self.viewer.add_labels(self.processor.nuclei, name='Segmentation updated')
            
            # Update chromosome count
            unique_chromosomes = len(np.unique(self.processor.nuclei)) - 1
            self.chromosome_counter.update_count(unique_chromosomes)
            
            # Set modified flag
            self.segmentation_modified = True
            
            # Clear the shapes layer
            shapes_layer.data = []
            
            show_info(f"Merged chromosomes. New total: {unique_chromosomes}")
            
        except Exception as e:
            show_info(f"Error merging chromosomes: {str(e)}")


    def _update_display(self, current_labels):
        """
        Helper method to update the display and processor data
        """
        self.processor.nuclei = current_labels.copy()
        self.processor.nuclei_split = current_labels.copy()
        
        if 'Segmentation updated' in self.viewer.layers:
            self.viewer.layers['Segmentation updated'].data = self.processor.nuclei
        else:
            self.viewer.add_labels(self.processor.nuclei, name='Segmentation updated')
        
        unique_chromosomes = len(np.unique(self.processor.nuclei)) - 1
        self.chromosome_counter.update_count(unique_chromosomes)
        
        self.viewer.layers['Shapes'].data = []

    def show_info(self, message):
        """
        Helper method to show information messages
        """
        print(message)
