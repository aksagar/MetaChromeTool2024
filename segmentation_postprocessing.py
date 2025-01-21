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

class SegmentationPostprocessing:
    def __init__(self, viewer, processor, chromosome_counter):
        """
        Initialize the SegmentationPostprocessing class.
        
        Parameters:
        -----------
        viewer : napari.Viewer
            The napari viewer instance
        processor : ImageProcessor
            The processor instance containing nuclei data
        chromosome_counter : ChromosomeCounter
            The counter instance for updating chromosome counts
        """
        self.viewer = viewer
        self.processor = processor
        self.chromosome_counter = chromosome_counter

    def split_chromosome(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                self.show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
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
                    
                else:  # Polygon or complex path
                    print("Processing polygon/complex path")
                    for i in range(len(shape_coords) - 1):
                        start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                        end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                        
                        rr, cc = draw.line(start_y, start_x, end_y, end_x)
                        valid_coords = (rr < current_labels.shape[0]) & (cc < current_labels.shape[1])
                        rr, cc = rr[valid_coords], cc[valid_coords]
                        mask[rr, cc] = True
                
                cross = np.array([[0,1,0],
                                [1,1,1],
                                [0,1,0]], dtype=bool)
                mask = morphology.dilation(mask, cross)
                
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
            
            self._update_display(current_labels)
            
        except Exception as e:
            self.show_info(f"Error splitting chromosomes: {str(e)}")
            raise e

    def remove_chromosome(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                self.show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
            for shape_coords in shapes_layer.data:
                mask = np.zeros_like(current_labels, dtype=bool)
                y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                mask[y, x] = True
                
                label_value = current_labels[y, x]
                
                if label_value > 0:
                    current_labels[current_labels == label_value] = 0
                    print(f"Removed chromosome with label: {label_value}")
            
            self._update_display(current_labels)
            
        except Exception as e:
            self.show_info(f"Error removing chromosome: {str(e)}")

    def merge_chromosomes(self):
        try:
            shapes_layer = self.viewer.layers['Shapes']
            
            if self.processor.nuclei is None:
                self.show_info("No segmented chromosomes found. Please segment chromosomes first.")
                return
                
            current_labels = self.processor.nuclei.copy()
            
            if len(shapes_layer.data) == 0:
                self.show_info("Please draw a line or select points to merge chromosomes.")
                return
            
            labels_to_merge = set()
            for shape_coords in shapes_layer.data:
                if len(shape_coords) == 1:  # Single point
                    y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                    label_value = current_labels[y, x]
                    if label_value > 0:
                        labels_to_merge.add(label_value)
                else:  # Line or multiple points
                    mask = np.zeros_like(current_labels, dtype=bool)
                    for i in range(len(shape_coords) - 1):
                        start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                        end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                        
                        rr, cc = draw.line(start_y, start_x, end_y, end_x)
                        valid_coords = (rr < current_labels.shape[0]) & (cc < current_labels.shape[1])
                        rr, cc = rr[valid_coords], cc[valid_coords]
                        mask[rr, cc] = True
                    
                    line_labels = set(current_labels[mask])
                    line_labels.discard(0)
                    labels_to_merge.update(line_labels)
            
            if len(labels_to_merge) < 2:
                self.show_info("Please select different chromosomes to merge.")
                return
            
            target_label = min(labels_to_merge)
            for label in labels_to_merge:
                current_labels[current_labels == label] = target_label
                
            print(f"Merged chromosomes with labels: {labels_to_merge} to label {target_label}")
            
            self._update_display(current_labels)
            
        except Exception as e:
            self.show_info(f"Error merging chromosomes: {str(e)}")

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