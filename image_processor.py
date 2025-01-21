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
import numpy as np
import pandas as pd
from skimage.io import imread
from skimage import morphology as morph
from scipy import ndimage as ndi
from cellpose import models

class ImageProcessor:
    def __init__(self):
        self.img_cenpc = None
        self.dna_fish_spots = None  # Placeholder for DNA-FISH spot locations
        self.cenpc_spots = None  # Placeholder for CENPC spot locations
        self.normalizedDNAFISH = None  # Placeholder for normalized DNA-FISH image
        self.normalizedCENPC = None  # Placeholder for normalized CENPC image
        self.df_centroid_dna_fish = None  # Placeholder for DNA-FISH centroids
        self.df_centroid_cenpc = None  # Placeholder for CENPC centroids
        self.nuclei = None  # Placeholder for segmented nuclei
        self.nuclei_merged = None  # Add this line
        self.spotLabelsDNAFISH = None  # Placeholder for DNA-FISH spot labels
        self.spotLabelsCENPC = None  # Placeholder for CENPC spot labels
        self.labels_dna_fish = None  # Placeholder for DNA-FISH labels from no-segmentation detection
        self.dna_fish_centroids = None
        self.cenpc_centroids = None



    def load_images(self, folder_path, dapi_id, dna_fish_id, cenpc_id, skip_segmentation):
        image_files = [f for f in os.listdir(folder_path) if f.endswith(('.tif', '.png', '.jpg'))]
        sorted_files = []
        
        # Check if we are skipping segmentation
        required_suffixes = (dna_fish_id, cenpc_id) if skip_segmentation else (dapi_id, dna_fish_id, cenpc_id)
        
        for suffix in required_suffixes:
            for f in image_files:
                name, ext = os.path.splitext(f)
                if name.endswith(suffix):
                    sorted_files.append(f)
        print(sorted_files)
        
        if (skip_segmentation and len(sorted_files) != 2) or (not skip_segmentation and len(sorted_files) != 3):
            print(f"Warning: {folder_path} does not contain the required images. Skipping this folder.")
            return None

        images = []
        for image_file in sorted_files:
            file_path = os.path.join(folder_path, image_file)
            image = imread(file_path)
            images.append(image)
            if cenpc_id in image_file:
                self.img_cenpc = image
        
        return images



    def load_images1(self, folder_path):
        image_files = [f for f in os.listdir(folder_path) if f.endswith(('.tif', '.png', '.jpg'))]
        sorted_files = []
        for suffix in ('435.tif', '525.tif', '679.tif'):
            sorted_files.extend([f for f in image_files if f.endswith(suffix)])
        
        images = []
        for image_file in sorted_files:
            file_path = os.path.join(folder_path, image_file)
            image = imread(file_path)
            images.append((image, file_path))
            if '679' in image_file:
                self.img679 = image
        
        return images

    def segment_image_BU(self, image):
        model = models.Cellpose(model_type='cyto')
        masks, _, _, _ = model.eval(image, diameter=None, channels=[0, 0])
        self.nuclei = masks
        print(masks.shape)
        return masks
    
    def segment_image(self, image, save_dir=None):
        """
        Segment the image using Cellpose and optionally save the results.
        
        Args:
            image: Input image to segment
            save_dir: Optional directory to save intermediate results
        
        Returns:
            The segmented image masks
        """
        trained_model_path = '/gpfs/gsfs10/users/sagarm2/cellpose_chr/newDataSet/trainingfiles/models/cellpose_1718127286.8010929'
        model = models.CellposeModel(gpu=True, pretrained_model=trained_model_path)

        masks, flows, styles = model.eval([image], diameter=None, channels=[0, 0])
        self.nuclei = masks[0]  # Use the first element of the list returned by eval
        print(f"Segmentation shape: {self.nuclei.shape}")

        # Save intermediate results if a directory is provided
        if save_dir is not None:
            try:
                intermediate_path = os.path.join(save_dir, "intermediate_results")
                os.makedirs(intermediate_path, exist_ok=True)
                
                # Save segmentation mask
                seg_file = os.path.join(intermediate_path, "segmentation.npy")
                np.save(seg_file, self.nuclei)
                print(f"Saved segmentation to: {seg_file}")
                
            except Exception as e:
                print(f"Error saving segmentation: {e}")

        return self.nuclei
    

    def segment_image_BU_original(self, image):
        trained_model_path = '/gpfs/gsfs10/users/sagarm2/cellpose_chr/newDataSet/trainingfiles/models/cellpose_1718127286.8010929'
        model = models.CellposeModel(gpu=True, pretrained_model=trained_model_path)

        masks, flows, styles = model.eval([image], diameter=None, channels=[0, 0])
        self.nuclei = masks[0]  # Use the first element of the list returned by eval
        print(self.nuclei.shape)
        return self.nuclei

        #results = model.eval([image_to_segment], channels=[0, 0])
        
        # Perform segmentation
        #save_mask_dir = '/gpfs/gsfs10/users/sagarm2/cellpose_chr'
        #mask_save_path = os.path.join(save_mask_dir, 'segmented_mask.tif')
        
        # Assuming results[0][0] contains the mask
        #segmented_mask = results[0][0]
        #self.nuclei = results
        # Save the segmented mask
        #tifffile.imwrite(mask_save_path, segmented_mask)
        
        #print(f"Segmentation completed. Mask saved to {mask_save_path}")
        #return segmented_mask

    def find_peaks(self, dataIn, n=5):
        neighborhood_size = (1,) * (dataIn.ndim - 2) + (n, n)
        data_max = ndi.maximum_filter(dataIn, neighborhood_size)
        data_min = ndi.minimum_filter(dataIn, neighborhood_size)
        peaks = data_max - data_min
        peaks[dataIn != data_max] = 0
        
        mask = np.ones(peaks.shape, dtype=bool)
        mask[..., n:-n, n:-n] = False
        peaks[mask] = 0
        
        return peaks

    def detect_spots1(self, image, channel, threshold=0.4):
        filtered = self.find_peaks(image, 5)
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold
        
        labeled_spots, _ = ndi.label(normalized2)
        
        normalizedIm = normalized2 > 0.25
        selem = morph.disk(4)
        imDilated = morph.dilation(normalizedIm, selem)
        
        if channel == 'DNA-FISH':
            self.normalizedDNAFISH = normalized
            self.spotLabelsDNAFISH = np.unique(self.nuclei[np.array(normalizedIm)])
        elif channel == 'CENPC':
            self.normalizedCENPC = normalized
            self.spotLabelsCENPC = np.unique(self.nuclei[np.array(normalizedIm)])
        
        return imDilated.astype(np.uint8)
    



    
    def detect_spots2(self, image, channel, threshold=0.4):
        print(f"Detecting spots for channel: {channel}")
        filtered = self.find_peaks(image, 5)
        print(f"Filtered image shape: {filtered.shape}")
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold
        
        labeled_spots, _ = ndi.label(normalized2)
        
        normalizedIm = normalized2 > 0.25
        selem = morph.disk(4)
        imDilated = morph.dilation(normalizedIm, selem)
        
        if channel == 'DNA-FISH':
            self.normalizedDNAFISH = normalized
            if self.nuclei is not None:
                self.spotLabelsDNAFISH = np.unique(self.nuclei[np.array(normalizedIm)])
            else:
                print("Error: nuclei is None for DNA-FISH")
        elif channel == 'CENPC':
            self.normalizedCENPC = normalized
            if self.nuclei is not None:
                self.spotLabelsCENPC = np.unique(self.nuclei[np.array(normalizedIm)])
            else:
                print("Error: nuclei is None for CENPC")
        
        return imDilated.astype(np.uint8)
    
    def detect_spots(self, image, channel, threshold=0.4):
        print(f"Detecting spots for channel: {channel}")
        filtered = self.find_peaks(image, 5)
        print(f"Filtered image shape: {filtered.shape}")
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold
        
        labeled_spots, _ = ndi.label(normalized2)
        
        normalizedIm = normalized2 > 0.25
        selem = morph.disk(4)
        imDilated = morph.dilation(normalizedIm, selem)
        
        if channel == 'DNA-FISH':
            self.normalizedDNAFISH = normalized
            self.dna_fish_spots = np.argwhere(imDilated)
            self.spotLabelsDNAFISH = np.unique(self.nuclei[np.array(normalizedIm)])
        elif channel == 'CENPC':
            self.normalizedCENPC = normalized
            self.cenpc_spots = np.argwhere(imDilated)
            self.spotLabelsCENPC = np.unique(self.nuclei[np.array(normalizedIm)])
        
        return imDilated.astype(np.uint8)
    

    from scipy import ndimage as ndi
    def detect_spots_cent(self, image, channel_type, threshold, save_dir=None):
        """
        Detect spots and save results if save_dir is provided.
        
        Args:
            image: Input image
            channel_type: Type of channel ('DNA-FISH' or 'CENPC')
            threshold: Detection threshold
            save_dir: Directory to save results (optional)
            
        Returns:
            numpy.ndarray: Array of centroids if spots are detected, None otherwise
        """
        try:
            if channel_type == 'DNA-FISH':
                labels, centroids = self.detect_spots_no_segmentation(image, threshold, channel='DNA-FISH')
                
                if labels is not None and centroids is not None and len(centroids) > 0:
                    self.labels_dna_fish = labels
                    self.dna_fish_centroids = centroids
                    
                    # Save results
                    if save_dir is not None:
                        try:
                            intermediate_path = os.path.join(save_dir, "intermediate_results")
                            os.makedirs(intermediate_path, exist_ok=True)
                            np.save(os.path.join(intermediate_path, "dna_fish_spots.npy"), self.labels_dna_fish)
                            np.save(os.path.join(intermediate_path, "dna_fish_centroids.npy"), self.dna_fish_centroids)
                            print(f"Saved DNA-FISH spots to: {intermediate_path}")
                        except Exception as save_error:
                            print(f"Error saving DNA-FISH results: {save_error}")
                    
                    print(f"Successfully detected {len(centroids)} DNA-FISH spots")
                    return self.dna_fish_centroids
                else:
                    print("No DNA-FISH spots detected or invalid results")
                    return None
                    
            elif channel_type == 'CENPC':
                labels, centroids = self.detect_spots_no_segmentation(image, threshold, channel='CENPC')
                
                if labels is not None and centroids is not None and len(centroids) > 0:
                    self.labels_cenpc = labels
                    self.cenpc_centroids = centroids
                    
                    # Save results
                    if save_dir is not None:
                        try:
                            intermediate_path = os.path.join(save_dir, "intermediate_results")
                            os.makedirs(intermediate_path, exist_ok=True)
                            np.save(os.path.join(intermediate_path, "cenpc_spots.npy"), self.labels_cenpc)
                            np.save(os.path.join(intermediate_path, "cenpc_centroids.npy"), self.cenpc_centroids)
                            print(f"Saved CENPC spots to: {intermediate_path}")
                        except Exception as save_error:
                            print(f"Error saving CENPC results: {save_error}")
                    
                    print(f"Successfully detected {len(centroids)} CENPC spots")
                    return self.cenpc_centroids
                else:
                    print("No CENPC spots detected or invalid results")
                    return None
            else:
                print(f"Unknown channel type: {channel_type}")
                return None
                        
        except Exception as e:
            print(f"Error in detect_spots_cent: {str(e)}")
            return None

    def detect_spots_cent_BU(self, image, channel, threshold=0.4):
        print(f"Detecting spots for channel: {channel}")
        filtered = self.find_peaks(image, 5)
        print(f"Filtered image shape: {filtered.shape}")
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold

        # Label the spots
        labeled_spots, _ = ndi.label(normalized2)
        
        # Calculate the centroids of each labeled spot
        centroids = ndi.center_of_mass(normalized2, labeled_spots, range(1, labeled_spots.max() + 1))
        centroids = np.array(centroids)  # Convert to numpy array for easier handling

        normalizedIm = normalized2 > 0.25
        selem = morph.disk(4)
        imDilated = morph.dilation(normalizedIm, selem)

        if channel == 'DNA-FISH':
            self.normalizedDNAFISH = normalized
            self.dna_fish_spots = np.argwhere(imDilated)
            self.dna_fish_centroids = centroids  # Store centroids
            self.spotLabelsDNAFISH = np.unique(self.nuclei[np.array(normalizedIm)])
        elif channel == 'CENPC':
            self.normalizedCENPC = normalized
            self.cenpc_spots = np.argwhere(imDilated)
            self.cenpc_centroids = centroids  # Store centroids
            self.spotLabelsCENPC = np.unique(self.nuclei[np.array(normalizedIm)])

        return imDilated.astype(np.uint8)
    

    
    def detect_spots_no_segmentation(self, image, threshold=0.4, channel=None):
        """
        Detect spots in an image without segmentation.
        
        Args:
            image: Input image
            threshold: Detection threshold (default: 0.4)
            channel: Channel type ('DNA-FISH' or 'CENPC')
            
        Returns:
            tuple: (labeled_spots, centroids) or (None, None) if no spots detected
        """
        try:
            print(f"Detecting spots independently with threshold: {threshold}")
            filtered = self.find_peaks(image, 5)
            print(f"Filtered image shape: {filtered.shape}")
            
            # Normalize the filtered image
            normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
            normalized2 = normalized > threshold
            
            # Label the spots
            labeled_spots, num_labels = ndi.label(normalized2)
            print(f"Number of spots detected: {num_labels}")
            
            if num_labels == 0:
                print(f"No spots detected for {channel}")
                return None, None
            
            # Calculate centroids
            centroids = ndi.center_of_mass(normalized2, labeled_spots, range(1, labeled_spots.max() + 1))
            centroids = np.array(centroids)  # Convert to numpy array for easier handling

            if len(centroids) > 0:
                # Store results based on channel type
                if channel == 'DNA-FISH':
                    self.labels_dna_fish = labeled_spots
                    self.dna_fish_centroids = centroids
                    print(f"DNA-FISH spots detected. Count: {len(centroids)}")
                elif channel == 'CENPC':
                    self.labels_cenpc = labeled_spots
                    self.cenpc_centroids = centroids
                    print(f"CENPC spots detected. Count: {len(centroids)}")
                else:
                    print(f"Warning: Unknown channel {channel}")
                    return None, None
                    
                return labeled_spots, centroids
            else:
                print(f"No valid centroids found for {channel}")
                return None, None
                
        except Exception as e:
            print(f"Error in detect_spots_no_segmentation: {str(e)}")
            return None, None


    def detect_spots_no_segmentation_BU(self, image, threshold=0.4):
        print(f"Detecting spots independently with threshold: {threshold}")
        filtered = self.find_peaks(image, 5)
        print(f"Filtered image shape: {filtered.shape}")
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold
        
        labeled_spots, num_labels = ndi.label(normalized2)
        
        # Calculate centroids using center_of_mass for better accuracy
        centroids = []
        for label in range(1, num_labels + 1):
            coords = np.where(labeled_spots == label)
            if len(coords[0]) > 0:
                y, x = ndi.center_of_mass(labeled_spots == label)
                centroids.append([x, y])  # Note: x, y order to match detect_spots
        
        if len(centroids) > 0:
            centroids = np.array(centroids)
            # Store centroids based on which channel is being processed
            if 'DNA-FISH' in str(image.name) if hasattr(image, 'name') else True:
                self.dna_fish_centroids = centroids
            else:
                self.cenpc_centroids = centroids
            print(f"Total spot count: {len(centroids)}")
        
        return normalized2.astype(np.uint8), labeled_spots
    
    def find_common(self, threshold_dna_fish, threshold_cenpc):
        """
        Find DNA-FISH spots in chromosomes that also contain CENPC spots and calculate CENPC intensity.
        
        Args:
            threshold_dna_fish: Threshold used for DNA-FISH detection
            threshold_cenpc: Threshold used for CENPC detection
            
        Returns:
            numpy.ndarray: Labeled image showing chromosomes with both types of spots
        """
        try:
            # Check if we have all required data
            if self.nuclei is None:
                print("No chromosome segmentation found. Please run segmentation first.")
                return None
                
            if self.dna_fish_centroids is None or self.cenpc_centroids is None:
                print("Spots not detected. Please run spot detection first.")
                return None

            # Initialize output image
            labelled_nuclei = np.zeros_like(self.nuclei, dtype=np.uint8)
            
            # Get unique chromosome labels
            chromosome_labels = np.unique(self.nuclei)
            chromosome_labels = chromosome_labels[chromosome_labels != 0]  # Remove background
            
            common_chromosomes = []
            dna_fish_locations = []
            cenpc_intensities = []
            
            # For each chromosome
            for label in chromosome_labels:
                # Create mask for this chromosome
                chromosome_mask = self.nuclei == label
                
                # Find DNA-FISH spots in this chromosome
                dna_fish_in_chromosome = []
                for y, x in self.dna_fish_centroids:
                    if chromosome_mask[int(y), int(x)]:
                        dna_fish_in_chromosome.append((y, x))
                        
                # Find CENPC spots in this chromosome
                cenpc_in_chromosome = []
                for y, x in self.cenpc_centroids:
                    if chromosome_mask[int(y), int(x)]:
                        cenpc_in_chromosome.append((y, x))
                
                # If chromosome has both types of spots
                if len(dna_fish_in_chromosome) > 0 and len(cenpc_in_chromosome) > 0:
                    common_chromosomes.append(label)
                    # Mark this chromosome in output
                    labelled_nuclei[chromosome_mask] = label
                    
                    # For each DNA-FISH spot in this chromosome
                    for dna_y, dna_x in dna_fish_in_chromosome:
                        dna_fish_locations.append((dna_y, dna_x))
                        # Get CENPC intensity at this location
                        if hasattr(self, 'img_cenpc') and self.img_cenpc is not None:
                            intensity = self.img_cenpc[int(dna_y), int(dna_x)]
                            cenpc_intensities.append(intensity)
            
            # Create DataFrames with results
            if len(dna_fish_locations) > 0:
                self.df_centroid_dna_fish = pd.DataFrame(dna_fish_locations, columns=['Y', 'X'])
                if cenpc_intensities:
                    self.df_centroid_dna_fish['CENPC_Intensity'] = cenpc_intensities
                print(f"Found {len(dna_fish_locations)} DNA-FISH spots in {len(common_chromosomes)} chromosomes with CENPC")
            else:
                print("No chromosomes found with both DNA-FISH and CENPC spots")
                self.df_centroid_dna_fish = pd.DataFrame(columns=['Y', 'X', 'CENPC_Intensity'])
                return np.zeros_like(self.nuclei, dtype=np.uint8)
                
            return labelled_nuclei
            
        except Exception as e:
            print(f"Error in find_common2: {str(e)}")
            return None
        
    def find_common_BU(self, threshold_dna_fish, threshold_cenpc):
        common_labels = np.intersect1d(self.spotLabelsDNAFISH, self.spotLabelsCENPC)
        self.df_centroid_cenpc = self.get_spot_location(self.normalizedCENPC, threshold_cenpc, common_labels)
        self.df_centroid_dna_fish = self.get_spot_location(self.normalizedDNAFISH, threshold_dna_fish, common_labels)
        labelled_nuclei = np.zeros_like(self.nuclei, dtype=np.uint8)
        for label in common_labels:
            labelled_nuclei[self.nuclei == label] = label
        return labelled_nuclei

    def find_common2(self, threshold_dna_fish, threshold_cenpc):
        if self.dna_fish_spots is None or self.cenpc_spots is None:
            raise ValueError("DNA-FISH spots or CENPC spots are not detected properly.")
        
        dna_fish_set = set(map(tuple, self.dna_fish_spots))
        cenpc_set = set(map(tuple, self.cenpc_spots))
        common_spots = np.array(list(dna_fish_set & cenpc_set))

        if common_spots.size == 0:
            self.df_centroid_dna_fish = pd.DataFrame(columns=['Y', 'X'])
            self.df_centroid_cenpc = pd.DataFrame(columns=['Y', 'X'])
        else:
            self.df_centroid_dna_fish = pd.DataFrame(common_spots, columns=['Y', 'X'])
            self.df_centroid_cenpc = pd.DataFrame(common_spots, columns=['Y', 'X'])

        labelled_nuclei = np.zeros_like(self.nuclei, dtype=np.uint8)
        for spot in common_spots:
            y, x = spot
            labelled_nuclei[y, x] = 1
        return labelled_nuclei
    # is the threhshold necessary here?


    
        #first argument is the spot image displayed in the viewer

    def get_spot_location(self, normIm, threshold, labels_to_get_centroid):
        spot_mask = normIm > threshold
        indices = self.nuclei[np.array(spot_mask)]
        spotLabels = labels_to_get_centroid
        spot_info = {'Label': [], 'Centroid': []}
        for i in range(1, len(spotLabels)):
            label = spotLabels[i]
            coords = np.argwhere((self.nuclei == label) & (spot_mask))
            if coords.size > 0:
                centroid = np.mean(coords, axis=0).astype(int)
                spot_info['Label'].append(label)
                spot_info['Centroid'].append(centroid.tolist())
        df = pd.DataFrame(spot_info)
        return df

    def gen_intensity_from_df(self, intensity_image, spots_df):
        """
        Generate intensity measurements from an image at specified spot locations.
        """
        try:
            if spots_df is None or spots_df.empty:
                print("No spots provided for intensity measurement")
                return None
                
            if intensity_image is None:
                print("No intensity image provided")
                return None
                
            # Create a copy of the input DataFrame
            result_df = spots_df.copy()
            
            # Ensure we have 'Y' and 'X' columns
            if 'Y' not in result_df.columns or 'X' not in result_df.columns:
                print("DataFrame must contain 'Y' and 'X' columns")
                return None
                
            # Calculate intensity at each spot location
            intensities = []
            for _, row in result_df.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < intensity_image.shape[0] and 0 <= x < intensity_image.shape[1]:
                    intensity = intensity_image[y, x]
                    intensities.append(intensity)
                else:
                    intensities.append(0)
                    
            # Add intensities to the DataFrame
            result_df['CENPC_Intensity'] = intensities
            
            print(f"Measured {len(intensities)} intensity values")
            return result_df
            
        except Exception as e:
            print(f"Error measuring intensities: {str(e)}")
            return None
        
        
    def calculate_intensity_all_dna_fish(self):
        """
        Calculate CENPC intensity at all DNA-FISH locations without segmentation.
        """
        try:
            if self.dna_fish_centroids is None:
                print("No DNA-FISH spots detected")
                return None
                
            if self.img_cenpc is None:
                print("CENPC image not found")
                return None
                
            # Create DataFrame with DNA-FISH centroids
            spots_df = pd.DataFrame(self.dna_fish_centroids, columns=['Y', 'X'])
            
            # Calculate CENPC intensity at each DNA-FISH location
            intensities = []
            for _, row in spots_df.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < self.img_cenpc.shape[0] and 0 <= x < self.img_cenpc.shape[1]:
                    intensity = self.img_cenpc[y, x]
                    intensities.append(intensity)
                else:
                    intensities.append(0)
            
            # Add intensities to DataFrame
            spots_df['CENPC_Intensity'] = intensities
            
            print(f"Calculated intensities for {len(spots_df)} DNA-FISH spots")
            return spots_df
            
        except Exception as e:
            print(f"Error calculating intensities: {str(e)}")
            return None

    def get_centroids_from_labels(self, labels):
        from scipy import ndimage
        
        # Get all unique labels except background (0)
        unique_labels = np.unique(labels)
        unique_labels = unique_labels[unique_labels != 0]
        
        # Calculate centroids for each label
        centroids = []
        for label in unique_labels:
            # Get coordinates of current label
            coords = np.where(labels == label)
            # Calculate centroid
            y, x = ndimage.center_of_mass(labels == label)
            centroids.append([x, y])
        
        return np.array(centroids)
    
    def merge_nuclei_with_line(self, line_coords):
        rr, cc = line_coords.T.astype(int)
        rr = np.clip(rr, 0, self.nuclei.shape[0] - 1)
        cc = np.clip(cc, 0, self.nuclei.shape[1] - 1)
        touched_labels = np.unique(self.nuclei[rr, cc])
        touched_labels = touched_labels[touched_labels > 0]  # Remove background
        if len(touched_labels) > 1:
            new_label = touched_labels[0]
            for label in touched_labels[1:]:
                self.nuclei[self.nuclei == label] = new_label
        return self.nuclei
    

    def remove_nuclei_with_line(self, line_coords):
        rr, cc = line_coords.T.astype(int)
        rr = np.clip(rr, 0, self.nuclei.shape[0] - 1)
        cc = np.clip(cc, 0, self.nuclei.shape[1] - 1)
        touched_labels = np.unique(self.nuclei[rr, cc])
        touched_labels = touched_labels[touched_labels > 0]  # Remove background
        for label in touched_labels:
            self.nuclei[self.nuclei == label] = 0  # Set the label to 0 (background)
        return self.nuclei
    

    def split_chromosome_with_line(self, line_coords):
        rr, cc = line_coords.T.astype(int)
        rr = np.clip(rr, 0, self.nuclei.shape[0] - 1)
        cc = np.clip(cc, 0, self.nuclei.shape[1] - 1)
        touched_labels = np.unique(self.nuclei[rr, cc])
        touched_labels = touched_labels[touched_labels > 0]  # Remove background

        if len(touched_labels) > 0:
            # Create a line mask based on the drawn coordinates
            line_mask = np.zeros_like(self.nuclei, dtype=bool)
            line_mask[rr, cc] = True

            # Create a copy of the original nuclei array to modify
            updated_nuclei = self.nuclei.copy()

            for label in touched_labels:
                mask = self.nuclei == label
                # Split the chromosome by labeling the connected components on each side of the line
                from scipy.ndimage import label
                left_labels, left_num_features = label(~line_mask & mask)
                right_labels, right_num_features = label(line_mask & mask)

                # Update the original nuclei array with the split components
                left_labels[left_labels > 0] += updated_nuclei.max()  # Ensure unique labels
                updated_nuclei[left_labels > 0] = left_labels[left_labels > 0]

                right_labels[right_labels > 0] += updated_nuclei.max()  # Ensure unique labels
                updated_nuclei[right_labels > 0] = right_labels[right_labels > 0]

            self.nuclei = updated_nuclei

        return self.nuclei

    def split_chromosome_with_line_BU(self, line_coords):
        rr, cc = line_coords.T.astype(int)
        rr = np.clip(rr, 0, self.nuclei.shape[0] - 1)
        cc = np.clip(cc, 0, self.nuclei.shape[1] - 1)
        touched_labels = np.unique(self.nuclei[rr, cc])
        touched_labels = touched_labels[touched_labels > 0]  # Remove background

        if len(touched_labels) == 1:
            label_to_split = touched_labels[0]
            mask = self.nuclei == label_to_split

            # Create a line mask based on the drawn coordinates
            line_mask = np.zeros_like(self.nuclei, dtype=bool)
            line_mask[rr, cc] = True

            # Remove pixels in the line from the original chromosome mask
            mask[line_mask] = False

            # Label the connected components in the mask, which effectively splits the chromosome
            from scipy.ndimage import label
            split_labels, num_features = label(mask)
            
            # Update the original nuclei array with split components
            split_labels[split_labels > 0] += self.nuclei.max()  # Ensure unique labels
            self.nuclei[split_labels > 0] = split_labels[split_labels > 0]

        return self.nuclei
