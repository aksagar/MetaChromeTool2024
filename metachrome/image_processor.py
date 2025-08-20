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


from skimage.draw import line
from scipy.ndimage import binary_dilation


class ImageProcessor:
    def __init__(self):
        self.img_cenpc = None
        self.img_dna_fish = None
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
            if dna_fish_id in image_file:  # Add this check
                self.img_dna_fish = image
        
        return images

    def get_spots_in_common_regions(self, df_spots, common_nuclei):
        """
        Filter spots to only include those in common regions.
        
        Parameters:
        -----------
        df_spots : pandas.DataFrame
            DataFrame containing spot coordinates with columns ['Y', 'X']
        common_nuclei : numpy.ndarray
            Label image of common regions
            
        Returns:
        --------
        pandas.DataFrame
            Filtered DataFrame containing only spots in common regions
        """
        try:
            if df_spots is None or df_spots.empty:
                return None
                
            # Create mask for spots in common regions
            spots_in_common = []
            for _, spot in df_spots.iterrows():
                y, x = int(spot['Y']), int(spot['X'])
                if y < common_nuclei.shape[0] and x < common_nuclei.shape[1]:
                    if common_nuclei[y, x] > 0:  # Check if spot is in a common region
                        spots_in_common.append(True)
                    else:
                        spots_in_common.append(False)
                else:
                    spots_in_common.append(False)
                    
            # Filter spots
            return df_spots[spots_in_common].copy()
            
        except Exception as e:
            print(f"Error filtering spots in common regions: {str(e)}")
            return None

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
        trained_model_path = 'cellpose_1718127286.8010929'
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
        
        if channel == 'Channel 1':
            self.normalizedDNAFISH = normalized
            self.dna_fish_spots = np.argwhere(imDilated)
            self.spotLabelsDNAFISH = np.unique(self.nuclei[np.array(normalizedIm)])
        elif channel == 'Channel 2':    

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
            if channel_type == 'Channel 1':
                labels, centroids = self.detect_spots_no_segmentation(image, threshold, channel='Channel 1')
                

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
                            print(f"Saved Channel 1 spots to: {intermediate_path}")
                        except Exception as save_error:
                            print(f"Error saving Channel 1 results: {save_error}")
                    

                    print(f"Successfully detected {len(centroids)} Channel 1 spots")
                    return self.dna_fish_centroids
                else:
                    print("No Channel 1 spots detected or invalid results")
                    return None
                    

            elif channel_type == 'Channel 2':
                labels, centroids = self.detect_spots_no_segmentation(image, threshold, channel='Channel 2')
                

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
                            print(f"Saved Channel 2 spots to: {intermediate_path}")
                        except Exception as save_error:
                            print(f"Error saving Channel 2 results: {save_error}")
                    

                    print(f"Successfully detected {len(centroids)} Channel 2 spots")
                    return self.cenpc_centroids
                else:

                    print("No Channel 2 spots detected or invalid results")
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
                if channel == 'Channel 1':
                    self.labels_dna_fish = labeled_spots
                    self.dna_fish_centroids = centroids
                    print(f"Channel 1 spots detected. Count: {len(centroids)}")

                elif channel == 'Channel 2':
                    self.labels_cenpc = labeled_spots
                    self.cenpc_centroids = centroids
                    print(f"Channel 2 spots detected. Count: {len(centroids)}")

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
    
    def find_common(self):
        """
        Find Channel 1 spots in chromosomes that also contain Channel 2 spots and calculate Channel 2 intensity.
        
        Returns:
            numpy.ndarray: Labeled image showing chromosomes with both types of spots
            None: If no spots are found or an error occurs
            
        Note:
            - Requires prior segmentation (self.nuclei must be set)
            - Requires prior spot detection for both channels
            - Updates self.common_nuclei with the labeled image
            - Updates self.df_centroid_dna_fish with spot locations and intensities
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
            
            common_nuclei = []
            channel1_locations = []
            channel2_intensities = []
            
            # For each chromosome
            for label in chromosome_labels:
                # Create mask for this chromosome
                chromosome_mask = self.nuclei == label
                
                # Find Channel 1 spots in this chromosome
                channel1_in_chromosome = []
                for y, x in self.dna_fish_centroids:
                    if chromosome_mask[int(y), int(x)]:
                        channel1_in_chromosome.append((y, x))
                        
                # Find Channel 2 spots in this chromosome
                channel2_in_chromosome = []
                for y, x in self.cenpc_centroids:
                    if chromosome_mask[int(y), int(x)]:
                        channel2_in_chromosome.append((y, x))
                
                # If chromosome has both types of spots
                if len(channel1_in_chromosome) > 0 and len(channel2_in_chromosome) > 0:
                    common_nuclei.append(label)
                    # Mark this chromosome in output
                    labelled_nuclei[chromosome_mask] = label
                    
                    # For each Channel 1 spot in this chromosome
                    for ch1_y, ch1_x in channel1_in_chromosome:
                        channel1_locations.append((ch1_y, ch1_x))
                        # Get Channel 2 intensity at this location
                        if hasattr(self, 'img_cenpc') and self.img_cenpc is not None:
                            intensity = self.img_cenpc[int(ch1_y), int(ch1_x)]
                            channel2_intensities.append(intensity)
            
            # Store the common nuclei
            self.common_nuclei = labelled_nuclei
            
            # Create DataFrames with results
            if len(channel1_locations) > 0:
                self.df_centroid_dna_fish = pd.DataFrame(channel1_locations, columns=['Y', 'X'])
                if channel2_intensities:
                    self.df_centroid_dna_fish['Channel2_Intensity'] = channel2_intensities
                print(f"Found {len(channel1_locations)} Channel 1 spots in {len(common_nuclei)} chromosomes with Channel 2")
            else:
                print("No chromosomes found with both Channel 1 and Channel 2 spots")
                self.df_centroid_dna_fish = pd.DataFrame(columns=['Y', 'X', 'Channel2_Intensity'])
                return np.zeros_like(self.nuclei, dtype=np.uint8)
                
            return labelled_nuclei
            
        except Exception as e:
            print(f"Error in find_common: {str(e)}")
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
            raise ValueError("Channel 1 spots or Channel 2 spots are not detected properly.")
        

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
    

    def measure_intensity_at_spots(self, intensity_image, spots_df, channel_name):
        """
        Measure intensity values at spot locations for a single channel.
        
        Args:
            intensity_image: Image to measure intensities from
            spots_df: DataFrame with spot locations (must have 'Y' and 'X' columns)
            channel_name: Name of the channel ('Channel1' or 'Channel2')
            
        Returns:
            DataFrame with spot locations and their corresponding intensities
        """
        try:
            if spots_df is None or spots_df.empty:
                print(f"No {channel_name} spots provided for intensity measurement")
                return None
                
            if intensity_image is None:
                print(f"No intensity image provided")
                return None
                
            # Calculate intensities at spot locations
            intensities = []
            for _, row in spots_df.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < intensity_image.shape[0] and 0 <= x < intensity_image.shape[1]:
                    # Calculate mean intensity in 5x5 region
                    y_min, y_max = max(0, y-2), min(intensity_image.shape[0], y+3)
                    x_min, x_max = max(0, x-2), min(intensity_image.shape[1], x+3)
                    intensity = np.mean(intensity_image[y_min:y_max, x_min:x_max])
                    intensities.append(intensity)
                else:
                    intensities.append(0)
                    
            # Create DataFrame with spots and intensities
            result_df = pd.DataFrame({
                'Y': spots_df['Y'],
                'X': spots_df['X'],
                'Intensity': intensities
            })
            
            print(f"Measured {len(intensities)} intensities at {channel_name} spots")
            return result_df
            
        except Exception as e:
            print(f"Error measuring intensities at {channel_name} spots: {str(e)}")
            return None
        

    def gen_intensity_from_df(self, intensity_image_ch2, spots_df_ch1, intensity_image_ch1=None, spots_df_ch2=None):
        """
        Generate intensity measurements for both channels at each other's spot locations.
        
        Args:
            intensity_image_ch2: Channel 2 intensity image
            spots_df_ch1: DataFrame with Channel 1 spot locations
            intensity_image_ch1: Channel 1 intensity image (optional)
            spots_df_ch2: DataFrame with Channel 2 spot locations (optional)
            
        Returns:
            DataFrame with both channels' spot locations and their corresponding intensities
        """
        try:
            if spots_df_ch1 is None or spots_df_ch1.empty:
                print("No Channel 1 spots provided for intensity measurement")
                return None
                
            if intensity_image_ch2 is None:
                print("No Channel 2 intensity image provided")
                return None
                
            # Calculate Channel 2 intensity at Channel 1 spot locations
            ch2_intensities = []
            for _, row in spots_df_ch1.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < intensity_image_ch2.shape[0] and 0 <= x < intensity_image_ch2.shape[1]:
                    # Calculate mean intensity in 5x5 region
                    y_min, y_max = max(0, y-2), min(intensity_image_ch2.shape[0], y+3)
                    x_min, x_max = max(0, x-2), min(intensity_image_ch2.shape[1], x+3)
                    intensity = np.mean(intensity_image_ch2[y_min:y_max, x_min:x_max])
                    ch2_intensities.append(intensity)
                else:
                    ch2_intensities.append(0)
                    
            # Create initial DataFrame with Channel 1 spots and Channel 2 intensities
            result_df = pd.DataFrame({
                'Channel1_Spot_Y': spots_df_ch1['Y'],
                'Channel1_Spot_X': spots_df_ch1['X'],
                'Channel2_Intensity': ch2_intensities
            })
            
            # If Channel 2 spots and Channel 1 intensity image are provided
            if spots_df_ch2 is not None and intensity_image_ch1 is not None and not spots_df_ch2.empty:
                # Calculate Channel 1 intensity at Channel 2 spot locations
                ch1_intensities = []
                ch2_spot_y = []
                ch2_spot_x = []
                
                for _, row in spots_df_ch2.iterrows():
                    y, x = int(row['Y']), int(row['X'])
                    ch2_spot_y.append(y)
                    ch2_spot_x.append(x)
                    
                    if 0 <= y < intensity_image_ch1.shape[0] and 0 <= x < intensity_image_ch1.shape[1]:
                        # Calculate mean intensity in 5x5 region
                        y_min, y_max = max(0, y-2), min(intensity_image_ch1.shape[0], y+3)
                        x_min, x_max = max(0, x-2), min(intensity_image_ch1.shape[1], x+3)
                        intensity = np.mean(intensity_image_ch1[y_min:y_max, x_min:x_max])
                        ch1_intensities.append(intensity)
                    else:
                        ch1_intensities.append(0)
                
                # Add Channel 2 spots and Channel 1 intensities
                result_df['Channel2_Spot_Y'] = pd.Series(ch2_spot_y)
                result_df['Channel2_Spot_X'] = pd.Series(ch2_spot_x)
                result_df['Channel1_Intensity'] = pd.Series(ch1_intensities)
                
                # Copy any additional columns from the original DataFrames
                for col in spots_df_ch1.columns:
                    if col not in ['Y', 'X']:
                        result_df[f'Channel1_{col}'] = spots_df_ch1[col]
                
                for col in spots_df_ch2.columns:
                    if col not in ['Y', 'X']:
                        result_df[f'Channel2_{col}'] = spots_df_ch2[col]
                
                print(f"Measured {len(ch2_intensities)} Channel 2 intensities and {len(ch1_intensities)} Channel 1 intensities")
            else:
                print(f"Measured {len(ch2_intensities)} Channel 2 intensity values")
            
            return result_df
            
        except Exception as e:
            print(f"Error measuring intensities: {str(e)}")
            return None       
    def calculate_intensity_all_cenpc(self, intensity_image_ch1=None):
        """
        Calculate Channel 1 intensity at all Channel 2 spot locations without segmentation.
        
        Args:
            intensity_image_ch1: Optional Channel 1 intensity image. If not provided, uses self.img_dna_fish
            
        Returns:
            pandas.DataFrame: DataFrame containing spot locations and corresponding intensities,
                            or None if spots or intensity image are missing
        """
        try:
            # Check for Channel 2 spots
            if self.cenpc_centroids is None or len(self.cenpc_centroids) == 0:
                print("No Channel 2 spots detected")
                return None
                
            # Use provided intensity image or fall back to stored image
            if intensity_image_ch1 is None:
                intensity_image_ch1 = self.img_dna_fish
                
            if intensity_image_ch1 is None:
                print("Channel 1 intensity image not found")
                return None
                
            # Create DataFrame with Channel 2 centroids
            spots_df = pd.DataFrame(self.cenpc_centroids, columns=['Y', 'X'])
            
            # Calculate Channel 1 intensity at each Channel 2 location
            intensities = []
            for _, row in spots_df.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < intensity_image_ch1.shape[0] and 0 <= x < intensity_image_ch1.shape[1]:
                    # Calculate mean intensity in 5x5 region
                    y_min, y_max = max(0, y-2), min(intensity_image_ch1.shape[0], y+3)
                    x_min, x_max = max(0, x-2), min(intensity_image_ch1.shape[1], x+3)
                    intensity = np.mean(intensity_image_ch1[y_min:y_max, x_min:x_max])
                    intensities.append(intensity)
                else:
                    intensities.append(0)
            
            # Add intensities to DataFrame
            spots_df['Channel1_Intensity'] = intensities
            
            print(f"Calculated intensities for {len(spots_df)} Channel 2 spots")
            return spots_df
            
        except Exception as e:
            print(f"Error calculating intensities: {str(e)}")
            return None
        
    def calculate_intensity_all_dna_fish(self, intensity_image_ch2=None):
        """
        Calculate Channel 2 intensity at all Channel 1 spot locations without segmentation.
        
        Args:
            intensity_image_ch2: Optional Channel 2 intensity image. If not provided, uses self.img_cenpc
            
        Returns:
            pandas.DataFrame: DataFrame containing spot locations and corresponding intensities,
                            or None if spots or intensity image are missing
        """
        try:
            # Check for Channel 1 spots
            if self.dna_fish_centroids is None or len(self.dna_fish_centroids) == 0:
                print("No Channel 1 spots detected")
                return None
                
            # Use provided intensity image or fall back to stored image
            if intensity_image_ch2 is None:
                intensity_image_ch2 = self.img_cenpc
                
            if intensity_image_ch2 is None:
                print("Channel 2 intensity image not found")
                return None
                
            # Create DataFrame with Channel 1 centroids
            spots_df = pd.DataFrame(self.dna_fish_centroids, columns=['Y', 'X'])
            
            # Calculate Channel 2 intensity at each Channel 1 location
            intensities = []
            for _, row in spots_df.iterrows():
                y, x = int(row['Y']), int(row['X'])
                if 0 <= y < intensity_image_ch2.shape[0] and 0 <= x < intensity_image_ch2.shape[1]:
                    # Calculate mean intensity in 5x5 region
                    y_min, y_max = max(0, y-2), min(intensity_image_ch2.shape[0], y+3)
                    x_min, x_max = max(0, x-2), min(intensity_image_ch2.shape[1], x+3)
                    intensity = np.mean(intensity_image_ch2[y_min:y_max, x_min:x_max])
                    intensities.append(intensity)
                else:
                    intensities.append(0)
            
            # Add intensities to DataFrame
            spots_df['Channel2_Intensity'] = intensities
            
            print(f"Calculated intensities for {len(spots_df)} Channel 1 spots")
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
    


    def delete_dna_fish_spots_with_line(self, viewer):
        """Delete DNA-FISH spots that intersect with the drawn line or points."""
        if self.dna_fish_centroids is None or len(self.dna_fish_centroids) == 0:
            return

        shapes_layer = viewer.layers['Shapes']
        
        # Get image shape from the first image if available, otherwise use default
        if hasattr(self, 'images') and self.images is not None and len(self.images) > 0:
            img_shape = self.images[0].shape
        else:
            img_shape = (1024, 1024)  # default shape
        
        # Create mask for all shapes
        line_mask = np.zeros(img_shape, dtype=bool)
        
        # Process all shapes (lines or points)
        for shape_coords in shapes_layer.data:
            if len(shape_coords) == 1:  # Single point
                y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                if 0 <= y < img_shape[0] and 0 <= x < img_shape[1]:
                    line_mask[y, x] = True
            else:  # Line
                for i in range(len(shape_coords) - 1):
                    start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                    end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                    
                    rr, cc = line(start_y, start_x, end_y, end_x)
                    valid_points = (rr >= 0) & (rr < img_shape[0]) & (cc >= 0) & (cc < img_shape[1])
                    rr, cc = rr[valid_points], cc[valid_points]
                    line_mask[rr, cc] = True
        
        # Buffer the mask
        line_mask = binary_dilation(line_mask, iterations=3)
        
        # Find spots that don't intersect with the mask
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
        for layer in viewer.layers:
            if 'Channel 1 Spots' in layer.name:
                viewer.layers.remove(layer)
        

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

    def delete_cenpc_spots_with_line(self, viewer):
        """Delete Channel 2 spots that intersect with the drawn line or points."""
        if self.cenpc_centroids is None or len(self.cenpc_centroids) == 0:
            return

        shapes_layer = viewer.layers['Shapes']
        
        # Get image shape from the first image if available, otherwise use default
        if hasattr(self, 'images') and self.images is not None and len(self.images) > 0:
            img_shape = self.images[0].shape
        else:
            img_shape = (1024, 1024)  # default shape
        
        # Create mask for all shapes
        line_mask = np.zeros(img_shape, dtype=bool)
        
        # Process all shapes (lines or points)
        for shape_coords in shapes_layer.data:
            if len(shape_coords) == 1:  # Single point
                y, x = int(shape_coords[0][0]), int(shape_coords[0][1])
                if 0 <= y < img_shape[0] and 0 <= x < img_shape[1]:
                    line_mask[y, x] = True
            else:  # Line
                for i in range(len(shape_coords) - 1):
                    start_y, start_x = int(shape_coords[i][0]), int(shape_coords[i][1])
                    end_y, end_x = int(shape_coords[i+1][0]), int(shape_coords[i+1][1])
                    
                    rr, cc = line(start_y, start_x, end_y, end_x)
                    valid_points = (rr >= 0) & (rr < img_shape[0]) & (cc >= 0) & (cc < img_shape[1])
                    rr, cc = rr[valid_points], cc[valid_points]
                    line_mask[rr, cc] = True
        
        # Buffer the mask
        line_mask = binary_dilation(line_mask, iterations=3)
        
        # Find spots that don't intersect with the mask
        kept_spots = []
        square_size = 5  # Half size of the square around each spot
        
        for spot in self.cenpc_centroids:
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
        self.cenpc_centroids = np.array(kept_spots) if kept_spots else np.array([])
        
        # Update the viewer
        for layer in viewer.layers:
            if 'Channel 2 Spots' in layer.name:
                viewer.layers.remove(layer)
        
        if len(self.cenpc_centroids) > 0:
            squares = [
                [[x - 5, y - 5], [x + 5, y - 5], [x + 5, y + 5], [x - 5, y + 5]]
                for x, y in self.cenpc_centroids
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
