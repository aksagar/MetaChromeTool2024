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
    
    def segment_image(self, image):
        trained_model_path = 'cellpose_1718127286.8010929'
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
    

    def detect_spots_no_segmentation(self, image, threshold=0.4):
        print(f"Detecting spots independently with threshold: {threshold}")
        filtered = self.find_peaks(image, 5)
        print(f"Filtered image shape: {filtered.shape}")
        normalized = (filtered - filtered.min()) / (filtered.max() - filtered.min())
        normalized2 = normalized > threshold
        
        labeled_spots, labels = ndi.label(normalized2)
        
        normalizedIm = normalized2 > 0.25
        selem = morph.disk(4)
        imDilated = morph.dilation(normalizedIm, selem)
        
        self.labels_dna_fish = labeled_spots
        print(labeled_spots)
        

        return imDilated.astype(np.uint8), labels
    
    def find_common(self, threshold_dna_fish, threshold_cenpc):
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

    def gen_intensity_from_df(self, in_image, df_to_look):
        mean_intensities = []
        
        for _, row in df_to_look.iterrows():
            y, x = row['Centroid']
            y1, y2 = max(y-2, 0), min(y+3, in_image.shape[0])
            x1, x2 = max(x-2, 0), min(x+3, in_image.shape[1])
            square = in_image[y1:y2, x1:x2]
            mean_intensity = np.mean(square)
            mean_intensities.append(mean_intensity)
        
        df_to_look['MeanIntensity'] = mean_intensities
        return df_to_look
    
    def calculate_intensity_all_dna_fish(self):
        if self.labels_dna_fish is None:
            raise ValueError("DNA-FISH spots are not detected properly.")

        # Calculate intensity of CENPC at all DNA-FISH locations
        spot_info = {'Centroid': [], 'MeanIntensity': []}

        labels_unique = np.unique(self.labels_dna_fish)
        print("labels here unique", labels_unique)

        print("here")
        for label in labels_unique:
            if label == 0:
                continue
            coords = np.argwhere(self.labels_dna_fish == label)
            centroid = np.mean(coords, axis=0).astype(int)

            y, x = centroid
            y1, y2 = max(y-2, 0), min(y+3, self.img_cenpc.shape[0])
            x1, x2 = max(x-2, 0), min(x+3, self.img_cenpc.shape[1])
            square = self.img_cenpc[y1:y2, x1:x2]
            mean_intensity = np.mean(square)

            spot_info['Centroid'].append(centroid.tolist())
            spot_info['MeanIntensity'].append(mean_intensity)

        df = pd.DataFrame(spot_info)
        return df

    
    
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