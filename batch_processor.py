import os
import pandas as pd
from image_processor import ImageProcessor

class BatchProcessor:
    def __init__(self, processor: ImageProcessor, control_widget_dna_fish, control_widget_cenpc):
        self.processor = processor
        self.control_widget_dna_fish = control_widget_dna_fish
        self.control_widget_cenpc = control_widget_cenpc

    def batch_processing(self, root_folder, dapi_id, dna_fish_id, cenpc_id, skip_segmentation):
        combined_df = pd.DataFrame()  # Initialize an empty DataFrame to store combined data

        for subdir, _, _ in os.walk(root_folder):
            images = self.processor.load_images(subdir, dapi_id, dna_fish_id, cenpc_id, skip_segmentation)
            if images is None:
                continue  # Skip processing if the required images are not loaded
            if all(img is not None for img in images):
                try:
                    threshold_dna_fish = self.control_widget_dna_fish.slider.value() / 100  # Convert slider value to 0-1 range
                    threshold_cenpc = self.control_widget_cenpc.slider.value() / 100  # Convert slider value to 0-1 range

                    if skip_segmentation:
                        spots_dna_fish, labels_dna_fish = self.processor.detect_spots_no_segmentation(images[0], threshold_dna_fish)
                        spots_cenpc, labels_cenpc = self.processor.detect_spots_no_segmentation(images[1], threshold_cenpc)
                        
                        df_with_cenpc_inten = self.processor.calculate_intensity_all_dna_fish()
                    else:
                        masks = self.processor.segment_image(images[0])
                        spots_dna_fish = self.processor.detect_spots(images[1], 'DNA-FISH', threshold_dna_fish)
                        spots_cenpc = self.processor.detect_spots(images[2], 'CENPC', threshold_cenpc)

                        common_nuclei = self.processor.find_common(threshold_dna_fish, threshold_cenpc)
                        if common_nuclei is None:
                            print(f"No common labels found in {subdir}")
                            continue

                        df_with_cenpc_inten = self.processor.gen_intensity_from_df(self.processor.img_cenpc, self.processor.df_centroid_dna_fish)

                    # Save the individual DataFrame in the current folder using the folder's name
                    folder_name = os.path.basename(subdir)
                    save_path = os.path.join(subdir, f"{folder_name}_intensity.csv")
                    df_with_cenpc_inten.to_csv(save_path, index=False)
                    print(f"DataFrame saved to {save_path}")

                    # Append the individual DataFrame to the combined DataFrame
                    df_with_cenpc_inten['Folder'] = folder_name  # Add a column for folder name
                    combined_df = pd.concat([combined_df, df_with_cenpc_inten], ignore_index=True)

                except Exception as e:
                    print(f"Error processing folder {subdir}: {e}")

        # Save the combined DataFrame in the root folder
        combined_save_path = os.path.join(root_folder, "combined_intensity.csv")
        combined_df.to_csv(combined_save_path, index=False)
        print(f"Combined DataFrame saved to {combined_save_path}")
