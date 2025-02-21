trained_model_path = '/home/sagarm2/cellpsoetrain/allfiles/models/cellpose_1715863264.6910446' # this one is functional

#cellpose_1717698700.292517 

# Load the trained model
model = models.CellposeModel(gpu=True, pretrained_model=trained_model_path)

# Load new images for segmentation
new_image_dir = '/home/sagarm2/cellpsoetrain/allfiles2/'
save_mask_dir = '/home/sagarm2/cellpsoetrain/predmask/'

new_image_files = sorted([f for f in os.listdir(new_image_dir) if f.endswith('.tif')])
new_images = [imread(os.path.join(new_image_dir, f)) for f in new_image_files]

results = model.eval(new_images, channels=[0, 0])