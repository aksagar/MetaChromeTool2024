# Metaphase Chromosome Analysis with Napari

This repository contains code for analyzing metaphase chromosomes using the [Napari](https://napari.org/stable/) platform. The code is tailored for visualization, segmentation, and quantitative analysis of chromosome images, with integrated tools for centromere detection and CENP-A level measurement.



## Description
This project uses Napari for analyzing metaphase chromosome images. It enables users to:

- Visualize and segment chromosome structures.
- Detect centromeres and measure CENP-A levels within metaphase chromosome regions.
- Perform quantitative analysis of chromosome features.

## Required Model Files

This repository includes pre-trained Cellpose models for chromosome segmentation. These large model files are managed using Git LFS:

- **cellpose_1718127286.8010929** (25 MB) - Required for Cellpose 2.0 segmentation
- **cellpose_1758654062.3588712** (581 MB) - Required for Cellpose SAM segmentation of fluorescence images of metaphase chromosomes
- **cpsam_GTG_nov7** (581 MB) - Required for G-band image segmentation

**Note:** If you clone this repository, Git LFS will automatically download these model files. If you installed via pip, the models will be included in the package.

## Installation Instructions

### Option 1: Install via pip (Recommended)

The easiest way to install MetaChrome is using pip:

```bash
pip install metachrome
```

After installation, you can run the application with:
```bash
metachrome
```

### Option 2: Install from source

To install from the source code, follow these steps:

### Step 1: Download the Repository

First, download all the files from the GitHub directory and navigate to the downloaded folder. You can do this by using the following commands:

#### Quick Installation (Recommended)

If you have Python installed, you can use the provided installation script:

```bash
python install.py
```

This will automatically install all dependencies and the package in development mode.

#### Manual Installation

#### Windows:
```powershell
# Navigate to the directory where you want to download the repository
cd path\to\your\directory

# Clone the repository
git clone <repository_url>

# Navigate to the downloaded directory
cd <repository_name>
```

#### Linux:
```bash
# Navigate to the directory where you want to download the repository
cd /path/to/your/directory

# Clone the repository
git clone <repository_url>

# Navigate to the downloaded directory
cd <repository_name>
```

Alternatively, if you downloaded the repository as a ZIP file, extract it and navigate to the extracted folder.

### Step 2: Install Anaconda

Before beginning, ensure that [Anaconda](https://www.anaconda.com/products/individual) is installed on your machine. You can download it from the official Anaconda website.

### Step 3: Create and Activate a Conda Environment

After launching Anaconda, follow the steps below to create a new environment and install the necessary packages for Napari.

1. Create a new environment with Python 3.10:
    ```bash
    conda create -y -n napari-env -c conda-forge python=3.10
    ```

2. Activate the newly created environment:
    ```bash
    conda activate napari-env
    ```

3. Install Napari and the required dependencies (PyQt):
    ```bash
    conda install -c conda-forge napari pyqt
    ```

4. Install additional libraries required for analysis:
    ```bash
    conda install -c conda-forge magicgui numpy scikit-image scipy matplotlib pandas qtpy
    ```

5. Install Cellpose using pip:
    ```bash
    pip install cellpose
    ```

### Step 4: Navigate to the Downloaded Folder and Run the Program

After installing the necessary dependencies, navigate to the folder where you downloaded the code and launch the program. To verify that you are in the correct folder, ensure that the directory contains the `main.py` file along with other required files.

#### Windows:
```powershell
cd path\to\downloaded\folder
python main.py
```

#### Linux:
```bash
cd /path/to/downloaded/folder
python main.py
```


    

For more detailed installation instructions, please refer to the [Napari Installation Guide](https://napari.org/stable/tutorials/fundamentals/installation.html).

## Usage

Once Napari is installed and running, you can load your metaphase chromosome images and use the provided tools for chromosome segmentation, centromere detection, and CENP-A level analysis. For further assistance, please refer to the documentation within the code or contact the author.

---

Feel free to contribute by submitting issues or pull requests. Your feedback is valuable for improving the analysis tools.
