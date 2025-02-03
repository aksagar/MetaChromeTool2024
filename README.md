# Metaphase Chromosome Analysis with Napari

This repository contains code for analyzing metaphase chromosomes using the [Napari](https://napari.org/stable/) platform. The code is tailored for visualization, segmentation, and quantitative analysis of chromosome images, with integrated tools for centromere detection and CENP-A level measurement.



## Description
This project uses Napari for analyzing metaphase chromosome images. It enables users to:

- Visualize and segment chromosome structures.
- Detect centromeres and measure CENP-A levels within metaphase chromosome regions.
- Perform quantitative analysis of chromosome features.

## Installation Instructions

To set up the environment and run Napari for chromosome analysis, follow these steps:

### Step 1: Install Anaconda

Before beginning, ensure that [Anaconda](https://www.anaconda.com/products/individual) is installed on your machine. You can download it from the official Anaconda website.

### Step 2: Create and Activate a Conda Environment

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

6. Navigate to the folder where you downloaded the code and launch the program (To verify you are in the right folder make sure this directory contains the main.py file along with the other files you downloaded)
    ```bash
    python main.py
    ```
    

For more detailed installation instructions, please refer to the [Napari Installation Guide](https://napari.org/stable/tutorials/fundamentals/installation.html).

## Usage

Once Napari is installed and running, you can load your metaphase chromosome images and use the provided tools for chromosome segmentation, centromere detection, and CENP-A level analysis. For further assistance, please refer to the documentation within the code or contact the author.

---

Feel free to contribute by submitting issues or pull requests. Your feedback is valuable for improving the analysis tools.
