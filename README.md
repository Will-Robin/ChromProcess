# ChromProcess

Python tools to deal with sets of chromatography data.

ChromProcess aims to provide a set of tools to aid in developing data processing pipelines for efficient and reproducible data analysis.

## Installation

### 1. Clone

Clone the repository from GitHub. If you're having trouble with cloning, just download the zip file and store it on your computer wherever you would like.

### 2. Create a virtual environment

Anaconda and Miniconda (conda) and pip are package managers for Python. I use Miniconda, but if you're in doubt what to choose, download [Anaconda](https://www.anaconda.com/products/individual-b#Downloads, 'Anaconda') for loads of science programming tools in one place.

#### Using conda:

Create a virtual environment with conda:

`conda create --name chromprocess-env`

Activate the virtual environment:

`conda activate chromprocess-env`

**Go to Install dependencies.**

#### Using pip:

Create a virtual environment:

`pip install virtualenv`

`virtualenv chromprocess-env`

Activate the virtual environment:

Mac:

`source chromprocess-env/bin/activate`

Windows:

`chromprocess-env\Scripts\activate`

**Go to Install dependencies.**

### Install dependencies

ChromProcess uses a few libraries of code which need to be installed in the virtual environment for it to work. These dependencies can be download using the conda or pip using commands similar to `conda install` or `pip install`. If you are using Anaconda/Miniconda, try to use `conda install etc.` if possible. If you're using pip, no conda for you...

The dependencies are given in environment.yml. For simplicity, you can install them by running the following (remember, the environment must be activated, see above):

conda:
- `conda install numpy`
- `conda install matplotlib`
- `conda install -c conda-forge netCDF4`
- `conda install -c anaconda scipy `

pip:
- `pip install numpy scipy matplotlib`
- `pip install netCDF4`

### Install ChromProcess

In command line/terminal, navigate to the folder (hint: use cd path/to/folder) containing the ChromProcess code, then type:

conda:
  - `conda develop ChromProcess`

  (you may have to run `conda install conda-build` first)

pip:
  - Install from the repository root using pip: `pip install .`,
  - Or in editable mode (so edits are immediately reflected): `pip install -e .`

## Acknowledgement

ChromProcess was created by William E. Robinson during free time and working in the group of Prof. Wilhelm T. S. Huck (Radboud University Nijmegen) supported by NWO, The Simons Foundation and the ERC.

