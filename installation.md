# Installation

## 1. Clone

`git clone https://github.com/Will-Robin/ChromProcess.git`

## 2. Create a virtual environment

### Using conda:

Create a virtual environment with conda:

`conda create --name chromprocess-env`

Activate the virtual environment:

`conda activate chromprocess-env`

**Go to Install dependencies.**

### Using pip:

Create a virtual environment:

`pip install virtualenv`

`virtualenv chromprocess-env`

Activate the virtual environment:

Mac:

`source chromprocess-env/bin/activate`

Windows:

`chromprocess-env\Scripts\activate`

**Go to Install dependencies.**

## Install dependencies

ChromProcess uses a few libraries of code which need to be installed in the
virtual environment for it to work. These dependencies can be download using the
conda or pip using commands similar to `conda install` or `pip install`. If you
are using Anaconda/Miniconda, try to use `conda install etc.` if possible.

The dependencies are given in environment.yml. For simplicity, you can install
them by running the following (activate the virtual environment before running
these commands):

conda:

- `conda install matplotlib`
- `conda install -c anaconda scipy `
- `conda install -c conda-forge netCDF4`

pip:

- `pip install scipy`
- `pip install numpy`
- `pip install netCDF4`
- `pip install matplotlib`

## Install ChromProcess

In command line/terminal, navigate to the folder (hint: use `cd path/to/folder`)
containing the ChromProcess code, then type:

conda:
  - `conda develop ChromProcess`

(you may have to run `conda install conda-build` first)

pip:
  - Install from the repository root using pip: `pip install .`,
  - Or in editable mode (so edits are immediately reflected): `pip install -e .`