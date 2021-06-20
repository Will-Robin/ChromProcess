# ChromProcess

A module for batch processing of chromatography data.

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

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

Go to Install dependencies.

#### Using pip:

Create a virtual environment:

`pip install virtualenv`

`virtualenv chromprocess-env`

Activate the virtual environment:

Mac:

`source chromprocess-env/bin/activate`

Windows:

`chromprocess-env\Scripts\activate`

Go to Install dependencies.

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

### Check install worked
Type `python` in the command line/terminal and then return. This opens a Python interpreted. Type:

`from ChromProcess import info_params` then the return key.

If no errors occur (i.e. nothing happens), the installation worked.

### Using the code in Scripts
You can import the ChromProcess code by putting `import ChromProcess` at the top of your scripts. You can import from specific files from ChromProcess using syntax like `from ChromProcess import info_params` or ` from ChromProcess.info_params import colour_assignments`. You can then use the code written in the files with similar names in the ChromProcess folder.

Take a look in the Scripts folder for some examples.

### Using ChromProcess in Jupyter Notebooks

Jupyter Notebooks and Jupyter Lab<sup>[1](jupyter-link)</sup> are excellent tools for executing blocks of code alongside notes and data plots.
For more information on which one to choose and how to install them, go to the annotated link.

There is a little to add to the setup of a Jupyter Notebook, since we have installed ChromProcess in a virtual environment. In short, an IPython kernel must be installed in the ChromProcess environment, and a file must be created so that the Jupyter Notebook can access the kernel.

Activate the ChromProcess environment (see above), then run:

`pip install --user ipykernel`

then:

`python -m ipykernel install --user --name=chromprocess-env`

Deactivate the environment (`conda deactivate` or `deactivate`).

The next time you open Jupyter Notebook, a ChromProcess kernel should be available to choose in the launcher, or in a dropdown list in the top right of the notebook window. Choose the ChromProcess kernel, and get scripting!

## Contributing

If you want to make an edit to this package:

1. Clone this repository from GitHub (or pull latest version if already cloned).
2. Create a new branch with a descriptive name for the change you want to make.
3. Edit the package until you are satisfied with the result.
4. Commit your changes to the branch you created (commit quick and often, and use descriptive commit messages)
5. Push everything to the Github repository
6. Create a pull request on Github, requesting to merge your branch with the main branch. Another person can then review the changes you have made, propose improvements, and accept them when everything is satisfactory.

The reviewing of pull requests is important to ensure that nothing is broken by the changes that you make. It also creates a nice opportunity to learn from somebody else how they would have approached the same edits.

## Acknowledgement

ChromProcess was created by William E. Robinson during free time and working in the group of Prof. Wilhelm T. S. Huck (Radboud University Nijmegen) supported by NWO and the Simons Foundation. The intention behind the development of this software is to provide a useful tool for many people to use. If you use lots of this code in your work, an acknowledgement in publications would be appreciated (and possibly a link to this repository if made public). Also consider citing/acknowledging the makers of ChromProcess's dependencies. The software is all free and the product of a huge amount of work by many people (and all we have to do is `pip install`!). Lots (all?) of them worked for free in the spirit of open source software.

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0
[jupyter-link]: https://jupyter.org
