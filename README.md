# VSFlow 

VSFlow is an open-source command-line tool built on top of the RDKit [1] for the ligand-based virtual screening 
of large compound libraries (databases). It includes a substructure-based, a fingerprint-based 
and a shape-based virtual screening tool. Additionally, it provides a tool to prepare databases for 
screening (molecule standardization, fingerprint and conformer generation). Screenings can be 
parallelized with Python's built-in multiprocessing package. Additionally, VSFlow accepts a wide 
range of input file formats. The screening results can be exported in various file formats, including Excel files.
As additional feature, VSFlow supports the visualization of the screening results as PDF file
and/or PyMOL file [2], allowing for a quick inspection of the results by the user. VSFlow is fully 
written in Python.

## Installation

First of all, you need a working installation of Anaconda (https://www.anaconda.com/products/individual) or Miniconda (https://conda.io/en/latest/miniconda.html). Both are available for all major platforms.  

Second, you need to clone the VSFlow GitHub repository to your system or download the zip file and unpack it (in the following called the repository folder).   

**All following instructions assume working with a bash shell!**  

Navigate into the repository folder.  

Now, you can install the required dependencies with the provided environment.yml file within the repository folder as follows:
```bash
conda env create --quiet --force --file environment.yml
conda activate vsflow
```
Alternatively, you can also create a new conda environment and install the dependencies manually:
```bash
conda create -n vsflow python=3.7
conda activate vsflow
conda install -c rdkit -c conda-forge -c viascience -c schrodinger rdkit xlrd xlsxwriter pdfrw fpdf pymol molvs matplotlib
```
The Python dependencies are:  
* Python = 3.7
* RDKit >= 2019.09.3
* FPDF >= 1.7.2
* PDFRW >= 0.4
* XlsxWriter >= 1.2.7
* Xlrd >= 1.2.0
* PyMOL >= 2.3.4
* Molvs >= 0.1.1
* Matplotlib >= 3.3.4

## General Usage
To run VSFlow as described below, you have to be within the repository folder and the conda environment has to be activated.
Now you can run VSFlow as follows:  
```bash
python vsflow.py {mode} {arguments}
```
For example, the following command will display all included modes (substructure, fpsim, shape, preparedb, managedb) and the general usage:
```bash
python vsflow.py -h
```
To display all possible arguments for a particular mode, type as follows:
```bash
python vsflow.py {mode} -h
```
For example, with the following command all arguments for mode substructure are shown:
```bash
python vsflow.py substructure -h
```
To execute VSFlow from outside the repository folder, the repository path has to be added to your PATH variable and 
the script has to be made executable. 
One possible way is described below:
```bash
export PATH="/path/to/repository:$PATH"
```
To permanently add the repository path to your PATH, add the above statement to your .bashrc or .bash_profile file.  
Next, make the script file executable:
```bash
chmod +x vsflow.py
```
Now you can run VSFlow also from outside the repository as follows:
```bash
vsflow.py {mode} {arguments}
```
## Example Usage

A detailed usage of VSFlow with many examples is provided in the GitHub Wiki:  
https://github.com/czodrowskilab/VSFlow/wiki

## References

[1] RDKit, Open-Source cheminformatics; http://www.rdkit.org.  
[2] The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.