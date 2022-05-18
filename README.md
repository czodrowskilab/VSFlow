# VSFlow - Virtual Screening Workflow

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
conda install -c conda-forge rdkit xlrd xlsxwriter pdfrw fpdf pymol-open-source molvs matplotlib 
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

Now, you can install VSFlow as follows:
```bash
pip install .
```

## General Usage
Always make sure the conda environment is activated.
Now you can run VSFlow as follows:  
```bash
vsflow {mode} {arguments}
```
For example, the following command will display all included modes (substructure, fpsim, shape, preparedb, managedb) and the general usage:
```bash
vsflow -h
```
To display all possible arguments for a particular mode, type as follows:
```bash
vsflow {mode} -h
```
For example, with the following command all arguments for mode substructure are shown:
```bash
vsflow substructure -h
```

## Example Usage

A detailed usage of VSFlow with many examples is provided in the GitHub Wiki:  
https://github.com/czodrowskilab/VSFlow/wiki

## References

[1] RDKit, Open-Source cheminformatics; http://www.rdkit.org.  
[2] The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.
