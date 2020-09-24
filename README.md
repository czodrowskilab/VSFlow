# VSFlow 

VSFlow is an open-source command line tool intended to facilitate virtual screening workflows. 
It is especially designed for screening large databases with numerous query molecules and for 
the preparation of large sets of molecules for subsequent steps such as, e.g., docking. The
program is written in Python and combines different tools and screening options. The main 
intention is the development of an all-in-one tool for different virtual screening purposes 
with a user-friendly access. Additionally, a wide range of input formats is accepted which 
reduces the need of tedious and time-consuming file conversions.  
To date, many of the functionalities of VSFlow are based on the RDKit [1] and therefore 
ready-to-use without additional licenses. However, VSFlow also serves as a wrapper for 
parts of the Openeye tools [2] which require additional licenses in order
to access these features.  
The functionalities implemented in VSFlow include substructure search, fingerprint-based 
similarity search, generation of 3D structures of molecules (all RDKit based), Omega, ROCS, 
combination of Omega and ROCS (all using Openeye tools) and a tool to search the PDB database 
for structures and ligands. [3] The results can be exported in various file formats, including Excel files.
As additional feature, VSFlow supports the direct visualization of the search results as PDF file
and/or PyMOL [4] file to allow for a quick inspection of the results by the user. 

## Prerequisites

The Python dependencies are:
* Python >= 3.7
* RDKit >= 2019.09.3
* Pandas >= 0.25
* FPDF >= 1.7.2
* PDFRW >= 0.4
* XlsxWriter >= 1.2.7
* Xlrd >= 1.2.0
* PyMOL >= 2.3.4

In order to use the functionalities of the vsflow modes omega, rocs and omrocs,
a license for the *OpenEye Applications 2020.0.4* is required.<sup>[2]</sup>
Please make sure *OpenEye Applications* are installed and are contained in your 'PATH'
variable. Additionally, the environment variable 'OE_LICENSE' (containing the path to
your *OpenEye* license file) has to be set.

### Installation

First, you need a working Miniconda/Anaconda installation. For download and installation
of Miniconda visit https://conda.io/en/latest/miniconda.html.

Now you can create an environment named "vsflow" with all needed dependencies using
the provided environment.yml file and activate it:
```bash
conda env create -f environment.yml
conda activate vsflow
```

You can also create a new environment by yourself and install all dependencies without the
environment.yml file:
```bash
conda create -n vsflow python=3.7
conda activate vsflow
conda install -c rdkit -c conda-forge -c viascience -c schrodinger rdkit pandas xlrd xlsxwriter pdfrw fpdf pymol 
```

### Usage
To use vsflow you have to be in the repository folder and your conda environment has to be activated.
#### Help
To display the available modes of vsflow enter the following:
```bash
python vsflow.py -h
```
* substructure: perform a substructure search in a database (provided as sdf file)
* omega: generate 3D conformers of your input molecules (needs OpenEye Applications)
* rocs: perform ROCS (shape-based) search in a database (provided as .oeb.gz file, needs OpenEye Applications)
* omrocs: combination of omega and rocs, perform ROCS search starting from 2D structures (needs OpenEye Applications)
* rd3d: generates 3D conformers of input molecules using RDKit 
* rd3d_rocs: combination of rd3d and rocs, starting from 2D input structures (needs OpenEye Applications)
* fpsim: perform fingerprint-based similarity search in database 
* pdb: perform search in pdb database for protein structures, ligands etc.

To display arguments for the respective mode enter the following:
```bash
python vsflow.py {mode} -h
```
For example:
```bash
python vsflow.py substructure -h
```

#### Databases
To use your compound libraries as databases in vsflow, please provide them as 2D (or 3D) sdf files
for substructure/fingerprint searches and as oeb.gz files to perform shaped-based 
ROCS searches. To generate .oeb.gz files from sdf files you can use the mode *omega* of 
vsflow as follows:
```bash
python vsflow.py omega -in {filename}.sdf -flag rocs -fo oeb.gz -out {filename}
```
To permanently integrate your database files in vsflow, add their path to the DATABASES.txt file as follows:  
name,path,identity  

name = chose a name to call the database in vsflow (free choice)  
path = path to the file  
identity = name of sdf tag which contains name/ID of compound (optional, if not provided,
the *_Name* tag is used by default)  

Example for a valid database file:   
name,path,identity 
synthesized_cpds,/lab/cpds/all_cpds.sdf,lab_id  
purchased_cpds,/lab/cpds/all_purchased.sdf,  
computed_3D,/insilico/qm/all_3D_cpds.oeb.gz,  

Please make sure to write every database in a separate row.  
To call one of the databases in vsflow, simply provide the name using the -db flag. For example:
```bash
python vsflow.py substructure -in {input_file} -out {name_of_output_file} -db synthesized_cpds
```
Alternatively, you can also directly enter the path to your database file using the -db flag
and specifiy the name tag using the -id flag (optional). For example:
```bash
python vsflow.py substructure -in {input_file} -out {name_of_output_file} -db /lab/cpds/all_cpds.sdf -id lab_id
```

#### Example Usage

Example 1:
```bash
python vsflow.py substructure -in query_mols.sdf -out sub_matches.sdf -db synthesized_cpds
```
Performs a substructure search for the compounds in query_mols.sdf in the database 
*synthesized_cpds* (which has been integrated for permanent use, see section databases above). 
The output is a sdf file with matching compounds from the database.

Example 2:
```bash
python vsflow.py substructure -in query_mols.sdf -out sub_matches.xlsx -db synthesized_cpds -pdf
```
Same as in Example 1, but the matching compounds are written as smiles to an excel file.
Additionally, a pdf file (named sub_matches.pdf) with a visualization of the compounds 
and highlighted matching atoms is generated (-pdf flag).  

Example 3:
```bash
python vsflow.py substructure -in query_mols.sdf -out sub_matches.xlsx -db synthesized_cpds -pdf -props
```
Same as Example 2, but additionally calculated molecular properties (MW, cLogP, TPSA,
HDon, HAcc, RotBonds, AromRings) of the matching compounds are written to the output file.

Example 4:
```bash
python vsflow.py substructure -in query_mols.sdf -out sub_matches.xlsx -db synthesized_cpds -pdf -props -np 12
```
Same as Example 3, but substructure search is run on multiple cores (useful for input files 
containing multiple compounds).

Example 5:  
Substructure search can also be performed by directly copying SMILES and SMARTS to the command line. 
Please make sure the SMILES or SMARTS is put in double quotes:
```bash
python vsflow.py substructure -smi "CN(C)C1=NC=C(C(O)=O)C=C1" -out sub_matches.sdf -db synthesized_cpds
```
```bash
python vsflow.py substructure -sma "c:1:c(:c(:c(:c:c:1)C)N)F" -out sub_matches.sdf -db synthesized_cpds
```
It is also possible to provide multiple SMARTS or SMILES at once:
```bash
python vsflow.py substructure -smi "CN(C)C1=NC=C(C(O)=O)C=C1" -smi "O=C(O)C1=CC=CC=C1" -smi "O=C(O)C1=CC=C(O)C=C1" -out sub_matches.sdf -db synthesized_cpds
```
Using the -mf flag, separate output files for every query molecule are generated:
```bash
python vsflow.py substructure -smi "CN(C)C1=NC=C(C(O)=O)C=C1" -smi "O=C(O)C1=CC=CC=C1" -smi "O=C(O)C1=CC=C(O)C=C1" -out sub_matches.sdf -db synthesized_cpds -mf
```
Example 6:
```bash
python vsflow.py fpsim -smi "CN(C)C1=NC=C(C(O)=O)C=C1" -out fp_sim_top10.sdf -db synthesized_cpds -pdf -props
```
Performs a fingerprint-based similarity search (default is feature-based Morgan Fingerprint with
radius 3 and 4096 bits, may be changed using the respective flags, enter python vsflow.py fpsim -h
to see all arguments) and by default outputs the 10 most similar compounds as sdf file. Additionally,
the compounds are visualized in as pdf file.

Additional examples and a dedicated documentation will be added shortly.

## References

[1] RDKit, Open-Source cheminformatics; http://www.rdkit.org.  
[2] *OMEGA* 4.0.0.4 and *ROCS* 3.4.0.4: OpenEye Scientific Software, Santa Fe, NM; http://www.eyesopen.com.  
[3] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242. http://www.rcsb.org  
[4] The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.