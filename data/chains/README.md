# Affimer
This program aims to design a specific Affimer of a specific antigen.

License : UPC

Authors:
* Bash and C implementation : Mohamed

Last revision : 2023-January-05

## Description
This program aims to design a specific Affimer of a specific antigen. 
This program was developed more specifically to allow modeling
Affimers whose purpose is to recognize specifically identified Flavivirus proteins, in particular NS1.
We will be using Patch Search [1, 2] for site comparison.

## Usage
A set of test data is necessary on the test folder

### Install the program
1. Download the programm on https://github.com/Moohmoo/Affimer by cloning him.
2. `cd src` 
3. `bash ./build.sh`
WARNING : Look at the parameters section.

### Parameters
```
usage: bash ./build.sh target chain_target affimer chain_affimer

This program aims to design a specific Affimer of a specific antigen. 

positional arguments:
    target                the name of target protein
    chain_target          the name of the target protein chain
    affimer               the name of affimer
    chain_affimer         the name of the affimer chain
```

### Program input
All files are in PDB format. The affimer and the target protein, they must be in the test folder and in the folder dedicated to it.
It is necessary to install an antigen antibody complex dataset and place it in the 'data/dataset' directory.

### Program output

Currently, the program only generates the modified sequence. To obtain the structure, there are sites such as Swiss Model to generate the structures by homologies.

### Examples
Building affimer
```
bash ./build.sh 5gs6 A 7ny8 C
```

### Important notes
All the pdb files contains antibody to be tested must be present in the `data/dataset` folder in order to run the program.

# Reference
*I. Rasolohery, G. Moroy, F. Guyon PatchSearch : a fast computational method for off-target detection,
Journal of chemical information and modeling, 2017*

*J. Rey, I Rasolohery, P. Tufféry, F. Guyon, G. Moroy PatchSearch : a web server for off-target protein
identification, Nucleic acids research, 2019*

# Contact
If you have any question or suggestion, please contact @Moohmoo
