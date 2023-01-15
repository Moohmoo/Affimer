# Affimer
Calculates the 2D/semi-3D folding of a protein using a REMC method described in the article by Thachuk C. et al.

License : UPC

Authors:
* Bash and C implementation : Mohamed

Last revision : 2023-January-05

## Description
Our program is based on an ab initio protein folding prediction method described by Thachuk C. et al.
It allows to determine a 2D and semi-3D folding by a simple evaluation of energy functions. It is a program
interpreted on the Python 3 version and is exploitable on any type of operating system.

## Usage
A set of test data is necessary on the test folder

### Install the program
1. Download the programm on https://github.com/Moohmoo/Affimer by cloning him.
2. `cd src` 
3. `bash ./execute.sh`

### Parameters
```
usage: parser.py [-h] [-o OUTPUT] [-g GRAPHIC]
             optimal n temp_min temp_max r file [file ...]

Calculates the 2D/semi-3D folding of a protein using a REMC method.

positional arguments:
    optimal               the optimal energy threshold
    n                     the number of local steps in a Monte Carlo search
    temp_min              the minimum temperature values
    temp_max              the maximum temperature values
    r                     the number of replicas to simulate
    file                  a protein fasta file

optional arguments:
    -h, --help            show this help message and exit
    -o OUTPUT, --output OUTPUT
                            directs the output to a name of your choice
    -g GRAPHIC, --graphic GRAPHIC
                            directs the output to a 2D or 3D graphic
```

### Program input
Input files can be in PDB format. The files can contain only one sequence.

### Program output
Depending on the parameters chosen, it is possible to have different outputs:
1. The first output is a summary of the current process on the terminal. It shows the energy of the
transient state
2. The second output that can be obtained is a 2D or semi-3D graph. This is the visualisation of
of the predicted structure by the method.
3. The second output can be a `.png` file with the visualisation according to the type of graph.
The file will be available in the `res` folder.

### Examples
REMC simulation
```
python3 parser.py -g 2d -o my_conf.png -10 500 220 250 2 ../test/P01013.fasta
```
Here we define an optimal energy at -10, a number of Monte Carlo iterations at 500, a minimum and maximum temperature of 220
and 250 respectively and a number of replicates of 2. The options -g define the type of representation and -o the output.

### Important notes
All the fasta files to be tested must be present in the `test` folder in order to run the program.
The file name of the output does not need an absolute path: the simple name is enough (see Example).

# Reference
*Thachuk C, Shmygelska A, Hoos HH. A replica exchange Monte Carlo
algorithm for protein folding in the HP model. BMC Bioinformatics. 2007 Sep
17;8:342. PubMed PMID: 17875212; PubMed Central PMCID: PMC2071922.*

# Contact
If you have any question or suggestion, please contact @Moohmoo
