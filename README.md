# The Gupta Lab Integrated Microbial Phylogeny and Supermatrix (GLIMPS) Pipeline
_A User-Friendly Pipeline for the production of Core Genome and Concatenated Protein based_
_Phylogenetic Trees and Protein based Comparative Genomic Analyses._

## Installation
* Clone repository to local directory.
* Extract PhyEco Marker HMM files
* Use "chmod +x" to make dependencies executable in Linux and OSX.
* Use "sudo apt-get install python-tk" to install the Tkinter module required for the GLIMPS GUI in Linux.
* Use "python GLIMPS_UI.py" to execute GLIMPS_UI.

## Command LinUsage
**GLIMPS_Pipeline.py** [*args]
-h, --help
Show this help message and exit
-i, --Input_Directory
The directory in which input genome files are located
-o, --Output_Directory
The directory in which Output files are to be placed
-d, --Protein_Distribution
The minimum proportion of the input genomes in which a protein family must be present
-t, --Threads
The minimum proportion of the input genomes in which a protein family must be present
-p, --Target_Proteins
The location of a file containing the protein targets for protein family identification (optional)
-m, --Marker_Proteins
The predefined marker protein family set to be used for analysis (optional)
--Single_Copy
Determines whether the pipeline utilizes only single copy homologs
--PAMatrix
Determines whether the pipeline produces a PA Matrix
--POCP
Determines whether the pipeline produces a POCP Matrix
--AAI
Determines whether the pipeline produces an AAI Matrix
--Fast_Cluster
Skips the HMM-based iterative clustering steps after CD-Hit during the core genome identification process
--Fast_Phylogeny
Skips the RAxML based tree building step after FastTree
--No_Tree
Skips all phylogenetic tree building steps
-f, --Alignment_Filtering
The method by which the alignments will be filtered
--cdhit
Path to the cd-hit executable
--jackhmmer
Path to the jackhmmer executable
--hmmbuild
Path to the hmmbuild executable
--hmmsearch
Path to the hmmsearch executable
--clustalo
Path to the clustalo executable
--trimal
Path to the trimal executable
--fasttree
Path to the fasttree executable
--raxml
Path to the raxml executable

## History
Version 0.5: Initial Public Release

## Credits
Mobolaji Adeolu
**McMaster University**
_Department of Biochemistry and Biomedical Sciences_

## License
**MIT License**

Copyright (c) 2016 Mobolaji Adeolu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.