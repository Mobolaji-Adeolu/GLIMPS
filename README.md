# The Gupta Lab Integrated Microbial Phylogeny and Supermatrix (GLIMPS) Pipeline
_A User-Friendly Pipeline for the production of Core Genome and Concatenated Protein based_
_Phylogenetic Trees and Protein based Comparative Genomic Analyses._

## Installation
* Clone repository to local directory.
* Extract PhyEco Marker HMM files
* Use "chmod +x" to make dependencies executable in Linux and OSX.
* Use "sudo apt-get install python-tk" to install the Tkinter module required for the GLIMPS GUI in Linux.
* Use "python GLIMPS_UI.py" to execute GLIMPS_UI.

## Command Line Usage
**GLIMPS_Pipeline.py** [*args]

-h, --help  _Show help message and exit_

-i, --Input_Directory  _The directory in which input genome files are located_

-o, --Output_Directory  _The directory in which Output files are to be placed_

-d, --Protein_Distribution  _The minimum proportion of the input genomes in which a protein family must be present_

-t, --Threads  _The minimum proportion of the input genomes in which a protein family must be present_

-p, --Target_Proteins  _The location of a file containing the protein targets for protein family identification (optional)_

-m, --Marker_Proteins  _The predefined marker protein family set to be used for analysis (optional)_

--Single_Copy  _Determines whether the pipeline utilizes only single copy homologs_

--PAMatrix  _Determines whether the pipeline produces a PA Matrix_

--POCP  _Determines whether the pipeline produces a POCP Matrix_

--AAI  _Determines whether the pipeline produces an AAI Matrix_

--Fast_Cluster  _Skips the HMM-based iterative clustering steps after CD-Hit during the core genome identification process_

--Fast_Phylogeny  _Skips the RAxML based tree building step after FastTree_

--No_Tree  _Skips all phylogenetic tree building steps_

-f, --Alignment_Filtering  _The method by which the alignments will be filtered_

--cdhit  _Path to the cd-hit executable_

--jackhmmer  _Path to the jackhmmer executable_

--hmmbuild  _Path to the hmmbuild executable_

--hmmsearch  _Path to the hmmsearch executable_

--clustalo  _Path to the clustalo executable_

--trimal  _Path to the trimal executable_

--fasttree  _Path to the fasttree executable_

--raxml  _Path to the raxml executable_

## History
Version 0.5: Initial Public Release

## Credits
Mobolaji Adeolu (Department of Biochemistry and Biomedical Sciences, McMaster University)

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