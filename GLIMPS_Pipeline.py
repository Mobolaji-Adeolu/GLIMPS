#!/usr/bin/env python

"""Gupta Lab Integrated Microbial Phylogenetic and Supermatrix Pipeline
Written by Mobolaji Adeolu (Adeolum@McMaster.ca)
Department of Biochemistry and Biomedical Sciences
McMaster University.
Copyright 2016."""

import os
import shutil
import sys
import subprocess
import multiprocessing
import random
import time
import argparse
import re
import tempfile
import tarfile


class GLIMPS_Writer:
    """Class to pipe stdout and sterr to parent process in asyncrounous threads"""
    def __init__(self, stdout_messenger, stderr_messenger):
        self.stdout_messenger = stdout_messenger
        self.stderr_messenger = stderr_messenger

    def write(self, s):
        if __name__ == '__main__':
            sys.stdout.write(s)
        else:
            self.stdout_messenger.put(s)

    def error(self, s):
        if __name__ == '__main__':
            sys.stderr.write(s)
        else:
            self.stderr_messenger.put(s)


def check_negative(value):
    """Checks for a negative interger value in command line argument"""
    int_value = int(value)
    if int_value < 0:
        raise argparse.ArgumentTypeError("%s is an invalid value. Use a positive interger." % value)
    return int_value


def check_arguments():
    """Reads command line arguments and generates help text"""
    arguments = argparse.ArgumentParser(
        description="Gupta Lab Integrated Microbial Phylogeny and Supermatrix Pipeline",
        epilog="Written by Mobolaji Adeolu (Adeolum@McMaster.ca), Department of Biochemistry and Biomedical Sciences, McMaster University. Copyright 2016.")
    Gene_Select = arguments.add_mutually_exclusive_group()
    arguments.add_argument(
        "-i", "--Input_Directory",
        action="store",
        type=str,
        default=os.path.join(os.path.dirname(sys.argv[0]), "Input"),
        required=False,
        metavar="Input Directory",
        help="The directory in which input genome files are located",
    )
    arguments.add_argument(
        "-o", "--Output_Directory",
        action="store",
        type=str,
        default=os.path.join(os.path.dirname(sys.argv[0]), "Output"),
        required=False,
        metavar="Output Directory",
        help="The directory in which Output files are to be placed",
    )
    arguments.add_argument(
        "-d", "--Protein_Distribution",
        action="store",
        type=float,
        default=0.8,
        required=False,
        choices=[x * 0.01 for x in range(0, 101)],
        metavar="Protein Distribution",
        help="The minimum proportion of the input genomes in which a protein family must be present",
    )
    arguments.add_argument(
        "-t", "--Threads",
        action="store",
        type=check_negative,
        default=0,
        metavar="Number of Threads",
        help="The minimum proportion of the input genomes in which a protein family must be present",
    )
    Gene_Select.add_argument(
        "-p", "--Target_Proteins",
        action="store",
        type=str,
        default="",
        required=False,
        metavar="Target Proteins",
        help="The location of a file containing the protein targets for protein family identification (optional)",
    )
    Gene_Select.add_argument(
        "-m", "--Marker_Proteins",
        action="store",
        type=str,
        default="",
        required=False,
        choices=["actinobacteria", "alphaproteobacteria", "archaea", "bacteria_and_archaea", "bacteria",
                 "bacteroidetes", "betaproteobacteria", "chlamydiae", "chloroflexi", "cyanobacteria",
                 "deinococcus-thermus", "deltaproteobacteria", "epsilonproteobacteria", "firmicutes",
                 "gammaproteobacteria", "spirochaetes", "thermotogae"],
        metavar="PhyEco Marker Protein Family",
        help="The predefined marker protein family set to be used for analysis (optional)",
    )
    arguments.add_argument(
        "--Single_Copy",
        action="store_true",
        default=False,
        required=False,
        help="Determines whether the pipeline utilizes only single copy homologs",
    )
    arguments.add_argument(
        "--PAMatrix",
        action="store_true",
        default=False,
        required=False,
        help="Determines whether the pipeline produces a PA Matrix",
    )
    arguments.add_argument(
        "--POCP",
        action="store_true",
        default=False,
        required=False,
        help="Determines whether the pipeline produces a POCP Matrix",
    )
    arguments.add_argument(
        "--AAI",
        action="store_true",
        default=False,
        required=False,
        help="Determines whether the pipeline produces an AAI Matrix",
    )
    arguments.add_argument(
        "--Fast_Cluster",
        action="store_true",
        default=False,
        required=False,
        help="Skips the HMM-based iterative clustering steps after CD-Hit during the core genome identification process",
    )
    arguments.add_argument(
        "--Fast_Phylogeny",
        action="store_true",
        default=False,
        required=False,
        help="Skips the RAxML based tree building step after FastTree",
    )
    arguments.add_argument(
        "--No_Tree",
        action="store_true",
        default=False,
        required=False,
        help="Skips all phylogenetic tree building steps",
    )
    arguments.add_argument(
        "-f", "--Alignment_Filtering",
        action="store",
        type=str,
        default="Trim",
        required=False,
        choices=["Trim", "Weight"],
        metavar="Alignment Filtering",
        help="The method by which the alignments will be filtered",
    )
    arguments.add_argument(
        "--cdhit",
        action="store",
        type=str,
        default="cd-hit",
        required=False,
        help="Path to the cd-hit executable",
    )
    arguments.add_argument(
        "--jackhmmer",
        action="store",
        type=str,
        default="jackhmmer",
        required=False,
        help="Path to the jackhmmer executable",
    )
    arguments.add_argument(
        "--hmmbuild",
        action="store",
        type=str,
        default="hmmbuild",
        required=False,
        help="Path to the hmmbuild executable",
    )
    arguments.add_argument(
        "--hmmsearch",
        action="store",
        type=str,
        default="hmmsearch",
        required=False,
        help="Path to the hmmsearch executable",
    )
    arguments.add_argument(
        "--clustalo",
        action="store",
        type=str,
        default="clustalo",
        required=False,
        help="Path to the clustalo executable",
    )
    arguments.add_argument(
        "--trimal",
        action="store",
        type=str,
        default="trimal",
        required=False,
        help="Path to the trimal executable",
    )
    arguments.add_argument(
        "--fasttree",
        action="store",
        type=str,
        default="fasttree",
        required=False,
        help="Path to the fasttree executable",
    )
    arguments.add_argument(
        "--raxml",
        action="store",
        type=str,
        default="raxml",
        required=False,
        help="Path to the raxml executable",
    )
    args = arguments.parse_args()
    return args


def make_dir(dir_name, Output):
    """Make/overwrite folders"""
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    else:
        shutil.rmtree(dir_name)
        count = 0
        while os.path.exists(dir_name) and count < 15:
            count += 1
            time.sleep(1)
        else:
            try:
                os.mkdir(dir_name)
            except OSError:
                Output.error("\nUnable to access directory: " + dir_name)
                sys.exit()


def Build_Output_Dirs(Output_Directory, stdout_messenger, stderr_messenger):
    """Builds output directory structure"""
    Output = GLIMPS_Writer(stdout_messenger, stderr_messenger)
    Genome_Dir = os.path.join(Output_Directory, "Data", "Genomes")
    Protein_Dir = os.path.join(Output_Directory, "Data", "Proteins")
    Alignment_Dir = os.path.join(Output_Directory, "Data", "Protein Alignments")
    Concatenated_Dir = os.path.join(Output_Directory, "Data", "Concatenated Alignments")
    Tree_Dir = os.path.join(Output_Directory, "Data", "Phylogenetic Trees")
    Log_Dir = os.path.join(Output_Directory, "Logs")
    GLIMPSe_Output_Dir = os.path.join(Output_Directory, "GLIMPSe Output")
    Dependency_Dir = os.path.join(os.path.dirname(sys.argv[0]), "Dependencies")
    Marker_Dir = os.path.join(os.path.dirname(sys.argv[0]), "PhyEco Marker Protein Families")
    if not os.path.exists(Output_Directory):
        os.mkdir(Output_Directory)
    make_dir(os.path.join(Output_Directory, "Data"), Output)
    make_dir(Genome_Dir, Output)
    make_dir(Protein_Dir, Output)
    make_dir(Alignment_Dir, Output)
    make_dir(Concatenated_Dir, Output)
    make_dir(Tree_Dir, Output)
    make_dir(Log_Dir, Output)
    make_dir(GLIMPSe_Output_Dir, Output)
    return Genome_Dir, Protein_Dir, Alignment_Dir, Concatenated_Dir, Tree_Dir, Log_Dir, GLIMPSe_Output_Dir, Dependency_Dir, Marker_Dir


def Prepare_Dependencies(Dependency_Dir, stdout_messenger, stderr_messenger):
    """Checks dependencies and sets dependency variables"""
    Output = GLIMPS_Writer(stdout_messenger, stderr_messenger)
    deps = check_arguments()
    Threads = deps.Threads
    with tempfile.TemporaryFile() as dump:
        try:
            subprocess.check_call(deps.cdhit, stdout=dump, stderr=dump)
            CDHIT = deps.cdhit
        except (subprocess.CalledProcessError, OSError):
            CDHIT = ""
        try:
            subprocess.check_call(deps.jackhmmer, stdout=dump, stderr=dump)
            JACKHMMER = deps.jackhmmer
        except (subprocess.CalledProcessError, OSError):
            JACKHMMER = ""
        try:
            subprocess.check_call(deps.hmmbuild, stdout=dump, stderr=dump)
            HMMBUILD = deps.hmmbuild
        except (subprocess.CalledProcessError, OSError):
            HMMBUILD = ""
        try:
            subprocess.check_call(deps.hmmsearch, stdout=dump, stderr=dump)
            HMMSEARCH = deps.hmmsearch
        except (subprocess.CalledProcessError, OSError):
            HMMSEARCH = ""
        try:
            subprocess.check_call(deps.clustalo, stdout=dump, stderr=dump)
            ClustalOmega = deps.clustalo
        except (subprocess.CalledProcessError, OSError):
            ClustalOmega = ""
        try:
            subprocess.check_call(deps.trimal, stdout=dump, stderr=dump)
            TrimAl = deps.trimal
        except (subprocess.CalledProcessError, OSError):
            TrimAl = ""
        try:
            subprocess.check_call([deps.fasttree, "-expert"], stdout=dump, stderr=dump)
            FastTree = deps.fasttree
        except (subprocess.CalledProcessError, OSError):
            FastTree = ""
        try:
            subprocess.check_call(deps.raxml, stdout=dump, stderr=dump)
            RAxML = deps.raxml
        except (subprocess.CalledProcessError, OSError):
            RAxML = ""
    if sys.platform.startswith("linux"):
        OS_Dir = "Linux"
        if CDHIT == "":
            CDHIT = os.path.join(Dependency_Dir, OS_Dir, "CD-HIT", "cd-hit")
        if JACKHMMER == "":
            JACKHMMER = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "jackhmmer")
        if HMMBUILD == "":
            HMMBUILD = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmbuild")
        if HMMSEARCH == "":
            HMMSEARCH = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmsearch")
        if ClustalOmega == "":
            ClustalOmega = os.path.join(Dependency_Dir, OS_Dir, "Clustal Omega", "clustalo")
        if TrimAl == "":
            TrimAl = os.path.join(Dependency_Dir, OS_Dir, "TrimAl", "trimal")
        if FastTree == "":
            FastTree = os.path.join(Dependency_Dir, OS_Dir, "FastTree", "fasttree")
        if RAxML == "":
            RAxML = os.path.join(Dependency_Dir, OS_Dir, "RAxML", "raxml")
        if Threads == 0:
            if isinstance(multiprocessing.cpu_count(), int):
                Threads = multiprocessing.cpu_count()
            else:
                Threads = 1
    elif sys.platform.startswith("darwin"):
        OS_Dir = "OSX"
        if CDHIT == "":
            CDHIT = os.path.join(Dependency_Dir, OS_Dir, "CD-HIT", "cd-hit")
        if JACKHMMER == "":
            JACKHMMER = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "jackhmmer")
        if HMMBUILD == "":
            HMMBUILD = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmbuild")
        if HMMSEARCH == "":
            HMMSEARCH = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmsearch")
        if ClustalOmega == "":
            ClustalOmega = os.path.join(Dependency_Dir, OS_Dir, "Clustal Omega", "clustalo")
        if TrimAl == "":
            TrimAl = os.path.join(Dependency_Dir, OS_Dir, "TrimAl", "trimal")
        if FastTree == "":
            FastTree = os.path.join(Dependency_Dir, OS_Dir, "FastTree", "fasttree")
        if RAxML == "":
            RAxML = os.path.join(Dependency_Dir, OS_Dir, "RAxML", "raxml")
        if Threads == 0:
            if isinstance(multiprocessing.cpu_count(), int):
                Threads = multiprocessing.cpu_count()
            else:
                Threads = 1
    elif sys.platform.startswith("win32"):
        OS_Dir = "Windows"
        if CDHIT == "":
            CDHIT = os.path.join(Dependency_Dir, OS_Dir, "CD-HIT", "cd-hit.exe")
        if JACKHMMER == "":
            JACKHMMER = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "jackhmmer.exe")
        if HMMBUILD == "":
            HMMBUILD = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmbuild.exe")
        if HMMSEARCH == "":
            HMMSEARCH = os.path.join(Dependency_Dir, OS_Dir, "HMMer", "hmmsearch.exe")
        if ClustalOmega == "":
            ClustalOmega = os.path.join(Dependency_Dir, OS_Dir, "Clustal Omega", "clustalo.exe")
        if TrimAl == "":
            TrimAl = os.path.join(Dependency_Dir, OS_Dir, "TrimAl", "trimal.exe")
        if FastTree == "":
            FastTree = os.path.join(Dependency_Dir, OS_Dir, "FastTree", "fasttree.exe")
        if RAxML == "":
            RAxML = os.path.join(Dependency_Dir, OS_Dir, "RAxML", "raxml.exe")
        if Threads == 0:
            if isinstance(multiprocessing.cpu_count(), int):
                Threads = multiprocessing.cpu_count()
            else:
                Threads = 1
    if CDHIT == "" or JACKHMMER == "" or HMMBUILD == "" or HMMSEARCH == "" or ClustalOmega == "" or TrimAl == "" or FastTree == "" or RAxML == "":
        Output.error("\nOperating system unrecognized.\n")
        sys.exit()
    with tempfile.TemporaryFile() as dump:
        try:
            subprocess.check_call(CDHIT, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute CD-Hit dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(JACKHMMER, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute JackHMMer dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(HMMBUILD, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute HMMBuild dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(HMMSEARCH, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute HMMSearch dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(ClustalOmega, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute Clustal Omega dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(TrimAl, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute TrimAl dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call([FastTree, "-expert"], stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute FastTree dependency. Check binaries and permissions.\n")
            sys.exit()
        try:
            subprocess.check_call(RAxML, stdout=dump, stderr=dump)
        except subprocess.CalledProcessError:
            pass
        except OSError:
            Output.error("\nUnable to execute RAxML dependency. Check binaries and permissions.\n")
            sys.exit()
    return CDHIT, JACKHMMER, HMMBUILD, HMMSEARCH, ClustalOmega, TrimAl, FastTree, RAxML, Threads


# Initial processing of genome files are handled by the following functions
def Replace_Chars(Text):
    """Replaces special characters in genome names"""
    Replacement_Dict = {" ": "_", "'": "", ":": "_", "/": "_", ",": "_", "(": "_", ")": "_"}
    for target, replacement in Replacement_Dict.iteritems():
        Text = Text.replace(target, replacement)
    return Text

# TODO Detect and tranlate nucleotide sequences
def Process_Genome_Files(Input_Dir, Genome_Dir, Log_Dir, Output):
    """Copies genome files to data directory and assigns IDs to proteins"""
    Genome_Dictionary = {}
    AllGenomeProts = []
    Genomes = []
    Accepted_Filetypes = [".fasta", ".fsa", ".faa", ".fas"]
    for fasta in os.listdir(Input_Dir):
        if os.path.splitext(fasta)[1].lower() in Accepted_Filetypes:
            Genomes.append(fasta)
        else:
            if len(os.path.splitext(fasta)[1]) > 0:
                Output.error(fasta + " is of an unsupported file type and will not be included in the analysis\n")

    def Assign_Protein_ID(Input_Genome, Output_Directory, GenID):
        """Local function to assign unique 10 digit ID to each protein"""
        ProteinID = 0
        with open(Input_Genome) as InGen:
            with open(os.path.join(Output_Directory, Replace_Chars(os.path.basename(Input_Genome))), "w+") as OutGen:
                Char_Count = 0
                for line in InGen:
                    if line.startswith(">") and line.endswith("\n"):
                        Char_Count = 0
                        ProteinID += 1
                        if ProteinID == 1:
                            OutGen.write(">" + str(GenID).zfill(4) + str(ProteinID).zfill(6) + "\n")
                        else:
                            OutGen.write("\n" + ">" + str(GenID).zfill(4) + str(ProteinID).zfill(6) + "\n")
                    elif line[0].isalpha() and ProteinID > 0:
                        line_content = []
                        for char in line:
                            if char in "OUBZJX":
                                char = "-"
                            if char.isalpha() and Char_Count < 60:
                                Char_Count += 1
                                line_content.append(char.upper())
                            elif char.isalpha() and Char_Count >= 60:
                                Char_Count = 1
                                line_content.append("\n")
                                line_content.append(char.upper())
                        else:
                            OutGen.write("".join(line_content))

    for GenomeID, Genome_File in enumerate(Genomes, 1):
        no_space_file_name = Replace_Chars(Genome_File)
        Assign_Protein_ID(os.path.join(Input_Dir, Genome_File), Genome_Dir, GenomeID)
        Genome_Dictionary[str(GenomeID).zfill(4)] = os.path.splitext(no_space_file_name)[0]
        with open(os.path.join(Genome_Dir, no_space_file_name), "r") as Gen:
            if GenomeID == 1:
                AllGenomeProts += Gen.readlines()
            else:
                AllGenomeProts.append("\n")
                AllGenomeProts += Gen.readlines()
        Count = GenomeID - 1
        if Count > 1 and Count % 10 == 0:
            Output.write("%s input files completed initial processing...\n" % Count)
    else:
        with open(os.path.join(Genome_Dir, "ConcatenatedGenomeFiles"), "w") as AllGenomes:
            AllGenomes.writelines(AllGenomeProts)
    with open(os.path.join(Log_Dir, "Genome_ID.tsv"), "w") as GenDict:
        GenDict.write("Genome ID\tGenome Name\n")
        for key in sorted(Genome_Dictionary.keys()):
            GenDict.write(key + "\t" + Genome_Dictionary[key] + "\n")
    return Genome_Dictionary


# Identification of protein families using CD-HIT and HMMer is handled by the following functions
def CDHIT_Subprocess(CDHIT, Log_Dir, Output, Input_File, Output_File, Threshold, Word_Size):
    """Runs subprocess for CDHIT"""
    CDHIT_Output = ""
    cmd = [CDHIT, "-i", Input_File, "-o", Output_File, "-c", Threshold, "-n", Word_Size, "-s", "0.5", "-M", "0", "-T",
           "0"]
    try:
        CDHIT_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as call_err:
        with open(os.path.join(Log_Dir, "CDHIT.txt"), "a") as CDHIT_Log:
            CDHIT_Log.write(CDHIT_Output)
            CDHIT_Log.write("\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output))
        Output.error("\n\nCD-HIT Error. Check Logs.\n")
        sys.exit()
    except OSError as os_err:
        with open(os.path.join(Log_Dir, "CDHIT.txt"), "a") as CDHIT_Log:
            CDHIT_Log.write(CDHIT_Output)
            CDHIT_Log.write("\n\n" + str(os_err.strerror))
        Output.error("\n\nCD-HIT Error. Check Logs.\n")
        sys.exit()
    with open(os.path.join(Log_Dir, "CDHIT.txt"), "a") as CDHIT_Log:
        CDHIT_Log.write(CDHIT_Output)


def Parse_Clusters(Clusters, Input_File):
    """Converts CDHIT output into a nested list of protein families"""
    with open(Input_File + ".clstr", "r") as Cluster_File:
        Clusters_Reps = [x[0] for x in Clusters]
        Primary_Cluster = ""
        Additional_Clusters = []
        for line in Cluster_File:
            if line[0].isdigit():
                if line.endswith("*\n"):
                    Primary_Cluster = re.search(">..........", line).group(0)
                else:
                    Additional_Clusters.append(re.search(">..........", line).group(0))
            else:
                if Primary_Cluster != "":
                    if Primary_Cluster not in Clusters_Reps:
                        Clusters.append([Primary_Cluster])
                        Clusters_Reps.append(Primary_Cluster)
                    if len(Additional_Clusters) > 0:
                        for clust in Additional_Clusters:
                            Secondary_Cluster = clust
                            if Primary_Cluster in Clusters_Reps and Secondary_Cluster in Clusters_Reps:
                                Primary_Index = Clusters_Reps.index(Primary_Cluster)
                                Secondary_Index = Clusters_Reps.index(Secondary_Cluster)
                                Clusters[Primary_Index] += Clusters[Secondary_Index]
                                del Clusters[Secondary_Index]
                                del Clusters_Reps[Secondary_Index]
                            else:
                                Primary_Index = Clusters_Reps.index(Primary_Cluster)
                                Clusters[Primary_Index].append(Secondary_Cluster)
                        else:
                            Additional_Clusters = []
    return Clusters


def Run_CDHIT(Genome_Dir, CDHIT, Log_Dir, ConcatenatedGenomeFile, Output):
    """Completes iterative CDHIT clustering"""
    Clusters = []
    Input_File = ConcatenatedGenomeFile
    Output_File = os.path.join(Genome_Dir, "CDHIT_Clusters_0.9")
    Threshold = "0.9"
    Word_Size = "5"
    CDHIT_Subprocess(CDHIT, Log_Dir, Output, Input_File, Output_File, Threshold, Word_Size)
    Wait_Count = 0
    while not os.path.exists(Output_File):
        time.sleep(0.1)
        Wait_Count += 1
        if Wait_Count > 30:
            Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
            sys.exit()
    Clusters = Parse_Clusters(Clusters, Output_File)
    Input_File = Output_File
    Output_File = os.path.join(Genome_Dir, "CDHIT_Clusters_0.8")
    Threshold = "0.8"
    Word_Size = "5"
    CDHIT_Subprocess(CDHIT, Log_Dir, Output, Input_File, Output_File, Threshold, Word_Size)
    Wait_Count = 0
    while not os.path.exists(Output_File):
        time.sleep(0.1)
        Wait_Count += 1
        if Wait_Count > 30:
            Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
            sys.exit()
    Clusters = Parse_Clusters(Clusters, Output_File)
    Input_File = Output_File
    Output_File = os.path.join(Genome_Dir, "CDHIT_Clusters_0.7")
    Threshold = "0.7"
    Word_Size = "4"
    CDHIT_Subprocess(CDHIT, Log_Dir, Output, Input_File, Output_File, Threshold, Word_Size)
    Wait_Count = 0
    while not os.path.exists(Output_File):
        time.sleep(0.1)
        Wait_Count += 1
        if Wait_Count > 30:
            Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
            sys.exit()
    Clusters = Parse_Clusters(Clusters, Output_File)
    Input_File = Output_File
    Output_File = os.path.join(Genome_Dir, "CDHIT_Clusters_0.6")
    Threshold = "0.6"
    Word_Size = "4"
    CDHIT_Subprocess(CDHIT, Log_Dir, Output, Input_File, Output_File, Threshold, Word_Size)
    Wait_Count = 0
    while not os.path.exists(Output_File):
        time.sleep(0.1)
        Wait_Count += 1
        if Wait_Count > 30:
            Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
            sys.exit()
    Clusters = Parse_Clusters(Clusters, Output_File)
    return Clusters, Output_File


def Parse_Fasta(fasta):
    """Converts Fasta file into a dictionary which uses fasta descriptions as keys"""
    Parsed = {}
    with open(fasta, "r") as Input:
        CurrentDesc = ""
        for line in Input:
            if line.startswith(">"):
                if line.count("|") > 0:
                    CurrentDesc = line.split("|")[-1].split("[")[0].strip()
                else:
                    CurrentDesc = line.rstrip()
                Parsed[CurrentDesc] = ""
            else:
                Parsed[CurrentDesc] += line.rstrip()
    return Parsed


def Parse_HMM(HMM):
    """Reads profile names from HMM profile files"""
    Names = []
    with open(HMM, "r") as Input:
        for line in Input:
            if line.startswith("NAME  "):
                Names.append(line[6:].strip())
    return Names


def Create_Prot_Files(Protein_Dir, Concat_Gen_Dict, Clusters):
    """Generates fasta file from fasta dictionary"""
    for index, item in enumerate(Clusters, 1):
        if not os.path.exists(os.path.join(Protein_Dir, "Protein_Family_" + str(index) + ".fasta")):
            with open(os.path.join(Protein_Dir, "Protein_Family_" + str(index) + ".fasta"), "w") as ProtFile:
                for sub_item in item:
                    ProtFile.write(sub_item + "\n")
                    ProtFile.write(Concat_Gen_Dict[sub_item] + "\n")


def Remove_Singletons(Nested_List):
    """Returns an output protein family list containing only multi-protein families"""
    copy_list = [x[:] for x in Nested_List]
    removal_list = []
    for item in copy_list:
        if len(item) < 2:
            removal_list.append(item)
    for item in removal_list:
        copy_list.remove(item)
    return copy_list


def Remove_Changed_Clusters(Directory, Old_Clusters, New_Clusters):
    """Compares two lists of proteins to identify differences.
    Deletes protein files from the data directory that differ between the two lists."""
    Rename_Dict = {}
    for ProtFam in Old_Clusters:
        Path = os.path.join(Directory, "Protein_Family_" + str(Old_Clusters.index(ProtFam) + 1))
        if ProtFam not in New_Clusters:
            if os.path.exists(Path + ".fasta"):
                os.remove(Path + ".fasta")
            if os.path.exists(Path):
                os.remove(Path)
        elif ProtFam in New_Clusters:
            if not Old_Clusters.index(ProtFam) == New_Clusters.index(ProtFam):
                if os.path.exists(Path + ".fasta"):
                    Rename_Dict[Old_Clusters.index(ProtFam)] = (os.path.basename(Path) + ".fasta", "Protein_Family_" + str(New_Clusters.index(ProtFam) + 1) + ".fasta")
                if os.path.exists(Path):
                    Rename_Dict[str(Old_Clusters.index(ProtFam))] = (os.path.basename(Path), "Protein_Family_" + str(New_Clusters.index(ProtFam) + 1))
    for key in sorted(Rename_Dict.keys()):
        if os.path.exists(os.path.join(Directory, Rename_Dict[key][1])):
            os.remove(os.path.join(Directory, Rename_Dict[key][1]))
        shutil.move(os.path.join(Directory, Rename_Dict[key][0]), os.path.join(Directory, Rename_Dict[key][1]))


def Run_HMMBUILD(Alignment_Dir, Protein_Alignment, HMMBUILD, Output):
    """Runs HMMBuild subprocess"""
    cmd = [HMMBUILD, "--cpu", "1", os.path.join(Alignment_Dir, os.path.splitext(Protein_Alignment)[0]),
           os.path.join(Alignment_Dir, Protein_Alignment)]
    HMMBUILD_Output = ""
    try:
        HMMBUILD_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as call_err:
        HMMBUILD_Output += "\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output)
        Output.error("\n\nHMMBUILD Error. Check Logs.\n")
    except OSError as os_err:
        HMMBUILD_Output += "\n\n" + str(os_err.strerror)
        Output.error("\n\nHMMBUILD Error. Check Logs.\n")
    return HMMBUILD_Output


def Run_HMMSEARCH(Alignment_Dir, Protein_Alignment, index, Genome_Dir, Representative_File, HMMSEARCH, Output):
    """Runs HMMSearch subprocess"""
    cmd = [HMMSEARCH, "-E", "1e-5", "--incE", "1e-20", "--noali", "--cpu", "1", os.path.join(Alignment_Dir, Protein_Alignment),
           os.path.join(Genome_Dir, Representative_File)]
    HMMSEARCH_Output = ""
    try:
        HMMSEARCH_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as call_err:
        HMMSEARCH_Output += "\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output)
        Output.error("\n\nHMMSEARCH Error. Check Logs.\n")
    except OSError as os_err:
        HMMSEARCH_Output += "\n\n" + str(os_err.strerror)
        Output.error("\n\nHMMSEARCH Error. Check Logs.\n")
    Return_Tuple = (index, HMMSEARCH_Output)
    return Return_Tuple


def HMM_Clustering(Genome_Dir, Alignment_Dir, HMMBUILD, HMMSEARCH, CDHIT_Clusters, CDHIT_Clusters_Only_Multi, Log_Dir, Representative_File, Threads, Output):
    """Runs HMMBuild and HMMSearch on all proteins in input directory and parses output"""
    Accepted_Filetypes = [".fasta", ".fsa", ".faa", ".fas"]
    Alignments = []
    Profiles = []
    HMM_Clusters = [x[:] for x in CDHIT_Clusters]
    HMM_Cluster_Reps = [x[0] for x in HMM_Clusters]
    HMMBuildPool = multiprocessing.Pool(Threads)
    HMMBUILDQueue = []
    for Protein_Alignment in os.listdir(Alignment_Dir):
        if os.path.splitext(Protein_Alignment)[1].lower() in Accepted_Filetypes and not os.path.exists(
                os.path.splitext(Protein_Alignment)[0]):
            HMMBUILDQueue.append(
                HMMBuildPool.apply_async(Run_HMMBUILD, args=(Alignment_Dir, Protein_Alignment, HMMBUILD, Output)))
            Alignments.append(Protein_Alignment)
            Profiles.append(os.path.splitext(Protein_Alignment)[0])
    HMMBuildOut = [result.get() for result in HMMBUILDQueue]
    HMMBuildPool.close()
    HMMBuildPool.join()
    with open(os.path.join(Log_Dir, "HMMBuild.txt"), "w") as HMMer_Log:
        HMMer_Log.writelines(HMMBuildOut)
    Alignments.sort()
    Profiles.sort()
    Output.write("Performing profile HMM searches using HMMSearch...\n")
    HMMSearchPool = multiprocessing.Pool(Threads)
    HMMSEARCHQueue = []
    for index, Protein_Alignment in enumerate(Profiles):
        HMMSEARCHQueue.append(HMMSearchPool.apply_async(Run_HMMSEARCH, args=(Alignment_Dir, Protein_Alignment, index, Genome_Dir, Representative_File, HMMSEARCH, Output)))
    HMMSearchOut = [result.get() for result in HMMSEARCHQueue]
    HMMSearchPool.close()
    HMMSearchPool.join()
    with open(os.path.join(Log_Dir, "HMMSearch.txt"), "w") as HMMer_Log:
        HMMer_Log.writelines([Out[1] for Out in HMMSearchOut])
    Exit_Terms = ["------ inclusion threshold ------", "[No hits detected that satisfy reporting thresholds]",
                  "Domain annotation for each sequence:"]
    Genome_Num = len([x for x in os.listdir(Genome_Dir) if os.path.splitext(x)[1].lower() in Accepted_Filetypes])
    for output in HMMSearchOut:
        Primary_Cluster = CDHIT_Clusters_Only_Multi[int(Profiles[output[0]].split("_")[-1]) - 1][0]
        Primary_Cluster_Whole = CDHIT_Clusters_Only_Multi[int(Profiles[output[0]].split("_")[-1]) - 1]
        try:
            Primary_Index = HMM_Clusters.index(Primary_Cluster_Whole)
        except (IndexError, ValueError):
            continue
        Primary_Genomes = [ID[:5] for ID in HMM_Clusters[Primary_Index]]
        for ID in Primary_Genomes:
            if Primary_Genomes.count(ID) > 1:
                Primary_Genomes.remove(ID)
        Primary_Len = len(Primary_Genomes)
        if Primary_Len < Genome_Num:
            for line in output[1].splitlines():
                if Exit_Terms[0] in line or Exit_Terms[1] in line or Exit_Terms[2] in line:
                    break
                elif line[60:90].strip().isdigit():
                    Secondary_Cluster = ">" + line[60:90].strip()
                    if Primary_Cluster != Secondary_Cluster:
                        try:
                            Secondary_Index = HMM_Cluster_Reps.index(Secondary_Cluster)
                            Secondary_Genomes = [ID[:5] for ID in HMM_Clusters[Secondary_Index]]
                            for ID in Secondary_Genomes:
                                if Secondary_Genomes.count(ID) > 1:
                                    Secondary_Genomes.remove(ID)
                            Secondary_Len = len(Secondary_Genomes)
                            if Primary_Len + Secondary_Len < Genome_Num:
                                HMM_Clusters[Primary_Index].extend(HMM_Clusters[Secondary_Index])
                                del HMM_Clusters[Secondary_Index]
                                del HMM_Cluster_Reps[Secondary_Index]
                            break
                        except (IndexError, ValueError):
                            continue
    return HMM_Clusters


def HMM_Representative_File(Genome_Dir, Concat_Gen_Dict, HMM_Clusters):
    """Generates a file with a single protein representing each protein family for accelerated searching"""
    with open(os.path.join(Genome_Dir, "HMM_Clusters"), "w") as Rep_File:
        for cluster in HMM_Clusters:
            Rep_File.write(cluster[0] + "\n")
            Rep_File.write(Concat_Gen_Dict[cluster[0]] + "\n")


def ClusterCoreProts(Genome_Dir, Protein_Dir, Alignment_Dir, CDHIT, ClustalOmega, HMMBUILD, HMMSEARCH, Log_Dir, Threads, Output, Fast_Cluster):
    """Automated core protein family identification"""
    ConcatenatedGenomeFile = os.path.join(Genome_Dir, "ConcatenatedGenomeFiles")
    Output.write("Iteratively clustering proteins...\n")
    CDHIT_Clusters, Cluster_File = Run_CDHIT(Genome_Dir, CDHIT, Log_Dir, ConcatenatedGenomeFile, Output)
    Output.write("Iterative clustering step complete.\n")
    Concat_Gen_Dict = Parse_Fasta(ConcatenatedGenomeFile)
    CDHIT_Clusters_Only_Multi = Remove_Singletons(CDHIT_Clusters)
    Create_Prot_Files(Protein_Dir, Concat_Gen_Dict, CDHIT_Clusters_Only_Multi)
    if Fast_Cluster:
        return CDHIT_Clusters_Only_Multi, Concat_Gen_Dict
    else:
        Output.write("Aligning " + str(len(os.listdir(Protein_Dir))) + " protein families...\n")
        ParallelAlignment(Protein_Dir, Alignment_Dir, ClustalOmega, Log_Dir, Threads, True, [], Output)
        Output.write("Protein family alignment step complete.\n")
        Output.write("Building HMM profiles for protein families...\n")
        HMM_Clusters_1 = HMM_Clustering(Genome_Dir, Alignment_Dir, HMMBUILD, HMMSEARCH, CDHIT_Clusters,
                                        CDHIT_Clusters_Only_Multi, Log_Dir, os.path.basename(Cluster_File), Threads,
                                        Output)
        HMM_Clusters_Only_Multi_1 = Remove_Singletons(HMM_Clusters_1)
        Remove_Changed_Clusters(Protein_Dir, CDHIT_Clusters_Only_Multi, HMM_Clusters_Only_Multi_1)
        Remove_Changed_Clusters(Alignment_Dir, CDHIT_Clusters_Only_Multi, HMM_Clusters_Only_Multi_1)
        Create_Prot_Files(Protein_Dir, Concat_Gen_Dict, HMM_Clusters_Only_Multi_1)
        Output.write("Aligning additionally identified protein families...\n")
        ParallelAlignment(Protein_Dir, Alignment_Dir, ClustalOmega, Log_Dir, Threads, True, [], Output)
        HMM_Representative_File(Genome_Dir, Concat_Gen_Dict, HMM_Clusters_1)
        Output.write("Building HMM profiles for additionally identified protein families...\n")
        HMM_Clusters_2 = HMM_Clustering(Genome_Dir, Alignment_Dir, HMMBUILD, HMMSEARCH, HMM_Clusters_1,
                                        HMM_Clusters_Only_Multi_1, Log_Dir, "HMM_Clusters", Threads, Output)
        HMM_Clusters_Only_Multi_2 = Remove_Singletons(HMM_Clusters_2)
        Remove_Changed_Clusters(Protein_Dir, HMM_Clusters_Only_Multi_1, HMM_Clusters_Only_Multi_2)
        Remove_Changed_Clusters(Alignment_Dir, HMM_Clusters_Only_Multi_1, HMM_Clusters_Only_Multi_2)
        Create_Prot_Files(Protein_Dir, Concat_Gen_Dict, HMM_Clusters_Only_Multi_2)
        return HMM_Clusters_Only_Multi_2, Concat_Gen_Dict


def Run_JACKHMMER(Input_Fasta_Dict, Target_Proteins, ConcatenatedGenomeFile, JACKHMMER, Log_Dir, Output):
    """Runs Jackhammer using input fasta file and parses output"""
    Protein_Fams = {}
    JACKHMMER_Output = ""
    cmd = [JACKHMMER, "-E", "1e-5", "--incE", "1e-20", "--noali", Target_Proteins, ConcatenatedGenomeFile]
    try:
        JACKHMMER_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as call_err:
        with open(os.path.join(Log_Dir, "JACKHMMER.txt"), "a") as HMMer_Log:
            HMMer_Log.write(JACKHMMER_Output)
            HMMer_Log.write("\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output))
        Output.error("\n\nJACKHMMER Error. Check Logs.\n")
    except OSError as os_err:
        with open(os.path.join(Log_Dir, "JACKHMMER.txt"), "a") as HMMer_Log:
            HMMer_Log.write(JACKHMMER_Output)
            HMMer_Log.write("\n\n" + str(os_err.strerror))
        Output.error("\n\nJACKHMMER Error. Check Logs.\n")
    with open(os.path.join(Log_Dir, "JACKHMMER.txt"), "a") as HMMer_Log:
        HMMer_Log.write(JACKHMMER_Output)
    Exit_Terms = ["------ inclusion threshold ------", "[No hits detected that satisfy reporting thresholds]"]
    Prot_Found = False
    Prot_Name = ""
    for line in JACKHMMER_Output.splitlines():
        if "Query:" in line or "Description:" in line:
            Prot_Found = True
            for key in Input_Fasta_Dict.keys():
                if key.lstrip(">") in line:
                    Prot_Name = Replace_Chars(key.lstrip(">"))
                    if Prot_Name in Protein_Fams.keys() and Protein_Fams[Prot_Name] == []:
                        count = 2
                        for x in range(10000):
                            if Prot_Name in Protein_Fams.keys() and Prot_Name + "_" + str(
                                    count) not in Protein_Fams.keys():
                                Prot_Name += "_" + str(count)
                                break
                            else:
                                count += 1
                    if Prot_Name not in Protein_Fams.keys():
                        Protein_Fams[Prot_Name] = []
                    break
        elif Exit_Terms[0] in line or Exit_Terms[1] in line:
            Prot_Found = False
        elif "@@ Round:" in line:
            Prot_Found = True
        elif line[60:90].strip().isdigit() and Prot_Found:
            for ID in Protein_Fams[Prot_Name]:
                if line[60:].startswith(ID[1:5]):
                    break
            else:
                Protein_Fams[Prot_Name].append(">" + line[60:90].strip())
    return Protein_Fams


def Find_Proteins(Genome_Dir, Protein_Dir, GLIMPSe_Output_Dir, Genome_Dictionary, Target_Proteins, JACKHMMER, PAMatrix, Log_Dir, Output):
    """Protein family identification based on input fasta files"""
    ConcatenatedGenomeFile = os.path.join(Genome_Dir, "ConcatenatedGenomeFiles")
    Concat_Gen_Dict = Parse_Fasta(ConcatenatedGenomeFile)
    Input_Fasta_Dict = Parse_Fasta(Target_Proteins)
    Output.write("Performing JACKHMMer serches for protein families...\n")
    Protein_Fams = Run_JACKHMMER(Input_Fasta_Dict, Target_Proteins, ConcatenatedGenomeFile, JACKHMMER, Log_Dir, Output)
    Create_Named_Prot_Files(Protein_Dir, Concat_Gen_Dict, Protein_Fams, Output)
    if PAMatrix:
        Output.write("Producing Presence Absence Matrix...\n")
        Build_PAMatrix(Protein_Fams, GLIMPSe_Output_Dir, Genome_Dictionary, Output)
    return Protein_Fams


def Marker_HMMSEARCH(Marker_Names, Marker_File, ConcatenatedGenomeFile, HMMSEARCH, Log_Dir, Output):
    """Runs HMMSearch using PhyEco Markers and parses output"""
    Protein_Fams = {}
    HMMSEARCH_Output = ""
    cmd = [HMMSEARCH, "-E", "1e-5", "--incE", "1e-20", "--noali", Marker_File, ConcatenatedGenomeFile]
    try:
        HMMSEARCH_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as call_err:
        with open(os.path.join(Log_Dir, "HMMSEARCH.txt"), "a") as HMMer_Log:
            HMMer_Log.write(HMMSEARCH_Output)
            HMMer_Log.write("\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output))
        Output.error("\n\nHMMSEARCH Error. Check Logs.\n")
    except OSError as os_err:
        with open(os.path.join(Log_Dir, "HMMSEARCH.txt"), "a") as HMMer_Log:
            HMMer_Log.write(HMMSEARCH_Output)
            HMMer_Log.write("\n\n" + str(os_err.strerror))
        Output.error("\n\nHMMSEARCH Error. Check Logs.\n")
    with open(os.path.join(Log_Dir, "HMMSEARCH.txt"), "a") as HMMer_Log:
        HMMer_Log.write(HMMSEARCH_Output)
    Exit_Terms = ["------ inclusion threshold ------", "[No hits detected that satisfy reporting thresholds]", "Domain annotation for each sequence:"]
    Prot_Found = False
    Prot_Name = ""
    for line in HMMSEARCH_Output.splitlines():
        if "Query:" in line:
            Prot_Found = True
            for Name in Marker_Names:
                if re.search(" \w.+? ", line).group(0).strip() in Name:
                    Prot_Name = Replace_Chars(Name)
                    if Prot_Name not in Protein_Fams.keys():
                        Protein_Fams[Prot_Name] = []
                    break
        elif Exit_Terms[0] in line or Exit_Terms[1] in line:
            Prot_Found = False
        elif line[60:90].strip().isdigit() and Prot_Found:
            for ID in Protein_Fams[Prot_Name]:
                if line[60:].startswith(ID[1:5]):
                    break
            else:
                Protein_Fams[Prot_Name].append(">" + line[60:90].strip())
    return Protein_Fams


def Find_Marker_Proteins(Genome_Dir, Protein_Dir, GLIMPSe_Output_Dir, Genome_Dictionary, HMMSEARCH, PAMatrix, Log_Dir, Marker_Dir, Marker_Proteins, Output):
    """Protein family identification using PhyEco marker sets"""
    Marker_File = os.path.join(Marker_Dir, Marker_Proteins + ".hmm")
    if not os.path.exists(Marker_File):
        Output.error("Cannot locate selected PhyEco marker set HMM profiles.\n")
        sys.exit()
    ConcatenatedGenomeFile = os.path.join(Genome_Dir, "ConcatenatedGenomeFiles")
    Concat_Gen_Dict = Parse_Fasta(ConcatenatedGenomeFile)
    Marker_Names = Parse_HMM(Marker_File)
    Output.write("Performing HMMSearch serches for protein families...\n")
    Protein_Fams = Marker_HMMSEARCH(Marker_Names, Marker_File, ConcatenatedGenomeFile, HMMSEARCH, Log_Dir, Output)
    Create_Named_Prot_Files(Protein_Dir, Concat_Gen_Dict, Protein_Fams, Output)
    if PAMatrix:
        Output.write("Producing Presence Absence Matrix...\n")
        Build_PAMatrix(Protein_Fams, GLIMPSe_Output_Dir, Genome_Dictionary, Output)
    return Protein_Fams


def Create_Named_Prot_Files(Protein_Dir, Concat_Gen_Dict, Protein_Fams, Output):
    """Creates output protein files for named protein families"""
    for item in Protein_Fams.keys():
        if len(Protein_Fams[item]) > 1:
            with open(os.path.join(Protein_Dir, item + ".fasta"), "w") as ProtFile:
                for sub_item in Protein_Fams[item]:
                    ProtFile.write(sub_item + "\n")
                    ProtFile.write(Concat_Gen_Dict[sub_item] + "\n")
        else:
            Output.error(
                "Less than 2 copies of " + item.lstrip(">") + " were identified in the input genomes. " + item.lstrip(
                    ">") + " will not be included in further analysis.\n")


def Build_PAMatrix(HMM_Clusters, GLIMPSe_Output_Dir, Genome_Dictionary, Output):
    """Generates PA Matrix"""
    ProtMatrix = {}
    IDs = sorted(Genome_Dictionary.keys())
    if isinstance(HMM_Clusters, list):
        Cluster_Dict = {}
        for num, item in enumerate(HMM_Clusters, 1):
            Cluster_Dict["Protein_Family_" + str(num)] = item
    elif isinstance(HMM_Clusters, dict):
        Cluster_Dict = HMM_Clusters.copy()
    else:
        Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
        sys.exit()
    with open(os.path.join(GLIMPSe_Output_Dir, "PA_Matrix.tsv"), "w") as PAMatrix:
        OutHeader = ""
        for ID in IDs:
            ProtMatrix[ID] = ""
            OutHeader = OutHeader + "\t" + Genome_Dictionary[ID]
        OutRows = []
        for Protein_Fam in Cluster_Dict.keys():
            OutRow = "\n" + Protein_Fam
            ProtDist = 0
            for ID in IDs:
                for Protein in Cluster_Dict[Protein_Fam]:
                    if Protein.startswith(">" + str(ID)):
                        OutRow = OutRow + "\t" + "1"
                        ProtDist += 1
                        break
                else:
                    OutRow = OutRow + "\t" + "0"
            else:
                OutRows.append((ProtDist, OutRow))
        OutRows.sort()
        OutRows.reverse()
        PAMatrix.write(OutHeader)
        for item in OutRows:
            for index, binary in enumerate(item[1].split("\t")):
                if index > 0:
                    ProtMatrix[IDs[index - 1]] = ProtMatrix[IDs[index - 1]] + binary
            if item[0] > 1:
                PAMatrix.write(item[1])
    return ProtMatrix, OutHeader


def Build_POCP(GLIMPSe_Output_Dir, Genome_Dictionary, ProtMatrix, OutHeader):
    """Generates Percent of Shared Proteins Matrix.
    Uses output from PA Matrix."""
    IDs = sorted(Genome_Dictionary.keys())
    with open(os.path.join(GLIMPSe_Output_Dir, "PoCP_Matrix.tsv"), "w") as PoCPMatrix:
        OutRows = []
        for ID_1 in IDs:
            OutRow = "\n" + Genome_Dictionary[ID_1]
            for ID_2 in IDs:
                Match = 0
                length_1 = ProtMatrix[ID_1].count("1")
                length_2 = ProtMatrix[ID_2].count("1")
                Average_lenght = (float(length_1) + float(length_2)) / 2
                for index, binary in enumerate(ProtMatrix[ID_1]):
                    if binary == "1" and ProtMatrix[ID_2][index] == "1":
                        Match += 1
                OutRow = OutRow + "\t" + str(float(Match) / float(Average_lenght))
            else:
                OutRows.append(OutRow)
        else:
            PoCPMatrix.write(OutHeader)
            PoCPMatrix.writelines(OutRows)


def Protein_Family_Filter(Protein_Dir, Alignment_Dir, Concat_Gen_Dict, HMM_Clusters, GLIMPSe_Output_Dir,
                          Genome_Dictionary, PAMatrix, POCP, Single_Copy, Output):
    """Generates protein matricies and filters paralogs/single copy"""
    Accepted_Clusters = [x[:] for x in HMM_Clusters]
    ProtMatrix, OutHeader = {}, ""
    if PAMatrix:
        Output.write("Producing Presence Absence Matrix...\n")
        ProtMatrix, OutHeader = Build_PAMatrix(HMM_Clusters, GLIMPSe_Output_Dir, Genome_Dictionary, Output)
    elif not PAMatrix and POCP:
        ProtMatrix, OutHeader = Build_PAMatrix(HMM_Clusters, GLIMPSe_Output_Dir, Genome_Dictionary, Output)
        os.remove(os.path.join(GLIMPSe_Output_Dir, "PA_Matrix.tsv"))
    if POCP:
        Output.write("Producing Percentage of Conserved Protein Family Matrix...\n")
        Build_POCP(GLIMPSe_Output_Dir, Genome_Dictionary, ProtMatrix, OutHeader)
    if Single_Copy:
        Removed_Clusters = []
        for item in Accepted_Clusters:
            IDs_Only = [x[1:5] for x in item]
            for ID in item:
                if IDs_Only.count(ID[1:5]) > 1:
                    Removed_Clusters.append(item)
                    break
        for item in Removed_Clusters:
            Accepted_Clusters[Accepted_Clusters.index(item)] = []
    elif not Single_Copy:
        for item in Accepted_Clusters:
            Removed_IDs = []
            Seen_IDs = []
            for ID in item:
                if ID[1:5] not in Seen_IDs:
                    Seen_IDs.append(ID[1:5])
                else:
                    Removed_IDs.append(ID)
            for ID in Removed_IDs:
                item.remove(ID)
        Accepted_Clusters = Remove_Singletons(Accepted_Clusters)
        Remove_Changed_Clusters(Protein_Dir, HMM_Clusters, Accepted_Clusters)
        Remove_Changed_Clusters(Alignment_Dir, HMM_Clusters, Accepted_Clusters)
        Create_Prot_Files(Protein_Dir, Concat_Gen_Dict, Accepted_Clusters)
    return Accepted_Clusters


def Determine_Protein_Distribution(Protein_Distribution, Protein_Clusters, Genome_Dictionary, Output):
    """Assesses the distribution of protein families among genomes"""
    Accepted_Proteins = []
    Genome_Number = len(Genome_Dictionary.keys())
    if isinstance(Protein_Clusters, list):
        Cluster_Dict = {}
        for num, item in enumerate(Protein_Clusters, 1):
            Cluster_Dict["Protein_Family_" + str(num)] = item
    elif isinstance(Protein_Clusters, dict):
        Cluster_Dict = Protein_Clusters.copy()
    else:
        Output.error("\n\nGLIMPSe pipeline Error. Check Logs.\n")
        sys.exit()
    Protein_Names = Cluster_Dict.keys()
    for Protein_Name in Protein_Names:
        if len(Cluster_Dict[Protein_Name]) >= Genome_Number * Protein_Distribution:
            Accepted_Proteins.append(str(Protein_Name) + ".fasta")
    if len(Accepted_Proteins) < 1:
        Output.error("\nNo identified proteins meet protein distribution threshold\n")
        sys.exit()
    else:
        Prot_Num = len(Accepted_Proteins)
        Output.write("GLIMPS pipeline has identified " + str(Prot_Num) + " protein(s) which meet selected criteria for core genome.\n")
    return Accepted_Proteins


def Calculate_AIs(Alignment_Dir, Alignment, AIs):
    """Calculates Aamino Acid Identity for a given amino acid alignment"""
    AI_Dict = {}
    Percent_ID = 0.0
    Aligned_Seqs = Parse_Fasta(os.path.join(Alignment_Dir, Alignment))
    Sorted_Seqs = sorted(Aligned_Seqs.keys())
    for index_1, key_1 in enumerate(Sorted_Seqs):
        for index_2 in range(index_1 + 1,len(Sorted_Seqs)):
            key_2 = Sorted_Seqs[index_2]
            for key in AIs.keys():
                if key_1[1:5] in key and key_2[1:5] in key:
                    identical = 0
                    length_1 = 0
                    length_2 = 0
                    for index, char in enumerate(Aligned_Seqs[key_1]):
                        if char == Aligned_Seqs[key_2][index]:
                            if char != "-":
                                identical += 1
                        if char != "-":
                            length_1 += 1
                        if Aligned_Seqs[key_2][index] != "-":
                            length_2 += 1
                    Percent_ID = float(identical) / float(min(length_1, length_2))
                    AI_Dict[key] = Percent_ID
                    break
    Return_List = [Alignment, AI_Dict]
    return Return_List


# Alignment and trimming of proteins is handled by the following functions
def Calculate_AAI(Alignment_Dir, Genome_Dictionary, Log_Dir, GLIMPSe_Output_Dir, Threads):
    """Generates AAI matrix"""
    IDs = sorted(Genome_Dictionary.keys())
    AIs = {}
    for ID_1 in IDs:
        for ID_2 in IDs:
            if IDs.index(ID_2) > IDs.index(ID_1):
                AIs[(ID_1, ID_2)] = []
    Accepted_Filetypes = [".fasta", ".fsa", ".faa", ".fas"]
    AIPool = multiprocessing.Pool(Threads)
    AIQueue = []
    for Alignment in os.listdir(Alignment_Dir):
        if os.path.splitext(Alignment)[1].lower() in Accepted_Filetypes:
            AIQueue.append(AIPool.apply_async(func=Calculate_AIs, args=(Alignment_Dir, Alignment, AIs)))
    AIOut = [result.get() for result in AIQueue]
    AIPool.close()
    AIPool.join()
    for result in AIOut:
        Align = result[0]
        for key in result[1].keys():
            AIs[key].append((Align, result[1][key]))
    AAIs = {}
    for key in sorted(AIs.keys()):
        AI_List = []
        total = float(0)
        if len(AIs[key]) < 1:
            AAIs[key] = "N/A"
        else:
            for item in AIs[key]:
                AI_List.append(item[1])
                total += item[1]
            Average = total / float(len(AI_List))
            AAIs[key] = Average
    AAI_Matrix = [""]
    for ID in IDs:
        AAI_Matrix[0] += "\t" + Genome_Dictionary[ID]
        AAI_Matrix.insert(IDs.index(ID) + 1, "\n" + Genome_Dictionary[ID])
        for key_1 in sorted(AAIs.keys()):
            if ID == key_1[1]:
                AAI_Matrix[IDs.index(ID) + 1] += "\t" + str(AAIs[key_1])
        else:
            AAI_Matrix[IDs.index(ID) + 1] += "\t" + "1"
            for key_2 in sorted(AAIs.keys()):
                if ID == key_2[0]:
                    AAI_Matrix[IDs.index(ID) + 1] += "\t" + str(AAIs[key_2])
    with open(os.path.join(Log_Dir, "Amino Acid Identities.tsv"), "w") as Log:
        Log.write("Genome 1\tGenome 2\tAmino Acid Identities\n")
        for key in sorted(AIs.keys()):
            Log.write(Genome_Dictionary[key[0]] + "\t" + Genome_Dictionary[key[1]])
            for item in AIs[key]:
                Log.write("\t" + os.path.splitext(item[0])[0] + ": " + str(item[1]))
            else:
                Log.write("\n")
    with open(os.path.join(GLIMPSe_Output_Dir, "AAI_Matrix.tsv"), "w") as Out:
        Out.writelines(AAI_Matrix)


def Align_Proteins(Protein_Dir, Protein_File, Alignment_Dir, ClustalOmega, Count, First_Run, Output):
    """Runs ClustalOmega"""
    ClustalOmega_Output = ""
    Accepted_Filetypes = [".fasta", ".fsa", ".faa", ".fas"]
    if os.path.splitext(Protein_File)[1].lower() in Accepted_Filetypes:
        cmd = [ClustalOmega, "-v", "--force", "--threads=1", "-i", os.path.join(Protein_Dir, Protein_File), "-o",
               os.path.join(Alignment_Dir, Protein_File)]
        try:
            ClustalOmega_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as call_err:
            ClustalOmega_Output += "\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output)
            Output.error("\n\nClustal Omega Error. Check Logs.\n")
            return ClustalOmega_Output
        except OSError as os_err:
            ClustalOmega_Output += "\n\n" + str(os_err.strerror)
            Output.error("\n\nClustal Omega Error. Check Logs.\n")
            return ClustalOmega_Output
    if not First_Run:
        if Count % 25 == 0:
            Output.write("Alignment of %d Protein families complete...\n" % Count)
    return ClustalOmega_Output


def ParallelAlignment(Protein_Dir, Alignment_Dir, ClustalOmega, Log_Dir, Threads, First_Run, Accepted_Proteins, Output):
    """Manages parallel instances of ClustalOmega"""
    AlignPool = multiprocessing.Pool(Threads)
    ClustalQueue = []
    if len(Accepted_Proteins) > 0:
        for Count, Accepted_Protein in enumerate(Accepted_Proteins, 1):
            ClustalQueue.append(AlignPool.apply_async(func=Align_Proteins, args=(Protein_Dir, Accepted_Protein, Alignment_Dir, ClustalOmega, Count, First_Run, Output)))
        ClustalOut = [result.get() for result in ClustalQueue]
        AlignPool.close()
        AlignPool.join()
        with open(os.path.join(Log_Dir, "ClustalOmega.txt"), "a") as Clustal_Log:
            Clustal_Log.writelines(ClustalOut)
    else:
        Accepted_Filetypes = [".fasta", ".fsa", ".faa", ".fas"]
        Accepted_Proteins_Files = []
        for Protein_File in os.listdir(Protein_Dir):
            if os.path.splitext(Protein_File)[1].lower() in Accepted_Filetypes and not os.path.exists(os.path.join(Alignment_Dir, Protein_File)):
                Accepted_Proteins_Files.append(Protein_File)
        for Count, Accepted_Proteins_File in enumerate(Accepted_Proteins_Files, 1):
            ClustalQueue.append(AlignPool.apply_async(func=Align_Proteins, args=(Protein_Dir, Accepted_Proteins_File, Alignment_Dir, ClustalOmega, Count, First_Run, Output)))
        ClustalOut = [result.get() for result in ClustalQueue]
        AlignPool.close()
        AlignPool.join()
        with open(os.path.join(Log_Dir, "ClustalOmega.txt"), "a") as Clustal_Log:
            Clustal_Log.writelines(ClustalOut)


# Not used. Based on the weighted TrimAl algorithm described in Chang, Di Tommaso, & Notredame (2014).
def Create_Weighted_Alignments(Alignment_Dir, Accepted_Proteins, TrimAl, Output):
    """Generates weighted alignments"""
    Alignment_Length = 0
    for Alignment_File in Accepted_Proteins:
        try:
            TrimAl_Output = subprocess.check_output([
                TrimAl,
                "-in",
                os.path.join(Alignment_Dir, Alignment_File),
                "-sgc"
            ])
        except subprocess.CalledProcessError as call_err:
            Output.error("\n\nTrimAl Error.\n")
            Output.error("\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output))
            sys.exit()
        except OSError as os_err:
            Output.error("\n\nTrimAl Error.\n")
            Output.error("\nError: " + str(os_err.strerror))
            sys.exit()
        seq_score = []
        for line in TrimAl_Output.splitlines():
            score_str = re.search("(\t\t).*\t(.*)", line)
            if score_str is not None:
                score_flt = float(score_str.group(2).strip())
                column_multiplier = round(10 * 0.9 * score_flt + 1)
                seq_score.append(column_multiplier)
        with open(os.path.join(Alignment_Dir, Alignment_File), "r") as Unweighted_Alignment:
            All_Sequences_In_Alignment = {}
            Sequence_Name = ""
            for line in Unweighted_Alignment:
                if line.startswith(">"):
                    Sequence_Name = line
                    All_Sequences_In_Alignment[Sequence_Name] = ""
                else:
                    All_Sequences_In_Alignment[Sequence_Name] += line.strip()
            with open(os.path.join(Alignment_Dir, "Weighted_" + Alignment_File), "w") as Weighted_Alignment:
                First_Line = True
                for Sequence_Key in All_Sequences_In_Alignment.keys():
                    if First_Line:
                        Weighted_Alignment.write(Sequence_Key)
                        First_Line = False
                        Alignment_Length += len(All_Sequences_In_Alignment[Sequence_Key])
                    else:
                        Weighted_Alignment.write("\n" + Sequence_Key)
                    Weighted_Sequence = []
                    for index, character in enumerate(All_Sequences_In_Alignment[Sequence_Key]):
                        count = 0
                        while count < seq_score[index]:
                            Weighted_Sequence.append(character)
                            count += 1
                    else:
                        Weighted_Sequence = "".join(Weighted_Sequence)
                        Weighted_Alignment.write(Weighted_Sequence)
    return Alignment_Length

# TODO Add option to filter all completely conserved sites in alignment
def Create_Trimmed_Alignments(Alignment_Dir, Accepted_Proteins, TrimAl, Output):
    """Generates trimmed alignments"""
    for Alignment_File in Accepted_Proteins:
        try:
            subprocess.call([
                TrimAl,
                "-in",
                os.path.join(Alignment_Dir, Alignment_File),
                "-out",
                os.path.join(Alignment_Dir, "Trimmed_" + Alignment_File),
                "-automated1"
            ])
        except subprocess.CalledProcessError as call_err:
            Output.error("\n\nTrimAl Error.\n")
            Output.error("\n\n" + "Command: " + str(call_err.cmd) + "\nError: " + str(call_err.output))
            sys.exit()
        except OSError as os_err:
            Output.error("\n\nTrimAl Error.\n")
            Output.error("\nError: " + str(os_err.strerror))
            sys.exit()


def Insert_Filler_Sequences(Accepted_Proteins, Alignment_Dir, Genome_Dictionary):
    """Inserts blank sequences from organisms missing in an alignment"""
    for alignment in Accepted_Proteins:
        if os.path.exists(os.path.join(Alignment_Dir, "Weighted_" + alignment)):
            Alignment_File = "Weighted_" + alignment
        else:
            Alignment_File = "Trimmed_" + alignment
        align_length = 0
        align_string = []
        desc_found = False
        for Genome_ID in Genome_Dictionary:
            key_found = False
            with open(os.path.join(Alignment_Dir, Alignment_File), "a+") as Align:
                Align.seek(0, 0)
                for line in Align:
                    if not desc_found and align_length == 0 and line.startswith(">"):
                        desc_found = True
                    elif desc_found and not line.startswith(">"):
                        align_string.append(line.strip())
                    elif desc_found and align_length == 0 and line.startswith(">"):
                        align_string = "".join(align_string)
                        align_length = len(align_string)
                        desc_found = False
                    if line.startswith(">" + Genome_ID):
                        key_found = True
                else:
                    if not key_found:
                        filler = "".join([">", Genome_ID, "\n", "-" * align_length, "\n"])
                        Align.write(filler)


def Concatenate_Alignments(Accepted_Proteins, Alignment_Dir, Concatenated_Dir, Genome_Dictionary, GLIMPSe_Output_Dir, Output):
    """Concatenates all trimmed/weighted alignments"""
    Insert_Filler_Sequences(Accepted_Proteins, Alignment_Dir, Genome_Dictionary)
    Alignment_Length = 0
    Concatenated_Sequence = {}
    Genome_IDs = Genome_Dictionary.keys()
    Current_ID = ""
    for alignment in Accepted_Proteins:
        if os.path.exists(os.path.join(Alignment_Dir, "Weighted_" + alignment)):
            Alignment_File = "Weighted_" + alignment
        else:
            Alignment_File = "Trimmed_" + alignment
        with open(os.path.join(Alignment_Dir, Alignment_File), "r") as align:
            for line in align:
                if line[1:5] in Genome_IDs:
                    Current_ID = line[1:5]
                    if Current_ID not in Concatenated_Sequence:
                        Concatenated_Sequence[Current_ID] = []
                elif not line.startswith(">"):
                    Concatenated_Sequence[Current_ID].append(line.rstrip())
    for Genome_ID in Genome_IDs:
        Concatenated_Sequence[Genome_ID] = "".join(Concatenated_Sequence[Genome_ID])
        if re.search("[a-zA-Z]", Concatenated_Sequence[Genome_ID]) is None:
            Output.write("Homologs for " + Genome_Dictionary[Genome_ID] + " were not identified for any protein. " +
                         Genome_Dictionary[Genome_ID] + " will be excluded from further analysis.\n")
        else:
            with open(os.path.join(Concatenated_Dir, "Concatenated_Alignment.fasta"), "a") as AlignOut:
                AlignOut.write(">" + Genome_Dictionary[Genome_ID] + "\n" + Concatenated_Sequence[Genome_ID] + "\n")
        if Alignment_Length == 0:
            Alignment_Length = len(Concatenated_Sequence[Genome_ID])
    shutil.copy(os.path.join(Concatenated_Dir, "Concatenated_Alignment.fasta"),
                os.path.join(GLIMPSe_Output_Dir, "Concatenated_Alignment.fasta"))
    return Alignment_Length


# Phylogenetic tree construction is handled by the following functions
# TODO Add option to build individual phylogenetic tree for each protein
def Run_FastTree(Tree_Dir, Concatenated_Dir, FastTree, Log_Dir, Output):
    """Builds a phylogenetic tree using FastTree"""
    Input_Alignment = os.path.join(Concatenated_Dir, "Concatenated_Alignment.fasta")
    Output_Tree = os.path.join(Tree_Dir, "FastTree.nwk")
    FastTree_Output = ""
    cmd = [FastTree,
           "-spr",
           "6",
           "-mlacc",
           "3",
           "-slownni",
           "-slow",
           "-nosupport",
           "-out",
           Output_Tree,
           Input_Alignment
           ]
    try:
        FastTree_Output = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        if os.path.getsize(Output_Tree) == 0:
            raise subprocess.CalledProcessError(1, str(cmd), str(FastTree_Output))
        else:
            with open(os.path.join(Log_Dir, "FastTree.txt"), "a") as FastTree_Log:
                FastTree_Log.write(FastTree_Output)
    except (subprocess.CalledProcessError, OSError):
        with open(os.path.join(Log_Dir, "FastTree.txt"), "a") as FastTree_Log:
            FastTree_Log.write(FastTree_Output)
        Output.error("\n\nFastTree Error. Check Logs.\n")
        sys.exit()


def Run_FastTree_Full(Tree_Dir, Concatenated_Dir, FastTree, Log_Dir, Output, GLIMPSe_Output_Dir):
    """Builds a phylogenetic tree using FastTree and generates statistical branch support"""
    Input_Alignment = os.path.join(Concatenated_Dir, "Concatenated_Alignment.fasta")
    Output_Tree = os.path.join(Tree_Dir, "FastTree.nwk")
    FastTree_Output = ""
    cmd = [FastTree,
           "-spr",
           "6",
           "-mlacc",
           "3",
           "-slownni",
           "-slow",
           "-out",
           Output_Tree,
           Input_Alignment
           ]
    try:
        FastTree_Output = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        if os.path.getsize(Output_Tree) == 0:
            raise subprocess.CalledProcessError(1, str(cmd), str(FastTree_Output))
        else:
            with open(os.path.join(Log_Dir, "FastTree.txt"), "a") as FastTree_Log:
                FastTree_Log.write(FastTree_Output)
    except (subprocess.CalledProcessError, OSError):
        with open(os.path.join(Log_Dir, "FastTree.txt"), "a") as FastTree_Log:
            FastTree_Log.write(FastTree_Output)
        Output.error("\n\nFastTree Error. Check Logs.\n")
        sys.exit()
    try:
        import dendropy
    except ImportError:
        dendropy = None
        shutil.copy(os.path.join(Tree_Dir, "FastTree.nwk"),
                    os.path.join(GLIMPSe_Output_Dir, "Final_Tree.nwk"))
    if dendropy:
        with open(os.path.join(Tree_Dir, "FastTree.nwk"), "r") as In_Tree:
            Tree = dendropy.Tree.get_from_string(In_Tree.readline(), "newick")
        Tree.ladderize(ascending=False)
        with open(os.path.join(GLIMPSe_Output_Dir, "Final_Tree.nwk"), "w") as Out_Tree:
            Out_Tree.write(Tree.as_string(schema='newick'))


def Run_RAxML(Tree_Dir, Concatenated_Dir, RAxML, Threads, Log_Dir, GLIMPSe_Output_Dir, Output):
    """Optimizes FastTree phylogeny using RAxML"""
    Input_Alignment = os.path.join(Concatenated_Dir, "Concatenated_Alignment.fasta")
    Input_FastTree = os.path.join(Tree_Dir, "FastTree.nwk")
    Input_RAxML = os.path.join(Tree_Dir, "RAxML_result.ML")
    RAxML_Output = ""
    cmd = [RAxML,
            "-f",
            "d",
            "-F",
            "-T",
            str(Threads),
            "-m",
            "PROTCATLG",
            "-n",
            "ML",
            "-p",
            str(random.randrange(1, 100000)),
            "-t",
            Input_FastTree,
            "-s",
            Input_Alignment,
            "-w",
            Tree_Dir
           ]
    try:
        RAxML_Output = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        count = 0
        while not os.path.exists(Input_RAxML) and count < 15:
                count += 1
                time.sleep(1)
        if not os.path.exists(Input_RAxML):
            raise subprocess.CalledProcessError(1, str(cmd), str(RAxML_Output))
    except (subprocess.CalledProcessError, OSError):
        with open(os.path.join(Log_Dir, "RAxML.txt"), "w") as RAxML_Log:
            RAxML_Log.write(RAxML_Output)
        Output.error("\n\nRAxML Error. Check Logs.\n")
        sys.exit()
    else:
        cmd = [RAxML,
               "-f",
               "J",
               "-T",
               str(Threads),
               "-m",
               "PROTCATLG",
               "-n",
               "SH",
               "-p",
               str(random.randrange(1, 100000)),
               "-t",
               Input_RAxML,
               "-s",
               Input_Alignment,
               "-w",
               Tree_Dir
               ]
        SH_Output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        with open(os.path.join(Log_Dir, "RAxML.txt"), "w") as RAxML_Log:
            RAxML_Log.write(RAxML_Output)
            RAxML_Log.write(SH_Output)
    try:
        import dendropy
    except ImportError:
        dendropy = None
        shutil.copy(os.path.join(Tree_Dir, "RAxML_fastTreeSH_Support.SH"),
                    os.path.join(GLIMPSe_Output_Dir, "Final_Tree.nwk"))
    if dendropy:
        with open(os.path.join(Tree_Dir, "RAxML_fastTreeSH_Support.SH"), "r") as In_Tree:
            Tree = dendropy.Tree.get_from_string(In_Tree.readline(), "newick")
        Tree.ladderize(ascending=False)
        with open(os.path.join(GLIMPSe_Output_Dir, "Final_Tree.nwk"), "w") as Out_Tree:
            Out_Tree.write(Tree.as_string(schema='newick'))


def GenomeID2GenomeName(Protein_Dir, Alignment_Dir, Genome_Dictionary):
    """Renames Sequences in Alignment Files"""
    for Protein_File in os.listdir(Protein_Dir):
        with open(os.path.join(Protein_Dir, Protein_File), "r") as ProtIn:
            Prot_Temp = ProtIn.readlines()
            for key in Genome_Dictionary:
                for item in Prot_Temp:
                    if key in item[1:5]:
                        Prot_Temp[Prot_Temp.index(item)] = ">" + Genome_Dictionary[key] + "\n"
        with open(os.path.join(Protein_Dir, Protein_File), "w") as ProtOut:
            ProtOut.writelines(Prot_Temp)
    for Alignment_File in os.listdir(Alignment_Dir):
        with open(os.path.join(Alignment_Dir, Alignment_File), "r") as AlignIn:
            Align_Temp = AlignIn.readlines()
            for key in Genome_Dictionary:
                for item in Align_Temp:
                    if key in item[1:5]:
                        Align_Temp[Align_Temp.index(item)] = ">" + Genome_Dictionary[key] + "\n"
        with open(os.path.join(Alignment_Dir, Alignment_File), "w") as AlignOut:
            AlignOut.writelines(Align_Temp)


def GLIMPSe_log(Alignment_Dir, Protein_Distribution, Alignment_Length, HMMer_Time, ClustalOmega_Time, FastTree_Time,
                RAxML_Time, Pipeline_Time, Alignment_Filtering, Log_Dir):
    """Generates log for pipeline"""
    with open(os.path.join(Log_Dir, "Pipeline Log.txt"), "w") as log:
        Accepted_Protein_Number = 0
        if Alignment_Filtering == "Trim":
            for alignment in os.listdir(Alignment_Dir):
                if alignment.startswith("Trimmed_"):
                    Accepted_Protein_Number += 1
            log.write("Number of accepted protein families in concatenated alignment = " + str(Accepted_Protein_Number))
            log.write(
                "\nMinimum proportion of organisms (organisms identified/total number of organisms) in accepted protein families = " + str(
                    Protein_Distribution * 100) + "%")
            log.write("\nLenght of trimmed concatenated alignment = " + str(Alignment_Length))
            log.write("\n\nThe following proteins were used in the concatenated alignment:")
            for alignment in os.listdir(Alignment_Dir):
                if alignment.startswith("Trimmed_"):
                    log.write("\n" + os.path.splitext(alignment)[0][8:])
        elif Alignment_Filtering == "Weight":
            for alignment in os.listdir(Alignment_Dir):
                if alignment.startswith("Weighted_"):
                    Accepted_Protein_Number += 1
            log.write("Number of accepted protein families in concatenated alignment = " + str(Accepted_Protein_Number))
            log.write(
                "\nMinimum proportion of organisms (organisms identified/total number of organisms) in accepted protein families = " + str(
                    Protein_Distribution * 100) + "%")
            log.write("\nLenght of unweighted concatenated alignment = " + str(Alignment_Length))
            log.write("\n\nThe following proteins were used in the concatenated alignment:")
            for alignment in os.listdir(Alignment_Dir):
                if alignment.startswith("Weighted_"):
                    log.write("\n" + os.path.splitext(alignment)[0][9:])
        log.write("\n\nProtein Family Identification Duration = " + str(round(HMMer_Time, 2)) + " seconds")
        log.write("\nDuration of Final ClustalOmega Alignments = " + str(round(ClustalOmega_Time, 2)) + " seconds")
        log.write("\nOperational time of FastTree = " + str(round(FastTree_Time, 2)) + " seconds")
        log.write("\nOperational time of RAxML = " + str(round(RAxML_Time, 2)) + " seconds")
        log.write("\nTotal operational time of pipeline = " + str(round(Pipeline_Time, 2)) + " seconds")


def Core_Pipeline(Input_Directory, Target_Proteins, Protein_Distribution, Alignment_Filtering, PAMatrix, POCP, AAI,
                  Single_Copy, Marker_Proteins, Fast_Cluster, Fast_Phylogeny, No_Tree, Genome_Dir, Protein_Dir,
                  Alignment_Dir, Concatenated_Dir, Tree_Dir, Log_Dir, GLIMPSe_Output_Dir, Marker_Dir, CDHIT, JACKHMMER,
                  HMMBUILD, HMMSEARCH, ClustalOmega, TrimAl, FastTree, RAxML, Threads, stdout_messenger,
                  stderr_messenger):
    """Main pipeline"""
    Output = GLIMPS_Writer(stdout_messenger, stderr_messenger)
    Pipeline_Start = time.time()
    Output.write(":::PIPELINE PREPERATION:::\n")
    if not os.path.exists(os.path.join(Marker_Dir, "bacteria.hmm")):
        Output.write("Extracting PhyEco Marker proteins...\n")
        try:
            with tarfile.open(os.path.join(Marker_Dir, "PhyEco Marker Protein Families.tar.bz2"), "r:bz2") as PhyEco:
                PhyEco.extractall(Marker_Dir)
            Output.write("PhyEco Marker proteins extracted.\n")
        except IOError:
            Output.error("Unable to detect or extract all PhyEco Markers.\n")
    Output.write("Initial processing of input files...\n")
    Genome_Dictionary = Process_Genome_Files(Input_Directory, Genome_Dir, Log_Dir, Output)
    Output.write("Initial processing of input files complete.\n")

    Output.write("\n:::PROTEIN FAMILY IDENTIFICATION:::\n")
    if os.path.exists(Target_Proteins):
        Output.write("Identifying protein families...\n")
        HMMer_Start = time.time()
        Protein_Clusters = Find_Proteins(Genome_Dir, Protein_Dir, GLIMPSe_Output_Dir, Genome_Dictionary, Target_Proteins, JACKHMMER, PAMatrix, Log_Dir, Output)
        HMMer_Time = time.time() - HMMer_Start
        Accepted_Proteins = Determine_Protein_Distribution(Protein_Distribution, Protein_Clusters, Genome_Dictionary, Output)
        Output.write("All protein families identiefied.\n")
    elif Marker_Proteins != "":
        Output.write("Identifying protein families...\n")
        HMMer_Start = time.time()
        Protein_Clusters = Find_Marker_Proteins(Genome_Dir, Protein_Dir, GLIMPSe_Output_Dir, Genome_Dictionary, HMMSEARCH, PAMatrix, Log_Dir, Marker_Dir, Marker_Proteins, Output)
        HMMer_Time = time.time() - HMMer_Start
        Accepted_Proteins = Determine_Protein_Distribution(Protein_Distribution, Protein_Clusters, Genome_Dictionary, Output)
        Output.write("All protein families identiefied.\n")
    else:
        Output.write("Identifying core protein families...\n")
        HMMer_Start = time.time()
        HMM_Clusters, Concat_Gen_Dict = ClusterCoreProts(Genome_Dir, Protein_Dir, Alignment_Dir, CDHIT, ClustalOmega, HMMBUILD, HMMSEARCH, Log_Dir, Threads, Output, Fast_Cluster)
        HMMer_Time = time.time() - HMMer_Start
        Filtered_Clusters = Protein_Family_Filter(Protein_Dir, Alignment_Dir, Concat_Gen_Dict, HMM_Clusters, GLIMPSe_Output_Dir, Genome_Dictionary, PAMatrix, POCP, Single_Copy, Output)
        Accepted_Proteins = Determine_Protein_Distribution(Protein_Distribution, Filtered_Clusters, Genome_Dictionary, Output)
        Output.write("All protein families identiefied.\n")

    Output.write("\n:::PROTEIN FAMILY ALIGNMENT AND TRIMMING:::\n")
    Output.write("Performing ClustalOmega alignments...\n")
    ClustalOmega_Start = time.time()
    if AAI:
        ParallelAlignment(Protein_Dir, Alignment_Dir, ClustalOmega, Log_Dir, Threads, False, [], Output)
        ClustalOmega_Time = time.time() - ClustalOmega_Start
        Output.write("All ClustalOmega alignments complete.\n")
        Output.write("Producing Average Amino Acid Identity Matrix...\n")
        Calculate_AAI(Alignment_Dir, Genome_Dictionary, Log_Dir, GLIMPSe_Output_Dir, Threads)
        Output.write("Average Amino Acid Identity Matrix Produced.\n")
    elif not AAI:
        ParallelAlignment(Protein_Dir, Alignment_Dir, ClustalOmega, Log_Dir, Threads, False, Accepted_Proteins, Output)
        ClustalOmega_Time = time.time() - ClustalOmega_Start
        Output.write("All ClustalOmega alignments complete.\n")
    else:
        Output.error("\n\nGLIMPSe pipeline error. Check Logs.\n")
        sys.exit()
    Alignment_Length = 0
    if Alignment_Filtering == "Trim":
        Output.write("Trimming alignments...\n")
        Create_Trimmed_Alignments(Alignment_Dir, Accepted_Proteins, TrimAl, Output)
        Output.write("All alignments trimmed.\n")
    elif Alignment_Filtering == "Weight":
        Output.write("Creating weighted alignments...\n")
        Alignment_Length = Create_Weighted_Alignments(Alignment_Dir, Accepted_Proteins, TrimAl, Output)
        Output.write("Weighted alignments created.\n")
    Output.write("Concatenating alignments...\n")
    if Alignment_Filtering == "Trim":
        Alignment_Length = Concatenate_Alignments(Accepted_Proteins, Alignment_Dir, Concatenated_Dir, Genome_Dictionary, GLIMPSe_Output_Dir, Output)
    elif Alignment_Filtering == "Weight":
        Concatenate_Alignments(Accepted_Proteins, Alignment_Dir, Concatenated_Dir, Genome_Dictionary, GLIMPSe_Output_Dir, Output)
    else:
        Output.error("\n\nGLIMPSe pipeline error. Check Logs.\n")
        sys.exit()
    Output.write("Alignments concatenated.\n")

    Output.write("\n:::PHYLOGENETIC TREE CONSTRUCTION:::\n")
    if No_Tree:
        FastTree_Time = 0
        RAxML_Time = 0
    else:
        if Fast_Phylogeny:
            Output.write("Building tree using FastTree...\n")
            FastTree_Start = time.time()
            Run_FastTree_Full(Tree_Dir, Concatenated_Dir, FastTree, Log_Dir, Output, GLIMPSe_Output_Dir)
            FastTree_Time = time.time() - FastTree_Start
            RAxML_Time = 0
            Output.write("Phylogenetic tree completed.\n")
        else:
            Output.write("Building initial tree using FastTree...\n")
            FastTree_Start = time.time()
            Run_FastTree(Tree_Dir, Concatenated_Dir, FastTree, Log_Dir, Output)
            FastTree_Time = time.time() - FastTree_Start
            Output.write("Using RaxML to optimize tree...\n")
            RAxML_Start = time.time()
            Run_RAxML(Tree_Dir, Concatenated_Dir, RAxML, Threads, Log_Dir, GLIMPSe_Output_Dir, Output)
            RAxML_Time = time.time() - RAxML_Start
            Output.write("Final phylogenetic tree completed.\n")

    Output.write("\nWriting output files...\n")
    GenomeID2GenomeName(Protein_Dir, Alignment_Dir, Genome_Dictionary)
    Pipeline_Time = time.time() - Pipeline_Start
    GLIMPSe_log(Alignment_Dir, Protein_Distribution, Alignment_Length, HMMer_Time, ClustalOmega_Time, FastTree_Time, RAxML_Time, Pipeline_Time, Alignment_Filtering, Log_Dir)
    Output.write("\nPipeline complete in %s seconds\nOutput files can be located in %s \n" % (str(round(Pipeline_Time, 2)), GLIMPSe_Output_Dir))


def main():
    args = check_arguments()
    stdout_messenger = multiprocessing.Manager().Queue()
    stderr_messenger = multiprocessing.Manager().Queue()
    Genome_Dir, Protein_Dir, Alignment_Dir, Concatenated_Dir, Tree_Dir, Log_Dir, GLIMPSe_Output_Dir, Dependency_Dir, Marker_Dir = Build_Output_Dirs(
        args.Output_Directory, stdout_messenger, stderr_messenger)
    CDHIT, JACKHMMER, HMMBUILD, HMMSEARCH, ClustalOmega, TrimAl, FastTree, RAxML, Threads = Prepare_Dependencies(Dependency_Dir, stdout_messenger, stderr_messenger)
    Core_Pipeline(args.Input_Directory, args.Target_Proteins, args.Protein_Distribution, args.Alignment_Filtering,
                  args.PAMatrix, args.POCP, args.AAI, args.Single_Copy, args.Marker_Proteins, args.Fast_Cluster,
                  args.Fast_Phylogeny, args.No_Tree, Genome_Dir, Protein_Dir, Alignment_Dir, Concatenated_Dir, Tree_Dir,
                  Log_Dir, GLIMPSe_Output_Dir, Marker_Dir, CDHIT, JACKHMMER, HMMBUILD, HMMSEARCH, ClustalOmega, TrimAl,
                  FastTree, RAxML, Threads, stdout_messenger, stderr_messenger)


if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
