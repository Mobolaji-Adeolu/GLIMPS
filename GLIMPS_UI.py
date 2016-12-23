#!/usr/bin/env python

"""Gupta Lab Integrated Microbial Phylogenetic and Supermatrix Pipeline
Written by Mobolaji Adeolu (Adeolum@McMaster.ca)
Department of Biochemistry and Biomedical Sciences
McMaster University.
Copyright 2016."""

import sys
try:
    import Tkinter
except ImportError:
    Tkinter = None
    sys.exit("Required python module, Tkinter, is unavailable. Please install TKinter to use the GLIMPS GUI.")
import os
import subprocess
import ScrolledText
import tkFileDialog
import tkMessageBox
import ttk
import multiprocessing


class GUI(ttk.Frame):
    """Controls the main UI, subclasses Themed Tk Frame"""
    def __init__(self, parent):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        self.Queue = multiprocessing.Manager().Queue()
        self.title = 'Error'
        self.text = 'The script returned an error. Check inputs/logs and try again.'

        # Variables for main UI screen
        self.inFolder = Tkinter.StringVar()
        self.inFolder.set('No Input Directory Selected')
        self.outFolder = Tkinter.StringVar()
        self.outFolder.set('No Output Directory Selected')
        self.TargetFasta = Tkinter.StringVar()
        self.TargetFasta.set('No Fasta File Selected')
        self.RadioVar = Tkinter.StringVar()
        self.RadioVar.set("core")
        self.MarkerSet = Tkinter.StringVar()
        self.MarkerSet.set("")
        self.hasIn = False
        self.hasOut = False
        self.LenOK = False
        self.hasFasta = False

        # Virables for STDOUT box
        self.progBarName = ttk.Label(self, text="Running...")
        self.progBar = ttk.Progressbar(self, mode='indeterminate')
        self.run = ttk.Button(self, text='Run Pipeline', command=self.check_inputs)
        self.settings = ttk.Button(self, text='Settings', command=self.SettingsUI)
        self.stdout = ScrolledText.ScrolledText(self, width=50, height=10, bg='#ffffff', font=("Helvetica", "12"), state=Tkinter.DISABLED)
        self.stdout.tag_config("normal", foreground="#000000")
        self.stdout.tag_config("error", foreground="#cc0000")
        sys.stdout = ScriptOutput(self.stdout, "normal")
        self.stdout_messenger = multiprocessing.Manager().Queue()
        sys.stderr = ScriptOutput(self.stdout, "error")
        self.stderr_messenger = multiprocessing.Manager().Queue()

        # Variables for settings screen
        self.Protein_Distribution = Tkinter.IntVar()
        self.Protein_Distribution.set(80)
        self.PAMatrix = Tkinter.IntVar()
        self.PAMatrix.set(0)
        self.POCP = Tkinter.IntVar()
        self.POCP.set(0)
        self.AAI = Tkinter.IntVar()
        self.AAI.set(0)
        self.Phylogeny = Tkinter.IntVar()
        self.Phylogeny.set(1)
        self.Single_Copy = Tkinter.IntVar()
        self.Single_Copy.set(0)
        self.Polymorphic = Tkinter.IntVar()
        self.Polymorphic.set(0)
        self.Fast_Clust = Tkinter.IntVar()
        self.Fast_Clust.set(0)
        self.Fast_Phylo = Tkinter.IntVar()
        self.Fast_Phylo.set(0)
        self.Threads = Tkinter.IntVar()
        if isinstance(multiprocessing.cpu_count(), int):
            self.Threads.set(multiprocessing.cpu_count())
        else:
            self.Threads.set(1)

        # Identify directories for dependencies
        import GLIMPS_Pipeline
        GLIMPS_CDHIT, GLIMPS_JACKHMMER, GLIMPS_HMMBUILD, GLIMPS_HMMSEARCH, GLIMPS_ClustalOmega, GLIMPS_TrimAl, GLIMPS_FastTree, GLIMPS_RAxML, Threads = GLIMPS_Pipeline.Prepare_Dependencies(os.path.join(os.path.dirname(sys.argv[0]), "Dependencies"), self.stdout_messenger, self.stderr_messenger)

        # Dependency variables
        self.CDHIT = Tkinter.StringVar()
        self.CDHIT.set(GLIMPS_CDHIT)
        self.JACKHMMER = Tkinter.StringVar()
        self.JACKHMMER.set(GLIMPS_JACKHMMER)
        self.HMMBUILD = Tkinter.StringVar()
        self.HMMBUILD.set(GLIMPS_HMMBUILD)
        self.HMMSEARCH = Tkinter.StringVar()
        self.HMMSEARCH.set(GLIMPS_HMMSEARCH)
        self.ClustalOmega = Tkinter.StringVar()
        self.ClustalOmega.set(GLIMPS_ClustalOmega)
        self.TrimAl = Tkinter.StringVar()
        self.TrimAl.set(GLIMPS_TrimAl)
        self.FastTree = Tkinter.StringVar()
        self.FastTree.set(GLIMPS_FastTree)
        self.RAxML = Tkinter.StringVar()
        self.RAxML.set(GLIMPS_RAxML)

        # Initiation of mainloop
        ttk.Style().configure('TLabelframe', padding=5)
        self.grid()
        self.createUI()
        self.parent.mainloop()

    def createUI(self):
        """Controls the main UI seen when GLIMPS is run"""
        # Input Folder Selection
        Input = ttk.LabelFrame(self, text="Select Input Directory")
        Input.grid(row=0, columnspan=4, sticky="NSWE", pady=5, padx=10)
        InputName = ttk.Entry(Input, textvariable=self.inFolder, state="readonly", width=50)
        InputName.grid(row=0, column=0, sticky="NSWE", padx=5)
        InputSelect = ttk.Button(Input, text='Browse...', command=self.getInputFolder)
        InputSelect.grid(row=0, column=1, sticky='NSWE', padx=5)

        # Output Folder Selection
        Output = ttk.LabelFrame(self, text='Select Output Directory')
        Output.grid(row=1, columnspan=4, sticky='NSWE', pady=5, padx=10)
        OutputName = ttk.Entry(Output, textvariable=self.outFolder, state="readonly", width=50)
        OutputName.grid(row=0, column=0, sticky='NSWE', padx=5)
        OutputSelect = ttk.Button(Output, text='Browse...', command=self.getOutputFolder)
        OutputSelect.grid(row=0, column=1, sticky='NSWE', padx=5)

        # Core Genome Identification Selection
        Core = ttk.LabelFrame(self, text='Identify Shared Protein Families in Core Genome')
        Core.grid(row=2, columnspan=4, sticky='NSWE', pady=5, padx=10)
        CoreMessage = "  Allow the GLIMPS pipeline to automatically identify all\n" \
                      + "  shared core protein families in the input genome files."

        # Target Fasta File Selection
        Fasta = ttk.LabelFrame(self, text='Select Target Protein Fasta File')
        Fasta.grid(row=3, columnspan=4, sticky='NSWE', pady=5, padx=10)
        FastaName = ttk.Entry(Fasta, state="disabled", textvariable=self.TargetFasta, width=46)
        FastaName.grid(row=0, column=1, sticky='NSWE', padx=5)
        FastaSelect = ttk.Button(Fasta, state="disabled", text='Browse...', command=self.getFastaFile)
        FastaSelect.grid(row=0, column=2, sticky='NSWE', padx=5)

        # PhyEco Marker Selection
        PyhEco = ttk.LabelFrame(self, text='Select PhyEco Marker Protein Family')
        PyhEco.grid(row=4, columnspan=4, sticky='NSWE', pady=5, padx=10)
        Markers = ("Bacteria and Archaea", "Archaea", "Bacteria", "Actinobacteria", "Alphaproteobacteria",
                   "Bacteroidetes", "Betaproteobacteria", "Chlamydiae", "Chloroflexi", "Cyanobacteria",
                   "Deinococcus-Thermus", "Deltaproteobacteria", "Epsilonproteobacteria", "Firmicutes",
                   "Gammaproteobacteria", "Spirochaetes", "Thermotogae")
        PyhEcoSelect = ttk.Combobox(PyhEco, state="disabled", values=Markers, textvariable=self.MarkerSet, width=46)
        PyhEcoSelect.grid(row=0, column=1, sticky='NSWE', padx=5)

        # Button Row
        self.run.grid(row=8, column=0, columnspan=2, pady=5)
        self.settings.grid(row=8, column=2, columnspan=2, pady=5)

        def RadioButtons():
            """Local function for Radio buttons"""
            if self.RadioVar.get() == "core":
                FastaName.configure(state="disabled")
                self.TargetFasta.set('No Fasta File Selected')
                FastaSelect.configure(state="disabled")
                PyhEcoSelect.configure(state="disabled")
                self.MarkerSet.set("")
            elif self.RadioVar.get() == "fasta":
                FastaName.configure(state="readonly")
                FastaSelect.configure(state="readonly")
                PyhEcoSelect.configure(state="disabled")
                self.MarkerSet.set("")
            elif self.RadioVar.get() == "marker":
                FastaName.configure(state="disabled")
                self.TargetFasta.set('No Fasta File Selected')
                FastaSelect.configure(state="disabled")
                PyhEcoSelect.configure(state="readonly")

        # RadioButtons
        R1 = ttk.Radiobutton(Core, variable=self.RadioVar, value="core", command=RadioButtons, text=CoreMessage)
        R1.grid(row=0, column=0, sticky='NSWE', padx=5)
        R2 = ttk.Radiobutton(Fasta, variable=self.RadioVar, value="fasta", command=RadioButtons)
        R2.grid(row=0, column=0, sticky='NSWE', padx=5)
        R3 = ttk.Radiobutton(PyhEco, variable=self.RadioVar, value="marker", command=RadioButtons)
        R3.grid(row=0, column=0, sticky='NSWE', padx=5)

    def SettingsUI(self):
        """Settings UI generated using Toplevel with widgets in a Themed Tk Frame"""
        SettingsWindow = Tkinter.Toplevel(self.parent)
        SettingsWindow.title("GLIMPS Settings")
        SettingsFrame = ttk.Frame(SettingsWindow)
        SettingsFrame.grid()

        # Output files produced by pipeline
        Outputs = ttk.LabelFrame(SettingsFrame, text="Outputs")
        Outputs.grid(row=0, columnspan=3, sticky="NSWE", padx=10, pady=5)
        PAmat = ttk.Checkbutton(Outputs, variable=self.PAMatrix, text='Presence Absence Matrix')
        PAmat.grid(row=0, column=0, sticky='NSW', padx=5, pady=5)
        POCPmat = ttk.Checkbutton(Outputs, variable=self.POCP,
                                  text='Percentage of Shared Protein Families Matrix (Only in automated clustering)')
        POCPmat.grid(row=1, column=0, sticky='NSW', padx=5, pady=5)
        AAImat = ttk.Checkbutton(Outputs, variable=self.AAI, text='Average Amino Acid Identity Matrix')
        AAImat.grid(row=2, column=0, sticky='NSW', padx=5, pady=5)
        Phylo = ttk.Checkbutton(Outputs, variable=self.Phylogeny, text='Phylogenetic Tree')
        Phylo.grid(row=3, column=0, sticky='NSW', padx=5, pady=5)

        # Pipeline algorithm settings
        Pipeline_Settings = ttk.LabelFrame(SettingsFrame, text="Pipeline Settings")
        Pipeline_Settings.grid(row=1, columnspan=3, sticky="NSWE", padx=10, pady=5)
        Poly = ttk.Checkbutton(Pipeline_Settings, variable=self.Polymorphic,
                             text='Use only polymorphic amino acid positions in phylogenetic analysis')
        Poly.grid(row=0, columnspan=2, sticky='NSW', padx=10, pady=5)
        SC = ttk.Checkbutton(Pipeline_Settings, variable=self.Single_Copy,
                             text='Use only single copy homologs (Only in automated clustering)')
        SC.grid(row=1, columnspan=2, sticky='NSW', padx=10, pady=5)
        FC = ttk.Checkbutton(Pipeline_Settings, variable=self.Fast_Clust,
                             text='Skip HMMer steps in automated core protein identification (Fast Clustering)')
        FC.grid(row=2, columnspan=2, sticky='NSW', padx=10, pady=5)
        FP = ttk.Checkbutton(Pipeline_Settings, variable=self.Fast_Phylo,
                             text='Skip RAxML step in phylogenetic tree construction (Fast Phylogeny)')
        FP.grid(row=3, columnspan=2, sticky='NSW', padx=10, pady=5)
        ThreadLabel = ttk.Label(Pipeline_Settings, text='Number of threads to utilize in pipeline')
        ThreadLabel.grid(row=4, column=0, sticky='NSE', padx=10, pady=5)
        Num_Thread = self.Threads.get()
        Thread = Tkinter.Spinbox(Pipeline_Settings, textvariable=self.Threads, from_=1, to=Num_Thread, width=3)
        Thread.grid(row=4, column=1, sticky='NSW', padx=10, pady=5)
        PDLabel = ttk.Label(Pipeline_Settings, text='Acceptnce threshold for distribution of core proteins (%)')
        PDLabel.grid(row=5, column=0, sticky='NSE', padx=10, pady=5)
        PD = Tkinter.Spinbox(Pipeline_Settings, textvariable=self.Protein_Distribution, from_=1, to=100, width=3)
        PD.grid(row=5, column=1, sticky='NSW', padx=10, pady=5)

        # Pipeline dependencies
        Dependencies = ttk.LabelFrame(SettingsFrame, text='Dependencies')
        Dependencies.grid(row=3, columnspan=3, sticky="NSWE", padx=10, pady=5)
        CDHITname = ttk.Label(Dependencies, text='CD-Hit')
        CDHITname.grid(row=0, column=0, sticky='E', padx=5)
        CDHITentry = ttk.Entry(Dependencies, textvariable=self.CDHIT, width=30)
        CDHITentry.grid(row=0, column=1, sticky='NSWE', padx=5)
        CDHITSelect = ttk.Button(Dependencies, text='Browse...', command=self.getCDHITFile)
        CDHITSelect.grid(row=0, column=2, sticky='E', padx=5)
        JACKHMMERname = ttk.Label(Dependencies, text='JackHMMer')
        JACKHMMERname.grid(row=1, column=0, sticky='E', padx=5)
        JACKHMMERentry = ttk.Entry(Dependencies, textvariable=self.JACKHMMER, width=30)
        JACKHMMERentry.grid(row=1, column=1, sticky='NSWE', padx=5)
        JACKHMMERSelect = ttk.Button(Dependencies, text='Browse...', command=self.getJACKHMMERFile)
        JACKHMMERSelect.grid(row=1, column=2, sticky='NSWE', padx=5)
        HMMbuildname = ttk.Label(Dependencies, text='HMMbuild')
        HMMbuildname.grid(row=2, column=0, sticky='E', padx=5)
        HMMbuildentry = ttk.Entry(Dependencies, textvariable=self.HMMBUILD, width=30)
        HMMbuildentry.grid(row=2, column=1, sticky='NSWE', padx=5)
        HMMbuildSelect = ttk.Button(Dependencies, text='Browse...', command=self.getHMMBUILDFile)
        HMMbuildSelect.grid(row=2, column=2, sticky='NSWE', padx=5)
        HMMsearchname = ttk.Label(Dependencies, text='HMMsearch')
        HMMsearchname.grid(row=3, column=0, sticky='E', padx=5)
        HMMsearchentry = ttk.Entry(Dependencies, textvariable=self.HMMSEARCH, width=30)
        HMMsearchentry.grid(row=3, column=1, sticky='NSWE', padx=5)
        HMMsearchSelect = ttk.Button(Dependencies, text='Browse...', command=self.getHMMSEARCHFile)
        HMMsearchSelect.grid(row=3, column=2, sticky='NSWE', padx=5)
        ClustalOname = ttk.Label(Dependencies, text='Clustal Omega')
        ClustalOname.grid(row=4, column=0, sticky='E', padx=5)
        ClustalOentry = ttk.Entry(Dependencies, textvariable=self.ClustalOmega, width=30)
        ClustalOentry.grid(row=4, column=1, sticky='NSWE', padx=5)
        ClustalOSelect = ttk.Button(Dependencies, text='Browse...', command=self.getClustalOmegaFile)
        ClustalOSelect.grid(row=4, column=2, sticky='NSWE', padx=5)
        TrimAlname = ttk.Label(Dependencies, text='TrimAl')
        TrimAlname.grid(row=5, column=0, sticky='E', padx=5)
        TrimAlentry = ttk.Entry(Dependencies, textvariable=self.TrimAl, width=30)
        TrimAlentry.grid(row=5, column=1, sticky='NSWE', padx=5)
        TrimAlSelect = ttk.Button(Dependencies, text='Browse...', command=self.getTrimAlFile)
        TrimAlSelect.grid(row=5, column=2, sticky='NSWE', padx=5)
        FastTreename = ttk.Label(Dependencies, text='FastTree')
        FastTreename.grid(row=6, column=0, sticky='E', padx=5)
        FastTreeentry = ttk.Entry(Dependencies, textvariable=self.FastTree, width=30)
        FastTreeentry.grid(row=6, column=1, sticky='NSWE', padx=5)
        FastTreeSelect = ttk.Button(Dependencies, text='Browse...', command=self.getFastTreeFile)
        FastTreeSelect.grid(row=6, column=2, sticky='NSWE', padx=5)
        RAxMLname = ttk.Label(Dependencies, text='RAxML')
        RAxMLname.grid(row=7, column=0, sticky='E', padx=5)
        RAxMLentry = ttk.Entry(Dependencies, textvariable=self.RAxML, width=30)
        RAxMLentry.grid(row=7, column=1, sticky='NSWE', padx=5)
        RAxMLSelect = ttk.Button(Dependencies, text='Browse...', command=self.getRAxMLFile)
        RAxMLSelect.grid(row=7, column=2, sticky='NSWE', padx=5)

        # Close settings button
        Close_Settings = ttk.Button(SettingsFrame, text='Close Settings', command=SettingsWindow.destroy)
        Close_Settings.grid(row=4, columnspan=3, pady=5)

    def progressBar_Start(self):
        """Initiates progress bar and STDOUT box"""
        self.progBarName.grid(row=5, columnspan=4)
        self.progBar.grid(row=6, columnspan=4, sticky='NSWE')
        self.progBar.start()
        self.stdout.grid(row=7, columnspan=4, sticky='NSWE')
        self.run.grid_remove()
        self.settings.grid_remove()
        Close = ttk.Button(self, text='Close GLIMPS', command=self.Quit)
        Close.grid(row=8, columnspan=4, pady=5)

    def progressBar_Stop(self):
        """Stops progress bar"""
        self.progBar.stop()
        self.progBarName.grid_remove()
        self.progBar.grid_remove()

    def Quit(self):
        """Closes/destroys UI and quits pipeline"""
        self.parent.destroy()
        sys.exit()

    def process_queue(self):
        """Checks queue for the output of the Pipeline, runs recursively until pipeline completes"""
        try:
            Result = self.Queue.get(0)
            self.progressBar_Stop()
            if Result == "Success":
                self.Messenger()
                sys.stdout.write("\nPipeline successfully completed.\n\n")
                self.title = "Success"
                self.text = 'Program completed successfully.\n\nOutput files located at ' + self.outFolder.get()
                if sys.platform == 'darwin':
                    subprocess.check_call(['open', '--', self.outFolder.get()])
                elif sys.platform == 'linux2':
                    subprocess.check_call(['xdg-open', '--', self.outFolder.get()])
                elif sys.platform == 'win32':
                    subprocess.check_call(['explorer', self.outFolder.get()])
            elif Result == "Fail":
                self.Messenger()
                self.title = "Error"
                self.text = "Pipeline failure.\nCheck Logs."
                sys.stderr.write("\nPipeline failure.\nCheck Logs.\n\n")
            self.Popup()
        except Exception:
            self.Messenger()
            self.parent.after(100, self.process_queue)

    def Messenger(self):
        """Local function which checks messenger queue for any STDOUT output from the pipeline"""
        for x in range(5):
            try:
                sys.stdout.write(self.stdout_messenger.get_nowait())
            except Exception:
                break
        for x in range(5):
            try:
                sys.stderr.write(self.stderr_messenger.get_nowait())
            except Exception:
                break

    def Popup(self):
        """Generates Info/Warning message boxes"""
        if self.title == "Error":
            tkMessageBox.showerror(self.title, self.text, parent=self)
        elif self.title == "Warning":
            if tkMessageBox.askokcancel(self.title, self.text, parent=self):
                self.LenOK = True
                self.check_inputs()
        else:
            tkMessageBox.showinfo(self.title, self.text, parent=self)

    def check_inputs(self):
        """Checks for input/output folder selection"""
        if not self.hasIn:
            self.title = 'Error'
            self.text = 'Please select an input folder.'
            self.Popup()
        elif not self.hasOut:
            self.title = 'Error'
            self.text = 'Please select an output folder.'
            self.Popup()
        elif not self.LenOK and sys.platform.startswith("win32") and len(self.outFolder.get()) > 100:
            self.title = 'Warning'
            self.text = 'Path to output directory is ' + str(len(self.outFolder.get())) + ' characters long. Paths in Windows are limited to 260 characters. This output directory may produce pipeline errors.'
            self.Popup()
        elif not self.hasFasta:
            self.TargetFasta.set('')
            self.runPipeline()
        else:
            self.runPipeline()

    def runPipeline(self):
        """Passes variables to AsyncTask class and initiates Pipeline"""
        MarkerSet = self.MarkerSet.get().lower().replace(" ", "_")
        AsyncTask(self.Queue, self.inFolder.get(), self.outFolder.get(), self.TargetFasta.get(), MarkerSet,
                  self.Threads.get(), self.Protein_Distribution.get(), self.PAMatrix.get(), self.POCP.get(),
                  self.AAI.get(), self.Phylogeny.get(), self.Polymorphic.get(),
                  self.Single_Copy.get(), self.Fast_Clust.get(), self.Fast_Phylo.get(), self.CDHIT.get(),
                  self.JACKHMMER.get(), self.HMMBUILD.get(), self.HMMSEARCH.get(), self.ClustalOmega.get(),
                  self.TrimAl.get(), self.FastTree.get(), self.RAxML.get(), self.stdout_messenger,
                  self.stderr_messenger).start()
        self.progressBar_Start()
        self.parent.after(100, self.process_queue)

    # File/Folder Selection functions
    def getInputFolder(self):
        fname = tkFileDialog.askdirectory(title='Select Input Folder', parent=self)
        if fname:
            self.inFolder.set(os.path.normpath(fname))
            self.hasIn = True

    def getOutputFolder(self):
        fname = tkFileDialog.askdirectory(title='Select Input Folder', parent=self)
        if fname:
            self.outFolder.set(os.path.normpath(fname))
            self.hasOut = True

    def getFastaFile(self):
        fname = tkFileDialog.askopenfilename(title='Select File', parent=self)
        if fname:
            self.TargetFasta.set(os.path.normpath(fname))
            self.hasFasta = True

    def getCDHITFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.CDHIT.set(os.path.normpath(fname))

    def getJACKHMMERFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.JACKHMMER.set(os.path.normpath(fname))

    def getHMMBUILDFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.HMMBUILD.set(os.path.normpath(fname))

    def getHMMSEARCHFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.HMMSEARCH.set(os.path.normpath(fname))

    def getClustalOmegaFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.ClustalOmega.set(os.path.normpath(fname))

    def getTrimAlFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.TrimAl.set(os.path.normpath(fname))

    def getFastTreeFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.FastTree.set(os.path.normpath(fname))

    def getRAxMLFile(self):
        fname = tkFileDialog.askopenfilename(title='Select Dependency File', parent=self)
        if fname:
            self.RAxML.set(os.path.normpath(fname))


class AsyncTask(multiprocessing.Process):
    """Class for running the Pipeline in a seperate thread, subclasses multiprocessing"""

    def __init__(self, Queue, Input_Directory, Output_Directory, Target_Proteins, MarkerSet, Threads,
                 Protein_Distribution, PAMatrix, POCP, AAI, Phylogeny, Polymorphic, Single_Copy, Fast_Clust, Fast_Phylo, CDHIT,
                 JACKHMMER, HMMBUILD, HMMSEARCH, ClustalOmega, TrimAl, FastTree, RAxML, stdout_messenger,
                 stderr_messenger):
        multiprocessing.Process.__init__(self)
        self.Queue = Queue
        self.Input_Directory = Input_Directory
        self.Output_Directory = Output_Directory
        self.Threads = Threads
        self.Target_Proteins = Target_Proteins
        self.Marker_Proteins = MarkerSet
        self.Protein_Distribution = Protein_Distribution * 0.01
        self.Alignment_Filtering = "Trim"
        if PAMatrix == 1:
            self.PAMatrix = True
        else:
            self.PAMatrix = False
        if POCP == 1:
            self.POCP = True
        else:
            self.POCP = False
        if AAI == 1:
            self.AAI = True
        else:
            self.AAI = False
        if Phylogeny == 1:
            self.Phylogeny = False
        else:
            self.Phylogeny = True
        if Polymorphic == 1:
            self.Polymorphic = True
        if Single_Copy == 1:
            self.Single_Copy = True
        else:
            self.Single_Copy = False
        if Fast_Clust == 1:
            self.Fast_Clust = True
        else:
            self.Fast_Clust = False
        if Fast_Phylo == 1:
            self.Fast_Phylo = True
        else:
            self.Fast_Phylo = False
        self.CDHIT = CDHIT
        self.JACKHMMER = JACKHMMER
        self.HMMBUILD = HMMBUILD
        self.HMMSEARCH = HMMSEARCH
        self.ClustalOmega = ClustalOmega
        self.TrimAl = TrimAl
        self.FastTree = FastTree
        self.RAxML = RAxML
        self.stdout_messenger = stdout_messenger
        self.stderr_messenger = stderr_messenger

    def run(self):
        """Runs pipeline in a seperate process and reports success or failure, overwrites run in multiprocessing"""
        try:
            import GLIMPS_Pipeline
            Genome_Dir, Protein_Dir, Alignment_Dir, Concatenated_Dir, Tree_Dir, Log_Dir, GLIMPSe_Output_Dir, Dependency_Dir, Marker_Dir = GLIMPS_Pipeline.Build_Output_Dirs(self.Output_Directory, self.stdout_messenger, self.stderr_messenger)
            GLIMPS_Pipeline.Core_Pipeline(self.Input_Directory, self.Target_Proteins, self.Protein_Distribution,
                                          self.Alignment_Filtering, self.PAMatrix, self.POCP, self.AAI, self.Single_Copy,
                                          self.Marker_Proteins, self.Fast_Clust, self.Fast_Phylo, self.Phylogeny, self.Polymorphic,
                                          Genome_Dir, Protein_Dir, Alignment_Dir, Concatenated_Dir, Tree_Dir, Log_Dir,
                                          GLIMPSe_Output_Dir, Marker_Dir, self.CDHIT, self.JACKHMMER, self.HMMBUILD,
                                          self.HMMSEARCH,
                                          self.ClustalOmega, self.TrimAl, self.FastTree, self.RAxML, self.Threads,
                                          self.stdout_messenger, self.stderr_messenger)
            self.Queue.put("Success")
        except Exception:
            self.Queue.put("Fail")


class ScriptOutput(object):
    """Class for redirecting STDOUT/STDERR to STDOUT box"""

    def __init__(self, Output_Window, tag):
        self.output = Output_Window
        self.tag = tag

    def write(self, string):
        if self.output.yview()[1] == 1:
            self.output.config(state=Tkinter.NORMAL)
            self.output.insert("end", string, self.tag)
            self.output.see("end")
            self.output.config(state=Tkinter.DISABLED)
        elif not self.output.yview()[1] == 1:
            self.output.config(state=Tkinter.NORMAL)
            self.output.insert("end", string, self.tag)
            self.output.config(state=Tkinter.DISABLED)

    def flush(self):
        pass


def main():
    root = Tkinter.Tk()
    root.wm_title("GLIMPS Pipeline")
    GUI(root)


if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
