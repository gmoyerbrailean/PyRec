# PyRec - Python-aided Recombineering

PyRec is a collection of Python scripts designed to aid the creation of gene knockouts through recombineering in bacteria. PyRec reads the gene of choice and creates five consecutive mutations, generating both a nonsense mutation and a unique restriction site. Also output are primers for a 1kb region centered around the mutation, to test for incorporation via restriction digest and PCR.

The scripts adapt to any organism, genes, and restriction enzymes used, so long as you provide the necessary information. In order to create 1kb PCR products for verification PyRec requires the full genome of the organism in fasta format. Any genes from that organism you wish to mutate you may place in the 'genes' directory (or anywhere else - just specify the path). Finally, provide the cut sites for any restriction enzymes you wish to use (the default being EcoRI, BamHI, and HindIII). 

Early versions of the program were developed for an organism with an open genome. As such, I built in the ability to use NCBI's command-line BLAST tools, blast+, to instead use a related genome. In the current version, however, this option is no longer fully supported, though the BLASTing code remains.

To launch PyRec, use the console:

	python /path/to/PyRec/File/PyRecGUI.py
	
(The program was originally a command line, and I never got around to implementing a way to launch the GUI without the console.)