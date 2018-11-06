*********************************
*    DNA Asymmetry browser.     *	
*				*				
*				*	
* Authors : BAL Matthieu	*
*           GHEZAIEL Morad	*
*				*
* Master 2 DLAD			*
*				*
*********************************


This python code can be used to detect and analyze DNA assymetry. The main code is subdivided into 2 pipelines :

1) GC-skew analysis (independantly runnable)

	- Performs DNA content, GC skew (WS : Window size is tunable) and statistical analysis
	- Plots for GC-skew (inter/intra specific) and statistical analysis
	- General informations for each sequences are displayed on the prompt

2) Genome signature analysis (need to launch the previous one)

	- Performs DN distribution analysis for each sequences using random generated sequences (these have to be obtained externally)
	- Perform inter and intra species statistical tests
	- Plot the ratio between random generated and biological sequences.


Working directory structure : 

- Species folders
- Random-seq folder
- Root:
	- GC skew profils (intraspecific)
	- GC skew profils (interspecific)
	- Genome signature profils (intraspecific)
	- Genome signature profils (interspecific)

Data format : FASTA

Dependencies that must be installed : pandas, scipy, numpy, matplotlib

Filenames are stored in a list (cf main). 

NB : This program was developped in a master project context with choosen sequence. If one wants to use it for differents sequences, 
     please do the following : 
	
1) Run the first pipeline 
2) Pull nucleotide distribution from the prompt
3) Go to http://users-birc.au.dk/biopv/php/fabox/random_sequence_generator.php 
4) Fill nucleotide distribution field with those pulled from the prompt
5) Generate n_random_sets(cf main) for each seqs
5) Save seqs as 10_AT_chr1.txt (1 set of 10 random sequences generated according to arabidopsis thaliana first chromosome nucleotide %)
6) Modify the n_random_sets value in the main
7) Modify the lists ontop of the main, according to changes
8) Launch it ! 








