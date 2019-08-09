# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

We provide an example Python script (written in version 2.7) to calculate the discordance score of genes given their expression readouts. With simplicity of the TRIAGE implementation, users may also want to write their own codes or simply incorporate TRIAGE as a part of their analysis.

"disc.py" is the main script to run the TRIAGE analysis. It requires (i) input expression matrix (e.g. "example_input.txt") and (ii) repressive tendency score (RTS) file (e.g. "human_rts.txt") both of which are plain text files. You can find the both files in the repository. Simply, the script reads in the both files and calculates the discordance score after value conversion (e.g. natural log-transformation) if specified. Please note that the script does NOT normalise the input data. While the normalisation across samples is not a requirement for TRIAGE, if one concerns, this should be done before running the script.

Input expression file is a matrix that defines expression values of genes across samples (where rows are genes and columns are samples). See "example_input.txt" below (which is also included in the repository). This toy example is RNA-seq data (RPKM normalised) for 3 selected Roadmap samples.  

<example_input.txt>

	Gene	Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium) 

	A1BG	3.065	1.983	4.45 

	A1CF	0.073	0.0	0.0 

	A2M	0.402	24.878	137.612 

	A2ML1	0.002	0.394	0.166 

	...

The output of TRIAGE is weighted expression readouts (i.e. discordance score) by the gene-specific RTS. 


<example_output.txt>. 

		Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium) 
  
	A1BG	0.008059178900858537	0.006280681756778532	0.009744106976712222 

	A1CF	0.00027847821321769683	0.0	0.0 

	A2M	7.800106276867326e-05	0.0007510159364827811	0.0011384327244910546 

	A2ML1	8.465243140668874e-07	0.00014073863700441038	6.506919860136177e-05 

	...

Implementation of TRIAGE is fast and highly scalable. Computational time required to analyse the example dataset (18,707 genes * 3 samples) is approximately 0.3 second (run on a personal Mac with 16 GB RAM). 

You may like to use an alternative RTS table for mouse datasets ("mouse_rts_mapped.txt"). Note that the mouse RTS values were obtained by directly mapping genes between human and mouse data. The mouse data currently only covers protein-coding genes.

Finally, users can modify a pseudo-count (default = 1) or whether to log-transform the expression value (default = True) by specifying parameters.

Parameters
  -i input file name (required)
  -o output file name (required)
  -p pseudo-count to be added to the input data, default = 1
  -l natural log-transformation of the input data, default = True    
  
  E.g. >>> python disc.py -i example_input.txt -o example_file.txt 
  
  --> Run TRIAGE on "example_input.txt" and output the discordance score as a text file "example_output.txt", with pseudo-count of 1 and natural log-transformation 

If any issues are found, please contact to w.shim@uq.edu.au
