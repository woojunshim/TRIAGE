# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

We provide an example Python script (written in version 2.7) to calculate the discordance score of genes given their expression readouts. With simplicity of the TRIAGE implementation, users may also want to write their own codes or simply incorporate TRIAGE as a part of their analysis.

"disc.py" is the main script to run the TRIAGE analysis. It requires (i) input expression matrix (e.g. "example_input.txt") and (ii) repressive tendency score (RTS) file (e.g. "human_rts.txt") both of which are plain text files. You can find the both files in the repository. 

Input expression file is a matrix that defines expression values of genes across samples (where rows are genes and columns are samples). See "example_input.txt" below (which is also included in the repository). 

<example_input.txt>

Gene	Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium). 
A1BG	3.065	1.983	4.45. 
A1CF	0.073	0.0	0.0. 
A2M	0.402	24.878	137.612. 
A2ML1	0.002	0.394	0.166. 
  
...  
...  


You can use an alternative RTS table ("mouse_rts_mapped.txt") for mouse datasets. Note that the mouse RTS values were obtained by directly mapping genes between human and mouse data. The mouse data currently only covers protein-coding genes.

Parameters
  -i input file name (required)
  -o output file name (required)
  -p pseudo-count to be added to the input data, default = 1
  -l natural log-transformation of the input data, default = True    
  
  E.g. >>> python disc.py -i example_input.txt -o output_file.txt 

If any issues are found, please contact to w.shim@uq.edu.au
