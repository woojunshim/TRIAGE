# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

We provide an example Python script (written in version 2.7) to calculate the discordance score of genes given their expression readout table. With simplicity of the TRIAGE implementation, users may want to write their own codes or simply incorporate TRIAGE as a part of their analysis.

"disc.py" is the main script to run the TRIAGE analysis. It requires (i) input expression matrix (e.g. "example_input.txt") and (ii) repressive tendency score (RTS) file (e.g. "human_rts.txt") both of which are plain text files. You can find the both in the repository. 

First the script reads input expression data in the text file format (tab-delimited). 
Input file is a matrix that defines expression values of genes across samples (where rows are genes and columns are samples).  
See "example_input.txt" below (which is also included in the repository) for an example input matrix. 

<example_input.txt>

Gene	Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium)
A1BG	3.065	1.983	4.45
A1CF	0.073	0.0	0.0
A2M	0.402	24.878	137.612
A2ML1	0.002	0.394	0.166
A4GALT	0.006	0.199	42.44
A4GNT	0.024	0.0	0.006
AAAS	12.397	36.851	12.696
AACS	2.594	3.098	2.428
AADAC	0.0	0.0	0.041
AADACL2	0.0	0.0	0.014
AADACL3	0.0	0.0	0.006
AADACL4	0.0	0.0	0.0
AADAT	0.029	18.335	1.234
AAED1	4.434	1.228	10.104
AAGAB	8.124	2.23	6.885
AAK1	62.15	13.05	4.018
AAMDC	2.795	7.934	42.186
AAMP	18.069	24.709	29.215
AANAT	0.233	0.072	0.017
AAR2	22.087	12.372	11.73
AARD	0.0	0.416	0.028
AARS	10.665	29.972	39.914
AARS2	9.836	6.63	5.326
AASDH	4.597	5.375	2.416

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
