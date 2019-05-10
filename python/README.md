# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

disc.py is the main script to run the analysis and human_rts.txt is the human repressive tendency score (RTS) table. These files must be located in the same directory. 

First the script reads input expression data in the text file format (tab-delimited).
Input file is a matrix that defines expression values of genes across samples where rows are genes and columns are samples.  
See 'example_input.txt' for an acceptable input matrix format.

You can use an alternative RTS file (e.g. mouse_rts_mapped.txt for mouse genes) if you like. Note that the mouse RTS values were obtained by directly mapping genes between human and mouse data. The mouse data currently only covers protein-coding genes.

Parameters
  -i input file name (required)
  -o output file name (required)
  -p pseudo-count to be added to the input data, default = 1
  -l natural log-transformation of the input data, default = True    
  
  E.g. >>> python disc.py -i example_input.txt -o output_file.txt 

If any issues are found, please contact to w.shim@uq.edu.au
