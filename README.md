# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

TRIAGE is a method that converts any expression readouts of genes (e.g. RNA-seq quantification, CAGE-seq or H3K36me3 tag density) to a new metric called the discordance score. In essense, the method introduces gene-specific weights to prioritize genes with cell-type specific regulatory function. It's fast and readily scalable to large datasets including multiple single-cell transcriptomes. Current implementation is available as a Python script (version 2.7). 

To run the script, first the source code and data file need to be downloaded. You can easily do this by cloning the git repository (git clone https://github.com/woojunshim/TRIAGE.git) 

disc.py is the main script to run the analysis and repressive_hg19.txt is the repressive tendency score (RTS) table. These files must be included in the same directory. 

To run the analysis, simply type, python disc.py -i (your input file) -o (output file), under the terminal. 

If any issues are found, please contact to w.shim@uq.edu.au
