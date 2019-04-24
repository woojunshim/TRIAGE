# TRIAGE
Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

TRIAGE is a method that converts any expression readouts of genes (e.g. RNA-seq, CAGE-seq or H3K36me3 tag density) to a new metric called the discordance score. In essense, the method introduces gene-specific weights, which indicate degree of association of the gene with broad H3K27me3 domains across diverse tissue and cell types, to prioritize genes with cell-type specific regulatory function. It's fast and scalable to large datasets including single-cell transcriptomes. Current implementation is available in Python (version 2.7) or R. We also provide a web interface where users can run ad-hoc analysis with various options for visualisation and filtering (http://bioinf.scmb.uq.edu.au/adhoc/) 

To download the source code, clone the git repository (git clone https://github.com/woojunshim/TRIAGE.git) 

More information is available under folders.

Python script was written by Woo Jun Shim (w.shim@uq.edu.au). R implementation and the web interface were written by Jun Xu (jun.xu@uq.edu.au).
