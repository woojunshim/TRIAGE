# TRIAGE
## Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

TRIAGE is a simple computational method that transforms any expression readouts of genes (e.g. RNA-seq, CAGE-seq or H3K36me3 tag density) into a new metric called the discordance score. In essense, the method introduces gene-specific weights, which indicate degree of association of the gene with broad H3K27me3 domains observed across diverse tissue and cell types, to prioritize genes with cell-type specific regulatory function. TRIAGE is fast and highly scalable to large datasets. For instance, it can process 1,000 single-cell transcriptomes in less than 2 minutes with a personal computer (2.4GHz, . Current implementation is available in Python (can be run in both versions 2 and 3) and R (version 3). 

Python script offers a bare-bones implementation of TRIAGE (suitable to run with commandlines) while R implementation provides more options for user interactions. We also provide a web interface where users can run ad-hoc analysis with various options for visualisation and filtering (http://bioinf.scmb.uq.edu.au/adhoc/). This web interface is equivalent to the implementation of TRIAGE in R.  

To download the source code, clone the git repository. 

> git clone https://github.com/woojunshim/TRIAGE.git 

Details of how to run an example dataset and the implementation are available under 'python' folder. 

Python script was written by Woo Jun Shim (w.shim@uq.edu.au). 
R implementation and the web interface were written by Jun Xu (jun.xu@uq.edu.au).
