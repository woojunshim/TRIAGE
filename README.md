## Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

### 1. Introduction

TRIAGE is a simple computational method that transforms any expression readouts of genes (e.g. RNA-seq, CAGE-seq or H3K36me3 tag density) into a new metric called the discordance score. In essense, the method introduces gene-specific weights, which indicate degree of association with broad H3K27me3 domains observed across diverse tissue and cell types, to prioritize genes with cell-type specific regulatory function. TRIAGE only requires the gene expression readout as the input. It does not require ChIP-seq data for the cell-type of interest or any external reference points for the analysis. 

### 2. Installation

We provide a Python script (can be run on both versions 2 and 3) as a bare-bones implementation of TRIAGE (suitable to run with commandlines), and a R script (version 3) for more visualisation options, which is also available [here](http://bioinf.scmb.uq.edu.au/adhoc/). 

TRIAGE is fast and highly scalable to a large number of datasets. For instance, it can process a thousand single-cell transcriptomes in less than 2 minutes with a personal computer (2.2 GHz Intel core i7, 16 GB RAM). 


To download the source code, clone the git repository. 

	> git clone https://github.com/woojunshim/TRIAGE.git 

### 3. Python implementation

We provide an example Python script to calculate the discordance score of genes given their expression readouts. With simplicity of the TRIAGE implementation, users may also want to write their own codes or incorporate TRIAGE as a part of their analysis.

Dependency: numpy, scipy. Please install the both before running the script. To check whether these have been installed, you can view all installed Python modules in your machine. 

	> pip list

To install a module (e.g. numpy)

	> pip install numpy

"disc.py" is the main script to run the TRIAGE analysis. It requires (i) input expression matrix (e.g. "example_input.txt") and (ii) repressive tendency score (RTS) file (e.g. "human_rts.txt") both of which are tab-delimited text files. You can find the both files in this repository. Simply, the script reads in the both files and calculates the discordance score after value conversion (e.g. natural log-transformation) if specified. Please note that the script does NOT normalise the input data. While the normalisation across samples is not a requirement for TRIAGE, if one concerns, this should be done before running the script.

Input expression file is a matrix that defines expression values of genes across samples (where rows are genes and columns are samples). See an extract from "example_input.txt" below. The first line should be "column names" to define each column uniquely. The first column is the gene symbol while the rest are values (i.e. gene expression) for samples, as shown below. Included example is RNA-seq data (RPKM normalised) for 3 selected Roadmap samples.  

<example_input.txt>

	Gene	Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium) 

	A1BG	3.065	1.983	4.45 

	A1CF	0.073	0.0	0.0 

	A2M	0.402	24.878	137.612 

	A2ML1	0.002	0.394	0.166 

	...

The output of TRIAGE is the discordance score. 


<example_output.txt> 

		Blood(T_helper_naive_cells)	Brain(Germinal_Matrix)	Heart(Right_Atrium) 
  
	A1BG	0.008059178900858537	0.006280681756778532	0.009744106976712222 

	A1CF	0.00027847821321769683	0.0	0.0 

	A2M	7.800106276867326e-05	0.0007510159364827811	0.0011384327244910546 

	A2ML1	8.465243140668874e-07	0.00014073863700441038	6.506919860136177e-05 

	...

Implementation of TRIAGE is fast and highly scalable. Computational time required to analyse the example dataset (18,707 genes * 3 samples) is approximately 0.3 second (run on a personal Mac with 2.2 GHz Intel core i7, 16 GB RAM). 

You may want to use an alternative RTS table for a different species (e.g. "mouse_rts_mapped.txt" for mouse). Note that the mouse RTS values were obtained by directly mapping genes between human and mouse data. The mouse data currently only covers protein-coding genes.

Finally, users can modify a pseudo-count (default = 1) or whether to log-transform the expression value (default = True) by specifying parameters. Statistical test for the discordance score is also available by specifying an additional -s option (i.e. the number of top genes for the statistical output). See below for adjustable parameters. 

Parameters

  -i input file name (required) ; input expression table file 
  
  -o output file name (required) ; output discordance score table file
  
  -p pseudo-count to be added to the input expresion data, default = 1 ; added to to avoid a negative output for log-transformation
  
  -l natural log-transformation of the input expression data, default = True 
  
  -r repressive tendency score file name, default = human_rts.txt ; specify this to use an alternative RTS file
  
  -f filtering of non-priority genes, default = False ; True if like to consider only the 1,359 prioirty genes
  
  -s number of top genes for the statistical output, default = 0 (i.e. no stat output) ; if set non-zero, a statistical test is performed for each gene
  
  -b bin size by the number of genes, default = 100 (only required for stat output) ; number of genes in a bin assumed to have similar RTS or experssion values, used to form the basis for the statistical background by permutations. 
  
  -u number of permutations, default = 10000 (only required for stat output) ; number of permutations to generate a distribution from a set of genes sharing similar RTS or expression values
  
  
  
  	Example1 > python3 disc.py -i input.txt -o output.txt -f True
  
  --> Run TRIAGE on "input.txt" and output the discordance score as a text file "output.txt", with pseudo-count of 1 and natural log-transformation, using the human RTS table. Only output priority genes (-f True). 
  
  	Example2 > python3 disc.py -i input.txt -o output.txt -r mouse_rts_mapped.txt -s 100 

  --> Run TRIAGE as above except (i) focus on all expressed genes in the input data, (ii) use a RTS for mouse genes, (iii) perform a permutation test for the statistical significance for top 100 genes by the discordance score.
  
### 4. Q & A 

#### My input file causea an error. Does it accept CSV file as the input? 

Answer: Please make sure that the text file is tab-delimited with headers in the first line and the gene symbol in the first column (See Python implementation). Current Python script only accepts a tab-delimited text file as the input. 

 	


### 5. R implementation

#### How to install:
    1) install R 3.5.0+
    2) install r packages: shiny, scales, zip, dplyr, ggplot2, reshape2, DESeq, FedData

#### How to use the tool:

##### 1. Run script
   *a) navigate to the downloaded folder ("adhoc" for custom loaded dataset, "cardiac" for cardiac scRNA temporal dataset)*
   
   *b) open the main R script ("app.r") in RStudio and click "Run App"*
   
   *c) follow the short guide at the beginning of the popup page*

##### 2. Analyse Data
   *a) Everytime the criteria are changed, please press "Refresh" to get results updated*
   
   *b) Use "Reset" + "Refresh" button to clear all criteria changes*
   
   *c) You need to navigate through tabs before downloading results (for adhoc tool)*


### 6. Contact

If any issues are found, please contact Woo Jun Shim (w.shim@uq.edu.au) for Python or Jun Xu (jun.xu@uq.edu.au) for R scripts. 

