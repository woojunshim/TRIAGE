# TRIAGE
# Transcriptional Regulatory Inference Analysis from Gene Expression (TRIAGE)

# We provide an example Python script (written in version 2.7) to calculate the discordance score of genes given their expression readouts. With simplicity of the TRIAGE implementation, users may also want to write their own codes or simply incorporate TRIAGE as a part of their analysis.

# "disc.py" is the main script to run the TRIAGE analysis. It requires (i) input expression matrix (e.g. "example_input.txt") and (ii) repressive tendency score (RTS) file (e.g. "human_rts.txt") both of which are plain text files. You can find the both files in the repository. 

# Input expression file is a matrix that defines expression values of genes across samples (where rows are genes and columns are samples). See "example_input.txt" below (which is also included in the repository). This toy example is RNA-seq data (RPKM normalised) for 3 selected Roadmap samples.  

# <example_input.txt>

#     Gene    Blood(T_helper_naive_cells) Brain(Germinal_Matrix)  Heart(Right_Atrium) 

#     A1BG    3.065   1.983   4.45 

#     A1CF    0.073   0.0 0.0 

#     A2M 0.402   24.878  137.612 

#     A2ML1   0.002   0.394   0.166 

#     ...

# The output of TRIAGE is weighted expression readouts (i.e. discordance score) by the gene-specific RTS. 


# <example_output.txt>. 

#         Blood(T_helper_naive_cells) Brain(Germinal_Matrix)  Heart(Right_Atrium) 
  
#     A1BG    0.008059178900858537    0.006280681756778532    0.009744106976712222 

#     A1CF    0.00027847821321769683  0.0 0.0 

#     A2M 7.800106276867326e-05   0.0007510159364827811   0.0011384327244910546 

#     A2ML1   8.465243140668874e-07   0.00014073863700441038  6.506919860136177e-05 

#     ...



# You may like to use an alternative RTS table for mouse datasets ("mouse_rts_mapped.txt"). Note that the mouse RTS values were obtained by directly mapping genes between human and mouse data. The mouse data currently only covers protein-coding genes.

# Finally, users can modify a pseudo-count (default = 1) or whether to log-transform the expression value (default = True) by specifying parameters.

# Parameters
#   -i input file name (required)
#   -o output file name (required)
#   -p pseudo-count to be added to the input data, default = 1
#   -l natural log-transformation of the input data, default = True    
  
#   E.g. >>> python disc.py -i example_input.txt -o example_file.txt 
  
#   --> Run TRIAGE on "example_input.txt" and output the discordance score as a text file "example_output.txt", with pseudo-count of 1 and natural log-transformation 

# If any issues are found, please contact to w.shim@uq.edu.au

import numpy as np
import sys
from optparse import OptionParser


def read_input(filename): 
    """ read an input file """
    temp = []
    results = {}        
    file_ = open(filename, 'r')
    for line in file_:
        line = line.strip().split()            
        temp.append(line)
    for t1 in temp[1:]:
        results[t1[0]] = {}
        for idx in range(1, len(temp[0])):
            t2 = temp[0][idx]
            results[t1[0]][t2] = float(t1[idx])
    return results

def read_ref(filename):
    """ read a repressive tendency score file """
    results = {}
    file_ = open(filename, 'r')
    for line in file_:
        line = line.strip().split()
        results[line[0]] = float(line[1])
    return results

def perform_analysis(exp, ref, pseudo_=1, log_conversion=True):  
    """ main function to calcualte the discordance score """    
    if pseudo_!=None:
        exp = add_pseudo(exp, count_=pseudo_)    
    if log_conversion==True:
        for g in exp:
            for c in exp[g]:
                exp[g][c] = np.log(exp[g][c])    
    results = calculate_scores(exp, ref)    
    return results


def calculate_scores(dic1, dic2):
    """ calculate discordance scores, dic1 = expression table, dic2 = RTS table """
    results = {}
    for g in dic1:
        if g in dic2:  # only considers genes with RTS
            results[g] = {}
            for c in dic1[g]:
                results[g][c] = dic1[g][c] * dic2[g]
    return results

def add_pseudo(input_data, count_):
    """ add a pseudo-count to the input data """
    for gene_ in input_data:
        for c in input_data[gene_]:
            input_data[gene_][c] += float(count_)
    return input_data


def write_file(input_data, output_file):
    """ write out an output as a text file"""
    output_ = open(output_file, 'w')
    rownames = []
    for row in input_data:
        colnames = list(input_data[row].keys())
        rownames.append(row)
    colnames.sort()
    rownames.sort()
    first_line = ''
    for col in colnames:
        first_line += '\t'+col
    output_.write(first_line+'\n')
    for row in rownames:
        line = str(row)
        for col in colnames:
            line += '\t'+str(input_data[row][col])
        output_.write(line+'\n')

if __name__ == '__main__':

    # Command line options   
    parser = OptionParser()
    parser.add_option('-i', '--i', dest='input_file', help='input filename')    
    parser.add_option('-o','--o', dest='output_file', help='output filename')
    parser.add_option('-p', '--p', dest='pseudo', help='pseudo count', default=1)
    parser.add_option('-l', '--l', dest='log_transform', help='log-transformation (True or False)', default=True)

    options = parser.parse_args()[0]  

    if options.log_transform == 'False':
        options.log_transform = False
    else:
        options.log_transform = True  

    if (options.input_file == None) or (options.output_file == None):
        sys.exit('Exiting: Input filename (-i) and output filename (-o) are required.')
    else:        
        exp = read_input(options.input_file)
        ref = read_ref('human_rts.txt')    # Specify a RTS file to be read        
        results = perform_analysis(exp, ref, pseudo_=options.pseudo, log_conversion=options.log_transform)        
        write_file(results, options.output_file)
        