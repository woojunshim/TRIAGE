
### Example Python script to calculate the discordance score ###

""" It requires input expression data in text file format. 
    The human repressive tendency score file (repressive_hg19.txt) must be also in the same directory where the script is run.

    Input file is a text file that defines expression values of genes across samples where rows are genes and columns are samples. 
    For example, see example_input.txt 

    Parameters
     -i input file name (required)
     -o output file name (required)
     -p pesudo-count to be added to the input data, default = 1
     -l log-transformation of the input data, default = True

     E.g. >>> python disc.py -i example_input.txt -o output_file.txt 
        Run the analysis with the human with the log-transformation of the input expression data, with a pseudocount of 1 """


import numpy as np
import sys
from optparse import OptionParser


def read_input(filename): 
    """ reads input files """
    temp = []
    results = {}        
    file_ = open(filename, 'r')
    for line in file_:
        line = line.strip().split()            
        temp.append(line)
    for t1 in temp[1:]:
        results[t1[0]] = {}
        for idx in range(len(temp[0])):
            t2 = temp[0][idx]
            results[t1[0]][t2] = float(t1[idx+1])
    return results

def read_ref(filename):
    """ read repressive tendency score file """
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
    for gene_ in input_data:   
        for c in input_data[gene_]:     
            input_data[gene_][c] += float(count_)
    return input_data


def write_file(input_data, output_file):
    """ writes out an output as a text file"""
    output_ = open(output_file, 'w')
    rownames = []
    for row in input_data:
        colnames = input_data[row].keys()
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
    parser.add_option('-p', '--p', dest='pseudo', help='pesudo count', default=1)
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
        ref = read_ref('repressive_hg19.txt')        
        results = perform_analysis(exp, ref, pseudo_=options.pseudo, log_conversion=options.log_transform)
        
        write_file(results, options.output_file)
   


