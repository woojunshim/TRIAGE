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
#   -r repressive tendency score file name, default = human_rts.txt
#   -f filtering of non-priority genes, default = False
  
#   E.g. >>> python disc.py -i example_input.txt -o example_file.txt -f True
  
#   --> Run TRIAGE on "example_input.txt" (-i example_input.txt) and output the discordance score as a text file "example_output.txt" (-o example_file.txt), with pseudo-count of 1 and natural log-transformation, using the human RTS table. Only output priority genes (-f True). 

# If any issues are found, please contact to w.shim@uq.edu.au
# Written by Woo Jun (Chris) Shim

import numpy as np
import sys
from optparse import OptionParser
import random
from random import shuffle

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

def read_ref(filename, col, numeric=True):
    """ read a repressive tendency score file """
    results = {}
    file_ = open(filename, 'r')
    for line in file_:
        line = line.strip().split()
        if numeric==True:
            results[line[0]] = float(line[col])
        else:
            results[line[0]] = line[col]
    return results

def read_file_as_list(filename, col=None, numeric_col=None):
    results = []
    temp = open(filename, 'r')
    for i in temp:
        i = i.replace(' ','_')
        i = i.strip().split()
        if i[0].startswith('#'):
            continue
        else:
            results.append([])
            if col==None:
                for no in range(len(i)):
                    c = i[no]
                    if numeric_col!=None:
                        if no in numeric_col:                       
                            c = float(i[no])
                    results[-1].extend([c])
            else:
                for c in col:
                    if numeric_col!=None:
                        if c in numeric_col:
                            value = float(i[c])
                        else:
                            value = i[c]
                    else:
                        value = i[c]
                    results[-1].extend([value])
    return results

def extract_column_from_table_as_dic(input_table, col, threshold=None):
    input_table = check_table(input_table)
    result = {}
    for g in input_table:
        if threshold==None:
            result[g] = input_table[g][col]
        else:
            if input_table[g][col] > threshold:
                result[g] = input_table[g][col]
    return result

def extract_column_from_table_as_list(input_table, col, threshold=0):
    input_table = check_table(input_table)
    result = [[g, float(input_table[g][col])] for g in input_table if float(input_table[g][col]) > threshold]
    return result

def order_list(input_list, col, reverse = False):
    return sorted(input_list, key=lambda x: x[col], reverse=reverse)

def check_table(filename, numeric=True):
    if type(filename)==str:
        filename = read_input(filename, numeric=numeric)
    return filename

def convert_dic_to_list(input_data):
    """ returns a list of lists (each with a key and an associated value) """
    results = [[g, input_data[g]] for g in input_data]

def perform_analysis(exp, ref, pseudo_=1, log_conversion=True, priority=False):  
    """ main function to calcualte the discordance score """    
    if pseudo_!=None:
        exp = add_pseudo(exp, count_=pseudo_)    
    if log_conversion==True:
        for g in exp:
            for c in exp[g]:
                exp[g][c] = np.log(exp[g][c])    
    results = calculate_scores(exp, ref, priority)    
    return results


def calculate_scores(dic1, dic2, priority):
    """ calculate discordance scores, dic1 = expression table, dic2 = RTS table """
    results = {}   
    for g in dic1:
        if priority==True:            
            if not g in priority_genes:                
                continue
            else:
                if g in dic2:  # only considers genes with RTS                
                    results[g] = {}
                    for c in dic1[g]:
                        results[g][c] = dic1[g][c] * dic2[g]
        else:
            if g in dic2:  # only considers genes with RTS                
                results[g] = {}
                for c in dic1[g]:
                    results[g][c] = dic1[g][c] * dic2[g]
    return results

def calculate_bin_rts(input_data, n=100):
    """ calculate average RTS for a set of consecutive genes. The bin size = n """
    values = [input_data[i][1] for i in range(n)]
    results = [np.mean(values)]
    total = len(input_data)-n+1
    for i in range(1, total-n+1):
        values.append(input_data[i+1][1])
        values.pop(0)
        results.append(np.mean(values))   
    return results

def find_closest_bin(input_data, value):
    """ calculates absolute difference between a value and each bin rts
        returns an index """
    dist = [np.abs(value-v) for v in input_data]
    return dist.index(np.min(dist))

def generate_permutations(input_data, iterations=10000):
    """ permutes values of input """
    results = []
    total = len(input_data)
    idx = [i for i in range(total)]
    for no in range(iterations):
        shuffle(idx)        
        results.append([input_data[idx[i]] for i in range(total)])
    return results  

def calculate_empirical_p(input_data, value, idx=None):
    """ returns an empirical p-value (one-sided) given a value """
    if idx == None:
        idx = random.randint(0, len(input_data[0])-1)
    cnt = [1 for i in range(len(input_data)) if input_data[i][idx] > value]
    return (np.sum(cnt)+1) / (len(input_data)+1)

def add_pseudo(input_data, count_):
    """ add a pseudo-count to the input data """
    for gene_ in input_data:
        for c in input_data[gene_]:
            input_data[gene_][c] += float(count_)
    return input_data

def write_file(input_data, output_file):
    """ write out an output as a text file"""
    output_ = open(output_file, 'w')
    if type(input_data)==dict:
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
    else:
        for i in input_data:
            line = i[0]
            for j in range(1,len(i)):
                line += '\t'+str(i[j])
            output_.write(line+'\n')

def get_colnames(input_data):
	return list(input_data[list(input_data.keys())[0]].keys())

if __name__ == '__main__':

    # Command line options   
    parser = OptionParser()
    parser.add_option('-i', '--i', dest='input_file', help='input filename')    
    parser.add_option('-o','--o', dest='output_file', help='output filename')
    parser.add_option('-p', '--p', dest='pseudo', help='pseudo count', default=1)
    parser.add_option('-l', '--l', dest='log_transform', help='log-transformation (True or False)', default=True)
    parser.add_option('-r', '--r', dest='ref_file', help='repressive tendency score (RTS) filename', default='human_rts.txt')
    parser.add_option('-f', '--f', dest='priority_focus', help='filtering: whether to output ONLY priority genes', default=False)
    parser.add_option('-s', '--s', dest='stats', help='Number of top genes for the statistical output', default=0)
    parser.add_option('-b', '--b', dest='bin_size', help='Bin size by the number of genes', default=100)
    parser.add_option('-u', '--u', dest='no_permutation', help='Number of permutations', default=10000)

    options = parser.parse_args()[0]  

    options.log_transform = bool(options.log_transform)
    options.priority_focus = bool(options.priority_focus)
    options.stats = int(options.stats)
    options.bin_size = int(options.bin_size)
    options.no_permutation = int(options.no_permutation)

    if (options.input_file == None) or (options.output_file == None):
        sys.exit('Exiting: At least input filename (-i) and output filename (-o) are required.')
    else:        
        exp = read_input(options.input_file)
        ref = read_ref(options.ref_file, col=1)    # Specify a RTS file to be read
        if options.priority_focus == True:
            temp = read_ref(options.ref_file, col=2, numeric=False)
            priority_genes = set([i for i in temp if temp[i]=='Y'])
        results = perform_analysis(exp, ref, pseudo_=options.pseudo, log_conversion=options.log_transform, priority=options.priority_focus)   
        if options.stats == 0:  # No statistics
            write_file(results, options.output_file)
        else:  # Permutation of expression values
            cols = get_colnames(results)            
            ref_list = read_file_as_list(options.ref_file, col=[0,1], numeric_col=[1])
            output_ = [['#gene','sample','DS','RTS','bin_RTS','exp','bin_exp','empirical_p']]
            for col in cols:            	
                exp_sample = extract_column_from_table_as_dic(exp, col=col, threshold=0)
                dis_sample = extract_column_from_table_as_list(results, col=col, threshold=0)                
                dis_sample = order_list(dis_sample, col=1, reverse=True)                
                genes = [i[0] for i in dis_sample[:options.stats]]                
                temp = [i for i in order_list(ref_list, col=1, reverse=True) if i[0] in exp_sample]                
                bin_rts = calculate_bin_rts(temp, n=options.bin_size)
                gene_bin_no = [find_closest_bin(bin_rts, float(ref[i[0]])) for i in dis_sample[:options.stats]]          
                for j in range(len(genes)):
                    gene = genes[j]
                    idx = gene_bin_no[j]
                    exp_set = [np.log(exp_sample[temp[i][0]]+1) for i in range(idx, idx+options.bin_size)]
                    perm_data = generate_permutations(exp_set, iterations=options.no_permutation)
                    p = calculate_empirical_p(perm_data, np.log(exp_sample[gene]+1), idx=None)
                    output_.append([gene, col, dis_sample[j][1], ref[gene], bin_rts[idx], np.log(exp_sample[gene]+1), np.mean(exp_set), p])
            write_file(output_, options.output_file)

        