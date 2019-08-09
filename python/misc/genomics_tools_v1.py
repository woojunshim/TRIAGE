''' A RANGE OF UTILITIES USEFUL FOR GENOMIC ANALYSIS.
DEVELOPED AND WRITTEN BY WOO JUN SHIM, 2016 '''

import numpy as np
import sys
import os
import scipy
from scipy import stats
import statistics
import collections
from collections import Counter
import random


###################
## GENERAL TOOLS ##
###################

# Designed and Written by Woo Jun Shim
# Nov. 2016

def read_csv(filename, delimiter=',', quote='"', fill_space=True):
    temp = open(filename, 'r')
    results = []
    for i in temp:
        i = i.strip()      
        results.append([])
        i = i.split(delimiter)         
        for m in i:
            m = m.strip(quote)            
            if fill_space==True:
                m = m.replace(' ', '_')           
            results[-1].extend([m])
    return results


def write_file_dic(input_, filename):
    ### EXPORT A DICTIONARY OF LISTS
    output_ = open(filename, 'w')
    for i in input_:
        line = str(i)
        for m in input_[i]:
            line += '\t'+str(m)
        output_.write(line+'\n')
    output_.close()


def read_file(filename, features, rowname='', header_=False):
    # GENERAL FILE LOADING FUNCTION #
    # can read any types of text files #
    ''' Reads a text file and create a data entry (OUTPUT: a list of lists if rowname not specified, a dictionary if rowname specified).
    filename = a string of input filename
    rowname = a string of column name (or column index, starting from 0) with which the data entry is sorted (if specified, output = a dictionary)
    numerical = a boolean to indicate whether the input data is numerical
    features = a list of strings of column names (or column indexes) indexes to extract
    NOTE: the first line of the input text file must be specified
    '''

    header = False
    features_idx = []
    rowname_idx = 0
    head = []
    if rowname != '':
        results = {}
    else:
        results = []
    data_ = open(filename, 'r')
    for line in data_:
        line = line.replace(' ', '_')
        line = line.strip().split()
        if len(line) != 0:
            if line[0].startswith('#'):
                head = line
            if header == False:
                cnt = 0
                header = {}
                for m in line:
                    header[m] = cnt
                    cnt += 1
                if str(features[0]).isdigit():
                    for m in features:
                        features_idx.append(int(m))
                else:
                    for m in features:
                        features_idx.append(int(header[m]))
                if rowname != '':
                    if rowname.isdigit():
                        rowname_idx = int(rowname)
                    else:
                        rowname_idx = int(header[rowname])
                header = True
            else:
                temp = []
                if rowname != '':
                    if line[rowname_idx] not in results:
                        results[line[rowname_idx]] = []
                    for m in features_idx:
                        temp.append(line[m])
                    results[line[rowname_idx]]=temp
                else:
                    for m in features_idx:
                        if int(m) < len(line):
                            temp.append(line[m])

                    results.extend(temp)

    if header_==False:
        return results
    else:
        return results, head

def read_file1(filename, fill_gap=True, cols=None, sort_by=None, reverse_=True, numerical=False):
    ### Reads in a whole text file and returns a list
    ### header = a boolean to indicate whether the first row is to be recoreded as the header (or the line startswith #)
    ### cols = a list of column indexes to take
    ### ignore_first = a boolean to indicate whether not to take the first line
    results = []
    input_ = open(filename, 'r')    
    for line in input_:
        if fill_gap==True:
            line = line.replace(' ','_')
        line = line.strip().split()        
        if len(line) != 0:            
            if not line[0].startswith('#'):
                if cols==None:
                    if sort_by != None:
                        line[sort_by] = float(line[sort_by])
                    results.append(line)
                else:
                    temp = []
                    for m in cols:   
                        pp = line[int(m)]
                        if numerical==False:
                            temp.append(pp)
                        else:
                            temp.append(float(pp))
                    results.append(temp)   
    if sort_by != None:
        results = sort_(results, idx=sort_by, reverse_=reverse_)
    return results

def read_file_dic(filename, col_idx=None, id_idx=0, numerical=True):
    ### READ A TEXT FILE AND RETURNS A DICTIONARY WITH A SINGLE VALUE DEFINED BY COL IDX
    temp = open(filename, 'r')
    results = {}
    for line in temp:
        line = line.strip().split()
        if '#' not in line[0]:
            if col_idx != None:
                if numerical == True:
                    value = float(line[col_idx])
                else:
                    value = line[col_idx]
            else:
                value = []
                for m in range(len(line)):
                    if m!=id_idx:
                        if numerical==False:
                            value.append(line[m])        
                        else:
                            value.append(float(line[m]))
            results[line[id_idx]] = value
    return results



def read_file_items(filename, col=None, numeric=False):
    ### reads in a text file which is composed of multiple line (each line with only one item)
    ### output is a list of items
    results = []
    input__ = open(filename, 'r')
    for line in input__:
        line = line.strip().split()        
        if len(line) != 0:
            if not '#' in line[0]:
                if col==None:
                    results.append(line[0])
                else:
                    results.append(line[col])
                if numeric==True:
                    results[-1] = float(results[-1])
    return results

def write_file_items(input_, filename):
    ### writes out a text file each line of which is an item in input_
    output_ = open(filename, 'w')
    for i in input_:
        output_.write(str(i)+'\n')
    output_.close()


def write_file(data_, output_file, tag='', header=False):
    ''' Export data entry(or entries) to a text file
    data_ = a list of data entry(s)
    output_file = a string of output filename
    header = a dictionary of descriptions for each key in data (for only dictionary type)
    tag = a string to be a prefix of the each line (This can be used to indicate different data entries as 'tag')
    '''

    output = open(output_file, 'w')
    if type(data_) == dict:
        order_ = []
        for item in data_:
            order_.append(item)
        order_.sort()
        for key in order_:            
            if header != False:
                output.write('#####' + str(header[key]) + '\n')

            if type(data_[key]) == list:
                for m in data_[key]:
                    output.write(tag + str(key))
                    for n in m:
                        output.write('\t'+str(n))
                    output.write('\n')


            else:
                output.write(str(key)+'\t'+str(data_[key])+'\n')
    elif type(data_) == list:
        for item in data_:
            if type(item) != list:
                output.write(tag+str(item))
            elif type(item) == list and len(item) > 1:
                output.write(tag + str(item[0]))
                for n in range(1, len(item)):
                    output.write('\t'+str(item[n]))
            elif type(item) == list and len(item) == 1:
                output.write(tag + str(item[0]))
            else:
                output.write(item+'\n')
            output.write('\n')



def intersection(data1, data2):
    ''' takes a pair of lists and returns a set of intersected items
    data1 or data2 = a list of elements
    '''
    a = collections.Counter(data1)
    b = collections.Counter(data2)
    return list((a & b).elements())


def tss(data_entry, chr_idx, position_list, strand_idx, id_idx):
    ''' Takes a data entry instance and convert into a TSS instance (dictionary)
    Note. Data entry for input MUST be a list (not dictionary).
    Note. Index starts from 0

    chr_idx = an integer indicating column number of chromosome ID
    position_idx = an list of integers indicating column number of a TSS position coordinate (i.e. [position1 (for +strand), position 2(for -strand)])
    strand_idx = an integer indicating column number of strand direction
    '''
    result = {}
    chr_idx = int(chr_idx)
    strand_idx = int(strand_idx)
    id_idx = int(id_idx)

    if type(data_entry) != list:
        sys.exit('Error: Input data entry must be a list.. Exiting..')
    for item in data_entry:
        if item[chr_idx] not in result:
            result[item[chr_idx]] = []
        if item[strand_idx] == '+':
            point = 0
        elif item[strand_idx] == '-':
            point = 1

        result[item[chr_idx]].append([int(item[position_list[0]]), int(item[position_list[1]]), int(item[position_list[point]]), item[strand_idx], item[id_idx]])

    for chr in result:
        result[chr] = sort_(result[chr], 2, reverse_=False)

    return result




def sort_(input_list, idx, reverse_ = False, numerical=True):
    ''' sort a list of entries based on value of items at a specified index (idx)
    and returns the sorted list'''
    if numerical==True:
        for no in range(len(input_list)):
            input_list[no][idx] = float(input_list[no][idx])

    return sorted(input_list, key=lambda x: x[idx], reverse=reverse_)

def union(a,b):
    return list(set(a) | set(b))

def intersect(a,b):
    return list(set(a) & set(b))


def assign_gene(data_entry, tss_entry, min_distance=0, max_distance=1000000, centre=False, width=False, gene_include=True):
    ''' Takes a data entry instance and a TSS instance. Assign each item to the closest TSS.
    Binary tree search-based method is used (recursively) until the number of TSS < 10. Then distance to each TSS is calculated.
    Returns a list of [unique ID (domain), assigned_gene, distance, TSSs within the domain, (width of domain)]
    Note. Input data entry MUST be in the form of [chr, start, end, unique ID]

    min_distance = an integer indicating the minimum distance required
    max_distance = an integer indicating the maximum distance allowed
    centre = a boolean to indicate whether the centre position of domain should be considered.
      If so, this will override the distance to both endpoints of domains
      Default = False, hence take the min. distance to either endpoint of the domain
    width = a boolean to indicate whether to include width of enquiry domains
    gene_include = a boolean to indicate whether to assign a domain to a gene that overlaps that regardless of the distance
    '''
    results = []
    if max_distance==None:
        max_distance=999999999999  # If max distance is not specified, there's virtually no limit
    if type(data_entry) != list:
        sys.exit('Error: Input data entry must be a list.. Exiting..')

    for item in data_entry:
        chrom = item[0]
        if chrom in tss_entry:
            temp_tss = __find_tss(item, tss_entry[chrom], centre)

            temp = __assign_gene(item, temp_tss, min_distance, max_distance, 9999999999, centre, width, gene_include)
            results.append(temp)


    return results

def __find_tss(item, tss_entries, centre):
    ''' recursively find a small set (i.e. ~10) of TSSs closest to a given interval as defined by user'''


    if len(tss_entries) > 10:
        mid_point = len(tss_entries) / 2
        mid_tss = tss_entries[mid_point]
        dist1 = int(item[1]) - int(mid_tss[2])  # d. btw left endpoint & TSS point
        dist2 = int(item[2]) - int(mid_tss[2])  # d. btw right endpoint & TSS point
        dist3 = int(item[1]) + (int(item[2]) - int(item[1])) / 2 - int(mid_tss[2])  # d. btw centre of the domain & TSS point

        if centre == False and dist1 >= 0:  # If the domain is located 'right' of the TSS
            new_entries = tss_entries[mid_point:]

        elif centre == False and dist1 < 0:  # 'left' of the TSS
            new_entries = tss_entries[:mid_point]

        elif centre == True and dist3 >= 0:
            new_entries = tss_entries[mid_point:]

        elif centre == True and dist3 < 0:
            new_entries = tss_entries[:mid_point]


        return __find_tss(item, new_entries, centre)

    else:

        return tss_entries

def __assign_gene(item, tss_entries, min_distance, max_distance, current_distance, centre, width=False, gene_include=True):
        gene_name = ''
        gene_included = []  # A list of genes that overlaps a given interval
        sign_ = 1
        for m in tss_entries:
            if gene_include == True:
                gene_start = int(m[0])
                gene_end = int(m[1])
                if (int(item[1]) >= gene_start and int(item[2]) <= gene_end) or (int(item[2]) >= gene_start and int(item[2]) <= gene_end) or (int(item[1]) >= gene_start and int(item[1]) <= gene_end) or (gene_end >= int(item[1]) and gene_end <= int(item[2])) or (gene_start >= int(item[1]) and gene_start <= int(item[2])):  # If the domain is overlapped by the gene

                    dist1 = int(item[1]) - int(m[2])
                    dist2 = int(item[2]) - int(m[2])
                    current_distance = min(abs(dist1), abs(dist2))
                    tt= min([int(item[1]) - int(m[2]), int(item[2]) - int(m[2])])
                    if tt < 0:
                        sign_ = -1
                    if tt >= 0:
                        sign_ = 1
                    gene_name = m[-1]

                    gene_included.append(gene_name)
                    break

            if centre == True:
                dist3 = abs(int(item[1]) + (int(item[2]) - int(item[1])) / 2 - m[2])
                if dist3 < current_distance and dist3 >= int(min_distance) and dist3 <= int(max_distance):
                    current_distance = dist3
                    tt= int(item[1]) + (int(item[2]) - int(item[1])) / 2 - m[2]
                    if tt < 0:
                        sign_ = -1
                    if tt >= 0:
                        sign_ = 1
                    gene_name = m[-1]


            else:
                dist1 = int(item[1]) - int(m[2])
                dist2 = int(item[2]) - int(m[2])
                dist_ = min(abs(dist1), abs(dist2))
                if dist_ < current_distance and dist_ >= int(min_distance) and dist_ <= int(max_distance):
                    current_distance = dist_
                    tt= min([int(item[1]) - int(m[2]),int(item[2]) - int(m[2])])
                    if tt < 0:
                        sign_ = -1
                    if tt >= 0:
                        sign_ = 1
                    gene_name = m[-1]



        current_distance = current_distance * sign_
        if width==False:
            return [item[-1], gene_name, current_distance, gene_included]
        else:
            width_ = int(item[2]) - int(item[1])
            return [item[-1], gene_name, current_distance, gene_included, width_]

def find_entry(data_entry, feature, idx_):
    ''' takes a data entry (in a list of lists format, see read_file) and returns elements with a feature
    feature = item(s) of interest (str, int or list)  **MODIFY**
    '''
    if type(data_entry) != list:
        sys.exit('Error: Input data entry must be a list.. Exiting..')
    results = []
    for item in data_entry:
        
        if feature in item[idx_]:
            results.append(item)
    return results

def extract_element(data_entry, feature_idx):
    ''' takes a data entry (in a list format) and extract only element of interest (defined by index)
    feature_idx = index (e.g. 1, 1:2 etc)
    '''
    if type(data_entry) != list:
        sys.exit('Error: Input data entry must be a list.. Exiting..')
    results = []
    for item in data_entry:
        results.extend([item[feature_idx]])
    return results

def extract_elements(data_entry, feature_idx):
    ''' same as extract_element except this can take multiple items
    '''
    results = []
    for item in data_entry:
        results.append([])
        for m in feature_idx:
            results[-1].extend([item[m]])
    return results


def combine_data(data_entries):
    ''' takes a list of data entries (in a list format) and combine them into a single data entry
    data_entries = a list of data entries
    '''
    results = []
    for item in data_entries:
        results.append(item)
    return results


def empirical_pvalue(input_value, data, lower=False):
    ''' calculate an empirical p-value given a set of data
    input_value = an query value
    data = a list of background numerical values
    lower = a boolean to indicate whether upper-tail or lower-tail (default=False, meaning upper-tail)
    '''
    cnt = 0
    total_items = len(data)
    for item in data:
        if lower==True and float(input_value) > float(item):
            cnt += 1
        elif lower==False and float(input_value) < float(item):
            cnt += 1
    return float(cnt+1) / float(total_items+1)

def width(query_id, data_entry):
    ''' calculates width of an BED interval (identified by 'query_id') found in data_entry
    query_id = a keyword to locate the interval of interest (unique ID for an entry) (e.g. Rank_1, ENSG**)
    Note. data_entry must be a list ([chr, start, end, gene_id])
    '''
    for item in data_entry:
        if query_id in item:
            return int(item[2]) - int(item[1])

def cbind(data_entry1, data_entry2):
    ''' joins two data entries (equivalent to cbind function in R).
    Size of the both data entries must be the same.
    '''
    for i in range(len(data_entry1)):
        if type(data_entry2[0]) == list:
            for j in range(len(data_entry2[i])):
                data_entry1[i].extend([data_entry2[i][j]])
        else:

            data_entry1[i].extend([data_entry2[i]])

    return data_entry1

def transpose(data_):
    ''' transpose a given dictionary (equivalent to 't' function in R)
    e.g. data2 = transpose(data1)
    data1[row][column] = ...
    data2[column][row] = ...
    '''
    results = {}
    for row in data_:
        for col in data_[row]:
            if col not in results:
                results[col] = {}
            if row not in results[col]:
                results[col][row] = data_[row][col]
    return results

def convert_id(enquiry, ref):
    ''' finds an enquiry ID and convert it to the corresponding one
    enquiry = a string of enquiry item
    ref = a conversion reference data in list format (e.g. [["gene1", "ENSGxxx"], [], ..,])
    '''
    for item in ref:
        if enquiry in item:
            idx = item.index(enquiry)
            if idx==0:
                idx_ = 1
            elif idx==1:
                idx_ = 0
            return item[idx_]



def write_table(input_data, output_file):
    ''' export a table (in a dictionary format) to a text file
    Note. this dictionary must have a key 'colnames' to define column names
    '''
    output = open(output_file, 'w')
    first_line= ''
    if 'colnames' in input_data:
        for item in input_data['colnames']:
            first_line += '\t'+item
        output.write(first_line+'\n')
    for row in input_data:

        if row =='colnames':
            continue
        else:
            line = row
            for m in input_data[row]:
                line += '\t'+str(m)
            output.write(line+'\n')

def write_table1(input_data, output_file):
    ''' writes out a table text file 
    takes a dictionary of dictionaries (e.g. input_data[row1][col1]=xx, input_data[row1][col2]=xx, ..) 
    '''
    output_ = open(output_file, 'w')
    rownames = []
    for row in input_data:
        colnames = input_data[row].keys()
        rownames.append(row)
    colnames.sort()
    rownames.sort()
    first_line = ''
    for col in colnames:
        first_line += '\t'+str(col)
    output_.write(first_line+'\n')
    for row in rownames:
        line = str(row)
        for col in colnames:
            line += '\t'+str(input_data[row][col])
        output_.write(line+'\n')



def read_table(input_file):
    ''' import a text file into a dictionary
    returns a dictionary with colnames (i.e. items in the first line) as the first keys &
    rownames as the second keys
    '''
    output = {}
    list_ = []
    input_ = open(input_file, 'r')
    header = False
    for line in input_:
        line = line.strip().split()
        if header == False:
            for item in line:
                output[item] = {}
                list_.append(item)
            header = True
        else:
            gene = line[0]
            for cell_no in range(len(list_)):
                cell = list_[cell_no]
                output[cell][gene] = line[cell_no+1]
    return output



def remove_column(data_, column_idx):
    ''' remove a user-defined column (in a list)
    '''
    results = []
    col = int(column_idx)
    for item in data_:
        temp = []
        for m in range(len(item)):
            if m == col:
                continue
            else:
                temp.append(item[m])
        results.append(temp)
    return results

def overlap_table(data_, row_idx, col_idx):
    ''' create a (overlap binary) table (as 1 overlapped and 0 non-overlapped)
    data_ = a list of items
    row_idx = an integer indicating index of elements in the list to become the row
    col_idx = an integer indicating index of elements in the list to become the column

    '''
    table = {}
    table['colnames'] = []
    rownames = set()
    for item in data_:
        if item[col_idx] not in table['colnames']:
            table['colnames'].append(item[col_idx])
        if item[row_idx] not in table:
            table[item[row_idx]] = []
            rownames.add(item[row_idx])
    for name in rownames:
        for i in range(len(table['colnames'])):
            table[name].append(0)
    for item in data_:
        col = item[col_idx]
        row = item[row_idx]
        idx = table['colnames'].index(col)
        table[row][idx] = 1
    return table

def select(data_, col_idx, min_, max_):
    ''' takes a list and return elements that satisfies conditions set by min_ and max_
    data_ = a list of elements
    col_idx = an integer indicating column index of interest
    min_ = a numerical value specifying the minimum value allowed (None = no limit)
    max_ = a numerical value specifying the maximum value allowed (None = no limit)
    '''
    results = []
    for item in data_:
        if min_==None:
            min_ = -99999999
        if max_==None:
            max_ = 99999999
        if (float(item[col_idx]) >= min_) and (float(item[col_idx]) <= max_):
            results.append(item)
    return results

def mean_(data_):
    ''' return a mean of input data
    data_ = a list of numerical values
    '''
    return np.mean(data_)

def standard1(input_, mean_, log_=True):
    ''' standardise the input value by dividing it by the population mean
    input_ = a numberical value of the input data
    mean_ = a numberical value of the mean
    log_ = a boolean value to indicate whehther the value should be log-converted
    '''
    result = float(input_) / float(mean_)
    if log_ == True:
        return np.log10(result)
    else:
        return result

def unique(input_):
    ''' takes a list of inputs and returns lists of unique elements and their smallest index numbers
    '''
    indexes = []
    elements = []
    for no in range(len(input_)):
        if input_[no] not in elements:
            elements.append(input_[no])
            indexes.append(no)
    return indexes, elements

def width_analysis(bed_file, tss_, inter_background, output_, inter_threshold=0.3, intra_threshold=0.3, remove_ncrna=False, dominant=False, list_all_genes=False):
    ''' perform all steps required for the width analysis
    bed_file = a BED text file (input data)
    tss_ = a RefSeq TSS text file downloaded from NCBI
    inter_background = a background (mean, variance) for inter-cell type analysis (e.g. 'summary_means_variances__.txt')
    output_ = a output filename
    inter- and intra-threshold = a value of converted width (log10(width/mean)) acting as a threshold for selection of domains
    remove_ncrna = a boolean to indicate whether to remove any ncRNA genes from intra-cell type analysis
    dominant = a boolean to indicate whether only dominant peaks should be included in the analysis
    list_all_genes = a boolean to indicate whether to export a list of all tagged genes. File name is output_+'all_genes.txt'

    '''
    tss_data = read_file(tss_,[1, 3, 4, 2, 5, 0])
    tss_human19 = tss(tss_data, chr_idx=0, position_list=[1, 2], strand_idx=3, id_idx=4)

    convert = {}
    convert['NM'] = set()
    convert['NR'] = set()
    convert_table = {}  # Convert gene symbols to refseq-mrna IDs
    for item in tss_data:
        if item[-1].startswith('NM'):
            convert['NM'].add(item[-2])
        elif item[-1].startswith('NR'):
            convert['NR'].add(item[-2])

        if item[-2] not in convert_table:
            convert_table[item[-2]] = item[-1]
        if item[-2] in convert_table and item[-1].startswith('NM'):  # If a gene symbol is annotated with multiple NM- and NR-, use NM-
            if convert_table[item[-2]].startswith('NR'):
                convert_table[item[-2]] = item[-1]


    bed_data = []
    ii = open(bed_file, 'r')
    max_ = 0.0
    min_ = 100000.0
    for line in ii:
        line = line.replace(' ', '_')
        line = line.strip().split()
        bed_data.append([line[0], line[1], line[2], line[3]])


    temp = assign_gene(bed_data, tss_human19, min_distance=0, max_distance=5000, centre=False, width=True, gene_include=True)

    ### NORMALISE WIDTHS WITHIN THE CELL TYPE
    for no in range(len(temp)):
        if float(temp[no][-1]) > max_:
            max_ = float(temp[no][-1])
        if float(temp[no][-1]) < min_:
            min_ = float(temp[no][-1])
    for no in range(len(temp)):
        if float(temp[no][-1]) - min_ == 0.0:
            temp[no][-1] += 1
        new_value = ((float(temp[no][-1]) - min_) / (max_ - min_))
        temp[no][-1] = new_value
    ###



    for no in range(len(temp)):
        if str(temp[no][0]) == str(bed_data[no][3]):
            bed_data[no].extend([temp[no][1], temp[no][-1]])

    if remove_ncrna == True:
        new_data = []
        for no in range(len(bed_data)):
            if bed_data[no][4] in convert['NM']:
                new_data.append(bed_data[no])
        bed_data = new_data

    # Inter-cell type
    inter_data = {}
    temp = read_file(inter_background, [0,1,2,3])
    for item in temp:
        if item[0] not in inter_data:
            inter_data[item[0]] = [item[1], item[2], item[3]]

    for no in range(len(bed_data)):
        gene = str(bed_data[no][4])
        if gene in inter_data:
            if float(inter_data[gene][0]) != 0.0:
                mean__ = float(inter_data[gene][-1])
        else:
            bed_data[no].extend([None])
            continue
        converted = float(bed_data[no][5]) / mean__
        converted = np.log10(converted)
        p = scipy.stats.norm.cdf(converted, loc=float(inter_data[gene][0]),
                                 scale=np.sqrt(float(inter_data[gene][1])))
        p = 1 - p
        bed_data[no].extend([converted, p])

    new_data = []
    for no in range(len(bed_data)):
        if bed_data[no][6] != None:
            new_data.append(bed_data[no])
    bed_data = new_data

    # Intra-cell type
    widths = []
    for item in bed_data:
        widths.append(float(item[5]))
    intra_raw_mean = mean_(widths)

    sum1 = 0.0
    for no in range(len(bed_data)):
        item = bed_data[no]
        converted = standard1(float(item[5]), intra_raw_mean, log_=True)
        bed_data[no].extend([converted])
        sum1 += converted
    intra_converted_mean = sum1 / len(bed_data)

    sum2 = 0.0
    for item in bed_data:
        sum2 += (float(item[-1]) - intra_converted_mean) * (float(item[-1]) - intra_converted_mean)
    intra_variance = sum2 / len(bed_data)

    for no in range(len(bed_data)):
        p = scipy.stats.norm.cdf(bed_data[no][-1], loc=float(intra_converted_mean),
                                 scale=np.sqrt(float(intra_variance)))
        p = 1 - p
        bed_data[no].extend([p])

    # Selection
    intra_selected = []
    inter_selected = []
    intersection_ = []
    gene_list = set()


    for item in bed_data:
        if float(item[6]) >= inter_threshold:
            inter_selected.append(item)
        if float(item[8]) >= intra_threshold:
            intra_selected.append(item)
        if float(item[6]) >= inter_threshold and float(item[8]) >= intra_threshold:
            intersection_.append(item)
        if list_all_genes == True:
            if str(item[4]) not in gene_list:
                gene_list.add(item[4])

    intra_selected = sort_(intra_selected, 5, reverse_=True)
    inter_selected = sort_(inter_selected, 6, reverse_=True)
    intersection_ = sort_(intersection_, 6, reverse_=True)

    if dominant == True:
        temp = []
        for m in intra_selected:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(intra_selected[int(m)])
        intra_selected = new_list

        temp = []
        for m in inter_selected:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(inter_selected[int(m)])
        inter_selected = new_list

        temp = []
        for m in intersection_:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(intersection_[int(m)])
        intersection_ = new_list

    intra_selected.insert(0,['#chr','start','end','domain','gene','raw_width','inter_log10(width/weighted_mean)','p-value(inter-)','intra_log10(width/mean)','p-value(intra-)'])
    inter_selected.insert(0,
        ['#chr', 'start', 'end', 'domain', 'gene', 'raw_width', 'inter_log10(width/weighted_mean)', 'p-value(inter-)',
         'intra_log10(width/mean)', 'p-value(intra-)'])
    intersection_.insert(0,
        ['#chr', 'start', 'end', 'domain', 'gene', 'raw_width', 'inter_log10(width/weighted_mean)', 'p-value(inter-)',
         'intra_log10(width/mean)', 'p-value(intra-)'])



    write_file(intra_selected, output_+'Intra-selected.txt')
    write_file(inter_selected, output_+'Inter-selected.txt')
    #write_file(intersection_, output_+'Intersection.txt')
    if list_all_genes == True:
        gene_list_ = []
        for item in gene_list:

            if item in convert_table:
                name = convert_table[item]
            else:
                name = ''
            gene_list_.append([item, name])
        write_file(gene_list_, output_+'all_genes.txt')


def width_analysis2(bed_file, tss_, inter_background, output_, inter_threshold=0.3, intra_threshold=0.3, remove_ncrna=False, dominant=False, list_all_genes=False):
    ''' perform all steps required for the width analysis
    NOTE: THIS FIRST CALCULATE Z-SCORES WITHIN CELL TYPES AND THEN CALCULATE MODIFIED INTER-CELL TYPE DS (i.e. x - mean)

    bed_file = a BED text file (input data)
    tss_ = a RefSeq TSS text file downloaded from NCBI
    inter_background = a background (mean, variance) for inter-cell type analysis (e.g. 'summary_means_variances__.txt')
    output_ = a output filename
    inter- and intra-threshold = a value of converted width (log10(width/mean)) acting as a threshold for selection of domains
    remove_ncrna = a boolean to indicate whether to remove any ncRNA genes from intra-cell type analysis
    dominant = a boolean to indicate whether only dominant peaks should be included in the analysis
    list_all_genes = a boolean to indicate whether to export a list of all tagged genes. File name is output_+'all_genes.txt'

    '''
    tss_data = read_file(tss_,[1, 3, 4, 2, 5, 0])
    tss_human19 = tss(tss_data, chr_idx=0, position_list=[1, 2], strand_idx=3, id_idx=4)

    convert = {}
    convert['NM'] = set()
    convert['NR'] = set()
    convert_table = {}  # Convert gene symbols to refseq-mrna IDs
    for item in tss_data:
        if item[-1].startswith('NM'):
            convert['NM'].add(item[-2])
        elif item[-1].startswith('NR'):
            convert['NR'].add(item[-2])

        if item[-2] not in convert_table:
            convert_table[item[-2]] = item[-1]
        if item[-2] in convert_table and item[-1].startswith('NM'):  # If a gene symbol is annotated with multiple NM- and NR-, use NM-
            if convert_table[item[-2]].startswith('NR'):
                convert_table[item[-2]] = item[-1]


    bed_data = []
    ii = open(bed_file, 'r')
    max_ = 0.0
    min_ = 100000.0
    for line in ii:
        line = line.replace(' ', '_')
        line = line.strip().split()
        bed_data.append([line[0], line[1], line[2], line[3]])


    temp = assign_gene(bed_data, tss_human19, min_distance=0, max_distance=5000, centre=False, width=True, gene_include=True)

    ### NORMALISE WIDTHS WITHIN THE CELL TYPE USING Z-SCORES
    array_ = []
    temp_mean = 0.0
    temp_sd = 0.0

    for no in range(len(temp)):
        temp_mean += float(temp[no][-1])
        array_.append(float(temp[no][-1]))
    temp_mean = temp_mean / float(len(temp))
    temp_sd = np.std(array_)

    for no in range(len(temp)):
        new_value = (float(temp[no][-1]) - temp_mean) / float(temp_sd)
        temp[no][-1] = new_value

    ###



    for no in range(len(temp)):
        if str(temp[no][0]) == str(bed_data[no][3]):
            bed_data[no].extend([temp[no][1], temp[no][-1]])

    if remove_ncrna == True:
        new_data = []
        for no in range(len(bed_data)):
            if bed_data[no][4] in convert['NM']:
                new_data.append(bed_data[no])
        bed_data = new_data

    # Inter-cell type
    inter_data = {}
    temp = read_file(inter_background, [0,1,2,3])
    for item in temp:
        if item[0] not in inter_data:
            inter_data[item[0]] = [item[1], item[2], item[3]]

    for no in range(len(bed_data)):
        gene = str(bed_data[no][4])
        if gene in inter_data:
            if float(inter_data[gene][0]) != 0.0:
                mean__ = float(inter_data[gene][-1])
        else:
            bed_data[no].extend([None])
            continue
        converted = float(bed_data[no][5]) / mean__
        converted = np.log10(converted)
        p = scipy.stats.norm.cdf(converted, loc=float(inter_data[gene][0]),
                                 scale=np.sqrt(float(inter_data[gene][1])))
        p = 1 - p
        bed_data[no].extend([converted, p])

    new_data = []
    for no in range(len(bed_data)):
        if bed_data[no][6] != None:
            new_data.append(bed_data[no])
    bed_data = new_data

    # Intra-cell type
    widths = []
    for item in bed_data:
        widths.append(float(item[5]))
    intra_raw_mean = mean_(widths)

    sum1 = 0.0
    for no in range(len(bed_data)):
        item = bed_data[no]
        converted = standard1(float(item[5]), intra_raw_mean, log_=True)
        bed_data[no].extend([converted])
        sum1 += converted
    intra_converted_mean = sum1 / len(bed_data)

    sum2 = 0.0
    for item in bed_data:
        sum2 += (float(item[-1]) - intra_converted_mean) * (float(item[-1]) - intra_converted_mean)
    intra_variance = sum2 / len(bed_data)

    for no in range(len(bed_data)):
        p = scipy.stats.norm.cdf(bed_data[no][-1], loc=float(intra_converted_mean),
                                 scale=np.sqrt(float(intra_variance)))
        p = 1 - p
        bed_data[no].extend([p])

    # Selection
    intra_selected = []
    inter_selected = []
    intersection_ = []
    gene_list = set()


    for item in bed_data:
        if float(item[6]) >= inter_threshold:
            inter_selected.append(item)
        if float(item[8]) >= intra_threshold:
            intra_selected.append(item)
        if float(item[6]) >= inter_threshold and float(item[8]) >= intra_threshold:
            intersection_.append(item)
        if list_all_genes == True:
            if str(item[4]) not in gene_list:
                gene_list.add(item[4])

    intra_selected = sort_(intra_selected, 5, reverse_=True)
    inter_selected = sort_(inter_selected, 6, reverse_=True)
    intersection_ = sort_(intersection_, 6, reverse_=True)

    if dominant == True:
        temp = []
        for m in intra_selected:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(intra_selected[int(m)])
        intra_selected = new_list

        temp = []
        for m in inter_selected:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(inter_selected[int(m)])
        inter_selected = new_list

        temp = []
        for m in intersection_:
            temp.append(m[4])
        new_list = []
        indexes, elements = unique(temp)
        for m in indexes:
            new_list.append(intersection_[int(m)])
        intersection_ = new_list

    intra_selected.insert(0,['#chr','start','end','domain','gene','raw_width','inter_log10(width/weighted_mean)','p-value(inter-)','intra_log10(width/mean)','p-value(intra-)'])
    inter_selected.insert(0,
        ['#chr', 'start', 'end', 'domain', 'gene', 'raw_width', 'inter_log10(width/weighted_mean)', 'p-value(inter-)',
         'intra_log10(width/mean)', 'p-value(intra-)'])
    intersection_.insert(0,
        ['#chr', 'start', 'end', 'domain', 'gene', 'raw_width', 'inter_log10(width/weighted_mean)', 'p-value(inter-)',
         'intra_log10(width/mean)', 'p-value(intra-)'])



    write_file(intra_selected, output_+'Intra-selected.txt')
    write_file(inter_selected, output_+'Inter-selected.txt')
    #write_file(intersection_, output_+'Intersection.txt')
    if list_all_genes == True:
        gene_list_ = []
        for item in gene_list:

            if item in convert_table:
                name = convert_table[item]
            else:
                name = ''
            gene_list_.append([item, name])
        write_file(gene_list_, output_+'all_genes.txt')

def width_analysis3(bed_file, tss_, max_width_, remove_ncrna=False, dominant=True, list_all_genes=False, centre_=False, only_width=True):
    ''' perform all steps required for the width analysis
    NOTE: THIS FIRST CALCULATE Z-SCORES WITHIN CELL TYPES AND THEN CALCULATE MODIFIED INTER-CELL TYPE DS (i.e. x - mean)

    bed_file = a BED text file (input data) or input data as a list of lists
    tss_ = a RefSeq TSS text file downloaded from NCBI
    inter_background = a background (mean, variance) for inter-cell type analysis (e.g. 'summary_means_variances__.txt')
    output_ = a output filename
    inter- and intra-threshold = a value of converted width (log10(width/mean)) acting as a threshold for selection of domains
    remove_ncrna = a boolean to indicate whether to remove any ncRNA genes from intra-cell type analysis
    dominant = a boolean to indicate whether only dominant peaks should be included in the analysis
    list_all_genes = a boolean to indicate whether to export a list of all tagged genes. File name is output_+'all_genes.txt'
    max_width_ = None (no limit) or a numerical value


    '''
    if type(tss_)==str:
        tss_data = read_file(tss_,[1, 3, 4, 2, 5, 0])   # for hg19
        #tss_data = read_file(tss_,[2, 4, 5, 3, 12, 1])  # for mm 10
        tss_human19 = tss(tss_data, chr_idx=0, position_list=[1, 2], strand_idx=3, id_idx=4)


    convert = {}
    convert['NM'] = set()
    convert['NR'] = set()
    convert_table = {}  # Convert gene symbols to refseq-mrna IDs
    for item in tss_data:
        if item[-1].startswith('NM'):
            convert['NM'].add(item[-2])
        elif item[-1].startswith('NR'):
            convert['NR'].add(item[-2])

        if item[-2] not in convert_table:
            convert_table[item[-2]] = item[-1]
        if item[-2] in convert_table and item[-1].startswith('NM'):  # If a gene symbol is annotated with multiple NM- and NR-, use NM-
            if convert_table[item[-2]].startswith('NR'):
                convert_table[item[-2]] = item[-1]


    bed_data = []
    max_ = 0.0
    min_ = 100000.0
    if type(bed_file) == str:
        ii = open(bed_file, 'r')

        for line in ii:
            line = line.replace(' ', '_')
            line = line.strip().split()
            bed_data.append([line[0], line[1], line[2], line[3]])
    else:
        bed_data = bed_file      

    temp = assign_gene(bed_data, tss_human19, min_distance=0, max_distance=max_width_, centre=centre_, width=True, gene_include=True)




    temp = sort_(temp, -1, reverse_=True)



    #for no in range(len(temp)):
    #    if str(temp[no][0]) == str(bed_data[no][3]):
    #        bed_data[no].extend([temp[no][1], temp[no][-1]])


    if dominant == True:
        sorted_list_= []
        genes = set()
        for item in temp:
            if item[1] != '':
                if item[1] not in genes:
                    genes.add(item[1])
                    sorted_list_.append(item)
    else:
        sorted_list_ = temp


    if remove_ncrna == True:
        mrna_list = set()
        temp__ = open('mRNA_genes.txt', 'r')
        for line in temp__:
            line = line.strip().split()
            if not line[0].startswith('#'):
                mrna_list.add(line[0])
        new_data = []
        for no in range(len(sorted_list_)):
            if sorted_list_[no][1] in mrna_list:
                new_data.append(sorted_list_[no])
    else:
        new_data = sorted_list_

    final = []
    if only_width==True:
        for item in new_data:
            final.append([item[1], item[-1]])        
        final.insert(0, ['width'])
        return final
    else:
        for item in new_data:
            final.append(item)
        return final


    #    temp = []
    #    for m in inter_selected:
    #        temp.append(m[4])
    #    new_list = []
    #    indexes, elements = unique(temp)
    #    for m in indexes:
    #        new_list.append(inter_selected[int(m)])
    #    inter_selected = new_list

    #    temp = []
    #    for m in intersection_:
    #        temp.append(m[4])
    #    new_list = []
    #    indexes, elements = unique(temp)
    #    for m in indexes:
    #        new_list.append(intersection_[int(m)])
    #    intersection_ = new_list

def intra_sd(data_):
    ### Take a list of widths and calculate width/SD (i.e. DS)
    ### Returns a list

    temp_sd = np.std(data_)
    results = []
    for item in data_:
        converted = float(item) / temp_sd
        results.append(converted)

    return results

def contribution_score(data_, log_=True):
    ### Take a list of DS scores and convert to contribution scores (i.e. DS / sum of DSs)
    ### Returns a list
    # log_ = a boolean to indicate whether to use log10 conversion
    temp_sum = 0.0
    results = []
    for item in data_:
        temp_sum += float(item)
    for item in data_:
        converted = float(item) /temp_sum
        if log_ == True:
            converted = np.log10(converted)
        results.append(converted)
    return results

def metacell_scores(input_,data_, ave=True):
    ### Calculate (max) 127 metacell scores for a given list of genes (input_)
    ### data_ = a dictionary of dictionary specifying a DS score for a given gene in a cell type
    ### e.g. data_[gene][cell_type] = 3.56
    # ave = a boolean to indicate whether to average
    results = {}
    total_no = len(input_)
    for gene in input_:
        if gene in data_:
            for cell in data_[gene]:
                if cell not in results:
                    results[cell] = 0.0
                results[cell] += float(data_[gene][cell])
    if ave==True:
        for cell in results:
            results[cell] = float(results[cell]) / total_no
    return results


def check_overlap(input_, data_):
    ### Take a list of genomic coordinate & a list of set of coordinates
    ### Return a number of overlapped elements.
    ### Note: It assumes both input and dataset are within the same chromosome.
    ### Note: Dataset must be first sorted in order
    ### e.g. check_overlap([1,10],[[1,2],[11,20],...]])
    ###      Returns 1
    cnt = 0
    for item in data_:
        if float(item[1]) == float(input_[0]):
            cnt += 1
        elif (float(item[0]) >= float(input_[0])) and (float(item[1]) <= float(input_[1])):
            cnt += 1
        elif float(item[1]) == float(input_[1]):
            cnt += 1
    return cnt


def ensg2symbol(filename):
    ### Take a NCBI gene info file and write put a text file that connect ENSG to gene symbol
    data_ = open(filename, 'r')
    results = {}
    for line in data_:
        line = line.strip().split()
        if not line[0].startswith('#'):
            symbol = line[2]
            if 'ENSG' in line[5]:
                temp = line[5].split('Ensembl:')

                temp_ = temp[1][0:14]

                ensg = temp_
                print ensg
                if ensg not in results:
                    results[ensg] = symbol
    print results
    #results.insert(0, ['#ENSG','Symbol'])
    write_file(results, '/Users/woojunshim/Research/Data/ensg2symbol.txt')

def replace(input_, output_, from_, to_):
    ### Read in a text file and substitute to_ for from_
    input_data_ = open(input_, 'r')
    output_data_ = open(output_, 'w')
    for line in input_data_:
        line = line.replace(from_, to_)
        output_data_.write(line)

def random_sample(input_list, n_samples=None, replacement=False):
    if replacement == False:
        if n_samples == None:
            return random.sample(input_list, len(input_list))
        else:
            return random.sample(input_list, n_samples)
    else:
        if n_samples == None:
            return np.random.choice(input_list, len(input_list))
        else:
            return np.random.choice(input_list, n_samples)



def main8(pathway, filename):
    ### PERFORMS THE STATISTICAL ANALYSIS ON RANKED GENE LIST USING CONTRIBUTION SCORE
    ### RUN THIS FOLLOWING 'main3' ###

    filename_ = pathway+filename
    dataset = []
    print 'Processing..',filename
    data_ = open(filename_, 'r')
    for line in data_:
        line = line.strip().split()
        if line[0].startswith('#'):
            continue
        else:
            if len(line) > 2:
                dataset.append([line[0],line[1],line[2]])
    temp = []
    for item in dataset:
        temp.append(float(item[2]))
    converted = intra_sd(temp)

    cont = contribution_score(converted, log_=True)

    for no in range(len(dataset)):
        dataset[no].extend([converted[no],cont[no]])

    background = read_file(
        '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/gene_contribution_table.txt', [1, 2],
        rowname='0')


    missing = 0
    for no in range(len(dataset)):
        item = dataset[no]
        value = float(item[-1])
        gene = str(item[1])
        if gene in background:
            if float(background[gene][0][1]) != float(0.0):
                z_ = statistics.z_score(value, mean_=float(background[gene][0][0]), sd_=float(background[gene][0][1]))
                p_ = statistics.p_value(z_)
                dataset[no].extend([p_])
            else:
                missing += 1
        else:
            missing += 1

    print 'Missing genes =', missing

    dataset.insert(0, ['#peak','gene','width','width/SD(DS)','contribution(log-converted)','p-value(contribution)'])
    write_file(dataset, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/RESULT_'+filename)

def run_file():
    ### Automate the analysis

    pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/'
    pathway1=''
    files = ['RESULT_CVP_K4.bed_assigned_genes.txt','RESULT_GSM2279997_7_2_H3K4me3_v2_peaks.narrowPeak_assigned_genes.txt']
    data_ = open('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb_ALL_standard.txt', 'r')
    data__ = {}
    for line in data_:
        line = line.strip().split()
        if line[0].startswith('#'):
            continue
        else:
            if line[0] not in data__:
                data__[line[0]] = {}
            if line[1] not in data__[line[0]]:
                data__[line[0]][line[1]] = line[3]
    result = []
    threshold = 2.0
    for file in files:
        input_ = []
        data_1 = []
        input__ = read_file(pathway+file, [0,1,3])
        for item in input__:
            if float(item[-1]) > float(threshold):
                input_.append(item)
                data_1.extend([item[1]])
        print data_1
        #print data__
        aa = metacell_scores(data_1, data__, ave=True)
        print aa
        for cell in aa:
            result.append([cell,aa[cell]])
        write_file(result,pathway+'METACELL_'+file)


def main9():

    pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/'
    files = ['RESULT_CVP_K4.bed_assigned_genes.txt','RESULT_GSM2279997_7_2_H3K4me3_v2_peaks.narrowPeak_assigned_genes.txt']
    data_ = open('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt', 'r')
    data__ = {}
    cell_list = []
    for line in data_:
        line = line.strip().split()
        if line[0].startswith('#'):
            continue
        else:
            if line[0] not in data__:
                data__[line[0]] = {}
            if line[1] not in data__[line[0]]:
                data__[line[0]][line[1]] = line[4]
            if line[1] not in cell_list:
                cell_list.append(line[1])
    results = {}
    results['colnames'] = cell_list
    for gene in data__:
        if gene not in results:
            results[gene] = []
        for cell in cell_list:
            if cell not in data__[gene]:
                results[gene].extend([0.0])
            else:
                results[gene].extend([data__[gene][cell]])
    write_table(results, pathway+'All_genes_score.txt')












def main2():
    ### Read in a data entry (as a list of lists)
    entry1 = read_file('CM_CP.txt', ['1','2'])
    print entry1

    ### Read in a data entry (as a dictionary sorted by '##Gene_ID')
    entry2 = read_file('CM_CP.txt', ['Symbol', 'H3K27ac-CP_Width'], rowname='##Gene_ID')
    print entry2

    ### Export entry1 & entry2 to an text file
    write_file([entry1, entry2], 'combined_entry.txt', tag='test_', header=False)

def main1():
    ### Create TSS data entry and assign genes to a set of domains
    tss_data = read_file('/Users/woojunshim/Research/Data/hg19_TSS.txt', [1,3,4,2,5])  # Numbers as column index in the input text file
    tss_human19 = tss(tss_data, chr_idx=0, position_list=[1,2] , strand_idx=3, id_idx=4)  # Create a TSS instance from 'tss_data' data entry
    write_file(tss_human19, 'tss.txt')

    ### Load a BED file and create a data entry
    E095_H3K4me3 = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/E095-H3K4me3.gappedPeak.bed', [0,1,2,3])

    ### Assign genes to the BED domains
    E095_H3K4me3_genes = assign_gene(E095_H3K4me3, tss_human19, min_distance=0, max_distance=1000000, centre=False, width=True, gene_include=True)


    ### Extract entries only with 'HAND2'
    E095_hand2 = find_entry(E095_H3K4me3_genes, 'HAND2-AS1')
    print E095_hand2

    ### Get domain IDs assigned to 'HAND2'
    #E095_hand2_ids = extract_element(E095_hand2,0)
    #print E095_hand2_ids

    ### Calculate width of domains assigned to 'HAND2'
    #E095_hand2_width =[]
    #for id in E095_hand2_ids:
    #    E095_hand2_width.append(width(id, E095_H3K4me3))
    #print E095_hand2_width

    ### Calculate empirical p-value for each HAND2-assigned domain
    #E095_hand2_pvalue = []
    #for item in E095_hand2_width:
    #    E095_hand2_pvalue.append(empirical_pvalue(item, E095_hand2_width, lower=False))
    #print E095_hand2_pvalue



def main():
    ### Assign genes to domains from multiple BED files in a directory
    pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/'
    output_pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/'
    ids_ = ['E015','E018','E010','E043','E029','E026','E052','E057','E053','E112','E070','E063','E100','E095','E111','E084']

    tss_data = read_file('/Users/woojunshim/Research/Data/hg19_TSS.txt',[1, 3, 4, 2, 5])  # Numbers as column index in the input text file
    tss_human19 = tss(tss_data, chr_idx=0, position_list=[1,2], strand_idx=3,id_idx=4)  # Create a TSS instance from 'tss_data' data entry

    genes_of_interest = ['HAND2', 'ISL1', 'NKX2-5', 'TNNI1', 'ATP2A2', 'COL1A1']

    for gene in genes_of_interest:
        output_pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/'+gene
        if not os.path.exists(output_pathway):
            os.makedirs(output_pathway)
        os.chdir(output_pathway)
        domains = {}
        widths = {}
        pvalues = {}
        ### Find domains IDs & width first
        for id in ids_:
            widths[id] = []  # widths['E095'] = [xxx,xxx,xxx,...,]
            domains[id] = []  # domains['E095'] = [[ID, width, pvalue],...,]
            H3K4me3 = read_file(pathway+id+'-H3K4me3.gappedPeak.bed', [0,1,2,3])
            assigned_genes = assign_gene(H3K4me3, tss_human19, min_distance=0, max_distance=50000, centre=False)
            entry = find_entry(assigned_genes, gene)
            entry_ids = extract_element(entry, 0)
            for m in entry_ids:
                widths[id].append(width(m, H3K4me3))
            for i in range(len(entry)):
                domains[id].append([entry[i][0], widths[id][i]])



        ### Calculate empirical p-values
        combined_widths = []
        for id in ids_:
            for m in widths[id]:
                combined_widths.append(float(m))
            pvalues[id] = []

        for id in ids_:
            for i in range(len(widths[id])):
                pvalues[id].append(empirical_pvalue(widths[id][i], combined_widths, lower=False))

            domains[id] = cbind(domains[id], pvalues[id])

            write_file([domains[id]], id+'_'+gene+'.txt', tag=id+'_')


def main3():
    ### USE THIS TO ASSIGN GENES ###


    pathway1 = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/'
    pathway2 = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/'
    pathway3 = '/Users/woojunshim/Research/Data/'
    user_pathway='/Users/woojunshim/Research/Data/'
    output_pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/'
    tissue_groups = {}
    #tissue_groups['IMR90'] = ['E017']
    #tissue_groups['ES_cell'] = ['E002', 'E008', 'E001', 'E015', 'E014', 'E016', 'E003', 'E024']
    #tissue_groups['iPSC'] = ['E020', 'E019', 'E018', 'E021', 'E022']
    #tissue_groups['ES_deriv.'] = ['E007', 'E009', 'E010', 'E013', 'E012', 'E011', 'E004', 'E005', 'E006']
    #tissue_groups['Blood&T_cell'] = ['E062', 'E034', 'E045', 'E033', 'E044', 'E043', 'E039', 'E041', 'E042',
    #                                 'E040', 'E037', 'E048', 'E038', 'E047']
    #tissue_groups['HSC&B_cell'] = ['E029', 'E031', 'E035', 'E051', 'E050', 'E036', 'E032', 'E046', 'E030']
    #tissue_groups['Mesench.'] = ['E026', 'E049', 'E025', 'E023']
    #tissue_groups['Myosat.'] = ['E052']
    #tissue_groups['Epithelial'] = ['E055', 'E056', 'E059', 'E061', 'E057', 'E058', 'E028', 'E027']
    #tissue_groups['Neurosph.'] = ['E054', 'E053']
    #tissue_groups['Thymus'] = ['E112', 'E093']
    #tissue_groups['Brain'] = ['E071', 'E074', 'E068', 'E069', 'E072', 'E067', 'E073', 'E070', 'E082', 'E081']
    #tissue_groups['Adipose'] = ['E063']
    #tissue_groups['Muscle'] = ['E100', 'E108', 'E107', 'E089', 'E090']
    #tissue_groups['Heart'] = ['E083', 'E104', 'E095', 'E105', 'E065']
    #tissue_groups['Smooth_muscle'] = ['E078', 'E076', 'E103', 'E111']
    #tissue_groups['Digestive'] = ['E092', 'E085', 'E084', 'E109', 'E106', 'E075', 'E101', 'E102', 'E110', 'E077',
    #                              'E079', 'E094']
    #tissue_groups['Other'] = ['E099', 'E086', 'E088', 'E097', 'E087', 'E080', 'E091', 'E066', 'E098', 'E096','E113']
    #tissue_groups['ENCODE2012'] = ['E114', 'E115', 'E116', 'E117', 'E118', 'E119', 'E120', 'E121', 'E122', 'E123','E124', 'E125', 'E126', 'E127', 'E128', 'E129']
    tissue_groups['user-defined'] = ['CVP_K4.bed','GSM2279997_7_2_H3K4me3_v2_peaks.narrowPeak.D32']
    tissue_groups['user-defined'] = ['GSE96290_ENCFF682XDY_peaks_hg19_fetal_cardiac_muscle_invitro.bed','GSE96290_ENCFF382BJH_replicated_peaks_hg19_fetal_cardiac_muscle_invitro.bed','GSM2279960_1_8_H3K4me3_v2_peaks.narrowPeak.D32']
    tissue_groups['user-defined'] = ['GSM2066619_H3K4me3_D11_fibroblast.peaks.txt']
    tissue_groups['user-defined'] = ['GSE63094_KK_RWPE1_H3K4me3_combinedreps_MACS14_peaks.bed']
    for group in tissue_groups:
        for epi in tissue_groups[group]:
            print epi
            if group != 'user-defined':
                bed_ = pathway2+epi+'-H3K4me3.gappedPeak.bed'
            else:
                bed_ = user_pathway+epi
            t = pathway3+'hg19_TSS_.txt'

            result_ = width_analysis3(bed_file=bed_, tss_=t, remove_ncrna=False, dominant=True, list_all_genes=True)
            write_file(result_, output_pathway+epi+'_assigned_genes.txt')



def main4():
    pathway3 = '/Users/woojunshim/Research/Data/'
    data_ = open(pathway3+'CVP_K4.txt', 'r')

    bed_data = []
    cnt = 1
    for line in data_:
        line = line.replace(' ','_')
        line = line.strip().split()
        if line[0].startswith('ENSG'):
            if len(line)!= 6:
                line.insert(2,'NA')
            bed_data.append(['chr'+line[3], line[4], line[5], 'domain_'+str(cnt)])
            cnt += 1
    write_file(bed_data, pathway3+'CVP_K4.bed')

def main5():
    ''' Create the table sorted by rank for each gene
    '''
    data_ = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt', [1,4], rowname='0')
    filter = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt', [8,9], rowname='0')

    for gene in data_:
        for i in range(len(data_[gene])):
            data_[gene][i][1] = float(data_[gene][i][1])

    results = {}
    total_cells = 127
    gene_list = []
    cell_types = {}
    for i in range(total_cells):
        results[str(i+1)] = []
        cell_types[str(i+1)] = []
    for gene in data_:
        if filter[gene][0][1] == 'True' and filter[gene][0][0] == 'False':  # CHANGE HERE
            gene_list.append(gene)
            temp = sort_(data_[gene], idx=1, reverse_=True)
            for i in range(total_cells):
                if len(temp) <= i:
                    results[str(i+1)].extend(['NA'])
                    cell_types[str(i+1)].extend(['NA'])
                else:
                    results[str(i+1)].extend([temp[i][1]])
                    cell_types[str(i+1)].extend([temp[i][0]])
    results['colnames'] = gene_list
    cell_types['colnames'] = gene_list
    write_table(results, output_file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_nontf_2.5kb.txt')
    write_table(cell_types, output_file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/cell_types_table_mrna_nontf_2.5kb.txt')

def main6():
    ''' Create the table sorted by rank for each gene (change of the value, i.e. slope)
    '''
    data_ = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt', [1,4], rowname='0')
    filter = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt', [7,8], rowname='0')

    for gene in data_:
        for i in range(len(data_[gene])):
            data_[gene][i][1] = float(data_[gene][i][1])

    results = {}
    total_cells = 99
    gene_list = []
    cell_types = {}
    for i in range(1, total_cells):
        results[str(i+1)] = []
        cell_types[str(i+1)] = []
    for gene in data_:

        gene_list.append(gene)
        temp = sort_(data_[gene], idx=1, reverse_=True)
        for i in range(1, total_cells):
            if len(temp) <= i:
                results[str(i+1)].extend(['NA'])
                cell_types[str(i+1)].extend(['NA'])
            else:
                results[str(i+1)].extend([float(temp[i][1])-float(temp[i-1][1])])
                cell_types[str(i+1)].extend([temp[i][0]])
    results['colnames'] = gene_list
    cell_types['colnames'] = gene_list
    write_table(results, output_file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_slope_5kb.txt')
    write_table(cell_types, output_file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/cell_types_table_only_non_tf_5kb.txt')

def main7():
    ''' Create the table of numbers of genes (TFs vs. non-TFs)
    '''
    data_ = read_file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt', [1,4,7,8], rowname='0')

    for gene in data_:
        for i in range(len(data_[gene])):
            data_[gene][i][1] = float(data_[gene][i][1])

    results = []
    for gene in data_:
        max_ = -2.0
        min_ = 2.0
        for item in data_[gene]:
            if float(item[1]) > max_:
                max_ = float(item[1])
            if float(item[1]) < min_:
                min_ = float(item[1])
        diff = max_ - min_
        results.append([gene, len(data_[gene]), str(max_), str(min_), str(diff), data_[gene][0][2], data_[gene][0][3]])

    write_file(results, output_file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/cell_type_distributions_5kb.txt')


def main10():
    ### TSS overlap analysis

    tissue_groups = {}
    tissue_groups['IMR90'] = ['E017']
    tissue_groups['ES_cell'] = ['E002', 'E008', 'E001', 'E015', 'E014', 'E016', 'E003', 'E024']
    tissue_groups['iPSC'] = ['E020', 'E019', 'E018', 'E021', 'E022']
    tissue_groups['ES_deriv.'] = ['E007', 'E009', 'E010', 'E013', 'E012', 'E011', 'E004', 'E005', 'E006']
    tissue_groups['Blood&T_cell'] = ['E062', 'E034', 'E045', 'E033', 'E044', 'E043', 'E039', 'E041', 'E042',
                                     'E040', 'E037', 'E048', 'E038', 'E047']
    tissue_groups['HSC&B_cell'] = ['E029', 'E031', 'E035', 'E051', 'E050', 'E036', 'E032', 'E046', 'E030']
    tissue_groups['Mesench.'] = ['E026', 'E049', 'E025', 'E023']
    tissue_groups['Myosat.'] = ['E052']
    tissue_groups['Epithelial'] = ['E055', 'E056', 'E059', 'E061', 'E057', 'E058', 'E028', 'E027']
    tissue_groups['Neurosph.'] = ['E054', 'E053']
    tissue_groups['Thymus'] = ['E112', 'E093']
    tissue_groups['Brain'] = ['E071', 'E074', 'E068', 'E069', 'E072', 'E067', 'E073', 'E070', 'E082', 'E081']
    tissue_groups['Adipose'] = ['E063']
    tissue_groups['Muscle'] = ['E100', 'E108', 'E107', 'E089', 'E090']
    tissue_groups['Heart'] = ['E083', 'E104', 'E095', 'E105', 'E065']
    tissue_groups['Smooth_muscle'] = ['E078', 'E076', 'E103', 'E111']
    tissue_groups['Digestive'] = ['E092', 'E085', 'E084', 'E109', 'E106', 'E075', 'E101', 'E102', 'E110', 'E077',
                                  'E079', 'E094']
    tissue_groups['Other'] = ['E099', 'E086', 'E088', 'E097', 'E087', 'E080', 'E091', 'E066', 'E098', 'E096', 'E113']
    tissue_groups['ENCODE2012'] = ['E114', 'E115', 'E116', 'E117', 'E118', 'E119', 'E120', 'E121', 'E122', 'E123',
                                   'E124', 'E125', 'E126', 'E127', 'E128', 'E129']
    threshold = 0.99 # as rank percentile of intra-DS
    pathway = '/Users/woojunshim/Research/Data/'
    tss1 = read_file(pathway+'hg19.cage_peak_phase1and2combined_coord.bed', [6,7], rowname='0')
    tss2 = read_file(pathway+'hg19_TSS_.txt', [3,4], rowname='1')
    pdf__ = open('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb_ALL_standard.txt', 'r')
    pdf_ = {}
    all_peaks = {}
    all_peaks_ = open('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/all_peaks_assigned.txt','r')
    for line in all_peaks_:
        line = line.strip().split()
        if len(line) == 6:
            if line[0] not in all_peaks:
                all_peaks[line[0]] = {}
            if line[2] not in all_peaks[line[0]]:
                all_peaks[line[0]][line[2]] = 0
            all_peaks[line[0]][line[2]] += 1


    for line in pdf__:
        line = line.strip().split()
        if not line[0].startswith('#'):
            if line[1] not in pdf_:
                pdf_[line[1]] = []
            pdf_[line[1]].append([line[0], line[2], float(line[3]), line[-2]])
    data_ = {}
    for group in tissue_groups:
        for epi in tissue_groups[group]:
            print 'Processing..', epi
            raw_data = {}
            data_[epi] = []
            temp = pdf_[epi]
            temp = sort_(temp, idx=2, reverse_=True)
            total_no = len(temp)
            cutoff = int(round(total_no*(1-threshold)))
            for no in range(cutoff):
                data_[epi].append([temp[no][1], temp[no][0], temp[no][-1]])  # peak ID
            raw_data_ = open('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/'+epi+'-H3K4me3.gappedPeak.bed', 'r')
            for line in raw_data_:
                line = line.strip().split()
                if line[3] not in raw_data:
                    raw_data[line[3]] = [line[0], line[1], line[2]]
            for no in range(len(data_[epi])):
                id = data_[epi][no][0]
                gene = data_[epi][no][1]
                dataset = tss1[raw_data[id][0]]
                dataset = sort_(dataset, idx=0, reverse_=False)
                cnt = check_overlap([int(raw_data[id][1]), int(raw_data[id][2])], dataset)
                width = int(raw_data[id][2]) - int(raw_data[id][1])
                data_[epi][no].extend([raw_data[id][0], raw_data[id][1], raw_data[id][2], width, cnt])
            for no in range(len(data_[epi])):
                id = data_[epi][no][0]
                gene = data_[epi][no][1]
                dataset = tss2[raw_data[id][0]]
                dataset = sort_(dataset, idx=0, reverse_=False)
                cnt = check_overlap([int(raw_data[id][1]), int(raw_data[id][2])], dataset)
                width = int(raw_data[id][2]) - int(raw_data[id][1])
                data_[epi][no].extend([cnt, all_peaks[epi][gene]])
    # data_['cell_type'] = [['peakID','gene','TF?','chr','start','end','width','TSS_counts(ByDominantPeak)','number_of_peaks']]
    write_file(data_, pathway+'TSS_H3K4me3_overlap_dominant.txt')























if __name__ == '__main__':
    #ensg2symbol('/Users/woojunshim/Research/Data/Homo_sapiens.gene_info')
    #main9()
    replace('/Users/woojunshim/Research/Data/Roadmap_IDs.txt','/Users/woojunshim/Research/Data/Roadmap_IDs_.txt', ' ', '_')

