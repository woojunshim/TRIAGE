source('/Users/woojunshim/Research/Scripts/R_scripts.R')

### Heatmap
data_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/Excluded_genes_table.txt')
library(gplots)
library(RColorBrewer)
hmcol = colorRampPalette(c("white","red"))(256)
heatmap.2(as.matrix(data_), col=hmcol, dendrogram='none',trace='none', main='Similarity of H3K4me3 domains in HAND2 locus', labCol=colnames(data_), labRow=rownames(data_))

### Histogram
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/'
genes = c('ATF1','ATP2A2','COL1A1','HAND2-AS1','ISL1','NKX2-5','TNNI1')
pdf('Histogram_converted_new.pdf')
for (i in genes){
  filename = paste(pathway,i,'.txt',sep='')
  data_ = read.table(filename)
  x_ = seq(1, length(as.numeric(data_[,2])))
  y_ = as.numeric(data_[,4])
  ## Transform raw values
  #mean_ = mean(y_)
  #y_ = log10(y_ / mean_)
  #n_ = paste('n=',length(y_),sep='')
  main_ = paste(i,n_,sep='   ')
  hist(y_, xlab='log10(width/mean)', ylab='Frequency', main=i, prob=TRUE)
  # lines(density(y_),col='blue', lwd=2)
  lines(density(y_, adjust=2),col='blue', lwd=2)
}
dev.off()




### Threshold analysis
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/'
data_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt')
data_TF = data_1[data_1$V8=='True',]
data_else = data_1[data_1$V8=='False',]

# Create 'count' table (incl. statistics by Fisher's Exact Test)
# Fisher's exact test (odd-ratios & p-value, one-sided)
# Null hypothesis = broad widths are not assiciated with known TFs  
x = seq(10.0, 0.0, by=-0.1)
count_table = matrix(nrow=4, ncol=length(x))
count_table = data.frame(count_table)
result_table = matrix(nrow=length(unique(data_TF$V2)), ncol=length(x))
result_table = data.frame(result_table)
colnames(count_table) = x
rownames(count_table) = c('TFs(domains)', 'non-TFs(domains)', 'Odd-ratio','p-value')
colnames(result_table) = x
rownames(result_table) = unique(data_TF$V2)

for (i in 1:length(x)){
  count1 = length(which(data_TF[,4]>x[i]))
  count2 = length(which(data_else[,4]>x[i]))
  count_table[1,i] = count1
  count_table[2,i] = count2
  matrix_ = matrix(c(count_table[1,i], length(data_TF[,1]) - count_table[1,i], count_table[2,i], length(data_else[,1]) - count_table[2,i]), nrow=2, byrow=TRUE)
  temp = fisher.test(matrix_, alternative='greater', or=1)
  count_table[3,i] = temp$estimate
  count_table[4,i] = temp$p.value
}

for (epi in rownames(result_table)){
  epi_ = data_1[data_1$V2==epi,]
  data_TF = epi_[epi_$V8=='True',]
  data_else = epi_[epi_$V8=='False',]
  for (i in 1:length(x)){
    count1 = length(which(data_TF[,4]>x[i]))
    count2 = length(which(data_else[,4]>x[i]))
    count_table[1,i] = count1
    count_table[2,i] = count2
    matrix_ = matrix(c(count_table[1,i], length(data_TF[,1]) - count_table[1,i], count_table[2,i], length(data_else[,1]) - count_table[2,i]), nrow=2, byrow=TRUE)
    temp = fisher.test(matrix_, alternative='greater', or=1)
    result_table[epi,i] = temp$p.value
  }
}
p_function = function(a){
  t = p.adjust(a, method='fdr')
  return (t)
}
corrected_result = apply(result_table, 2, p_function)
write.table(t(corrected_result), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/threshold_FET_pvalue_standard_ALL_FDR.txt', quote=FALSE, sep='\t')
boxplot.matrix(as.matrix(-log10(corrected_result)), use.cols = TRUE, xlab='Width/SD',ylab='-log10(FDR))',main="TF enrichment analysis (Fisher's exact test)", outpch='.', col='light blue', outcol='light blue')

### Threshold analysis 2 by standard score (z-score)
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/'
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt')
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/intra/excluding_ncrnas/E_combined_mRNAs.txt')

data_ = data_[complete.cases(data_$V4),]  # remove NaN
data_ = data_[which(data_$V9=='True'),]  # only mRNAs
data_TF = data_[which(data_$V8=='True'),]
data_else = data_[which(data_$V8=='False'),]

# Create 'count' table (incl. statistics by Fisher's Exact Test)
# Fisher's exact test (odd-ratios & p-value, one-sided)
# Null hypothesis = broad widths are not assiciated with known TFs  
x = seq(1.5, -2.0, by=-0.1)
count_table = matrix(nrow=4, ncol=length(x))
count_table = data.frame(count_table)
colnames(count_table) = x
rownames(count_table) = c('TFs(domains)', 'non-TFs(domains)', 'Odd-ratio','p-value')
for (i in 1:length(x)){
  count1 = length(which(data_TF[,5]>x[i]))
  count2 = length(which(data_else[,5]>x[i]))
  count_table[1,i] = count1
  count_table[2,i] = count2
  matrix_ = matrix(c(count_table[1,i], length(data_TF[,1]) - count_table[1,i], count_table[2,i], length(data_else[,1]) - count_table[2,i]), nrow=2, byrow=TRUE)
  temp = fisher.test(matrix_, alternative='greater', or=1)
  count_table[3,i] = temp$estimate
  count_table[4,i] = temp$p.value
}
write.table(t(count_table), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/thrshold_FET_inter_5kb.txt', quote=FALSE, sep='\t')

# Draw side-by-side bar cumulative plots for distribution 
total_TF = nrow(data_TF)
total_else = nrow(data_else)
table_ = count_table
vec = vector()
for (i in 1:ncol(count_table)){
  table_[1,i] = count_table[1,i] / total_TF
  table_[2,i] = count_table[2,i] / total_else
  vec = c(vec,table_[1,i], table_[2,i])
}
vec = rev(vec)
x_label = seq(0.0, 4.0, 0.1)
main_ = c('Cumulative distribution of H3K4me3 peaks\n','(mRNA genes only)')
barplot(matrix(vec,nr=2), beside=T, 
        col=c("aquamarine3","coral"), 
        names.arg=x_label, main=main_,
        xlab='z-score', ylab='proportion of peaks')
legend("topright", c("NonTF","TF"), pch=15, 
       col=c("aquamarine3","coral"), 
       bty="n")

# Draw side-by-side bar histogram plots for distribution 

total_TF = nrow(data_TF)
total_else = nrow(data_else)
table_ = count_table
vec = vector()
prev1 = 0
prev2 = 0
for (i in 1:ncol(count_table)){
  table_[1,i] = (count_table[1,i] - prev1) / total_TF
  table_[2,i] = (count_table[2,i] - prev2) / total_else
  prev1 = count_table[1,i]
  prev2 = count_table[2,i]
  vec = c(vec,table_[1,i], table_[2,i])
}
vec = rev(vec)
x_label = seq(-1.5, 1.7, 0.1)
main_ = c('Distribution of H3K4me3 peaks\n','(mRNA genes only, intra-cell types)')

pdf("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/histogram_intra_width_BIG.pdf", width = 20, height = 20)
barplot(matrix(vec,nr=2), beside=T, 
        col=c("aquamarine3","coral"), 
        names.arg=x_label, main=main_,
        xlab='log10(width/mean)', ylab='proportion of peaks')
legend("topright", c("NonTF","TF"), pch=15, 
       col=c("aquamarine3","coral"), 
       bty="n")
dev.off()

# count numbers of TFs captured by a given threshold 
# calculate effectiveness of capturing (see 'metadata.txt')
threshold = seq(1.5, -2.0, by=-0.1) 
tf_table = matrix(nrow=5, ncol=length(threshold))
tf_table = data.frame(tf_table)
colnames(tf_table) = threshold
rownames(tf_table) = c('TFs','effectiveness(TFs)','non_TFs','effectiveness(Non-TFs)', 'TF_ratio/non_TF_ratio')
total_tf = unique(data_TF[,1])
total_non_tf = unique(data_else[,1])
total_tf_widths = length(data_TF[,1])
total_non_tf_widths = length(data_else[,1])
for (i in 1:length(threshold)){
  thre = threshold[i]
  tf_list=unique(data_TF[which(data_TF[,4]>thre),][,1])
  non_tf_list=unique(data_else[which(data_else[,4]>thre),][,1])
  tf.effect = (length(tf_list) / length(total_tf))  / (count_table[1,i] / total_tf_widths)
  non.tf.effect = (length(non_tf_list) / length(total_non_tf)) / (count_table[2,i] / total_non_tf_widths)
  tf_table[1,i] = length(tf_list)
  tf_table[3,i] = length(non_tf_list)
  tf_table[2,i] = tf.effect
  tf_table[4,i] = non.tf.effect
  tf_table[5,i] = (length(tf_list) / length(total_tf)) / (length(non_tf_list) / length(total_non_tf) )
}
write.table(t(tf_table), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/count_effectiveness2_width_mRNAs_intra.txt', quote=FALSE, sep='\t')

# calculate average change rates(tangent line for each section) above a certain threshold
# for each TF or gene above a threshold 

tf_list = data_TF[which(data_TF[,4]>0.5),]
non_tf_list = data_TF[which(data_else[,4]>0.5),]
data_TF$rate_change = rep('-', nrow(data_TF))
data_else$rate_change = rep('-', nrow(data_else))

for (tf in unique(tf_list[,1])){
  temp = tf_list[which(tf_list[,1]==tf),]
  temp = temp[order(temp$V4, decreasing = TRUE),]
  for (m in 1:nrow(temp)){
    if (m != nrow(temp)){
      rate = temp[m,4] - temp[m+1,4]
      data_TF[rownames(temp)[m],6] = rate
    } 
  }
}


tf_list = data_TF[which(data_TF[,4]>0.5),]
colnames(tf_list) = c('TF','cell_type','domain','Log10(Width/Mean)','p-value(AUC)','Rate_change')
write.table(tf_list, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/rate_change_TFs.txt', quote=FALSE, sep='\t')

######### RUN AGAIN
list_ = unique(non_tf_list[,1])
list_ = list_[1:843]
for (tf in list_){
  temp = non_tf_list[which(non_tf_list[,1]==tf),]
  temp = temp[order(temp$V4, decreasing = TRUE),]
  for (m in 1:nrow(temp)){
    if (m != nrow(temp)){
      rate = temp[m,4] - temp[m+1,4]
      data_else[rownames(temp)[m],6] = rate
    } 
  }
}

non_tf_list = data_else[which(data_else[,4]>=0.5),]
colnames(non_tf_list) = c('gene','cell_type','domain','Log10(Width/Mean)','p-value(AUC)','Rate_change')
write.table(non_tf_list, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/rate_change_non_TFs.txt', quote=FALSE, sep='\t')

### Histogram for widths in a cell type
### WITHIN A CELL TYPE ANALYSIS
### Fisher's exact test

data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/assigned_genes.txt')
col_ = unique(data_$V1)

x = seq(1.0, 0.0, by=-0.1)
x = round(x, digits=2)
or_table = matrix(nrow=length(x), ncol=length(col_))
or_table = data.frame(or_table)
colnames(or_table) = col_
rownames(or_table) = x
p_table = or_table

for (i in 1:length(x)){
  for (j in col_){
    temp = data_[which(data_[,1]==j),]
    temp1 = temp[which(temp[,5]>x[i]),]
    if (nrow(temp1)>0){
      count1 = length(which(temp1[,6]=='True'))
      count2 = nrow(temp1) - count1
      matrix_ = matrix(c(count1, length(which(temp[,6]=='True')) - count1, count2, nrow(temp)-length(which(temp[,6]=='True'))-count2), nrow=2, byrow=TRUE)
      temp_ = fisher.test(matrix_, alternative='greater', or=1)
      or_table[i, j] = temp_$estimate
      p_table[i, j] = temp_$p.value
    }
  }
}

write.table(t(p_table), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/p_table.txt', quote=FALSE, sep='\t')

### Heatmap 
library(gplots)
m = seq(1.0, 0.5, by=-0.1)
for (m_ in m){
  data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/selected_above_width_0.0.txt')
  data_ = data_[which(data_$V5>=m_),]
  col_ = unique(data_$V1)
  row_ = unique(data_$V2)
  col_ = c('E017','E002','E008','E001','E015','E014','E016','E003','E024','E020'
           ,'E019','E018','E021','E022','E007','E009','E010','E013','E012','E011'
           ,'E004','E005','E006','E062','E034','E045','E033','E044'
           ,'E043','E039','E041','E042','E040','E037','E048','E038','E047','E029'
           ,'E031','E035','E051','E050','E036','E032','E046','E030','E026','E049'
           ,'E025','E023','E052','E055','E056','E059','E061','E057','E058','E028'
           ,'E027','E054','E053','E112','E093','E071','E074','E068','E069','E072'
           ,'E067','E073','E070','E082','E081','E063','E100','E108','E107','E089'
           ,'E090','E083','E104','E095','E105','E065','E078','E076','E103','E111'
           ,'E092','E085','E084','E109','E106','E075','E101','E102','E110','E077'
           ,'E079','E094')
  data_table = matrix(nrow=length(row_), ncol=length(col_))
  data_table = data.frame(data_table)
  colnames(data_table) = col_
  rownames(data_table) = row_
  for (i in 1:nrow(data_)){
    epi = as.character(data_$V1[i])
    gene = as.character(data_$V2[i])
    value = as.numeric(data_$V5[i])
    data_table[gene,epi] = value
  }
  file = ''
  file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/genes_table_width_',m_,'.txt',sep='')
  write.table(data_table, file=file_, quote=FALSE, sep='\t')
  
  tf_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/data_TF.txt')
  tf_ = unique(tf_$V1)
  data_table_tf = data_table[which(rownames(data_table) %in% tf_),] 
  
  file_ = ''
  file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/genes_table_width_',m_,'.pdf',sep='')
  pdf(file_, width = 20, height = 100)
  title_ = paste('TFs (log10(width/mean)>',m,sep='')
  heatmap.2(as.matrix(data_table_tf), main=title_, dendrogram='none', trace = 'none', Rowv=FALSE, Colv=FALSE, col=hmcol, na.color='black',labRow = rownames(data_table_tf), labCol=colnames(data_table_tf))
  dev.off()
}

data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/selected_above_width_0.0.txt')
data_ = data_[which(data_$V5>=1.0),]
col_ = unique(data_$V1)
row_ = unique(data_$V2)
col_ = c('E017','E002','E008','E001','E015','E014','E016','E003','E024','E020'
         ,'E019','E018','E021','E022','E007','E009','E010','E013','E012','E011'
         ,'E004','E005','E006','E062','E034','E045','E033','E044'
         ,'E043','E039','E041','E042','E040','E037','E048','E038','E047','E029'
         ,'E031','E035','E051','E050','E036','E032','E046','E030','E026','E049'
         ,'E025','E023','E052','E055','E056','E059','E061','E057','E058','E028'
         ,'E027','E054','E053','E112','E093','E071','E074','E068','E069','E072'
         ,'E067','E073','E070','E082','E081','E063','E100','E108','E107','E089'
         ,'E090','E083','E104','E095','E105','E065','E078','E076','E103','E111'
         ,'E092','E085','E084','E109','E106','E075','E101','E102','E110','E077'
         ,'E079','E094')
data_table = matrix(nrow=length(row_), ncol=length(col_))
data_table = data.frame(data_table)
colnames(data_table) = col_
rownames(data_table) = row_
for (i in 1:nrow(data_)){
  epi = as.character(data_$V1[i])
  gene = as.character(data_$V2[i])
  value = as.numeric(data_$V5[i])
  data_table[gene,epi] = value
}
write.table(data_table, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/genes_table_width_1.0.txt', quote=FALSE, sep='\t')

tf_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/data_TF.txt')
tf_ = unique(tf_$V1)
data_table_tf = data_table[which(rownames(data_table) %in% tf_),] 


library(RColorBrewer)
hmcol = colorRampPalette(c("white","red"))(256)
pdf("selected_TF_width_1.0.pdf", width = 20, height = 100)
heatmap.2(as.matrix(data_table_tf), main='TFs (log10(width/mean)>1.0)', dendrogram='none', trace = 'none', Rowv=FALSE, Colv=FALSE, col=hmcol, na.color='black',labRow = rownames(data_table_tf), labCol=colnames(data_table_tf))
dev.off()

### Miscellenous analysis
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_group_specific_table_0.3_5kb_.txt')
data_ = data_[which(rownames(data_) %in% tf_),] 
pdf("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_specific_above0.5.pdf", width = 20, height = 100)
heatmap(as.matrix(data_), col=hmcol)
pdf("Tissue_group_specific_0.3_5kb_big_.pdf", width = 20, height = 50)
heatmap.2(as.matrix(data_),  trace = 'none', col=hmcol, na.color='black')
dev.off()

##### Calculate gain of TFs per an added domain to deal with
##### This is to find the most efficient threshold to cover most TFs with minimum
##### number of domains. This also indicates a threshold where it significantly 
##### start to fill genes 'already' marked by a wider domain. 
##### We need both 'count_table' & 'tf_table' for this.
count_table1 = data.frame(t(count_table))
tf_table1 = data.frame(t(tf_table))
result_table = data.frame(matrix(nrow=nrow(count_table1), ncol=2))
rownames(result_table) = rownames(count_table1)
colnames(result_table) = c('TF','NonTF')

prev1_tf = 0
prev1_non_tf = 0
prev2_tf = 0
prev2_non_tf = 0
for (i in 1:nrow(count_table1)){
  diff1_tf = count_table1[i,1] - prev1_tf
  diff1_non_tf = count_table1[i,2] - prev1_non_tf
  diff2_tf = tf_table1[i,1] - prev2_tf
  diff2_non_tf = tf_table1[i,3] - prev2_non_tf
  prev1_tf = count_table1[i,1]
  prev1_non_tf = count_table1[i,2]
  prev2_tf = tf_table1[i,1]
  prev2_non_tf = tf_table1[i,3]
  result_table[i,1] = diff2_tf / diff1_tf
  result_table[i,2] = diff2_non_tf / diff1_non_tf
}

write.table(result_table, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/coverage_inter.txt', quote=FALSE, sep='\t')

# Plot 
x_range = range(as.numeric(rownames(result_table)))
y_range = range(complete.cases(result_table))
x_value = as.numeric(rownames(result_table))
plot(x_range, y_range, type='n', xlab='log10(width/weighted mean)', ylab='Coverage score', xlim=rev(x_range))
colors = c('red', 'darkblue')
for (i in 1:ncol(result_table)){
  lines(x_value, result_table[,i], col=colors[i], type='b',
        lwd=1.5, lty=1)
}
title('Coverage score (Inter-cell type)')
legend('topright', cex=0.8, col=colors, lty=1, colnames(result_table))

### Analyse dynamics for a given gene across cell types
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb_ALL.txt')

# Example1: HAND2                
name= 'HAND2'
epi = c('E083','E095','E065','E104','E105')
example1 = data_[data_$V1==name,]
example1 = example1[order(example1$V4, decreasing=TRUE),]

#If want to use -log10(p-values)
example1$V6 = -log10(example1$V6)

plot(x=rev(seq(1, nrow(example1), 1)), y=rev(example1$V5), xlab='Rank position',ylab='log10(width/mean)', main=name) # by raw width

x_ = c(1, nrow(example1))
y_ = c(example1$V5[1], example1$V5[nrow(example1)])
lines(x_, y_, type='l', col='red')

# Or linear regression line
abline(lm(example1$V5 ~ seq(1, nrow(example1), 1)), col='red')

# labels
labels_ = example1$V2[1:10]
col_ = ifelse(labels_ %in% epi, 'red', 'black')
a=4
text(x=seq(1+a,length(labels_)+a,1), y=example1$V5[1:length(labels_)], labels=labels_, col=col_, cex=0.6)


### Analyse dynamics for genes within a cell type
data_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant.txt')

# Example1: HAND2                
name= 'E083'
epi = c('E083','E095','E065','E104','E105')
example1 = data_[data_$V2==name,]
example1 = example1[order(example1$V4, decreasing=TRUE),]

#If want to use -log10(p-values)
example1$V6 = -log10(example1$V6)

plot(x=rev(seq(1, nrow(example1), 1)), y=rev(example1$V4), xlab='Rank (z-score)',ylab='raw width', main=name) # by raw width
x_ = c(1, nrow(example1))
y_ = c(example1$V4[1], example1$V4[nrow(example1)])
lines(x_, y_, type='l', col='red')

# labels
labels_ = example1$V2[1:10]
col_ = ifelse(labels_ %in% epi, 'red', 'black')
a=4
text(x=seq(1+a,length(labels_)+a,1), y=example1$V5[1:length(labels_)], labels=labels_, col=col_, cex=0.6)

### Randomly select widths given a variance (mean is always 0)
### and plot 
data_1 = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_means_variances_dominant.txt")

name = 'HAND2'
example2 = data_1[data_1$V1==name,]

no_points = nrow(example1)
random_background = rnorm(no_points, 0, example2$V3)
random_background = sort(random_background, decreasing=T)

name_ = paste(name, '_randomly selected')
plot(x=rev(seq(1, length(random_background), 1)), y=rev(random_background), xlab='Rank',ylab='log10(width/mean)', main=name_) # by raw width

abline(lm(random_background ~ seq(1, length(random_background), 1)), col='red')

### AUTOMATE THE DYNAMIC ANALYSIS for a given set of genes (specified)

data_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb.txt')
data_1 = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_means_variances_dominant.txt")
genes = c('HAND2','NKX2-5','TBX5','GATA4','GATA6','MEIS2','ISL1','MYL4','COL1A1','CKM','TNNI1','ATP2A2')
data__ = data_
var_ = var(data__$V5)
pdf("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/results/Dynamics_plots_all.pdf")
for (name in genes){
  example1 = data__[data__$V1==name,]
  example1 = example1[order(example1$V5, decreasing=TRUE),]
  example2 = data_1[data_1$V1==name,]
  no_points = nrow(example1)
  random_background = rnorm(no_points, 0, sqrt(var_))
  random_background = sort(random_background, decreasing=T)
  x_ = rev(seq(1, nrow(example1), 1))
  y_lim = range(max(example1$V5), min(example1$V5), max(random_background), min(random_background))
  plot(x=x_, y=rev(example1$V5), ylim=y_lim, xlim=range(x_), xlab='Rank position',ylab='log10(width/mean)', main=name, col='red') # by raw width
  points(x=x_, y=rev(random_background), col='blue')
  abline(lm(example1$V5 ~ seq(1, nrow(example1), 1)), col='red')
  abline(lm(random_background ~ seq(1, length(random_background), 1)), col='blue')
  legend('topright',c('data','bg'), col=c('red','blue'), pch=1, cex=0.8)
}
dev.off()

### Plot background bands
data_ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_2.5kb.txt")
data_$position = rownames(data_)
data_ = data_[order(as.numeric(data_$position)),]
data_$position = NULL

data__ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_tf_2.5kb.txt")
data__$position = rownames(data__)
data__ = data__[order(as.numeric(data__$position)),]
data__$position = NULL

data___ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_nontf_2.5kb.txt")
data___$position = rownames(data___)
data___ = data___[order(as.numeric(data___$position)),]
data___$position = NULL
#rownames(data_) = paste('Rank_',seq(0,101), sep='')
#data__ = t(data_[,1:ncol(data_)])
tf_ = c('HAND2','TBX5','ISL1','GATA4','GATA6','NKX2.5')
str_= c('MYL4','CKM','TNNI1','MYH6','PLN','MYL3','MYL2')
groups = c(tf_, str_)

file_name = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/','Inter-dynamics_structural_cardiac.pdf',sep='')
group = tf_
#pdf(file_name)
data_ = log10(data_)
boxplot.matrix(as.matrix(data_), use.cols = FALSE, xlab='Rank',ylab='log10(width/SD)',main='Cardiac lineage regulators', outpch='.')
cols_ = rainbow(length(group))
xx = seq(1,nrow(data_))
for (i in 1:length(group)){
  yy = data_[,group[i]]
  lines(xx,yy,col=cols_[i], lwd=2)
}
legend(x=115, y=1.8, group, col = cols_, lty=1, lwd=2, cex=0.5, bty='n')
#dev.off()



boxplot.matrix(as.matrix(data__), use.cols = FALSE, xlab='Rank',ylab='log10(width/mean)',main='TFs vs. non-TFs')
cols_ = rainbow(ncol(data__))
xx = seq(1,nrow(data__))
for (i in 1:ncol(data__)){
  yy = data__[,i]
  lines(xx,yy,col=cols_[i])
}

### Plot ranges 
tt_factor = rep('TF',ncol(data__))
tt_factor = c(tt_factor, rep('non-TF',ncol(data___)))
tt = cbind(data__, data___)

### Hierarchical clustering of selected top xxx genes in a cell type
data_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant__.txt')
# If slope of the change is to be used, use this. 
data_ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_slope.txt")

sample_data_1 = data_1[data_1$V2=='E095',]
sample_data_1 = sample_data_1[order(sample_data_1$V5, decreasing = TRUE),]
sample_data_1_selected = sample_data_1[1:200,]
gene_names = sample_data_1_selected$V1
new_data = data_[,as.character(gene_names)]

# If slope of the change is to be used, use this. 

genes.cor = cor(new_data, use='pairwise.complete.obs', method='pearson')
genes.cor.dist = as.dist(1-genes.cor)
genes.tree = hclust(genes.cor.dist, method='average')

library(gplots)
hmcol = colorRampPalette(c("white","red"))(256)
heatmap.2(as.matrix(new_data), Rowv=FALSE, Colv=TRUE, col=hmcol, dendrogram='none',trace='none', labCol=colnames(new_data), labRow=rownames(new_data))

### Plot side-by-side boxplots (TFs vs. non-TFs)
# If only mRNA genes are to be considered..
data_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt')
data___ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_nontf_2.5kb.txt")
data___$position = rownames(data___)
data___ = data___[order(as.numeric(data___$position)),]
data___$position = NULL
tf_list = data_1[data_1$V9=='True',]
tf_list = unique(tf_list$V1)

data__ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/sorted_table_mrna_tf_2.5kb.txt")
data__$position = rownames(data__)
data__ = data__[order(as.numeric(data__$position)),]
data__$position = NULL


y_range_ = range(max(data__, na.rm=TRUE), min(data__, na.rm=TRUE), max(data___, na.rm=TRUE), min(data___, na.rm=TRUE))
boxplot.matrix(as.matrix(data__), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:100-0.3, col='green',ann=FALSE,xlab='Rank', ylab='Log10(width/mean)', outpch='.', outcol='green')
boxplot.matrix(as.matrix(data___), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:100, col='purple',ann=FALSE,xaxt='n', yaxt='n',add=TRUE, outpch='.', outcol='purple')
legend(x=79, y=1.5, cex=0.8, c('TFs(n=1,347)','non-TFs(n=16,505)'), col=c('green','purple'), pch='.', bty='n')

### Plot histograms number of cell types (TFs vs. non-TFs)
data_ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/cell_type_distributions_5kb.txt")
data_ = data_[data_$V7=='True',]  #mRNA only


data_TF = data_[data_$V6=='True',]
data_nonTF = data_[data_$V6=='False',]

data_ = data_TF
t_factor = as.factor(data_$V6)

## Clustering analysis (number of cell types & difference in log10 as inputs)
t = data.frame(cell_types=data_$V2, diff=data_$V5)
t.std = scale(t, center=TRUE, scale=TRUE)
rownames(t.std) = data_$V1
diff.cor = cor(t(t.std), use='pairwise.complete.obs', method='pearson')
diff.cor.dist = as.dist(1-diff.cor)
diff.tree = hclust(diff.cor.dist, method='average')

features.cor = cor(t.std, use='pairwise.complete.obs', method='spearman')
features.cor.dist = as.dist(1-features.cor)
features.tree = hclust(features.cor.dist, method='average')

library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
colours_ = ifelse(t_factor=='True', 'red', 'blue')
heatmap(t(t.std), col=col, ColSideColors = colours_)

#hist(data_TF$V2, freq=FALSE)

h_tf = hist(data_TF$V2, breaks=100, freq=FALSE)
h_nontf = hist(data_nonTF$V2, breaks=100,freq=FALSE)
results = vector()
temp_tf = 0
temp_nontf = 0
for (i in 1:99){
  temp_tf = temp_tf + h_tf$density[i] 
  temp_nontf = temp_nontf + h_nontf$density[i]
  temp = temp_tf / temp_nontf
  results = c(results, temp)
}
vec = data.frame(tf= as.numeric(h_tf$density), nontf = as.numeric(h_nontf$density))
x_label = seq(1,99)
main_ = 'Histograms'
barplot(matrix(vec, nr=2), beside=T, 
        col=c("aquamarine3","coral"), 
        names.arg=x_label, main=main_,
        xlab='Rank Bins', ylab='Proportion')
legend("topright", c("TF","nonTF"), pch=15, 
       col=c("aquamarine3","coral"), 
       bty="n")

#### Fisher's exact test
no_rank = 100
x = seq(1, no_rank, by=1)
count_table = matrix(nrow=4, ncol=length(x))
count_table = data.frame(count_table)
colnames(count_table) = x
rownames(count_table) = c('TFs', 'non-TFs', 'Odd-ratio','p-value')
for (i in 1:length(x)){
  count1 = length(which(data_TF$V2<=x[i]))
  count2 = length(which(data_nonTF$V2<=x[i]))
  count_table[1,i] = count1
  count_table[2,i] = count2
  matrix_ = matrix(c(count_table[2,i], nrow(data_nonTF) - count_table[2,i], count_table[1,i], nrow(data_TF) - count_table[1,i]), nrow=2, byrow=TRUE)
  temp = fisher.test(matrix_, alternative='greater', or=1)
  count_table[3,i] = temp$estimate
  count_table[4,i] = temp$p.value
}

write.table(t(count_table), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Dynamics_table.txt', quote=FALSE, sep='\t')

plot(x, -log10(count_table[4,]), xlab='Rank', ylab='-log10(p-value)', main='Enrichment of non-TF genes')

#### Calculate dynamics in three sections (0~10%, 10~90%, 90~100%)
#### change rate = (max. - min.) / number of cell types in the section 
data_TF = data_1[data_1$V8=='True',]
data_nonTF = data_1[data_1$V8=='False',]
data_nonTF = data_nonTF[data_nonTF$V9=='True',]
data_ = rbind(data_TF, data_nonTF)
results_ = c('name','section1','section2','section3','no_cell_types','Dynamics_difference','Max_dynamics','Min_dynamics','Section1_mean','Section2_mean','Section3_mean','TF')
for (gene in unique(data_$V1)){
  temp = data_[data_$V1==gene,]
  temp = temp[order(temp$V5, decreasing='True'),]
  no_terms = nrow(temp)
  if (no_terms > 3){
    interval = ceiling(no_terms * 0.05) 
    top_ = as.numeric((temp$V5[1] - temp$V5[interval+1]) / (interval+1))
    top_mean = as.numeric((sum(temp$V5[1:(interval+1)])) / (interval+1))
    middle_ = as.numeric((temp$V5[interval+1] - temp$V5[no_terms-interval]) / (no_terms - 2*interval+1))
    middle_mean = as.numeric((sum(temp$V5[(interval+1):(no_terms-interval)])) / (no_terms - 2*interval+1))
    bottom_ = as.numeric((temp$V5[no_terms-interval] - temp$V5[no_terms]) / (interval+1))
    bottom_mean = as.numeric((sum(temp$V5[(no_terms-interval):no_terms])) / (interval+1))
    results_ = rbind(results_, c(gene, top_, middle_, bottom_, no_terms, max(temp$V5) - min(temp$V5), max(temp$V5), min(temp$V5), top_mean, middle_mean, bottom_mean, as.character(temp$V8[1])))
  }
}

write.table(results__, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Dynamics_total_5kb_new.txt', quote=FALSE, sep='\t')

results_nonTF = c('name','section1','section2','section3','no_cell_types')
for (gene in unique(data_nonTF$V1)){
  temp = data_nonTF[data_nonTF$V1==gene,]
  temp = temp[order(temp$V5, decreasing='True'),]
  no_terms = nrow(temp)
  if (no_terms > 3){
    interval = ceiling(no_terms * 0.1) 
    top_ = (temp$V5[1] - temp$V5[interval+1]) / interval
    middle_ = (temp$V5[interval+1] - temp$V5[no_terms-interval]) / (no_terms - 2*interval - 1)
    bottom_ = (temp$V5[no_terms-interval] - temp$V5[no_terms]) / interval
    results_nonTF = rbind(results_nonTF, c(gene, top_, middle_, bottom_, no_terms))
  }
}

write.table(results_nonTF, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Dynamics_3sections_nonTF.txt', quote=FALSE, sep='\t')


### Divide values by the section 2 (in each gene)
results_TF = data.frame(results_TF)
results_TF = results_TF[2:nrow(results_TF),]
results_nonTF = data.frame(results_nonTF[2:nrow(results_nonTF),])

results_TF$X2 = as.numeric(as.character(results_TF$X2))
results_TF$X3 = as.numeric(as.character(results_TF$X3))
results_TF$X4 = as.numeric(as.character(results_TF$X4))
results_TF$X5 = as.numeric(as.character(results_TF$X5))
results_nonTF$X2 = as.numeric(as.character(results_nonTF$X2))
results_nonTF$X3 = as.numeric(as.character(results_nonTF$X3))
results_nonTF$X4 = as.numeric(as.character(results_nonTF$X4))
results_nonTF$X5 = as.numeric(as.character(results_nonTF$X5))
results_TF$X6 = rep('TF', nrow(results_TF))
results_nonTF$X6 = rep('nonTF', nrow(results_nonTF))

results_TF$X2 = as.numeric(as.character(results_TF$X2)) / as.numeric(as.character(results_TF$X3))

results_TF$X4 = as.numeric(as.character(results_TF$X4)) / as.numeric(as.character(results_TF$X3))

results_TF$X3 = as.numeric(as.character(results_TF$X3)) / as.numeric(as.character(results_TF$X3))

results_TF$X5 = as.numeric(as.character(results_TF$X5))




results_nonTF$X2 = as.numeric(as.character(results_nonTF$X2)) / as.numeric(as.character(results_nonTF$X3))

results_nonTF$X4 = as.numeric(as.character(results_nonTF$X4)) / as.numeric(as.character(results_nonTF$X3))

results_nonTF$X3 = as.numeric(as.character(results_nonTF$X3)) / as.numeric(as.character(results_nonTF$X3))

results_nonTF$X5 = as.numeric(as.character(results_nonTF$X5))




final_table = matrix(nrow=nrow(results_TF)+nrow(results_nonTF), ncol=3)
colnames(final_table) = c('section1', 'section2', 'section3')
rownames(final_table) = c(as.character(results_TF$X1), as.character(results_nonTF$X1))
final_table_factor = as.factor(c(results_TF$X6, results_nonTF$X6))
final_table[,1] = c(results_TF$X2, results_nonTF$X2)
final_table[,2] = c(results_TF$X3, results_nonTF$X3)
final_table[,3] = c(results_TF$X4, results_nonTF$X4)

final_table = matrix(nrow=nrow(results_TF), ncol=2)
colnames(final_table) = c('section1', 'section3')
rownames(final_table) = as.character(results_TF$X1)
final_table[,1] = results_TF$X2
#final_table[,2] = results_TF$X3
final_table[,2] = results_TF$X4
heatmap(as.matrix(t(final_table)), col=col)

library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
colours_ = ifelse(final_table_factor=='TF', 'red', 'blue')
heatmap(as.matrix(t(final_table)), col=col, ColSideColors = colours_)

# Further analysis  (Fisher's exact test for enrichment of TFs by dynamics difference)
data_2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/cell_type_distributions.txt')
data_3 = data_2[data_2$V7=='True',]
data_3_TF = data_3[data_3$V6=='True',]
data_3_nonTF = data_3[data_3$V6=='False',]

x = seq(2.5, 0.0, by=-0.1)
count_table = matrix(nrow=4, ncol=length(x))
count_table = data.frame(count_table)
colnames(count_table) = x
rownames(count_table) = c('TFs(domains)', 'non-TFs(domains)', 'Odd-ratio','p-value')
for (i in 1:length(x)){
  count1 = length(which(data_3_TF[,5]>=x[i]))
  count2 = length(which(data_3_nonTF[,5]>=x[i]))
  count_table[1,i] = count1
  count_table[2,i] = count2
  matrix_ = matrix(c(count_table[1,i], length(data_3_TF[,1]) - count_table[1,i], count_table[2,i], length(data_3_nonTF[,1]) - count_table[2,i]), nrow=2, byrow=TRUE)
  temp = fisher.test(matrix_, alternative='greater', or=1)
  count_table[3,i] = temp$estimate
  count_table[4,i] = temp$p.value
}

write.table(t(count_table), file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Dynamics_difference_over_80_cell_types_FET.txt', quote=FALSE, sep='\t')

### Analysis with 'Dynamics_total.txt'
results__ = results_[2:nrow(results_),]
results__ = data.frame(results__)
results__ = results


results__$X2 = as.numeric(as.character(results__$X2))
results__$X3 = as.numeric(as.character(results__$X3))
results__$X4 = as.numeric(as.character(results__$X4))
results__$X5 = as.numeric(as.character(results__$X5))
results__$X6 = as.numeric(as.character(results__$X6))
results__$X7 = as.numeric(as.character(results__$X7))
results__$X8 = as.numeric(as.character(results__$X8))
results__$X9 = as.numeric(as.character(results__$X9))
results__$X10 = as.numeric(as.character(results__$X10))
results__$X11 = as.numeric(as.character(results__$X11))
# colnames(results__) = as.character(results[1,])

write.table(temp, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Selected_genes.txt', quote=FALSE, sep='\t')


### Read in data table
data_file = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb_tissue.txt')
epis = unique(data_file$V2)
for (epi in epis){
  temp = data_file[data_file$V2==epi,]
  temp = temp[order(as.numeric(temp$V5), decreasing=TRUE),]
  temp = temp[temp$V11=='True',]
  name = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Tissue_group_clustering/intersection/',epi,'.txt')
  write.table(temp, file=name, sep='\t', quote = FALSE)
}

data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Fold_enrich_0.3_5kb_.txt')
data__ = data_[complete.cases(data_),]
heatmap(t(as.matrix(data1)), ColSideColors = col__)

### Read in 'Tissue_Enrichment_Scores.txt' and plot heatmap
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_Enrichment_Scores_pvalues.txt')
labels_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_list_training.txt')
cols = rainbow(length(unique(labels_[,2])))

# For p-values
data1 = -log10(data1)
data1.std = data1
data1.std= scale(data1, scale=TRUE, center=TRUE)

library(gplots)
library(RColorBrewer)
hmcol = colorRampPalette(c("white","red"))(256)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/plots/Heatmap_tissue_enrichment_scores_trainingsets_big.pdf',width = 30, height = 20 )
t=heatmap.2(as.matrix(t(data1.std)), col=hmcol, dendrogram='none',trace='none', main='Tissue Enrichment Scores (Training sets)', labCol=rownames(data1.std), labRow=colnames(data1.std), cexRow=0.6)
col_idx = rownames(data1)[t$colInd]
dev.off()

heatmap(as.numeric(data1), ColSideColors = cols)

col__ = vector()
for (i in col_idx){
  if (labels_[which(labels_[,1]==i), 2] == 'Heart'){
    col__ = c(col__, cols[1])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Brain'){
    col__ = c(col__, cols[2])
  }
  if (labels_[which(labels_[,1]==i), 2] == 'Epithelial'){
    col__ = c(col__, cols[3])
  }
  if (labels_[which(labels_[,1]==i), 2] == 'Adipose'){
    col__ = c(col__, cols[4])
  }
  if (labels_[which(labels_[,1]==i), 2] == 'Blood&T_cell'){
    col__ = c(col__, cols[5])
  }
  if (labels_[which(labels_[,1]==i), 2] == 'ES_deriv.'){
    col__ = c(col__, cols[6])
  }  
  if (labels_[which(labels_[,1]==i), 2] == 'ES_cell'){
    col__ = c(col__, cols[7])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'IMR90'){
    col__ = c(col__, cols[8])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'iPSC'){
    col__ = c(col__, cols[9])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Muscle'){
    col__ = c(col__, cols[10])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Digestive'){
    col__ = c(col__, cols[11])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Thymus'){
    col__ = c(col__, cols[12])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'HSC&B_cell'){
    col__ = c(col__, cols[13])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Mesench.'){
    col__ = c(col__, cols[14])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Smooth_muscle'){
    col__ = c(col__, cols[15])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Neurosph.'){
    col__ = c(col__, cols[16])
  } 
  if (labels_[which(labels_[,1]==i), 2] == 'Myosat.'){
    col__ = c(col__, cols[17])
  } 
}

heatmap.2(as.matrix(t(data1)), col=hmcol, trace='none', main='Tissue Enrichment Scores (Training sets)', labCol=rownames(data1),  labRow=colnames(data1), cexRow=0.6, ColSideColors = col__)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")
cl.row <- hclustfunc(distfunc(t(data1)))
cl.col <- hclustfunc(distfunc(data1))
gr.row <- cutree(cl.row, 17)
gr.col <- cutree(cl.col, 17)
col1 <- rainbow(17)
col2 <- rainbow(17)
col3 = rainbow(17)
heatmap.2(as.matrix(t(data1)), col=hmcol, hclustfun=hclustfunc, distfun=distfunc, trace='none',RowSideColors=col1[gr.col], ColSideColors=col2[gr.row])

### Calculate distance matrix between training set and test set
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_Enrichment_Scores.txt')
data2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Tissue_Enrichment_Scores_TestSet.txt')
results = matrix(nrow=nrow(data2), ncol=nrow(data1))
results = data.frame(results)
rownames(results) = rownames(data2)
colnames(results) = rownames(data1)
data1.std = scale(data1, scale=TRUE, center=TRUE)
data2.std = scale(data2, scale=TRUE, center=TRUE)
for (i in 1:nrow(data2)){
  for (j in 1:nrow(data1)){
    results[i,j] = dist(rbind(data2.std[i,],data1.std[j,]), method='euclidean')
  }
}
a=as.matrix(results)
temp = heatmap(a)
temp_labels = labels_[temp$colInd,2]
mycol = vector(length=ncol(a))
mycol = cols[as.numeric(temp_labels)]
hmcol = colorRampPalette(c("white","red"))(256)
heatmap.2(as.matrix(t(results)), col=hmcol, dendrogram='none',trace='none', main='Euclidean distance', labCol=rownames(results), labRow=colnames(results), cexRow=0.6)

### Clustering analysis of genes based on fold enrichment 
data_4 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/Results/Fold_enrich_0.3_5kb_.txt')
data_4 = data_4[complete.cases(data_4),]
data_4 = scale(data_4, scale=TRUE, center=TRUE)

hmcol = colorRampPalette(c("white","red"))(256)
heatmap.2(as.matrix(data_4), col=hmcol, dendrogram='none',trace='none', main='Hierarchical clustering (genes vs. tissues) ', labCol=colnames(data_4), cexRow=0.6)

heatmap(as.matrix(data_4))

### Plot probability plots (p(c=TF | x > m)) where e.g. x is dynamics difference, m is a threshold)
data_5 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/Dynamics_total_5kb_.txt')
data_5 = data.frame(data_5)

results__ = results_
results__ = data.frame(results__)
col.name = results__[1,]

results__$X2 = as.numeric(as.character(results__$X2))
results__$X3 = as.numeric(as.character(results__$X3))
results__$X4 = as.numeric(as.character(results__$X4))
results__$X5 = as.numeric(as.character(results__$X5))
results__$X6 = as.numeric(as.character(results__$X6))
results__$X7 = as.numeric(as.character(results__$X7))
results__$X8 = as.numeric(as.character(results__$X8))
results__$X9 = as.numeric(as.character(results__$X9))
results__$X10 = as.numeric(as.character(results__$X10))
results__$X11 = as.numeric(as.character(results__$X11))
results__ = results__[2:nrow(results__),]

results__$X13 = results__$X2 / results__$X3
results__$X14 = results__$X4 / results__$X3

col.name = c('name','section1','section2','section3','no_cell_types','Dynamics_difference','Max_dynamics','Min_dynamics','Section1_mean','Section2_mean','Section3_mean','Section1/Section2','Section3/Section2')
#colnames(results__) = col.name
# 1. Dynamics difference between TF vs non-TF
# p = proportion of a set with > threshold
p_tf = vector()
p_non_tf = vector()
for (i in seq(0.0, 2.5, by=0.1)){
  temp_tf = results__[results__$X12=='True',]
  temp_non_tf = results__[results__$X12=='False',]
  p1 = length(which(temp_tf$X6>i)) / length(which(results__$X6>i))
  p2 = length(which(temp_non_tf$X6>i)) / length(which(results__$X6>i))
  p_tf = c(p_tf, p1)
  p_non_tf = c(p_non_tf, p2)
}
y_lim = range(p_tf, p_non_tf)
plot(x=seq(0.0, 2.5, by=0.1), y=p_tf, ylim=y_lim, col='dark green', type='o', xlab='Dynamics difference', ylab='Probability')
lines(x=seq(0.0, 2.5, by=0.1), y=p_non_tf, type='o', col='purple')

# 2. Dynamics difference between TF vs non-TF (Cumulative)
# p = proportion of a set with < threshold
p_tf = vector()
p_non_tf = vector()
for (i in seq(0.1, 3.0, by=0.1)){
  temp_tf = results__[results__$X12=='True',]
  temp_non_tf = results__[results__$X12=='False',]
  p1 = length(which(temp_tf$X6<i)) / nrow(temp_tf)
  p2 = length(which(temp_non_tf$X6<i)) / nrow(temp_non_tf)
  p_tf = c(p_tf, p1)
  p_non_tf = c(p_non_tf, p2)
}
y_lim = range(p_tf, p_non_tf)
plot(x=seq(0.1, 3.0, by=0.1), y=p_tf, ylim=y_lim, col='dark green', type='o', xlab='Dynamics difference', ylab='Cumulative proportion')
lines(x=seq(0.1, 3.0, by=0.1), y=p_non_tf, type='o', col='purple')
legend(x=2.3, y=0.2, cex=0.7, c('TFs(n=1,291)','Non-TFs(n=15,270)'), col=c('dark green','purple'), lty=1)

# 3. Plot histograms comparing dynamics difference (TFs vs. non-TFs) 
# and calculate 95% CI (assuming gaussian density) based on sample mean +/- SE & others
library(lattice)
histogram(~X14 | X12, data=results__, xlab='Dynamics Difference', type='density')

feature = 14
mean_tf = mean(temp_tf[,feature])
mean_non_tf = mean(temp_non_tf[,feature])
se_tf = sqrt(var(temp_tf[,feature])) / sqrt(nrow(temp_tf))
ci_tf = c(mean_tf-1.96*se_tf, mean_tf+1.96*se_tf)

se_non_tf = sqrt(var(temp_non_tf[,feature])) / sqrt(nrow(temp_non_tf))
ci_non_tf = c(mean_non_tf-1.96*se_non_tf, mean_non_tf+1.96*se_non_tf)

t.test(temp_tf[,feature], temp_non_tf[,feature])
wilcox.test(temp_tf[,feature], temp_non_tf[,feature])

# 4. No.cell type plot (cumulative)
# p = proportion of a set with < threshold
p_tf = vector()
p_non_tf = vector()
for (i in seq(1, 100, by=1)){
  temp_tf = results__[results__$X12=='True',]
  temp_non_tf = results__[results__$X12=='False',]
  p1 = length(which(temp_tf$X5<=i)) / nrow(temp_tf)
  p2 = length(which(temp_non_tf$X5<=i)) / nrow(temp_non_tf)
  p_tf = c(p_tf, p1)
  p_non_tf = c(p_non_tf, p2)
}
y_lim = range(p_tf, p_non_tf)
plot(x=seq(1, 100, by=1), y=p_tf, ylim=y_lim, col='dark green', type='o', xlab='Number of cell types', ylab='Cumulative proportion')
lines(x=seq(1, 100, by=1), y=p_non_tf, type='o', col='purple')
legend(x=0.0, y=1.0, cex=0.7, c('TFs(n=1,291)','Non-TFs(n=15,270)'), col=c('dark green','purple'), lty=1)

# 5. Logistic Regression
# using cross-validation (5-fold)
k = 2
interval_tf = round(nrow(temp_tf) / k)
interval_non_tf = round(nrow(temp_non_tf) / k)
testset = rbind(temp_tf[1:interval_tf,], temp_non_tf[1:interval_tf,])
trainingset = rbind(temp_tf[(interval_tf+1):nrow(temp_tf),], temp_non_tf[(interval_tf+1):nrow(temp_tf),])
trainingset$X12 = ifelse(trainingset$X12=='True',1,0)
testset$X12 = ifelse(testset$X12=='True',1,0)
test.glm = glm(X12 ~ X5+X6, family=binomial(link='logit'), data=trainingset)
summary(test.glm)

prediction = predict(test.glm, newdata=subset(testset, select=c(5,6)), type='response')
fitted.results = ifelse(prediction > 0.5, 1, 0)
error = mean(fitted.results != testset$X12)

#### NEW WAY OF ANALYSIS 
#### WITHOUT CORRECTING BY INCORPORATING WEIGHTS FOR TISSUE GROUPS
#### USING 127 ROADMAP CELL TYPES
data_new = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_5kb_ALL_standard.txt')
epi_list = unique(data_new$V2)
gene_list = unique(data_new$V1)
threshold = 0.3  # DS threshold

## TOO SLOW TO CREATE THE TABLE IN R
## BETTER TO PROCESS ALL DATA USING PYTHON AND READ IN THE TEXT FILE
### THUS, SKIP THIS
summary_table = data.frame(matrix(nrow=length(gene_list), ncol=length(epi_list)))
colnames(summary_table) = epi_list
rownames(summary_table) = gene_list
for (i in 1:nrow(data_new)){  
  epi = as.character(data_new[i,2])
  gene = as.character(data_new[i,1])
  summary_table[epi,gene] = data_new[i,5]
}
##############

spearman_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_spearmans_cor_standard.txt')
pdf(file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/plots/Spearman_between_cell_types1_normalised.pdf',width = 20, height = 20)
heatmap(as.matrix(spearman_))
dev.off()

library(gplots)
library(RColorBrewer)
pdf(file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/plots/Spearman_between_cell_types2_standard.pdf',width = 20, height = 20)
hmcol = colorRampPalette(c("blue","white","red"))(256)
heatmap.2(as.matrix(spearman_), col=hmcol, dendrogram='none',trace='none', main='Spearmans correlation (rho)', cexRow=0.6)
dev.off()

################################
### EXPRESSION DATA ANALYSIS ###
################################
exp_data = read.table('/Users/woojunshim/Nathan/Data/57epigenomes.RPKM.pc', header=T)
gene_list = exp_data$gene_id
exp_data = exp_data[,2:ncol(exp_data)] 
cell_list = colnames(exp_data)
cell_list = cell_list[2:length(cell_list)]
rownames(exp_data) = gene_list
exp_table = data.frame(matrix(nrow=nrow(exp_data), ncol=ncol(exp_data)))
colnames(exp_table) = colnames(exp_data)

# unit score, a=observed value, b=min, c=sd
us = function(a,b,c){
  return ((a-b) / c)
}

# Convert raw values to unit scores
mins = apply(exp_data, 2, min)
sds = apply(exp_data, 2, sd)
exp_table = us(t(exp_data), mins, sds)
exp_table = t(exp_table)
rownames(exp_table) = gene_list

# dynamics score, a=observed value, b=mean 
ds = function(a,b){
  return (log10(a/b))  
}

# Calcuate expression dynamics across cell types
means = apply(exp_table, 1, mean)
exp_table_ = ds(exp_table, means)

# Export the outcome
write.table(exp_data, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/rpkm.txt', sep='\t', quote=FALSE)
write.table(exp_table, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/expression_unit_score.txt', quote=FALSE)
write.table(exp_table_, file='/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/expression_dynamics_score.txt', quote=FALSE)

#### Mapping analysis
h3k4me3 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/mapping_analysis_top10_tf_h3k4me3_ds.txt')
exp = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/mapping_analysis_top10_tf_exp_ds.txt')

#### H3K4me3 ranked list (with TF categories)
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/intra/h3k4me3_rank_TABLE_final.txt')
data2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/intra/expression_rank_TABLE_final.txt')
heatmap(t(data1), Rowv = F, Colv = F)
data1.t =t(data1)
write.table(data1.t, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/intra/h3k4me3_rank_TABLE_us1.txt', sep='\t', quote=FALSE)

cat1 = vector()
cat2 = vector()
cat3 = vector()
for (i in 1:nrow(data1)){
  for (j in 1:ncol(data1)){
    if (data1[i,j] == 0){
      cat1 = c(cat1, j)
    }
    if (data1[i,j] == 1){
      cat2 = c(cat2, j)
    }
    if (data1[i,j] == 2){
      cat3 = c(cat3, j)
    }
  }
}


#### PLOT TRIPLE BOX PLOTS FOR GENE CATAGORIES 
epi_factor = factor(rownames(data1))
x_range = c(0.5, length(epi_factor)*2+0.5)
y_range = c(0,100)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/intra/boxplots_expression.pdf', width=50, height=50)
temp = data1[1,]
temp = as.numeric(temp)
temp = temp[!is.na(temp)]
total_ = length(temp)
a = (1-which(temp=='0')/total_) * 100
b = (1-which(temp=='1')/total_) * 100
c = (1-which(temp=='2')/total_) * 100
boxplot(c, xlim=x_range, ylim=y_range, boxwex=0.25, at=0.75, col='red',xaxt='n')
boxplot(a, xlim=x_range, ylim=y_range, boxwex=0.25, at=1.05, col='blue',xaxt='n', add=TRUE)
boxplot(b, xlim=x_range, ylim=y_range, boxwex=0.25, at=1.35, col='orange',xaxt='n',add=TRUE)
for (i in 1:length(epi_factor)-1){
  temp = data1[i,]
  temp = as.numeric(temp)
  temp = temp[!is.na(temp)]
  total_ = length(temp)
  a = (1-which(temp=='0')/total_) * 100
  b = (1-which(temp=='1')/total_) * 100
  c = (1-which(temp=='2')/total_) * 100
  boxplot(c, xlim=x_range, ylim=y_range, boxwex=0.25, at=2*i+0.75, col='red',xaxt='n', add=TRUE)
  boxplot(a, xlim=x_range, ylim=y_range, boxwex=0.25, at=2*i+1.05, col='blue',xaxt='n', add=TRUE)
  boxplot(b, xlim=x_range, ylim=y_range, boxwex=0.25, at=2*i+1.35, col='orange',xaxt='n',add=TRUE)
}
axis(1, at=seq(1,112,by=2)+0.05, labels=epi_factor, cex.axis=1.2)
dev.off()

### Plot side-by-side boxplots (conditional probability)
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/cp_expressed_given_TF.txt')
data2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/bg_expressed_given_h3k4me3.txt')
y_range_ = c(0.0,1.0)
boxplot.matrix(as.matrix(data1), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:nrow(data1)-0.3, col='green',ann=FALSE,xlab='Rank percentile(DS)', ylab='Probability', outcol='green')
boxplot.matrix(as.matrix(data2), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:nrow(data2), col='purple',ann=FALSE,xaxt='n', yaxt='n',add=TRUE, outcol='purple')
legend(x=12, y=0.3, cex=0.8, c('P(expressed | TF, above threshold)','P(expressed | above threshold)'), col=c('green','purple'), pch=1)

data3 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/bg_tf_given_h3k4me3.txt')
data4 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/bg_tf_given_exp.txt')
y_range_ = range(data3, data4)
boxplot.matrix(as.matrix(data3), ylim=y_range_, use.cols=FALSE, boxwex=0.25, at=1:nrow(data3)-0.3, col='green',ann=FALSE,xlab='Rank percentile', ylab='Probability', outcol='green')
boxplot.matrix(as.matrix(data4), ylim=y_range_, use.cols=FALSE, boxwex=0.25, at=1:nrow(data4), col='orange',ann=FALSE,xaxt='n', yaxt='n',add=TRUE, outcol='orange')
legend(x=13, y=0.35, cex=0.8, c('P(expressed TF | H3K4me3)','P(expressed TF | Expression)'), col=c('green','orange'), pch=1)

### Plot side-by-side boxplots comparing proportional overlaps between TFs and non-TFs
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/h3k4me3_tf_exp_table.txt')
data2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Expression/h3k4me3_nontf_exp_table.txt')
y_range_ = c(0.0,1.0)
boxplot.matrix(as.matrix(data1), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:nrow(data1)-0.3, col='green',ann=FALSE,xlab='Rank percentile(DS)', ylab='Prop. overlap', outcol='green')
boxplot.matrix(as.matrix(data2), use.cols=FALSE, ylim=y_range_, boxwex=0.25, at=1:nrow(data2), col='purple',ann=FALSE,xaxt='n', yaxt='n',add=TRUE, outcol='purple')
legend(x=15, y=0.3, cex=0.8, c('TF(expressed)','non-TF(expressed)'), col=c('green','purple'), pch=1)

### COMBINED TRIPLE-BOXPLOTS
a=vector()
b=vector()
c=vector()
for (i in 1:nrow(data1)){
  temp = data1[i,]
  temp = as.numeric(temp)
  temp = temp[!is.na(temp)]
  total_ = length(temp)
  a_ = (1-which(temp=='0')/total_) * 100
  b_ = (1-which(temp=='1')/total_) * 100
  c_ = (1-which(temp=='2')/total_) * 100  
  a = c(a,a_)
  b = c(b,b_)
  c = c(c,c_)
}
x_range = c(0,4)
y_range = c(0,100)
boxplot(c, xlim=x_range, ylim=y_range, boxwex=0.25, at=1, col='red',xaxt='n')
boxplot(a, xlim=x_range, ylim=y_range, boxwex=0.25, at=2, col='blue',xaxt='n', add=TRUE)
boxplot(b, xlim=x_range, ylim=y_range, boxwex=0.25, at=3, col='orange',xaxt='n',add=TRUE)

#### Poisson distribution analysis (pilot)
gene = 'HAND2'
values = vector()
temp1 = data_1[data_1$V1==gene,]
epi_list = unique(data_1$V2)
for (epi in epi_list){
  temp2 = data_1[data_1$V2==epi,]
  value = temp2[which(temp2$V1==gene),4]
  total_u = sum(temp2$V4)
  u = value / total_u
  total_v = sum(temp1$V4)
  v = value / total_v
  values = c(values, u*v)
}
converted_values = log10(values/mean(values))

#### Histogram (TFs & non-TFs, pooled dataset) by width/sd
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Threshold/threshold_FET_pvalue_standard.txt')
x_lab_ = seq(10,0,-0.1)
y_range_ = c(0,100)
t1_ = data_$TFs.domains./data_[101,1]
t2_ = data_$non.TFs.domains./data_[101,2]
plot(t1_,x=rev(x_lab), xaxt='n', xlab='Width/SD', main='Cumulative plot (pooled)', ylab='Proportion', col='green', type='l')
axis(1, at=seq(0,10,0.1), labels=x_lab)
lines(t2_, x=rev(x_lab), col='purple')
legend(x=0, y=0.9, col=c('green','purple'), c('TFs','non-TFs'), lwd='1', cex=0.8)

#### Test analysis
data_new = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard_percentile.txt')
temp = data_new[data_new$V2=='E083',]
temp1 = data_new[data_new$V2=='E095',]
temp2 = data_new[data_new$V2=='E073',]
gene = 'HAND2'
con1 = temp[which(temp$V1==gene),]$V4 / sum(temp$V4)
con2 = temp1[which(temp1$V1==gene),]$V4 / sum(temp1$V4)
con3 = temp2[which(temp2$V1==gene),]$V4 / sum(temp2$V4)

#### TSS overlap analysis
data_ = read.table('/Users/woojunshim/Research/Data/TSS_H3K4me3_overlap_dominant_UCSC.txt')
temp = data_[data_$V1=='E083',]
TF = temp[temp$V4=='True',]
nonTF = temp[temp$V4=='False',]
mean(TF$V10)
mean(nonTF$V10)
sd(TF$V10)
sd(nonTF$V10)
wilcox.test(TF$V10, nonTF$V10)
wilcox.test(TF$V9, nonTF$V9)
mean(TF$V9)
mean(nonTF$V9)

#### Distance analysis
data_1 = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/Distance_Analysis.txt")
data_2 = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/ANALYSIS_CVP_K4.bed_assigned_genes.txt")

#### Mapping between GO terms and gene symbols
## Beware that these are non-transitively extracted terms
library(topGO)
library(org.Hs.eg.db)

library(GO.db)  # If need to include parent GO terms
a = as.list(GOBPPARENTS)
a_idx = list(names(a))[[1]]
b = as.list(GOBPCHILDREN)
b_idx = list(names(b))[[1]]

mapping_bp = annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
mapping_mf = annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "symbol")
names_bp = names(mapping_bp)
names_mf = names(mapping_mf)
data_ = read.table("/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/ANALYSIS_ANALYSIS_E095-H3K4me3.gappedPeak.bed_assigned_genes.txt")
benayoun = data_$V2[1:nrow(data_)]
rp = data_$V2[order(data_$V14)]

## 1. Plot GO term saturation curve (transitive)
# Get all the data
go_term1 = 'GO:0007507'  # Heart development
go_term2 = 'GO:0003007'  # Heart morphogenesis
go_term1_idx = which(a_idx==go_term1)
go_terms1 = c(go_term1, list(a[go_term1_idx][[1]])[[1]])
positive_genes1 = vector()
for (term in go_terms1){
  temp = list(mapping_bp[which(names_bp==term)])[[1]][[1]]
  positive_genes1 = c(positive_genes1, temp)
}
positive_genes1 = unique(positive_genes1)


positive_genes1 = list(mapping_bp[which(names_bp==go_term1)])[[1]][[1]]
positive_genes2 = list(mapping_bp[which(names_bp==go_term2)])[[1]][[1]]
combined = union(positive_genes1, positive_genes2)
total_benayoun1 = length(intersect(benayoun, positive_genes1))
total_rp1 = length(intersect(rp, positive_genes1))
total_benayoun2 = length(intersect(benayoun, positive_genes2))
total_rp2 = length(intersect(rp, positive_genes2))
total_combined = length(intersect(benayoun, combined))
vector_benayoun1 = vector()
vector_benayoun2 = vector()
vector_benayoun = vector()
vector_rp1 = vector()
vector_rp2 = vector()
vector_rp = vector()
for (i in 1:length(rp)){
  vector_benayoun1 = c(vector_benayoun1, length(intersect(benayoun[1:i],positive_genes1))/total_benayoun1)
  vector_benayoun2 = c(vector_benayoun2, length(intersect(benayoun[1:i],positive_genes2))/total_benayoun2)
  vector_benayoun = c(vector_benayoun, length(intersect(benayoun[1:i],combined))/total_combined)
  vector_rp1 = c(vector_rp1, length(intersect(rp[1:i],positive_genes1))/total_rp1)
  vector_rp2 = c(vector_rp2, length(intersect(rp[1:i],positive_genes2))/total_rp2)
  vector_rp = c(vector_rp, length(intersect(rp[1:i],combined))/total_combined)
}
# Plot cumulative plots
# 1. genes associated with GO:0007507(heart development) 
plot(x=seq(1,length(rp),1), y=predict(loess(vector_rp1~seq(1,length(rp),1))), col='green', type='l')
lines(x=seq(1,length(benayoun),1), y=predict(loess(vector_benayoun1~seq(1,length(benayoun),1))), col='purple', type='l')

# Only for selected cardiac developmental genes (within the top list)
# For GO:0007507
cardiac_genes = c('WNT5A', 'TBX20', 'SOX4', 'TPM1', 'PTEN', 'ZIC3', 'TGFB2', 'ZFP36L1', 'HEY1', 'HAND1', 'HAND2', 'GATA4', 'TMEM100', 'TWIST1', 'BMP4', 'BMP2', 'TBX3', 'SMAD7', 'NKX2-6', 'SMAD6', 'FZD2', 'ISL1', 'TENM4', 'ADM', 'ID2', 'GSK3B', 'ZMIZ1', 'SALL1', 'VEGFA', 'PDGFRA', 'PRDM1', 'KDM6B','NKX2-5')
threshold = 100  # top 100 genes

cardiac_genes_rp = c('MEF2A', 'KDM6A', 'NRP1', 'EFNA1', 'LMO4', 'TBX20', 'GJA1', 'GLI3', 'PTEN', 'ZIC3', 'TGFB1', 'CTNNB1', 'TGFB2', 'GATA2', 'MYOCD', 'HEY1', 'ROBO1', 'GATA4', 'MKKS', 'SHC1', 'YAP1', 'TMEM100', 'TWIST1', 'NKX2-6', 'EFNB2', 'PROX1', 'HES1', 'HHEX', 'ACVR2B', 'CRKL', 'ADM', 'JUN', 'ZMIZ1', 'SIX1', 'VEGFA', 'PDGFRA', 'DSP', 'PRDM1', 'KDM6B', 'WNT5A', 'FGFR2', 'FGFR1', 'BMPR2', 'SOX4', 'ZBTB14', 'TPM1', 'SRF', 'EPHB4', 'MSX2', 'ZFP36L1', 'APLNR', 'CHD7', 'HAND1', 'HAND2', 'HEXIM1', 'FAT4', 'CACYBP', 'CAMK2D', 'PLXND1', 'AXIN2', 'BCOR', 'NKX2-5', 'IFT140', 'BMP4', 'BMP2', 'TBX3', 'SMAD7', 'RBM20', 'SMAD6', 'FZD1', 'SMAD3', 'FZD2', 'SNAI2', 'ISL1', 'SNAI1', 'TAB2', 'FOXP1', 'SOD2', 'TENM4', 'SALL4', 'ID2', 'C3ORF58', 'GSK3A', 'ID1', 'GSK3B', 'SALL1', 'LRP6', 'RBPJ', 'LRP2', 'BMP5', 'BMPR1A')
cardiac_genes_benayoun = c('MEF2A', 'KDM6A', 'NRP1', 'EFNA1', 'LMO4', 'TBX20', 'GJA1', 'REST', 'PTEN', 'ZIC3', 'TGFB1', 'CTNNB1', 'TGFB2', 'GATA2', 'MYOCD', 'HEY1', 'ROBO1', 'GATA4', 'MKKS', 'SHC1', 'YAP1', 'TMEM100', 'TWIST1', 'CDK1', 'NKX2-6', 'EFNB2', 'PROX1', 'HES1', 'HHEX', 'ADAMTS6', 'ACVR2B', 'CRKL', 'ADM', 'JUN', 'ZMIZ1', 'SIX1', 'VEGFA', 'PDGFRA', 'DSP', 'PRDM1', 'KDM6B', 'WNT5A', 'FGFR2', 'FGFR1', 'BMPR2', 'SOX4', 'TPM1', 'SRF', 'EPHB4', 'ZFP36L1', 'MSX2', 'APLNR', 'CHD7', 'HAND1', 'HAND2', 'FAT4', 'HEXIM1', 'CAMK2D', 'AXIN2', 'BCOR', 'PCSK5', 'NKX2-5', 'IFT140', 'BMP4', 'BMP2', 'TBX3', 'SMAD7', 'SMAD6', 'FZD1', 'SMAD3', 'TEAD2', 'FZD2', 'SNAI2', 'ISL1', 'PPP1R13L', 'FOXP1', 'SOD2', 'TENM4', 'SALL4', 'ID2', 'C3ORF58', 'GSK3A', 'ID1', 'GSK3B', 'SALL1', 'PRKAR1A', 'RBPJ', 'LRP2', 'BMP5', 'BMPR1A')
cardiac_genes = union(cardiac_genes_rp, cardiac_genes_benayoun)
threshold = 1000

key_regulators = c('TBX3','TBX20','GATA4','WNT2','TBX2','NKX2-5','HAND2','TBX5','ISL1','HAND1','FOXC1','FOXC2','WNT11','WNT5A','BMP2','SMARCD3','MEF2C','PDGFRA','GATA6','MESP1','MESP2')
key_regulators = c(key_regulators, 'MEIS2','MEIS1','SALL1','TWIST1','JARID2')
cardiac_genes = key_regulators
threshold = length(rp)

total_benayoun = length(intersect(benayoun[1:threshold], cardiac_genes))
total_rp = length(intersect(rp[1:threshold], cardiac_genes))
#total = 511 # Total number of genes in the dataset with this term
random_rate = 511/16792  # Background (i.e. expected ratio of genes with this term )
random_rate = 28/9629
random_rate = length(intersect(rp[1:length(rp)], cardiac_genes)) / length(rp)
vector_benayoun = vector()
vector_rp = vector()
vector_random = vector()
vector_expr = vector()
for (i in 1:threshold){
  vector_benayoun = c(vector_benayoun, length(intersect(benayoun[1:i],cardiac_genes)))
  vector_rp = c(vector_rp, length(intersect(rp[1:i],cardiac_genes)))
  vector_expr = c(vector_expr, length(intersect(symbol_[1:i],cardiac_genes)))
  vector_random = c(vector_random, random_rate*i)
}
plot(x=seq(1,threshold,1), y=predict(loess(vector_rp~seq(1,threshold,1))), col='green', type='l', xlab='Rank position', ylab='Number of genes', main='Heart Development GO:0007507')
lines(x=seq(1,threshold,1), y=predict(loess(vector_benayoun~seq(1,threshold,1))), col='purple', type='l')
#lines(x=seq(1,threshold,1), y=vector_random, col='black', type='l')
lines(x=c(0,threshold), y=c(1,vector_random[threshold]), col='black', type='l')
legend(x=0, y=30, c('Width alone','Width+DS+Distance','Random'), col=c('purple','green','black'), lty=1, cex=0.8)
legend(x=7800, y=4, c('Width alone','Width+DS+Distance','Random'), col=c('purple','green','black'), lty=1, cex=0.7)
legend(x=75, y=3, c('Width alone','Width+DS+Distance'), col=c('purple','green'), lty=1, cex=0.7)
legend(x=10000, y=7, c('H3K4me3 width','RPKM','Random'), col=c('purple','gold','black'), lty=c(1,1,3), cex=0.7)


plot(x=seq(1,threshold,1), y=vector_rp, col='purple', type='l', xlab='Rank position', ylab='Number of genes', main=c('28 Cardiac lineage regulators', 'Adult left ventricle (E095)'))
lines(x=seq(1,threshold,1), y=vector_benayoun, col='purple', type='l')
lines(x=seq(1, threshold,1), y=vector_expr, col='gold', type='l')
lines(x=seq(1,threshold,1), y=vector_random, col='black', type='l', lty=3)
lines(x=c(1,threshold), y=c(0,vector_random[threshold]), col='black', type='l')

expr = read.table('/Users/woojunshim/Nathan/Data/57epigenomes.RPKM.pc',header=T)

library(biomaRt)
listMarts(host="grch37.ensembl.org")  # To get a list of available biomart databases
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
attr<-listAttributes(ensembl)
filters<-listFilters(ensembl)
ensg2symbol<-getBM(attributes=c("ensembl_gene_id","ensembl_peptide_id","entrezgene","hgnc_symbol","gene_biotype","wikigene_description"),value=T,mart=ensembl)
idx = which(ensg2symbol$ensembl_gene_id==expr$gene_id)

expr = expr[order(expr$E095, decreasing=T),]
expr = expr[order(expr$E065, decreasing=T),]
expr_ = expr$gene_id
symbol_ = vector()

for (i in 1:nrow(expr)){
  symbol_ = c(symbol_, ensg2symbol$hgnc_symbol[which(ensg2symbol$ensembl_gene_id==expr_[i])])
}

#### Extract the expression rank position 
temp = data_new[data_new$V2=='E083',]
genes = as.character(temp$V1)
results = vector()
for (i in 1:length(genes)){
  results = c(results, which(symbol_==genes[i]))
}

## Plot a heatmap for enrichment
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/enrichment_counts_normalised_FPKM.txt')
colnames(file_) = seq(95,0,-5)
file_$TF = c('Eomes','T','Pou5f1','CTCF')
plot(x=seq(1,ncol(file_)), y=file_[2,], xlab='bin(size=50)', ylab='Prop. positive genes', main='ChIP-seq enrichment (NR5A2)')

## Heatmap using ggplot2
library(ggplot2)
library(reshape)
table = melt(file_, id='TF')
p = ggplot(table, aes(x=variable, y=TF, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=(max(table$value)+min(table$value))/2, mid='grey70') + ggtitle('Normalised TF ChIP-seq counts per 1 kbp') + labs(x='Percentile Rank', y='', subtitle='Ranked by the FPKM value')
p + geom_tile()

## (continued) Plot line plots (for percentile rank)
file_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/enrichment_counts_normalised_H3K4me3.txt')
colnames(file_1) = seq(95,0,-5)
file_1$TF = c('Eomes','T','Pou5f1','CTCF')
file_1$method = 'H3K4me3'
file_2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/enrichment_counts_normalised_FPKM.txt')
colnames(file_2) = seq(95,0,-5)
file_2$TF = c('Eomes','T','Pou5f1','CTCF')
file_2$method = 'FPKM'

file_ = rbind(file_1, file_2)
file_ = file_[file_$TF=='CTCF',]
table = melt(file_, id=c('TF','method'))
ggplot(table, aes(x=variable, y=value, colour=method, group=method)) + geom_point() + geom_line() + ggtitle('CTCF') + labs(x='Percentile Rank', y='Normalised ChIP-seq counts per 1 kbp')

## (continued) Plot line plots with regression lines (for 100-regions bin)
file_1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/enrichment_counts_normalised_H3K4me3.txt')
colnames(file_1) = seq(1,ncol(file_1),1)
file_1$TF = c('Eomes','T','Pou5f1','CTCF')
file_1$method = 'H3K4me3'
file_2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/enrichment_counts_normalised_FPKM.txt')
colnames(file_2) = seq(1,ncol(file_2),1)
file_2$TF = c('Eomes','T','Pou5f1','CTCF')
file_2$method = 'FPKM'

file_ = rbind(file_1, file_2)
file_ = file_[file_$TF=='Pou5f1',]
table = melt(file_, id=c('TF','method'))
p = ggplot(table, aes(x=variable, y=value, colour=method, group=method))
x_scale = c(1,seq(10,120,10))
p + geom_point() + stat_smooth(method=lm) + scale_x_discrete(breaks=x_scale) + ggtitle('Pou5f1') + labs(x='Bin number(size=100)', y='Normalised ChIP-seq counts per 1 kbp')

### Bar plots (GO terms)using ggplot2
# All genes with a TF ChIP-seq binding
tf = 'Pou5f1'
file0 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/functional_analysis/short/',tf,'_GO_all_genes_short.txt',sep='')
file_ = read.table(file0, header=T)
file_$no = seq(1, nrow(file_), 1)
subtable_ = file_[1:10, c(2,4,5)]  # Top 10 terms
df_ = subtable_[order(subtable_$FDR, decreasing = T),]
df_$name = factor(df_$Term, levels=df_$Term)
title_ = paste(tf,'(All Genes)')
ggplot(df_, aes(x=factor(name), y=-log10(FDR))) + geom_bar(stat='identity', fill='light blue', colour='black') + labs(x='', title=title_) + theme(axis.text.y = element_text(face='bold', color='#993333', size=7), axis.text.x = element_text(face='bold', color='black', size=12)) + coord_flip() 
save0 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/plots/','GO_enrichment_',tf,'_all_genes.png',sep='')
ggsave(save0,width=10,height=5)

# H3K4me3 top 100 genes with the ChIP-seq bindings
file1 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/functional_analysis/short/',tf,'_GO_top100_short_H3K4me3.txt',sep='')
file_1 = read.table(file1, header=T)
file_1$no = seq(1, nrow(file_1), 1)
subtable_1 = file_1[1:10, c(2,4,5)]  # Top 10 terms
df_1 = subtable_1[order(subtable_1$FDR, decreasing = T),]
df_1$name = factor(df_1$Term, levels=df_1$Term)
title_ = paste(tf,'(H3K4me3)')
ggplot(df_1, aes(x=factor(name), y=-log10(FDR))) + geom_bar(stat='identity', fill='light blue', colour='black') + labs(x='', title=title_) + theme(axis.text.y = element_text(face='bold', color='#993333', size=10), axis.text.x = element_text(face='bold', color='black', size=12)) + coord_flip() 
save1 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/plots/','GO_enrichment_',tf,'_top100_H3K4me3.png',sep='')
ggsave(save1,width=10,height=5)

# FPKM top 100 genes with the ChIP-seq bindings
file2 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/functional_analysis/short/',tf,'_GO_top100_short_FPKM.txt',sep='')
file_2 = read.table(file2, header=T)
file_2$no = seq(1, nrow(file_2), 1)
subtable_2 = file_2[1:10, c(2,4,5)]  # Top 10 terms
df_2 = subtable_2[order(subtable_2$FDR, decreasing = T),]
df_2$name = factor(df_2$Term, levels=df_2$Term)
title_ = paste(tf,'(FPKM)')
ggplot(df_2, aes(x=factor(name), y=-log10(FDR))) + geom_bar(stat='identity', fill='light blue', colour='black') + labs(x='', title=title_) + theme(axis.text.y = element_text(face='bold', color='#993333', size=10), axis.text.x = element_text(face='bold', color='black', size=12)) + coord_flip() 
save2 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/plots/','GO_enrichment_',tf,'_top100_FPKM.png',sep='')
ggsave(save2,width=10,height=5)

## Side-by-side bar plots (for the 3 groups) for selected GO terms
## Fold enrichment used (as numbers of genes are not the same)
library(reshape2)
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/High_dynamics/top100/'
name = 'GSE96290_fetal_cardiac_muscle'
f1 = paste(pathway, 'GO_',name,'_top100_short.txt', sep='')
f2 = paste(pathway, 'GO_',name,'_filtered_short.txt', sep='')
f3 = paste(pathway, 'GO_',name,'_filtered_random_short.txt', sep='')
file_ = read.table(f1, header=T)
file_1 = read.table(f2, header=T)
file_2 = read.table(f3, header=T)
tf = 'CTCF'
terms_ = c('heart_morphogenesis(GO:0003007)','cardiac_muscle_tissue_development(GO:0048738)','cardiocyte_differentiation(GO:0035051)','cardiac_chamber_morphogenesis(GO:0003206)','cardiac_septum_morphogenesis(GO:0060411)')
file_$group = rep('top 100',nrow(file_))
file_1$group = rep('filtered',nrow(file_1))
file_2$group = rep('reduced',nrow(file_2))
combined = rbind(file_, file_1, file_2)
combined_ = combined[combined$Term %in% terms_,]
combined_melt = combined_[1:15,c(2,3,5)]
df.long = melt(combined_melt)
df_2 = combined_melt[order(combined_melt$Fold, decreasing = T),]
df_2$name = factor(df_2$Term, levels=df_2$Term)
name_ = paste(name, sep=' ')
ggplot(df.long, aes(Term, value, fill=group)) + geom_bar(stat='identity',position='dodge') + labs(y='Fold enrichment', x='', title=name_) + coord_flip() 
save_ = paste(pathway,'GO_enrichment_comparison_',name,'.png',sep='')
ggsave(save_,width=10,height=5)

## Heatmap for selected cardiac regulators + structural genes (Paige et al.)
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard_percentile.txt')
structural = c('MYL4','CKM','TNNI1','ATP2A2','PLN','MYL3','MYL2')
regulatory = c('HAND2','TBX5','ISL1','GATA4','GATA6','NKX2-5')
s_table = data_[data_$V1 %in% structural, c(1,2,11)]
r_table = data_[data_$V1 %in% regulatory, c(1,2,11)]
s_table$group = 'S'
r_table$group = 'R'
table_ = melt(rbind(s_table, r_table))
ggplot(table_, aes(x=V1, y=V2, fill=value)) + geom_tile() + scale_fill_gradient2(midpoint=50, mid='grey70')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/plots/Cardiac_genes_heatmap.png', width=10, height=10)

## Converting Ensembl IDs to gene symbols
data1 = read.table('/Users/woojunshim/Research/Data/TF/TF_animalTFDB_ENSG_short.txt', stringsAsFactors = F)
data2 = read.table('/Users/woojunshim/Research/Data/TF/hs.tf.ass', stringsAsFactors = F)
symbol1 = vector()
symbol1_ = vector()
symbol2 = vector()
symbol2_ = vector()
for (i in 1:nrow(data1)){
  symbol1 = c(symbol1, ensg2symbol$hgnc_symbol[which(ensg2symbol$ensembl_gene_id==data1[i,1])])
}
for (i in 1:nrow(data2)){
  symbol2 = c(symbol2, ensg2symbol$hgnc_symbol[which(ensg2symbol$ensembl_peptide_id==data2[i,2])])
}
write(symbol1, '/Users/woojunshim/Research/Data/TF/TF_animalTFDB_short.txt', sep='\n')
write(symbol2, '/Users/woojunshim/Research/Data/TF/TF_DBD_short.txt', sep='\n')

## Extract all ENSG & gene symbol IDs 
table_ = data.frame(Ensembl=ensg2symbol$ensembl_gene_id, gene_symbol=ensg2symbol$hgnc_symbol)
write.table(table_, '/Users/woojunshim/Research/Data/Ensembl_gene_symbols_conversion.txt', quote=F, sep='\t')

## TF characterisation
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test/ANALYSIS_GSE96290_ENCFF682XDY_peaks_hg19_fetal_cardiac_muscle_invitro.bed_assigned_genes.txt')
tf = data_[data_$V8=='True',]
non_tf = data_[data_$V8=='False',]
wilcox.test(tf$V10, non_tf$V10)
mean(tf$V10)
mean(non_tf$V10)
plot(density(non_tf$V10), col='green')
lines(density(tf$V10), col='red')
temp = data_[data_$V10>3000, ]

## Dynamics analysis
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt')
name = 'ATP2A2'
cut_off = 4.20801191985
temp = data_[data_$V1==name,]
temp_ = temp[order(temp$V5, decreasing = T),]
plot(x=seq(1,nrow(temp_)), y=temp_$V5, main=name, xlab='Rank',ylab='Width/SD')
abline(h=cut_off, col='red',lty=2)
t1 = abs(sum(temp_$V5-2.3145439)) / nrow(temp_)
t2 = abs(sum(temp_$V5-2.3145439)) / nrow(temp_)

## Dynamics heatmap (for those selected by the high dynamics range)
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/tf_dynamics_range_table_GROUP.txt')
annotation = read.table('/Users/woojunshim/Research/Data/Roadmap_IDs_.txt')
for (i in 1:ncol(data_)){
  idx = which(annotation$V1==colnames(data_)[i])
  colnames(data_)[i] = as.character(annotation$V2[idx])
}
hc.rows = hclust(dist(data_))
hc.cols = hclust(dist(t(data_)))
temp = data_[rownames(data_) %in% cardiac_genes,]
#pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Dynamics/plots/DS/cardiac_structural_heatmap(high_dynamics)_.pdf', width=50, height=50)
heatmap(as.matrix(temp), scale='none', cexCol=0.4)
heatmap(as.matrix(temp), scale='none', col=colorRampPalette(c('white','red'))(256), cexCol=0.7)
#dev.off()
cardiac_genes = c('MYLK3','MYH6','MYL4','PLN','MYH7','ATP2A2','MYL2')
cardiac_genes = c('TBX3','TBX20','GATA6','GATA4','HAND2','HAND1','MEIS2','MEIS1','TBX5','NKX2-5','ISL1')

### Tissue type enrichment
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/tf_tissue_enrichment_table_pvalues.txt')
data_ = -log10(data_)  
data_ = ifelse(data_ <= 0.05, 1, 0)
library(gplots)
heatmap.2(as.matrix(data_), key=FALSE, main='TFs(-log10(p-value))', labRow='',trace='none', scale='none',col=colorRampPalette(c('white','red'))(256), cexCol=0.6)

### Expression data (57epigenomes.RPKM)
data_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.pc')

### Cut-off_list analysis
data1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/tf_cut_off_list_1.0_.txt')
data2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/non_tf_cut_off_list_1.0_.txt')
data1 = data1[which(data1$V2<9999),]
data2 = data2[which(data2$V2<9999),]
plot(density(data2$V2), col='green', xlab='DS(Width/SD)', main='Thresholds (DS) for the non-linear section')
lines(density(data1$V2), col='red')
legend('topright', col=c('red','green'), c('TFs','Non-TFs'), lty=1)

### JUST READING IN A DATA
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard_1.5_combined.txt')
true_table = data_[data_$V12!=0,]
false_table = data_[data_$V12==0,]
a=length(which(true_table$V11=='True'))
b=nrow(true_table)-a
c=length(which(false_table$V11=='True'))
d=nrow(false_table)-c
mat_ = matrix(c(a,c,b,d), nrow=2)
fisher.test(mat_)  

data__ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Test1/summary_results_pdf_dominant_2.5kb_standard_1.5_combined_.txt')
data_ = data__[data__$V2=='E083',]
data_ = data_[order(data_$V8, decreasing = T),]
data_ = data__
plot(x=seq(1,nrow(data_)), y=data_$V8)

### Box plots comparing filtering specificity between members & non-members
top_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/filtering_specificity_top100.txt')
library(ggplot2)
library(reshape2)
colnames(top_) = c('Non-members','Members','Non-members/Members')
top_ = melt(top_[,1:2])
colnames(top_) = c('Membership','Filtering_specificity')
ggsave('/Users/woojunshim/Research/Data/highly_expressed_genes/Top_100_TFs.png',width=10,height=5)

combined_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/filtering_specificity_broadpeak.txt')
combined = melt(combined_)
ggplot(combined, aes(x=variable,y=value)) + geom_boxplot() + ggtitle('Specificity ratios (non-members/members)') + ylab('Filtering specificity') + coord_cartesian(ylim = c(0, 3)) 
ggsave('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/Specificity_ratios.png',width=5,height=5)

### Fisher's exact test line plot
library(ggplot2)
library(reshape2)

group.colors = c('red','blue')
for (epi in rownames(top_)){
  filename = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/',epi,'_performance_fet_top150_broadpeak_percentile.txt',sep='')
  data_ = read.table(filename)
  data_ = -log10(data_)
  data_ = t(data_)
  colnames(data_) = c('Filtering','Simple extraction')
  rownames(data_) = seq(nrow(data_)-1,0, -1)
  data_ = melt(data_)
  colnames(data_) = c('rank','method','value')
  ggplot(data=data_, aes(x=rank, y=value, group=method, colour=method)) + geom_line() + ggtitle('Left ventricle(E095)')+ ylab('-log10(p-value)') +xlab('Percentile Rank') + scale_colour_manual(values=c("red","dark green","dark blue")) + scale_x_reverse()
  filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/',epi,'_performance_fet_top150_broadpeak_percentile.png',sep='')
  ggsave(filename_,width=10,height=5,dpi=300)
}


#### ENRICHMENT PLOT (FISHER'S EXACT TEST)
epis = c('H3K4me3','H3K27me3','H3K4me1','H3K36me3','H3K9me3')
#epis = c('H3K9me3')
for (epi in epis){
  filename = paste('/Users/woojunshim/Nathan/Scripts/E095-',epi,'.broadPeakbenayoun_sliding_fet_E095-',epi,'.broadPeak.txt', sep='')
  data_ = read.table(filename)
  data_ = -log10(data_)
  data_[1,] = rev(data_[1,])
  data_[2,] = rev(data_[2,])
  data_ = t(data_)
  colnames(data_) = c('GO:0007507 (Heart development)','GO:0003677 (DNA binding)')
  rownames(data_) = seq(nrow(data_)-1,0, -1)
  data_ = melt(data_)
  data_$pp = rev(data_$Var1)+1
  colnames(data_) = c('rank1','GO_term','value','rank2')
  ggplot(data=data_, aes(x=rank2, y=value, group=GO_term, colour=GO_term)) + geom_line() + ggtitle(epi)+ ylab('-log10(p-value)') +xlab('Rank position') + scale_colour_manual(values=c("red","dark green","dark blue")) + coord_cartesian(ylim = c(0, 50)) 
  filename_ = paste('/Users/woojunshim/Nathan/Results/57_epigenomes/E095_',epi,'_performance_fet.png',sep='')
  ggsave(filename_,width=10,height=5,dpi=300)
}

filename = paste('/Users/woojunshim/Nathan/Scripts/E095-',epi,'.broadPeakbenayoun_sliding_fet_E095-',epi,'.broadPeak.txt', sep='')
data_ = read.table(filename)
data_ = -log10(data_)
data_[1,] = rev(data_[1,])
data_[2,] = rev(data_[2,])
data_ = t(data_)
colnames(data_) = c('GO:0007507 (Heart development)','GO:0003677 (DNA binding)')
rownames(data_) = seq(nrow(data_)-1,0, -1)
data_ = melt(data_)
data_$pp = rev(data_$Var1)+1
colnames(data_) = c('rank1','GO_term','value','rank2')
ggplot(data=data_, aes(x=rank2, y=value, group=GO_term, colour=GO_term)) + geom_line() + ggtitle('Left ventricle (E095)')+ ylab('-log10(p-value)') +xlab('Rank position') + scale_colour_manual(values=c("red","dark green","dark blue")) 
filename_ = paste('/Users/woojunshim/Nathan/Results/57_epigenomes/E095_',epi,'_performance_fet.png',sep='')
ggsave(filename_,width=10,height=5,dpi=300)


epi='E095'
filename = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/',epi,'_sliding_fet_top150.txt',sep='')
filename = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/Assigned_genes_dMS__.txt_fet.txt'
data_ = read.table(filename)
data_ = -log10(data_)
data_ = t(data_)
colnames(data_) = c('Filtering','Bottom-up')
rownames(data_) = seq(1, nrow(data_))
data_ = melt(data_)
colnames(data_) = c('rank','group','value')
ggplot(data=data_, aes(x=rank, y=value, group=group, colour=group)) + geom_line() + ylab('-log10(p-value)')
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/plots/dMS.png',sep='')
ggsave(filename_,width=10,height=5)

### Density heat map
file_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/density_ratio_expression_top150_.txt')
file__ = file_[,2:ncol(file_)]
rownames(file__) = file_[,1]
colnames(file__) = seq(95,0,by=-5)
file__ = t(file__)
library(ggplot2)
library(reshape)
table = melt(file__)
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=1, mid='grey70', limits=c(0,8.0)) + ggtitle('Odd ratios of expressed TFs (among all genes with RPKM>1.0)') + labs(x='Percentile Rank (RPKM)', y='Cell type') + scale_x_reverse(breaks = seq(0,95,by=5)) 
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/plots/odd_ratios_top150_TFs_expression_.png',sep='')
ggsave(filename_,width=5,height=5)

# Boxplots 
file_1 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/new/density_ratio_genes_h3k4me3_a.txt')
file_2 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/new/density_ratio_genes_rpkm_a.txt')
file__1 = file_1[,2:ncol(file_1)]
file__2 = file_2[,2:ncol(file_2)]
rownames(file__1) = file_1[,1]
colnames(file__1) = seq(95,0,by=-5)
rownames(file__2) = file_2[,1]
colnames(file__2) = seq(95,0,by=-5)
file__1 = t(file__1)
file__2 = t(file__2)
table1 = melt(file__1)
table2 = melt(file__2)
table1$method = rep('H3K4me3 width', nrow(table1))
table2$method = rep('RPKM (Expression)', nrow(table2))
table = rbind(table1, table2)
table$lev = factor(table$Var1, levels=seq(95,0,by=-5))
p = ggplot(table, aes(lev, value, fill=method))
p + geom_boxplot() + ggtitle('Differentially expressed genes') + xlab('Percentile Rank') + ylab('Enrichment score')
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/new/enrichment_DE_genes_a.png',sep='')
ggsave(filename_,width=10,height=5)


### best p-values boxplot
file_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/best_p-values.txt')
colnames(file_) = c('cell_type','filtering','bottom_up')
table = melt(file_)
table$value = -log10(table$value)
p = ggplot(table, aes(variable, value)) 
p + geom_boxplot() + ggtitle('Enrichment of cell type specific top 150 TFs') + xlab('Percentile Rank') + ylab('Odd ratio')
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/plots/best_p-values.png',sep='')
ggsave(filename_,width=10,height=5)

# Side-by-side for top30, 50, 150 and expressed
file1 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top30/best_pvalues_top30_broadpeak.txt')
file2 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/best_pvalues_top150_broadpeak.txt')
file3 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/expressed/best_pvalues_expressed_broadpeak.txt')
file4 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/best_p-values_expressed.txt')
file1$group = rep('top30',nrow(file1))
file2$group = rep('top150',nrow(file2))
file3$group = rep('expressed',nrow(file3))
table_ = rbind(file1,file2,file3)
colnames(table_) = c('cell type','filtering','simple extraction','group')
table = melt(table_)
table$value = -log10(table$value)
table$lev = factor(table$group, levels=c('top30','top150','expressed'))
colnames(table) = c('cell type','group','method','value','lev')
p = ggplot(table, aes(lev, value, fill=method)) 
p + geom_boxplot() + ggtitle('Most significant p-values (FET) ') +  ylab('-log10(p-value)') + xlab('')
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/comparison_best_p-values.png',sep='')
ggsave(filename_,width=5,height=5)

### PLOT CORRELATION PLOTS
library(corrplot)
data_ = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/Rank_table_H3K4me3_.txt')
cor_ = round(cor(data_, use='complete.obs'), digits=2)
corrplot(cor_, shade.col=NA, tl.col='black', tl.srt=45)

### BOXPLOTS FOR CORE 5 HISTONE MARKS (FROM THE LAST YEAR'S ANALYSIS)
file1 = read.table('/Users/woojunshim/Nathan/Results/57_epigenomes/H3K4me1_best_positions_percentile_.txt')
file2 = read.table('/Users/woojunshim/Nathan/Results/57_epigenomes/H3K4me3_best_positions_percentile.txt')
file3 = read.table('/Users/woojunshim/Nathan/Results/57_epigenomes/H3K9me3_best_positions_percentile.txt')
file4 = read.table('/Users/woojunshim/Nathan/Results/57_epigenomes/H3K27me3_best_positions_percentile.txt')
file5 = read.table('/Users/woojunshim/Nathan/Results/57_epigenomes/H3K36me3_best_positions_percentile.txt')
temp1_1 = file1[,1:2]
temp2_1 = file2[,1:2]
temp3_1 = file3[,1:2]
temp4_1 = file4[,1:2]
temp5_1 = file5[,1:2]
temp1_2 = file1[,c(1,3)]
temp2_2 = file2[,c(1,3)]
temp3_2 = file3[,c(1,3)]
temp4_2 = file4[,c(1,3)]
temp5_2 = file5[,c(1,3)]
temp1_1$method = rep('Broadest peak',nrow(temp1_1))
temp2_1$method = rep('Broadest peak',nrow(temp2_1))
temp3_1$method = rep('Broadest peak',nrow(temp3_1))
temp4_1$method = rep('Broadest peak',nrow(temp4_1))
temp5_1$method = rep('Broadest peak',nrow(temp5_1))
temp1_2$method = rep('Sum of peaks', nrow(temp1_2))
temp2_2$method = rep('Sum of peaks', nrow(temp2_2))
temp3_2$method = rep('Sum of peaks', nrow(temp3_2))
temp4_2$method = rep('Sum of peaks', nrow(temp4_2))
temp5_2$method = rep('Sum of peaks', nrow(temp5_2))
temp1_1$mark = rep('H3K4me1',nrow(temp1_1))
temp2_1$mark = rep('H3K4me3',nrow(temp2_1))
temp3_1$mark = rep('H3K9me3',nrow(temp3_1))
temp4_1$mark = rep('H3K27me3',nrow(temp4_1))
temp5_1$mark = rep('H3K36me3',nrow(temp5_1))
temp1_2$mark = rep('H3K4me1',nrow(temp1_2))
temp2_2$mark = rep('H3K4me3',nrow(temp2_2))
temp3_2$mark = rep('H3K9me3',nrow(temp3_2))
temp4_2$mark = rep('H3K27me3',nrow(temp4_2))
temp5_2$mark = rep('H3K36me3',nrow(temp5_2))
colnames(temp1_1) = c('cell_type','value','method','mark')
colnames(temp2_1) = c('cell_type','value','method','mark')
colnames(temp3_1) = c('cell_type','value','method','mark')
colnames(temp4_1) = c('cell_type','value','method','mark')
colnames(temp5_1) = c('cell_type','value','method','mark')
colnames(temp1_2) = c('cell_type','value','method','mark')
colnames(temp2_2) = c('cell_type','value','method','mark')
colnames(temp3_2) = c('cell_type','value','method','mark')
colnames(temp4_2) = c('cell_type','value','method','mark')
colnames(temp5_2) = c('cell_type','value','method','mark')
table_ = rbind(temp1_1,temp2_1,temp3_1,temp4_1,temp5_1,temp1_2,temp2_2,temp3_2,temp4_2,temp5_2)
table = melt(table_)
p = ggplot(table_, aes(mark, -log10(value), colour=method)) 
p + geom_boxplot() + ggtitle('Best p-values') + xlab('Histone mark') + ylab('-log10(p-value)')
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/plots/best_pvalues_5marks.png',sep='')
ggsave(filename_,width=10,height=5)

### Compare TFs vs nonTFs 
# 1. number of H3K4me3 peak counts
tf_ = read.csv('/Users/woojunshim/Research/Data/TF/TF_combined.txt')
mrna_ = read.csv('/Users/woojunshim/Research/Data/mRNA_genes.txt')
data__ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_pdf_dominant_2.5kb_standard.txt_1.0_.txt')
deg_list = read.csv('/Users/woojunshim/Research/Data/highly_expressed_genes/DEG_list_top30.txt')
mrna_idx = which(rownames(table_) %in% mrna_[,1])
table_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/Peak_count_table.txt')
table_ = table_[mrna_idx,]
tf_idx = which(rownames(table_) %in% tf_[,1])
nontf_idx = which(!(rownames(table_) %in% tf_[,1]))
tf = table_[tf_idx,]
nontf = table_[nontf_idx,]
tf_list = apply(as.matrix(tf), 1, FUN=mean)
nontf_list = apply(as.matrix(nontf), 1, FUN=mean)

# 2. Compare 
### Statistics for assigned peaks
data_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/ranksum_peak_count_mrna.txt')

## Side-by-side bar plots (for the 3 groups) for selected GO terms
## Fold enrichment used (as numbers of genes are not the same)
library(reshape2)
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/GO/'
name = 'E083_mrna_top5%'
f1 = paste(pathway, 'GO_',name,'_dominant_short.txt', sep='')
f2 = paste(pathway, 'GO_',name,'_average_short.txt', sep='')
f3 = paste(pathway, 'GO_',name,'_sum_short.txt', sep='')
file_ = read.table(f1, header=T)
file_1 = read.table(f2, header=T)
file_2 = read.table(f3, header=T)
tf = 'Fetal heart (E083)'
terms_ = c('heart_development(GO:0007507)','cardiac_muscle_tissue_development(GO:0048738)','circulatory_system_development(GO:0072359)','cardiovascular_system_development(GO:0072358)','heart_morphogenesis(GO:0003007)')
file_$group = rep('Broadest',nrow(file_))
file_1$group = rep('Mean',nrow(file_1))
file_2$group = rep('Sum',nrow(file_2))
combined = rbind(file_, file_1, file_2)
combined_ = combined[combined$Term %in% terms_,]
combined_melt = combined_[1:15,c(2,4,5)]
df.long = melt(combined_melt)
df_2 = combined_melt[order(combined_melt$FDR, decreasing = T),]
df_2$name = factor(df_2$Term, levels=df_2$Term)
name_ = paste(name, sep='')
ggplot(df.long, aes(Term, -log10(value), fill=group)) + geom_bar(stat='identity',position='dodge') + labs(y='-log10(FDR)', x='', title=tf) + coord_flip() 
save_ = paste(pathway,'GO_enrichment_comparison_mrna_top5percent_E083_.png',sep='')
ggsave(save_,width=10,height=5,dpi=300)

# Side-by-side boxplots
file1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/Dominant_best_p-values_top150.txt')
file2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/Average_best_p-values_top150.txt')
file3 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/Sum_best_p-values_top150.txt')
file1$group = rep('Broadest',nrow(file1))
file2$group = rep('Mean',nrow(file2))
file3$group = rep('Sum',nrow(file3))
table_ = rbind(file1,file2,file3)
colnames(table_) = c('cell_type','value','group')
table = melt(table_)
table$value = -log10(table$value)
table$lev = factor(table$group, levels=c('Broadest','Mean','Sum'))
colnames(table) = c('cell type','group','method','value','lev')
p = ggplot(table, aes(lev, value)) 
p + geom_boxplot() + ggtitle('') +  ylab('Percentile Rank') + xlab('')
filename_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/comparison_methods_best_positions.png',sep='')
ggsave(filename_,width=5,height=5, dpi=300)

### Wilcoxon test for overlap (dominant vs smaller peaks)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Assigned_genes/all_peaks_2.5kb/more/variance_table.txt')
wilcox.test(V12 ~ V13, data=file_)
c=file_$V12[which(file_$V13=='True')]
d=file_$V12[which(file_$V13=='False')]


# Side-by-side boxplots
file_ = read.table('/Users/woojunshim/Research/Data/rp/Roadmap/results/best_pvalues_top30.txt')
file1 = file_[,c(1,2)]
file2 = file_[,c(1,3)]
file3 = file_[,c(1,4)]
file1$group = rep('Broadest',nrow(file1))
file2$group = rep('RP',nrow(file2))
file3$group = rep('Difference_from_threshold',nrow(file3))
colnames(file1) = c('cell_type','value','group')
colnames(file2) = c('cell_type','value','group')
colnames(file3) = c('cell_type','value','group')
table_ = rbind(file1,file2,file3)
colnames(table_) = c('cell_type','value','group')
table = melt(table_)
table$value = -log10(table$value)
table$lev = factor(table$group, levels=c('Broadest','RP','Difference_from_threshold'))
colnames(table) = c('cell type','group','method','value','lev')
p = ggplot(table, aes(lev, value)) 
p + geom_boxplot() + ggtitle('') +  ylab('-log10(p-value)') + xlab('')
filename_ = paste('/Users/woojunshim/Research/Data/rp/Roadmap/plots/comparison_best_pvalues_top30.png',sep='')
ggsave(filename_,width=5,height=5, dpi=300)

### Partition around medois (PAM) for clustering
library("cluster")
file_=read.table('/Users/woojunshim/Research/Data/Paige/1/analysis/ds_table_.txt')
file_=read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/broadpeak_tissue_enrichment_pvalues_sig.txt')
file_ = file_[which(rownames(file_) %in% mrna_$X.gene),]  # only mRNA genes
file_ = -log10(file_)
file__ = ifelse(file_<0.05, 1, 0)
heatmap(file__)
data_ = file_[,2:6]
rownames(data_) = file_$V1
data_ = data_[complete.cases(data_),]
colnames(data_) = c('day0','day2','day5','day9','day14')
# Finding the best k (cluster) number
l = daisy(data_)
asw = numeric(10)
for (k in 2:10){
  asw[k] = pam(l, k)$silinfo$avg.width
}
k.best = which.max(asw)
k.best = 4
clu = pam(l, k=k.best)
data_$cluster = clu$clustering
data1 = data_[clu$cluster==1,]
yrange = range(data1)
plot(c(1,5), yrange, type='n')
for (n in 1:nrow(data1)){
  lines(x=seq(1,5,1), y=data1[n,1:5], type='l', lty=1, col='blue')
}
si = silhouette(clu$clustering, dist(data_, "canberra"))
plot(clu, col=c("red", "green", "blue", "purple"))
clusplot(clu, clu$cluster, color=T, shade=T, labels=2, lines=0)

# Hierarchical clustering
data_ = file_
genes.cor = cor(t(data_), method='pearson')
genes.cor.dist = dist(data_, method='euclidian')
genes.cor.dist = as.dist(1- genes.cor)
genes.tree = hclust(genes.cor.dist, method='ward.D')
times.cor = cor(data_, method='spearman')
#times.cor.dist = dist(t(data_), method='euclidian')
times.cor.dist = as.dist(1-times.cor)
times.tree = hclust(times.cor.dist, method='ward.D')
library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
cut_ = cutree(genes.tree, k=7)
my.color = rep('green', times=nrow(data_))
my.color[cut_==2] = 'blue'
my.color[cut_==3] = 'red'
my.color[cut_==4] = 'orange'
my.color[cut_==5] = 'brown'
my.color[cut_==6] = 'purple'
my.color[cut_==7] = 'cyan'
#my.color
heatmap(as.matrix(data_), scale='row', col=col, Rowv=as.dendrogram(genes.tree), Colv=as.dendrogram(times.tree), RowSideColors = my.color, labRow = F)
names_green = names(cut_)[cut_=='1']
names_blue = names(cut_)[cut_=='2']
names_red = names(cut_)[cut_=='3']
names_orange = names(cut_)[cut_=='4']
names_brown = names(cut_)[cut_=='5']
names_purple = names(cut_)[cut_=='6']
names_cyan = names(cut_)[cut_=='7']
results = matrix(nrow=nrow(data_), ncol=2)
results[,1] = c(names_green,names_blue,names_red,names_orange,names_brown,names_purple,names_cyan)
results[,2] = c(rep('green',length(names_green)),rep('blue',length(names_blue)),rep('red',length(names_red)),rep('orange',length(names_orange)),rep('brown',length(names_brown)),rep('purple',length(names_purple)),rep('cyan',length(names_cyan)))
write.table(results, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7cluster_table_.txt', sep='\t', quote=F)
summary_table = matrix(nrow=length(unique(my.color))+1, ncol=2)
summary_table[,1] = c('green','blue','red','orange','brown','purple','cyan','total')
summary_table[,2] = c(length(names_green),length(names_blue),length(names_red),length(names_orange),length(names_brown),length(names_purple),length(names_cyan),length(my.color))
write.table(summary_table, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7cluster_summary_.txt', sep='\t', quote=F)
# ggplot plotting 
library(ggplot2)
library(reshape2)
data_$cluster = my.color
color_ = 'orange'
data1 = as.matrix(data_[data_$cluster==color_,1:5])
table_ = matrix(nrow=5, ncol=6)
rownames(table_) = colnames(data1)
colnames(table_) = c('mean_','max_','min_','time','upper.q','lower.q')
for (i in 1:nrow(table_)){
  table_[i,1] = mean(data1[,i])
  table_[i,2] = max(data1[,i])
  table_[i,3] = min(data1[,i])
  table_[i,4] = rownames(table_)[i]
  table_[i,5] = quantile(data1[,i], 0.75, type=1)
  table_[i,6] = quantile(data1[,i], 0.25, type=1)
}
yrange = range(data_[,1:5])
plot(c(1,5), c(1,5), type='n', xaxt='n', xlab='', ylab='DS')
axis(1, at=1:5, labels=colnames(data1))
lines(x=seq(1,5,1), y=table_[,1], col=color_, lwd=3)
#lines(x=seq(1,5,1), y=table_[,2], col=color_, lwd=1, lty=3)
#lines(x=seq(1,5,1), y=table_[,3], col=color_, lwd=1, lty=3)
lines(x=seq(1,5,1), y=table_[,5], col=color_, lwd=2, lty=6)
lines(x=seq(1,5,1), y=table_[,6], col=color_, lwd=2, lty=6)

### DISTANCE ANALYSIS
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hotspots_distances.txt')
mrna_ = read.csv('/Users/woojunshim/Research/Data/mRNA_genes.txt')
tf_ = read.csv('/Users/woojunshim/Research/Data/TF/TF_combined.txt')
file__ = file_[file_$V4!=0, ]
file__ = file__[file__$V6 %in% mrna_$X.gene,]
file__ = file__[file__$V6 %in% tf_$X.TF_combined,]
plot(x=file__$V7, y=file__$V4, type='h', xlab='Distance to TSS', ylab='Overlap count')

### GENE ONTOLOGY ANALYSIS using TopGO
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
goterm = as.matrix(Term(GOTERM))
write.table(goterm, '/Users/woojunshim/Research/Data/GO_Terms/GO_table.txt', sep='\t', quote=F)
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
source("/Users/woojunshim/Nathan/Scripts/fisherTopGO.R")
source("/Users/woojunshim/Nathan/Scripts/GOenrichTopTerms.R")
background_file = read.table('/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7cluster_table_.txt')
mrna_ = read.table('/Users/woojunshim/Research/Data/mRNA_genes.txt')
all_genes = mrna_$V1
groups = c('blue','brown','cyan','red','orange','purple','green')
top_no = 100 # top 100 GO terms by FDR
for (group in groups){
  genelist = background_file[background_file$V2==group, 1]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### PLOT HORIZONTAL BAR GRAPHS (GO TERMS)
library(ggplot2)
for (group in groups){
  file__ = paste('/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/',group,'_GO_BP.txt',sep='')
  file_ = read.table(file__)
  table_ = data.frame(tt=gsub('_', ' ', file_$V3[1:5]), fdr=-log10(file_$V4[1:5]))
  table_$tt = reorder(table_$tt, table_$fdr)
  ggplot(table_, aes(x=tt, y=fdr)) + geom_bar(stat='identity', fill=group, colour='black') + coord_flip() + ylab('-log10(FDR)') + xlab('') + theme(axis.text=element_text(size=12, face='bold')) + scale_y_continuous(limits=c(0, 15.5))
  output_= paste('/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/',group,'_FDR.png',sep='')
  ggsave(output_, width=10, height=5)
}

### PLOT LINE GRAPHS (7 CLUSTERS)
yrange = range(data_[,1:5])
pdf('/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7clusters_comparison.pdf', width=10, height=5)
plot(c(1,5), c(1,3), type='n', xaxt='n', xlab='', ylab='DS')
axis(1, at=1:5, labels=colnames(data_)[1:5])
for (group in groups){
  data1 = data_[data_$cluster==group,]
  lines(x=seq(1,5,1), y=apply(as.matrix(data1[1:5]), 2, mean), col=group, lwd=3)
}
dev.off()

### PLOT LINE GRAPHS FOR FREQUENCIES OF ELEMENTARY INTERVALS
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hotspots_category_counts.txt')
idx = which((file_$bin>=-100) & (file_$bin<=100))
file_ = file_[idx,]
colnames(file_) = c('bin','x>50','10<x<=50','x<=10')
library(ggplot2)
table_1 = data.frame(bin=file_$bin, cat=file_$`x>50`, group=rep('x>50'))
table_2 = data.frame(bin=file_$bin, cat=file_$`10<x<=50`, group=rep('10<x<=50'))
table_3 = data.frame(bin=file_$bin, cat=file_$`x<=10`, group=rep('x<=10'))
table_ = rbind(table_1, table_2, table_3)
ggplot(table_, aes(x=bin, y=cat, colour=group, group=group))  + geom_line() +  labs(x='Distance to TSS (kb)', y='Frequency')
ggsave('/Users/woojunshim/Research/Data/broadPeaks/overlap/elementary_intervals_distribution.png', width=10, height=5)

### HEATMAP FOR SIGNIFICANT GENES (P<0.05) ACROSS CELL TYPES ANALYSIS
file_=read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/gene_entropy_enrichment_pvalues_.txt')
file_ = file_[which(rownames(file_) %in% mrna_$X.gene),]  # only mRNA genes
file__ = ifelse(file_<0.05, 1, 0)
idx = vetor()
for (no in 1:nrow(file__)){
  if (sum(file__[no,])!=0){
    idx = c(idx, no)
  }
}
file__ = file__[idx,]
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/across_tissue_types_fet_significant.pdf', width=5, height=5)
no_labs = rep('', ncol(file__))
heatmap(as.matrix(file__),col=c('white','red'), labRow=no_labs, cexCol=0.7)
dev.off()

# GO analysis for each groups of significant genes
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
groups = colnames(file__)
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (no in 1:ncol(file__)){
  group = colnames(file__)[no]
  idx = which(file__[,no]==1)
  genelist = rownames(file__)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/GO/',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### HEATMAP FOR SIGNIFICANT GO TERMS ACROSS DEFINED GROUPS
### RUN THIS FOR PLOTTING ANY PREDEFINED GROUPS 
groups = colnames(file__)
groups = c('orange','blue','cyan','purple','green','brown','red')
pathway_ = '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/Tissues/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/GO/'
top_no = 5  # Number of GO terms to be included
terms = vector()
for (group in groups){
  filename_ = paste(pathway_,group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (i in 1:top_no){
    terms = c(terms, as.character(file_$V3[i]))
  }
}
terms= unique(terms)
table_ = data.frame(matrix(data=rep(1,length(groups)*length(terms)), ncol=length(groups), nrow=length(terms)))
colnames(table_) = groups
rownames(table_) = terms
for (no in 1:ncol(table_)){
  group = colnames(table_)[no]
  filename_ = paste(pathway_,group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (m in 1:nrow(table_)){
    t = rownames(table_)[m]
    idx = which(file_$V3==t)
    if (length(idx)!=0){
      table_[t, group] = file_[idx,4] 
    }
  }
}
table_ = -log10(table_)
library(ggplot2)
library(reshape)
table_$tt = gsub('_', ' ', rownames(table_))
table_$tt <- factor(table_$tt, levels = table_$tt)
table = melt(table_)

p = ggplot(table, aes(y=tt, x=variable, fill=value)) 
p + geom_tile()+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
filename_ = paste(pathway_,'GO_enrichment_plot_tissue.png',sep='')
ggsave(filename_,width=10,height=10)

# Percentile box plots
# Boxplots 
file_1 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/fet_filtering_percentile.txt')
file_2 = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/fet_extraction_percentile.txt')
file__1 = file_1[,2:ncol(file_1)]
file__2 = file_2[,2:ncol(file_2)]
rownames(file__1) = file_1[,1]
colnames(file__1) = seq(99,0,by=-1)
rownames(file__2) = file_2[,1]
colnames(file__2) = seq(99,0,by=-1)
file__1 = t(file__1)
file__2 = t(file__2)
file__1 = -log10(file__1)
file__2 = -log10(file__2)
table1 = melt(file__1)
table2 = melt(file__2)
table1$method = rep('Filtering', nrow(table1))
table2$method = rep('Simple Extraction', nrow(table2))
table = rbind(table1, table2)
table$lev = factor(table$Var1, levels=seq(99,0,by=-1))
p = ggplot(table, aes(lev, value, fill=method))
p + geom_boxplot() + ggtitle('Top 150 TFs') + xlab('Percentile Rank') + ylab('-log10(p-value)') + scale_x_discrete(breaks=seq(95,0,-5))
filename_ = paste('/Users/woojunshim/Research/Data/highly_expressed_genes/broadpeak/top150/top150_filtering_extraction_fet.png',sep='')
ggsave(filename_,width=10,height=5)


## TopGo 
# GO analysis for each groups of significant genes
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
groups = c('day0_day2','day2_day5','day5_day9','day9_day14')
for (group in groups){
  filename_ = paste('/Users/woojunshim/Research/Data/Paige/1/',group,'_.txt',sep='')
  file_ = read.table(filename_)
  idx = which(file_$V5<0.01)
  genelist = file_[idx,1]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/Paige/1/',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

#Plots 
pathway_ = '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/Tissues/'
pathway_ = '/Users/woojunshim/Research/Data/Paige/1/'
top_no = 10  # Number of GO terms to be included
terms = vector()
new_groups = c('day2','day5','day9','day14')
for (group in groups){
  filename_ = paste(pathway_,group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (i in 1:top_no){
    terms = c(terms, as.character(file_$V3[i]))
  }
}
terms= unique(terms)
table_ = data.frame(matrix(data=rep(1,length(groups)*length(terms)), ncol=length(groups), nrow=length(terms)))
colnames(table_) = groups
rownames(table_) = terms
for (no in 1:ncol(table_)){
  group = colnames(table_)[no]
  filename_ = paste(pathway_,group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (m in 1:nrow(table_)){
    t = rownames(table_)[m]
    idx = which(file_$V3==t)
    if (length(idx)!=0){
      table_[t, group] = file_[idx,4] 
    }
  }
}
table_ = -log10(table_)
library(ggplot2)
library(reshape)
table_$tt = gsub('_', ' ', rownames(table_))
table_$tt <- factor(table_$tt, levels = table_$tt)
#colnames(table_) = new_groups
table = melt(table_)

p = ggplot(table, aes(y=tt, x=variable, fill=value)) 
p + geom_tile()+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(breaks=groups,labels=new_groups)
filename_ = paste(pathway_,'GO_enrichment_plot_temporal.png',sep='')
ggsave(filename_,width=10,height=10)

### GO analysis for Tsankov data sets
### Also calculates Wilcoxon rank sum tests on normalised ChIP-seq binding events
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/functional_analysis/'
input_data = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Analysis_results/Tsankov/results.txt')
all_genes = mrna_$X.gene
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
gene = 'Eomes'

## Calculate Wilcox test for the top 5 % gene
data_ = input_data[order(input_data$V3, decreasing = TRUE),]
no_ = round(nrow(data_) * 0.05)
a = data_$V10[1:no_] / (data_$V3[1:no_]/1000)
data_ = input_data[order(input_data$V7, decreasing = TRUE),]
no_ = round(nrow(data_) * 0.05)
b = data_$V10[1:no_] / (data_$V3[1:no_]/1000)
wilcox.test(a,b)

sub_data = input_data[which(input_data$V8!=0),]  # Remove genes without ChIP-seq bindings
groups = c('Width','FPKM')
colnames(sub_data) = c('peak_id','gene','Width','chr','start','end','FPKM','Eomes','T','Pou5f1','CTCF')
for (group in groups){
  data_ = sub_data[order(sub_data[,which(colnames(sub_data)==group)], decreasing=TRUE),]
  genelist = data_[1:100,2]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste(pathway_ ,gene,'_',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}


top_no = 10  # Number of GO terms to be included
terms = vector()
new_groups = c('H3K4me3','FPKM')
for (group in groups){
  filename_ = paste(pathway_,gene,'_',group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (i in 1:top_no){
    terms = c(terms, as.character(file_$V3[i]))
  }
}
terms= unique(terms)
table_ = data.frame(matrix(data=rep(1,length(groups)*length(terms)), ncol=length(groups), nrow=length(terms)))
colnames(table_) = groups
rownames(table_) = terms
for (no in 1:ncol(table_)){
  group = colnames(table_)[no]
  filename_ = paste(pathway_,gene,'_',group,'_GO_BP.txt', sep='')
  file_ = read.table(filename_)
  for (m in 1:nrow(table_)){
    t = rownames(table_)[m]
    idx = which(file_$V3==t)
    if (length(idx)!=0){
      table_[t, group] = file_[idx,4] 
    }
  }
}
table_ = -log10(table_)
library(ggplot2)
library(reshape)
table_$tt = gsub('_', ' ', rownames(table_))
table_$tt <- factor(table_$tt, levels = table_$tt)
#colnames(table_) = new_groups
table = melt(table_)

p = ggplot(table, aes(y=tt, x=variable, fill=value)) 
p + geom_tile()+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(breaks=groups,labels=new_groups)
filename_ = paste(pathway_,gene,'_GO_enrichment_plot.png',sep='')
ggsave(filename_,width=10,height=10, dpi=300)

data_ = read.table('/Users/woojunshim/Research/Data/Paige/1/h3k4me3_day14_assigned.txt')

### CORRELATION ANALYSIS ON 'DS_change_table.txt'
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_ds/DS_change_table.txt')
cor_ = cor(file_, method='pearson')

### DO I NEED TO CONSIDER NORMALISATION AGAIN? 
### THEY SEEM TO BE HIGHLY CORRELATED ANYWAY
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_dominant_1.0_broadpeak_mrna_.txt')
epigenomes = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/epigenomes_list.txt', stringsAsFactors = F)
for (epi in epigenomes$V1){
  temp = file_[file_$V2==epi,]
  med_ = median(temp$V4)
  mad_ = mad(temp$V4)
  temp$V4 = (temp$V4-med_) / mad_
  output_file = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_normalisation/',epi,'.txt',sep='')
  write.table(temp, output_file, quote=F, sep='\t')
}

### PCA FOR THE CONTRIBUTION ANALYSIS (ROADMAP DATA)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_ds/mle/epigenomes/Roadmap_summary.txt')
pca_ = prcomp(file_, center=T, scale.=T)
file_$Tissue = rep('', nrow(file_))
for (no in 1:nrow(file_)){
  epi_name = rownames(file_)[no]
  idx = which(epigenomes$V1==epi_name)
  t = as.character(epigenomes$V2[idx])
  file_$Tissue[no] = t
}
library(ggfortify)
autoplot(pca_, data=file_, colour='Tissue')

### CLUSTERING OF 'ALL_DS.TXT' DS FOR ALL GENES ACROSS 111 EPIGENOMES
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/broadpeak_combined_foreground_counts_intra.txt.txt')
data_ = file_
genes.cor = cor(t(data_), method='pearson')
genes.cor.dist = dist(data_, method='euclidian')
genes.cor.dist = as.dist(1- genes.cor)
genes.tree = hclust(genes.cor.dist, method='ward.D')
times.cor = cor(data_, method='spearman')
#times.cor.dist = dist(t(data_), method='euclidian')
times.cor.dist = as.dist(1-times.cor)
times.tree = hclust(times.cor.dist, method='ward.D')
library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
cut_ = cutree(genes.tree, k=10)
my.color = rep('green', times=nrow(data_))
my.color[cut_==2] = 'blue'
my.color[cut_==3] = 'red'
my.color[cut_==4] = 'orange'
my.color[cut_==5] = 'brown'
my.color[cut_==6] = 'purple'
my.color[cut_==7] = 'cyan'
my.color[cut_==8] = 'grey'
my.color[cut_==9] = 'pink'
my.color[cut_==10] = 'black'
#my.color
col_labels =vector()
for (no in 1:ncol(file_)){
  epi_name = colnames(file_)[no]
  idx = which(epigenomes$V1==epi_name)
  t = as.character(epigenomes$V2[idx])
  col_labels = c(col_labels, t)
}
pol = rainbow(length(unique((col_labels))))
col__ = pol[as.numeric(as.factor(col_labels))]
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust_all_genes_pearson_binary_10clusters_.pdf',width=10,height=10)
heatmap(as.matrix(data_), scale='none', col=c('white','red'), Rowv=as.dendrogram(genes.tree), Colv=as.dendrogram(times.tree), RowSideColor=my.color, labRow = F)
dev.off()
file_$Colour = col__
names_green = names(cut_)[cut_=='1']
names_blue = names(cut_)[cut_=='2']
names_red = names(cut_)[cut_=='3']
names_orange = names(cut_)[cut_=='4']
names_brown = names(cut_)[cut_=='5']
names_purple = names(cut_)[cut_=='6']
names_cyan = names(cut_)[cut_=='7']
results = matrix(nrow=nrow(data_), ncol=2)
results[,1] = c(names_green,names_blue,names_red,names_orange,names_brown,names_purple,names_cyan)
results[,2] = c(rep('green',length(names_green)),rep('blue',length(names_blue)),rep('red',length(names_red)),rep('orange',length(names_orange)),rep('brown',length(names_brown)),rep('purple',length(names_purple)),rep('cyan',length(names_cyan)))
write.table(results, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7cluster_table_.txt', sep='\t', quote=F)
summary_table = matrix(nrow=length(unique(my.color))+1, ncol=2)
summary_table[,1] = c('green','blue','red','orange','brown','purple','cyan','total')
summary_table[,2] = c(length(names_green),length(names_blue),length(names_red),length(names_orange),length(names_brown),length(names_purple),length(names_cyan),length(my.color))
write.table(summary_table, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7cluster_summary_.txt', sep='\t', quote=F)

### PLOT LINES WITH TISSUE GROUP COLOURS
a=unique(col__)
b=unique(col_labels)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/colour_legend_tissue_groups_new_final.pdf', width=10, height=10)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
#cc = as.numeric(as.factor(col_labels))
legend('center', b, col=a, lty=1, lwd=2, cex=0.8)
dev.off()

### IDENTIFY NUMBER OF CLUSTERS 
require(cluster)
library(factoextra)
fviz_nbclust(file_, hcut, method = "silhouette",hc_method = "ward.D")  # Gives 2 as the optimal 

library('NbClust')
results1 = NbClust(data_, diss=genes.cor.dist, distance=NULL, min.nc=2, max.nc=20, index='silhouette', method='ward.D')
results1$Best.partition[which(names(results1$Best.partition)=='C15orf41')]
# Let's test
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_dominant_1.0_broadpeak_mrna_.txt')
temp = file_[file_$V2=='E071',]
temp = temp[order(temp$V5, decreasing=T),]
ref_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_clusters.txt')
ll = vector()
for (n in 1:500){
  gene_ = temp[n,1]
  cluster = results1$Best.partition[which(names(results1$Best.partition)==gene_)]
  ll = c(ll, cluster)
}



bb = vector()
for (n in 1:nrow(temp)){
  gene_ = temp[n,1]
  cluster = results1$Best.partition[which(names(results1$Best.partition)==gene_)]
  bb = c(bb, cluster)
}

# GO analysis
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
ref_list = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_clusters.txt')
groups = unique(as.character(ref_list$V2))
for (group in groups){
  idx = which(ref_list$V2==group)
  genelist = ref_list$V1[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

group ='E071'
cluster_ = 8
genes = names(results1$Best.partition)[results1$Best.partition==cluster_]
geneList<-factor(as.integer(all_genes %in% genes))
names(geneList)<-all_genes
TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genes,nodeSize=5,annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
aa = score(resultFisher)
aa = p.adjust(aa, method='fdr')  # FDR
aa = sort(aa, decreasing=FALSE)  # sort by FDR
bb = as.character(go_table[names(aa),1])
results = matrix(ncol=3,nrow=top_no)  # result table 
results[,1] = names(aa)[1:top_no]
results[,2] = bb[1:top_no]
results[,3] = as.numeric(aa)[1:top_no]
colnames(results) = c('#GO_ID','GO_Term','FDR')
file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/',group,'_cluster',cluster_,'.txt', sep='')
write.table(results, file_, sep='\t', quote=F)


names_green = names(cut_)[cut_=='1']
names_blue = names(cut_)[cut_=='2']
names_red = names(cut_)[cut_=='3']
names_orange = names(cut_)[cut_=='4']
names_brown = names(cut_)[cut_=='5']
names_purple = names(cut_)[cut_=='6']
names_cyan = names(cut_)[cut_=='7']
names_grey = names(cut_)[cut_=='8']
names_pink = names(cut_)[cut_=='9']
names_black = names(cut_)[cut_=='10']
results = matrix(nrow=nrow(data_), ncol=2)
results[,1] = c(names_green,names_blue,names_red,names_orange,names_brown,names_purple,names_cyan,names_grey,names_pink,names_black)
results[,2] = c(rep('green',length(names_green)),rep('blue',length(names_blue)),rep('red',length(names_red)),rep('orange',length(names_orange)),rep('brown',length(names_brown)),rep('purple',length(names_purple)),rep('cyan',length(names_cyan)),rep('grey',length(names_grey)),rep('pink',length(names_darkgreen)),rep('black',length(names_black)))
write.table(results, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_clusters.txt', sep='\t', quote=F)
summary_table = matrix(nrow=length(unique(my.color))+1, ncol=2)
summary_table[,1] = c('green','blue','red','orange','brown','purple','cyan','grey','pink','black','total')
summary_table[,2] = c(length(names_green),length(names_blue),length(names_red),length(names_orange),length(names_brown),length(names_purple),length(names_cyan),length(names_grey),length(names_pink),length(names_black),length(my.color))
write.table(summary_table, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_cluster_summary.txt', sep='\t', quote=F)

#### ANALYSIS FOR ROADMAP DATA USING NEWLY DEVELOPED CLUSTERING APPROACH
library(topGO)
library(GO.db)
cutoff_table = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/Roadmap_cutoff_positions.txt', stringsAsFactors=F)
ref_table = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_clusters.txt',stringsAsFactors=F)
file__ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_dominant_1.0_broadpeak_mrna_.txt', stringsAsFactors=F)
mrna_ = read.csv('/Users/woojunshim/Research/Data/mRNA_genes.txt')
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
all_genes = mrna_$X.gene
epigenomes = cutoff_table$V1
for (epi in epigenomes){
  temp = file__[file__$V2==epi,]
  temp = temp[order(temp$V5, decreasing=T),]
  background = vector()
  foreground = vector()
  cutoff = cutoff_table[cutoff_table$V1==epi, 2]
  for (no in 1:nrow(temp)){
    gene_ = temp[no,1]
    background = c(background, ref_table[ref_table$V2==gene_,3])
  }
  for (no in 1:cutoff){
    gene_ = temp[no,1]
    foreground = c(foreground, ref_table[ref_table$V2==gene_,3])
  }
  all_clusters = unique(foreground)
  t1 = table(foreground)
  t2 = table(background)
  results = matrix(nrow=length(all_clusters), ncol=4)
  colnames(results) = c('#cluster','FET','Odds','No.genes_in_the_cluster')
  results[,1] = all_clusters
  sig_clusters = vector()
  for (no in 1:length(all_clusters)){
    cluster = all_clusters[no]
    a = t1[cluster]
    b = sum(t1) - a
    c = t2[cluster]
    d = sum(t2) - c
    con_t = matrix(c(a,b,c,d), nrow=2)
    p_ = fisher.test(con_t)$p.value
    odd_ = fisher.test(con_t)$estimate
    results[no,2] = p_
    results[no,3] = odd_
    results[no,4] = a
    if (p_<0.05 && odd_>1.0){
      sig_clusters = c(sig_clusters, cluster)
    }
  }
  output_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/',epi,'_clusters.txt',sep='')
  results = results[order(results[,3], decreasing=T),]
  write.table(results, output_, sep='\t', quote=F)
  
  # GO analysis comparing between selected genes vs. simple extraction from Benayoun's
  selected= vector()
  benayoun= vector()
  cnt = 0
  for (no in 1:cutoff){
    gene_ = temp[no,1]
    colour_ = ref_table[which(ref_table$V2==gene_),3]
    if (colour_ %in% sig_clusters){
      selected = c(selected, gene_)
      cnt = cnt +1
    }
  }
  if (cnt!=0){
    for (no in 1:cnt){
      gene_ = temp[no,1]
      benayoun = c(benayoun, gene_)
    }
  }
  if (cnt!=0){
    geneList<-factor(as.integer(all_genes %in% selected))
    names(geneList)<-all_genes
    TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genes,nodeSize=5,annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
    resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
    aa = score(resultFisher)
    aa = p.adjust(aa, method='fdr')  # FDR
    aa = sort(aa, decreasing=FALSE)  # sort by FDR
    bb = as.character(go_table[names(aa),1])
    results = matrix(ncol=3,nrow=top_no)  # result table 
    results[,1] = names(aa)[1:top_no]
    results[,2] = bb[1:top_no]
    results[,3] = as.numeric(aa)[1:top_no]
    colnames(results) = c('#GO_ID','GO_Term','FDR')
    file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/epigenomes/',epi,'_selected_genes_by_clusters.txt', sep='')
    write.table(results, file_, sep='\t', quote=F)
    
    geneList<-factor(as.integer(all_genes %in% benayoun))
    names(geneList)<-all_genes
    TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genes,nodeSize=5,annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
    resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
    aa = score(resultFisher)
    aa = p.adjust(aa, method='fdr')  # FDR
    aa = sort(aa, decreasing=FALSE)  # sort by FDR
    bb = as.character(go_table[names(aa),1])
    results = matrix(ncol=3,nrow=top_no)  # result table 
    results[,1] = names(aa)[1:top_no]
    results[,2] = bb[1:top_no]
    results[,3] = as.numeric(aa)[1:top_no]
    colnames(results) = c('#GO_ID','GO_Term','FDR')
    file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/epigenomes/',epi,'_simple_extraction.txt', sep='')
    write.table(results, file_, sep='\t', quote=F)
  }
}

### Generate grid table for conditional probabilities
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_conditional_probabilities_odds_ratio.txt',stringsAsFactors=F)
file_ = t(round(file_, digits = 3))
library(gridExtra)
library(grid)
grid.table(file_)

library(ggplot2)
library(reshape2)
table = melt(file_)
colnames(table) = c('Var1','Var2','probability')
ggplot(table, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=probability)) + geom_text(aes(label=probability)) + scale_fill_gradient(low='white', high='red') + xlab('Condition B') + ylab('Condition A') 
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/conditional_probability.pdf', width=10, height=10)

### PLOT GENERAL PROBABILITY OF OVER-REPRESENTATION OF EACH CLUSTER (BAR GRAPH)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/Roadmap_odds_ratios.txt',stringsAsFactors=F)
results = vector()
names_ = vector()
for (no in 1:ncol(file_)){
  total_cnt =0
  cnt = 0
  names_ = c(names_, colnames(file_)[no])
  for (n in 1:nrow(file_)){
    total_cnt = total_cnt +1
    if (file_[n,no] > 1.0){
      cnt = cnt +1
    }
  }
  results = c(results, cnt/total_cnt)
}
names(results) = names_
table = t(data.frame(results))
table_ = melt(table)
table_$Var2 = reorder(table_$Var2, table_$value)
colnames(table_) = c('results','Cluster','Probability')
ggplot(table_, aes(x=Cluster, y=Probability)) + geom_bar(stat='identity') + geom_text(aes(label=round(Probability,3)), vjust=1.5, colour="white") 
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_probability.pdf', width=10, height=10)

table = data.frame(cluster=results)
rownames(table) = names_
write.table(table, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_probability.txt', sep='\t', quote=F)

### CALCULATE DIFFERENCE BETWEEN ALL PAIRS OF CONDITIONAL PROBABILITIES
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_conditional_probabilities_odds_ratio_.txt',stringsAsFactors=F)
gp = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_marginal_probabilities.txt',stringsAsFactors=F)
results = data.frame(matrix(nrow=9, ncol=9))
colnames(results) = colnames(file_)
rownames(results) = rownames(file_)
all = vector()
for (a in 1:nrow(file_)){
  name_a = rownames(file_)[a]
  for (b in 1:ncol(file_)){
    name_b = colnames(file_)[b]
    if (a==b){
      results[a,b] = 0.0
    }
    p_a = gp[gp$V1==name_a,2]
    p_b = gp[gp$V1==name_b,2]
    diff_ = file_[a,b] - p_a
    all = c(all, diff_)
    results[a,b] = diff_
  }
}
table = t(t(data.frame(results)))
table_ = melt(table)
table_$Var2 = reorder(table_$Var2, as.numeric(table_$Var2))
colnames(table_) = c('A','B','Difference')
ggplot(table_, aes(x=B, y=A, reorder(table_$B, as.numeric(table_$B)))) + geom_tile(aes(fill=Difference)) + geom_text(aes(label=round(Difference,3))) + scale_fill_gradient2(midpoint=0, high="red",mid="white",low="blue") + xlab('Condition B') + ylab('Condition A') + ggtitle('p(A|B) - p(A)')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_dependence.pdf', width=10, height=10)
write.table(results, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_dependence.txt', sep='\t', quote=F)

z_table = (as.matrix(results) - mean(as.matrix(results))) / sd(as.matrix(results))
write.table(t(z_table), '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/joint_probability_difference_z.txt', sep='\t', quote=F)
rt = melt(results)
colnames(rt) = c('cluster','difference')
ggplot(rt, aes(x=difference)) + geom_histogram(bins=20) + geom_density(col=2)
hist(as.matrix(results), prob=T, col='grey', breaks=10)
lines(density(as.matrix(results)), col='red', lwd=2)

### CORRELATION BETWEEN CLUSTERS (ODDS RATIOS)
file_ = read.file('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/Roadmap_odds_ratios.txt')
cor_= cor(file_)
table = t(data.frame(cor_))
table_ = melt(table)
table_$Var2 = reorder(table_$Var2, table_$value)
colnames(table_) = c('Cluster1','Cluster2','r')
ggplot(table_, aes(x=Cluster1, y=Cluster2)) + geom_tile(aes(fill=r)) + geom_text(aes(label=round(r,3))) + scale_fill_gradient2(midpoint=0, high="red",mid="white",low="blue") +xlab('') +ylab('')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_pearsons.pdf', width=10, height=10)

### PLOT NETWORK GRAPH 
library(igraph)
set.seed(555)
file_p = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/positive_clusters.txt', stringsAsFactors = F)
file_n = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/negative_clusters.txt', stringsAsFactors = F)
labels_p = vector()
labels_n = vector()
for (i in 1:nrow(file_p)){
  labels_p = c(labels_p, file_p[i,1], file_p[i,2])
}
for (i in 1:nrow(file_n)){
  labels_n = c(labels_n, file_n[i,1], file_n[i,2])
}
labels = labels_n
gd <- graph(labels, directed=T)
col_ = unique(labels)
labels__ = rep('', length(col_))
plot(gd, vertex.size=20, vertex.color=col_, edge.arrow.size=0.5)

### CREATE A TABLE OF DIFFERENCE (CONDITIONAL PROBABILITY - MARGINAL PROBABILITY)
cond = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_conditional_probabilities_odds_ratio.txt')
marg = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_probability.txt')
for (i in 1:nrow(cond)){
  cluster = rownames(cond)[i]
  marg_ = marg[marg$V1==cluster,2]
  for (j in 1:ncol(cond)){
    if (i==j){
      cond[i,j] = 0.0
    } else {
      cond[i,j] = cond[i,j] - marg_
    }
  }
}
write.table(cond, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/cluster_probability_difference.txt', sep='\t', quote=F)

### CREATE DATA POINTS FOR IGRAPH
pos_pairs = vector()
pos_values = vector()
pos_colours = vector()
neg_pairs = vector()
neg_values = vector()
neg_colours = vector()
threshold_ = 0.05
for (i in 1:ncol(cond)){
  start = colnames(cond)[i]
  for (j in 1:nrow(cond)){
    end = rownames(cond)[j]
    value_ = cond[j,i]
    if (value_ > +threshold_){
      pos_pairs = c(pos_pairs, start, end)
      pos_values = c(pos_values, value_)
      pos_colours = c(pos_colours, start, end)
    }
    if (value_ <= -threshold_){
      neg_pairs = c(neg_pairs, start, end)
      neg_values = c(neg_values, -value_)
      neg_colours = c(neg_colours, start, end)
    }
  }
}
pos_values = pos_values / max(pos_values)
neg_values = neg_values / max(neg_values)

library(igraph)
set.seed(555)
labels = neg_pairs
gd <- graph(labels, directed=T)
col_ = unique(labels)
labels__ = rep('', length(col_))
plot(gd, vertex.size=20, vertex.color=col_, edge.arrow.size=neg_values)

### GENES CONDITIONAL PROBABILITIES
cond = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_ds/h3k4me3_day14_assigned_.txt_cond_prob')
cond = read.table('/Users/woojunshim/Research/Data/Paige/1/h3k4me3_day14_assigned_.txt_cond_prob')
cond = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_ds/E095_new_ds.txt_cond_prob')
cond['NKX2-5','ATP2A2']
cond['HAND2','ATP2A2']
length(which(cond['ISL1',]<0.05))
colnames(cond)[which(cond['HAND2',]<0.05)]

### GO ANALYSIS FOR GENES IN THE SELECTED TOP CLUSTER
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_dominant_1.0_broadpeak_mrna_.txt',stringsAsFactors=F )
file_ = file_[which(rownames(file_) %in% mrna_$X.gene), ]
file_ = data_
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = rownames(file_)
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
groups = colnames(file_)
for (group in groups){
  idx = which(file_[,group]==1)
  genelist = rownames(file_)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file__ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/1/',group,'_GO_BP_1.0.txt', sep='')
  write.table(results, file__, sep='\t', quote=F)
}

### CLUSTERING OF GENES WITHIN THE INTRA-THRESHOLD ACROSS 111 EPIGENOMES
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/broadpeak_combined_foreground_counts_intra.txt')
lo = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new_broadpeak__tf_cut_off_list.txt')
genes_ = lo[which(lo$V2!=999999),1]
data_ = file_[which(rownames(file_) %in% genes_),]
genes.cor = cor(t(data_), method='pearson')
#genes.cor.dist = dist(data_, method='euclidian')
genes.cor.dist = as.dist(1- genes.cor)
genes.tree = hclust(genes.cor.dist, method='ward.D')
times.cor = cor(data_, method='spearman')
#times.cor.dist = dist(t(data_), method='euclidian')
times.cor.dist = as.dist(1-times.cor)
times.tree = hclust(times.cor.dist, method='ward.D')
library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
cut_ = cutree(genes.tree, k=10)
my.color = rep('green', times=nrow(data_))
my.color[cut_==2] = 'blue'
my.color[cut_==3] = 'red'
my.color[cut_==4] = 'orange'
my.color[cut_==5] = 'brown'
my.color[cut_==6] = 'purple'
my.color[cut_==7] = 'cyan'
my.color[cut_==8] = 'grey'
my.color[cut_==9] = 'pink'
my.color[cut_==10] = 'black'
#my.color
epigenomes = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/epigenomes_list_.txt')
col_labels =vector()
for (no in 1:ncol(file_)){
  epi_name = colnames(file_)[no]
  idx = which(epigenomes$V1==epi_name)
  t = as.character(epigenomes$V2[idx])
  col_labels = c(col_labels, t)
}
pol = rainbow(length(unique((col_labels))))
col__ = pol[as.numeric(as.factor(col_labels))]
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/new/hclust_15407genes_2.0.pdf',width=10,height=10)
heatmap(as.matrix(data_), scale='none', col=c('white','red'), Rowv=as.dendrogram(genes.tree), Colv=as.dendrogram(times.tree), ColSideColor=col__, labRow = F)
dev.off()
file_$Colour = col__
names_green = names(cut_)[cut_=='1']
names_blue = names(cut_)[cut_=='2']
names_red = names(cut_)[cut_=='3']
names_orange = names(cut_)[cut_=='4']
names_brown = names(cut_)[cut_=='5']
names_purple = names(cut_)[cut_=='6']
names_cyan = names(cut_)[cut_=='7']
results = matrix(nrow=nrow(data_), ncol=2)
results[,1] = c(names_green,names_blue,names_red,names_orange,names_brown,names_purple,names_cyan)
results[,2] = c(rep('green',length(names_green)),rep('blue',length(names_blue)),rep('red',length(names_red)),rep('orange',length(names_orange)),rep('brown',length(names_brown)),rep('purple',length(names_purple)),rep('cyan',length(names_cyan)))
write.table(results, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7cluster_table_.txt', sep='\t', quote=F)
summary_table = matrix(nrow=length(unique(my.color))+1, ncol=2)
summary_table[,1] = c('green','blue','red','orange','brown','purple','cyan','total')
summary_table[,2] = c(length(names_green),length(names_blue),length(names_red),length(names_orange),length(names_brown),length(names_purple),length(names_cyan),length(my.color))
write.table(summary_table, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7cluster_summary_.txt', sep='\t', quote=F)

### GO ANALYSIS FOR GENES SELECTED BY THRESHOLD
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/summary_results_dominant_1.0_broadpeak_mrna_.txt',stringsAsFactors=F )
gene_cluster = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/gene_clusters.txt')
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
cut_off_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/Roadmap_cutoff_positions.txt')
groups = epigenomes$V1
for (group in groups){
  c_file_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/',group,'_clusters.txt',sep='')
  c_file = read.table(c_file_)
  i = 1
  cluster_ = as.character(c_file$V2[i])
  cut_off = cut_off_[cut_off_$V1==group, 2]
  temp = file_[file_$V2==group,]
  temp = temp[order(temp$V5, decreasing=T),]
  temp = temp[1:cut_off,]
  ref_genes = gene_cluster$V2[which(gene_cluster$V3==cluster_)]
  idx = which(temp$V1 %in% ref_genes)
  genelist = temp$V1[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file__ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/GO/',group,'_GO.txt', sep='')
  write.table(results, file__, sep='\t', quote=F)
}

groups = unique(as.character(table_$A))
groups = c(groups, 'purple')
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/Tissues/'
pathway = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/FET/clusters/GO/'
groups = c('E082','E107','E105','E112','E003','E020','E013','E062','E029','E049','E052','E058')
top_no = 5  # Number of GO terms to be included
terms = vector()
for (group in groups){
  filename_ = paste(pathway,group,'_GO.txt', sep='')
  file_ = read.table(filename_)
  for (i in 1:top_no){
    terms = c(terms, as.character(file_$V3[i]))
  }
}
terms= unique(terms)
table_ = data.frame(matrix(data=rep(1,length(groups)*length(terms)), ncol=length(groups), nrow=length(terms)))
colnames(table_) = groups
rownames(table_) = terms
for (no in 1:ncol(table_)){
  group = colnames(table_)[no]
  filename_ = paste(pathway,group,'_GO.txt', sep='')
  file_ = read.table(filename_)
  for (m in 1:nrow(table_)){
    t = rownames(table_)[m]
    idx = which(file_$V3==t)
    if (length(idx)!=0){
      table_[t, group] = file_[idx,4] 
    }
  }
}
# If want to change to the description of epigenomes, add this
po = read.table('/Users/woojunshim/Research/Data/Roadmap_IDs_.txt', stringsAsFactors = F)
for (i in 1:ncol(table_)){
  epi = colnames(table_)[i]
  colnames(table_)[i] = paste(gsub('_', ' ', po[which(po$V1==epi), 2]), '(',epi,')',sep='')
}

table_ = -log10(table_)
library(ggplot2)
library(reshape)
table_$tt = gsub('_', ' ', rownames(table_))
table_$tt <- factor(table_$tt, levels = table_$tt)
table = melt(table_)

p = ggplot(table, aes(y=tt, x=variable, fill=value)) 
p + geom_tile()+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
filename_ = paste(pathway,'Roadmap_GO_enrichment.png',sep='')
ggsave(filename_,width=10,height=10)

### GO analysis for Paige data 
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
cut_off_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/hclust/Roadmap_cutoff_positions.txt')
groups = c('day0','day2','day5','day9','day14')
for (group in groups){
  c_file_ = paste('/Users/woojunshim/Research/Data/Paige/1/',group,'_results.txt',sep='')
  c_file = read.table(c_file_, stringsAsFactors = F)
  cut_off = nrow(c_file)
  b_file = paste('/Users/woojunshim/Research/Data/Paige/1/h3k4me3_',group,'_assigned_.txt',sep='')
  temp = read.table(b_file, stringsAsFactors = F)
  temp = temp[order(temp$V3, decreasing = T),]
  genelist = temp$V2[1:cut_off]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file__ = paste('/Users/woojunshim/Research/Data/Paige/1/',group,'_GO_benayoun.txt', sep='')
  write.table(results, file__, sep='\t', quote=F)
  
  genelist = c_file$V1[1:cut_off]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file__ = paste('/Users/woojunshim/Research/Data/Paige/1/',group,'_GO_our.txt', sep='')
  write.table(results, file__, sep='\t', quote=F)
}

### ENTROPY-BASED APPROACH
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_specificity_table.txt', stringsAsFactors = F)
name = 'NFIA'
temp = file_[name,]
#temp = 1/temp
temp = as.numeric(temp[order(temp, decreasing=T)])
plot(y=as.numeric(temp), x=seq(1,length(temp)))
shapiro.test(temp)
cnt = 0
for (no in 1:nrow(file_)){
  temp = file_[no,]
  temp = as.numeric(temp[order(temp, decreasing=T)])
  t = shapiro.test(temp)
  if (t$p.value > 0.05){
    cnt = cnt + 1
  }
}

### CLUSTERING OF ENTROPY GENE TABLE
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/significant_cell_types_table_sorted1.txt')
data_ = file_
genes.cor = cor(t(data_), method='pearson')
#genes.cor.dist = dist(data_, method='euclidian')
genes.cor.dist = as.dist(1- genes.cor)
genes.tree = hclust(genes.cor.dist, method='ward.D')
times.cor = cor(data_, method='spearman')
#times.cor.dist = dist(t(data_), method='euclidian')
times.cor.dist = as.dist(1-times.cor)
times.tree = hclust(times.cor.dist, method='ward.D')
library('RColorBrewer')
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256)
col <- rev(col)
cut_ = cutree(genes.tree, k=8)
my.color = rep('green', times=nrow(data_))
my.color[cut_==2] = 'blue'
my.color[cut_==3] = 'red'
my.color[cut_==4] = 'orange'
my.color[cut_==5] = 'brown'
my.color[cut_==6] = 'purple'
my.color[cut_==7] = 'cyan'
my.color[cut_==8] = 'grey'
#my.color[cut_==9] = 'pink'
#my.color[cut_==10] = 'black'
#my.color
col_labels =vector()
for (no in 1:ncol(data_)){
  epi_name = colnames(data_)[no]
  idx = which(epigenomes$V1==epi_name)
  t = as.character(epigenomes$V2[idx])
  col_labels = c(col_labels, t)
}
pol = rainbow(length(unique((col_labels))))
col__ = pol[as.numeric(as.factor(col_labels))]
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/significant_selected_cell_types_heatmap_top5_final.pdf',width=10,height=10)
heatmap(as.matrix(data_), scale='none', col=c('white','red'), Rowv=as.dendrogram(genes.tree), Colv=as.dendrogram(times.tree), ColSideColors = col__, labRow = F)
dev.off()
file_$Colour = col__
names_green = names(cut_)[cut_=='1']
names_blue = names(cut_)[cut_=='2']
names_red = names(cut_)[cut_=='3']
names_orange = names(cut_)[cut_=='4']
names_brown = names(cut_)[cut_=='5']
names_purple = names(cut_)[cut_=='6']
names_cyan = names(cut_)[cut_=='7']
results = matrix(nrow=nrow(data_), ncol=2)
results[,1] = c(names_green,names_blue,names_red,names_orange,names_brown,names_purple,names_cyan)
results[,2] = c(rep('green',length(names_green)),rep('blue',length(names_blue)),rep('red',length(names_red)),rep('orange',length(names_orange)),rep('brown',length(names_brown)),rep('purple',length(names_purple)),rep('cyan',length(names_cyan)))
write.table(results, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7cluster_table_.txt', sep='\t', quote=F)
summary_table = matrix(nrow=length(unique(my.color))+1, ncol=2)
summary_table[,1] = c('green','blue','red','orange','brown','purple','cyan','total')
summary_table[,2] = c(length(names_green),length(names_blue),length(names_red),length(names_orange),length(names_brown),length(names_purple),length(names_cyan),length(my.color))
write.table(summary_table, '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/7cluster_summary_.txt', sep='\t', quote=F)


pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/Colour_legend_tissue_groups_final.pdf', width=10, height=10)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
cc = as.numeric(as.factor(col_labels))
legend('center', unique(as.character(epigenomes_$V2)), col=unique(col__), lty=1, lwd=2, cex=0.8)
dev.off()

### PEARSON'S CORRELATION BETWEEN GENES
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/significant_cell_types_table_sorted1.txt', stringsAsFactors = F)
file_ = t(file_)
cor_ = cor(file_)
cor_ = round(cor_, 4)
write.table(cor_, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/gene_pearson_short_sorted1.txt', quote=F, sep='\t')

file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_pearson_short.txt', stringsAsFactors = F)
gene = 'SKIDA1'
idx = which(file_[gene,])
mean(as.matrix(file_[gene,]))
hist(as.matrix(file_[gene,]))
file_[gene,'GATA4']
'TTLL10' %in% colnames(file_)[which(file_[gene,]>0.4)]

### EXTRACT ONLY TOP 5% GENES BY H3K4ME3 WIDTH TO PLOT HEAT MAP
gene_list = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/top5_genes_h3k4me3.txt', stringsAsFactors = F)
data_ = file_[which(rownames(file_) %in% gene_list$V1), ]
epigenomes_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/selected_epigenomes.txt', stringsAsFactors = F)
data_ = data_[, which(colnames(data_) %in% epigenomes_$V1)]
idx = vector()
for (no in 1:nrow(data_)){
  if (sum(data_[no,]) != 0){
    idx = c(idx, no)
  }
}
data_ = data_[idx,]

write.table(data_, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/selected_tissues_counts.txt', sep='\t', quote=F)

### CORRELATION PLOT FOR CARDIAC GENES 
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cardiac_genes_table.txt', stringsAsFactors = F)
table = t(data.frame(file_))
table_ = melt(table)
table_$Var2 = reorder(table_$Var2, unique(table_$Var2))
colnames(table_) = c('Gene1','Gene2','r')
ggplot(table_, aes(x=Gene2, y=Gene1)) + geom_tile(aes(fill=r)) + geom_text(aes(label=round(r,3))) + scale_fill_gradient2(midpoint=0, high="red",mid="white",low="blue") +xlab('') +ylab('')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cardiac_genes_table_cor.pdf', width=15, height=15)

### SUMMARY TABLE PLOT FOR CARDIAC GENES
genes_ = read.table('/Users/woojunshim/Research/Data/cardiac_regulators_modified_.txt', stringsAsFactors = F)
genes = genes_$V1
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cardiac_genes_E083.txt', stringsAsFactors = F)
file_ = file_[which(file_$V1 %in% genes),]
colnames(file_) = c('Gene','Mean(candidates)','Mean(background)','t.statistics','p-value')
file_ = file_[order(file_$Gene, decreasing=F),]
library(gridExtra)
library(grid)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/cardiac_genes_summary_E083.pdf',width=10, height=10)
grid.table(file_)
dev.off()

### CALCULATE PEARSON'S CORRELATION AND SUBTRACT DIFFERENCE FOR EACH GENE 
files = seq(18,111,by=3)
results = data.frame(matrix(ncol=(length(files)-1), nrow=16588))
colnames(results) = seq(1,(length(files)-1))
for (no in 1:(length(files)-1)){
  file1 = files[no]
  file2 = files[no+1]
  filename1 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/',file1,'_z_binary_.txt', sep='')
  filename2 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/',file2,'_z_binary_.txt', sep='')
  file1_ = read.table(filename1, stringsAsFactors = F)
  file2_ = read.table(filename2, stringsAsFactors = F)
  cor1 = cor(t(file1_))
  cor2 = cor(t(file2_))
  sum_ = cor2-cor1
  for (i in 1:nrow(cor1)){
    results[i,no] = length(which(abs(sum_[i,])>0.1))
  }
}
rownames(results) = rownames(cor1)
write.table(results, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/cor_difference_0.1.txt', sep='\t', quote=F)

for (no in 1:(length(files)-1)){
  file1 = files[no]
  file2 = files[no+1]
  filename1 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/',file1,'_z_binary_.txt', sep='')
  filename2 = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/',file2,'_z_binary_.txt', sep='')
  file1_ = read.table(filename1, stringsAsFactors = F)
  file2_ = read.table(filename2, stringsAsFactors = F)
  cor1 = cor(t(file1_))
  cor2 = cor(t(file2_))
  sum_ = cor2-cor1
  for (i in 1:nrow(cor1)){
    results[i,no] = length(which(abs(sum_[i,])>0.05))
  }
}
rownames(results) = rownames(cor1)
write.table(results, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/cor_difference_0.05.txt', sep='\t', quote=F)

### PLOT SATURATION CURVE 
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/cor_difference_0.1.txt', stringsAsFactors = F)
file_ = file_/nrow(file_)
file__ = apply(file_, 2, FUN=mean)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/saturation/saturation_analysis_0.1.pdf', width=10, height=10)
plot(x=seq(1,31), y=file__, type='o', col='dark blue', xlab='number of sets added', ylab='p (|r2-r1|>0.1)')
dev.off()

### PAIGE EXPRESSION DATA, using oligo
pathway = '/Users/woojunshim/Research/Data/Paige/expression/day2/'
setwd(pathway)
library(oligo)
library(pd.huex.1.0.st.v2)
library(hugene10sttranscriptcluster.db)
celfiles = list.celfiles()
affyraw1 = read.celfiles(celfiles[1])
affyraw2 = read.celfiles(celfiles[2])
affyraw = read.celfiles(celfiles)
eset = rma(affyraw, target='probeset')
aa = getNetAffx(eset, 'probeset')
table_ = data.frame(matrix(nrow=length(aa@data$geneassignment), ncol=2))
for (no in 1:length(aa@data$geneassignment)){
  table_[no,1] = aa@data$probesetid[no]
  t = aa@data$geneassignment[no]
  t = strsplit(t, split='//')
  table_[no,2] = t[[1]][2]
}
write.table(table_, '/Users/woojunshim/Research/Data/Paige/expression/affy_annotation_table_.txt', sep='\t', quote=F)

write.exprs(eset,file="expression_data_h7_esc.txt")

my_frame = data.frame(exprs(eset))
x = hugene10sttranscriptclusterSYMBOL
mapped_probes = mappedkeys(x)
xx = as.list(x[mapped_probes])
yy = data.frame(xx)
annot = matrix(ncol=2, nrow=length(yy))
annot[,1] = gsub('X','', names(yy))
for (i in 1:length(yy)){
  annot[i,2] = as.character(yy[[i]])
}
write.table(annot, '/Users/woojunshim/Research/Data/Paige/expression/affy_gene_annotation.txt', sep='\t', quote=F)

# DIFFERENTIAL EXPRESSION ANALYSIS 
library(DEGseq)
file_c = read.table('/Users/woojunshim/Research/Data/Paige/expression/h7_esc/expression_data_h7_esc_id.txt', stringsAsFactors = F)
file_i = read.table('/Users/woojunshim/Research/Data/Paige/expression/expression_data_id.txt', stringsAsFactors = F)
#write.table(file_, '/Users/woojunshim/Research/Data/Paige/expression/h7_esc/expression_data_h7_esc_id.txt', sep='\t', quote=F)
#matrix1 = readGeneExp(file='/Users/woojunshim/Research/Data/Paige/expression/expression_data_id.txt', geneCol=9, valCol=c(1:2))
#matrix2 = readGeneExp(file='/Users/woojunshim/Research/Data/Paige/expression/expression_data_id.txt', geneCol=9, valCol=c(3:8))
DEGexp(geneExpMatrix1 = file_i, geneCol1=9, expCol1=c(7:8), geneExpMatrix2 = file_i, geneCol2=9, expCol2=c(1:6), groupLabel1 = 'day2', groupLabel2 = 'rest', outputDir='/Users/woojunshim/Research/Data/Paige/expression/DEG/day2')
file_ = read.table('/Users/woojunshim/Research/Data/Paige/expression/h7_hsc_13133.txt', stringsAsFactors = F)

widths = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/All_widths.txt', stringsAsFactors = F)
spec = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_specificity_table.txt', stringsAsFactors = F)
background = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_specificity_.txt', stringsAsFactors = F)
back_width = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_widths_.txt', stringsAsFactors = F)
gene = 'OXT'
temp4 = back_width[gene,]
temp3 = background[gene,]
temp2 = spec[gene,]
temp1 = widths[gene,]
temp4 = temp4[order(temp4)]
temp3 = temp3[order(temp3)]
temp2 = temp2[order(temp2)]
temp1 = temp1[order(temp1)]
temp1 = temp1[which(temp1!=0)]
temp4 = temp4[which(temp4!=0)]
filename_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/', gene,'_width_.pdf', sep='')
y_range_ = range(temp1, temp4, na.rm=T)
pdf(filename_, width=10, height=10)
plot(x=seq(1,length(temp1)), y=temp1, ylim=y_range_, main=gene, xlab='Rank', ylab='H3K4me3 width (bp)', col='red')
points(x=seq(1,length(temp4)), y=temp4, col='dark blue')
legend("topleft", c('actual data','linear dynamics'), col=c('red','dark blue'), pch=c(1,1))
dev.off()

filename_ = paste('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/', gene,'_Q_.pdf', sep='')
y_range_ = range(temp2, temp3, na.rm=T)
pdf(filename_, width=10, height=10)
plot(x=seq(1,length(temp2)), y=temp2, ylim=y_range_, main=gene, xlab='Rank', ylab='Q', col='red')
points(x=seq(1,length(temp3)), y=temp3, col='dark blue')
legend("topleft", c('actual data','linear dynamics'), col=c('red','dark blue'), pch=c(1,1))
dev.off()



background = seq(min(temp2, na.rm=T), max(temp2, na.rm=T), length.out=length(which(!is.na(temp2))))
hist(as.numeric(temp2))
shapiro.test(as.numeric(temp2))

# DIFFERENTIAL EXPRESSION GENES FOR ROADMAP
library(DEGseq)
gg = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/filtering_specificity_combined.txt', stringsAsFactors = F)
epigenomes = rownames(gg)
#epigenomes = c(epigenomes, 'E000')
file_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
file_ = file_[1:nrow(file_),epigenomes]
file_$id = rownames(file_)
pathway = '/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/between_them/'
for (epi in epigenomes){
  col_idx = which(colnames(file_)==epi)
  other_idx = which(colnames(file_)!=epi)
  #other_idx = 47
  output__ = paste(pathway,epi,sep='')
  DEGexp(geneExpMatrix1 = file_, geneCol1=47, expCol1=col_idx, geneExpMatrix2 = file_, geneCol2=47, expCol2=other_idx[1:45], groupLabel1 = epi, groupLabel2 = 'others', outputDir=output__)
}

### DIFFERENTIALLY EXPRESSED GENES FOR PAIGE DATA
library(DEGseq)
gg = read.table('/Users/woojunshim/Research/Data/highly_expressed_genes/filtering_specificity_combined.txt', stringsAsFactors = F)
epigenomes = rownames(gg)
file_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/Paige/expression/expression_data_ave_symbol.txt', stringsAsFactors = F)
file_ = file_[1:nrow(file_),epigenomes]
file_$id = rownames(file_)
pathway = '/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/'
for (epi in epigenomes){
  col_idx = which(colnames(file_)==epi)
  other_idx = which(colnames(file_)!=epi)
  output__ = paste(pathway,epi,sep='')
  DEGexp(geneExpMatrix1 = file_, geneCol1=47, expCol1=col_idx, geneExpMatrix2 = file_, geneCol2=47, expCol2=other_idx[1:45], groupLabel1 = epi, groupLabel2 = 'rest', outputDir=output__)
}

### HISTOGRAM FOR ENTROPY
library(ggplot2)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_entropy.txt', stringsAsFactors = F)
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf = tf_list$V1
file_$tf = ifelse(file_[,1] %in% tf, 1, 0)
ggplot(file_, aes(x=V2)) + geom_histogram(data=subset(file_, tf=='1'), binwidth=bwidth, fill='red', alpha=0.2) + geom_histogram(data=subset(file_, tf=='0'), binwidth=bwidth, fill='blue', alpha=0.2) + geom_histogram(data=file_[,1:2], binwidth=bwidth, fill='green', alpha=0.2) + xlab('Entropy')
length(which(file_[,2]>6.7))
data_ = as.numeric(file_[,2])
breaks = pretty(range(data_), n = nclass.FD(data_), min.n = 1)
bwidth = breaks[2] - breaks[1]
hist(as.numeric(file_[,2]), breaks=50)
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/entropy_dist.png', width=10, height=10)

### HISTOGRAM FOR CELL TYPE COUNTS
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/cell_type_counts.txt', stringsAsFactors = F)
hist(as.numeric(file_[,2]), breaks=111, xlab='No.cell types',ylab='No.genes', main='')
file_[file_$V1=='OXT',]
length(which(file_$V2==5))

### CHECK 'BACKGROUND_SPECIFICITY.TXT' FILE
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_specificity_.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_specificity_table.txt', stringsAsFactors = F)
count_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/cell_type_counts.txt', stringsAsFactors = F)
genes = count_[which(count_[,2]==88),1]
temp = file_[which(rownames(file_) %in% genes),]
min_v1 = vector()
max_v1 = vector()
for (no in 1:nrow(temp)){
  tt = range(temp[no,], na.rm=T)
  min_v1 = c(min_v1, tt[1])
  max_v1 = c(max_v1, tt[2])
}
hist(min_v_)
hist(max_v_)
r101_min = min_v
r101_max = max_v
t = temp[1, order(temp[1,])]

### CHECKING ALL_WIDTHS.TXT AND BACKGROUND_WIDTHS.TXT FILES
file1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/All_widths.txt', stringsAsFactors = F)
file2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_widths_.txt', stringsAsFactors = F)
temp1 = file1['UBE2Q1',]
temp2 = file2['UBE2Q1',]

### DRAW WIDTHS & SPECIFICITY DISTRIBUTIONS 
file1 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/All_widths.txt', stringsAsFactors = F)
file2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_widths_.txt', stringsAsFactors = F)
temp1 = file1['GRIN1',]
temp2 = file2['HAND2',]
temp1 = temp1[order(temp1, decreasing=F)]
temp2 = temp2[order(temp2, decreasing=F)]
temp1 = temp1[which(temp1 != 0)]
temp2 = temp2[which(temp2 != 0)]
range(temp1)
range(temp2)

file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_pearson_short_new.txt', stringsAsFactors = F)
gene = 'TNNI3'
epi_ ='E065'
file_[gene,'MYH6']
hist(as.numeric(file_[gene,]))

file2 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/specificity_p_table.txt', stringsAsFactors = F)
file2[gene, epi_]

file3 = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/significant_cell_types_table_.txt', stringsAsFactors = F)
file3[gene, epi_]

### PLOT BOX PLOTS FOR Q VALUES (111 CELL TYPES)
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/cell_type_counts.txt', stringsAsFactors = F)
genes = file_[which(file_$V2==111),1]

file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/background_specificity_.txt', stringsAsFactors = F)
temp = file_[which(rownames(file_) %in% genes),]
temp1 = apply(temp, 1, sort)
spec = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/gene_specificity_table.txt', stringsAsFactors = F)
spec1 = spec[which(rownames(spec) %in% c('ATP5A1','GTF2B','GATA6','HAND2')),]
spec1 = apply(spec1, 1, sort)
y_range_ = range(spec1, temp1)
pdf('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cell_counts_111.pdf', width=10, height=5)
boxplot(temp1, use.cols = F, xlab='Rank', ylab='Q', ylim=y_range_)
cols_ = rainbow(ncol(spec1))
xx = seq(1,nrow(spec1))
for (i in 1:ncol(spec1)){
  yy = spec1[,i]
  lines(xx,yy,col=cols_[i], lwd=2)
}
legend(x=1, y=22, colnames(spec1), col = cols_, lty=1, lwd=2, cex=0.8, bty='n')
dev.off()

### PLOT 
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cardiac_genes_table.txt', stringsAsFactors = F)
table = t(data.frame(file_))
table_ = melt(table)
table_$Var2 = reorder(table_$Var2, unique(table_$Var2))
colnames(table_) = c('Gene1','Gene2','r')
ggplot(table_, aes(x=Gene2, y=Gene1)) + geom_tile(aes(fill=r)) + geom_text(aes(label=round(r,3))) + scale_fill_gradient2(midpoint=0, high="red",mid="white",low="blue") +xlab('') +ylab('')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/example_cardiac_genes_table_cor.pdf', width=15, height=15)

### ANALYSE GENES WITH LESS THAN 111 CELL TYPES
cell_counts = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/cell_type_counts.txt', stringsAsFactors = F)
widths = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/All_widths.txt', stringsAsFactors = F)
genes = cell_counts[which(cell_counts$V2==100),1]
temp = widths[which(rownames(widths) %in% genes),]
temp = apply(temp, 1, sort)
boxplot(temp[,1:20], use.cols = T)
colnames(temp)[8]

length(which(cell_counts$V2<30))

cnt = vector()
for (no in 1:111)
  
### CUMULATIVE PLOT FOR CELL TYPE COUNTS
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/cell_type_counts.txt', stringsAsFactors = F)
counts = vector()
for (no in 1:111){
  no = length(which(file_$V2<=no))/nrow(file_)
  counts = c(counts, no)
}
plot(x=seq(1,111), y=counts, type='l', col='dark blue', ylab='cumulative proportion', xlab='no.cell types')

### PLOT FOR CARDIAC VS BRAIN 
file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/cardiac_genes_one_side_p_nontf.txt', stringsAsFactors = F)
file_ = -log10(file_)
table = t(data.frame(file_))
table_ = melt(table)
y=factor(table_$Var1, levels=unique(c('E065','E095','E104','E105','E070','E071','E082')))
table_$Var1 = factor(table_$Var1, levels=unique(c('E065','E095','E104','E105','E070','E071','E082')))
colnames(table_) = c('cell_type','gene','value')
ggplot(table_, aes(x=cell_type, y=gene)) + geom_tile(aes(fill=value)) + geom_text(aes(label=round(value,3))) + scale_fill_gradient2(midpoint=0, high="red",mid="white",low="blue") +xlab('') +ylab('') + labs(fill='-log10(p-value)')
ggsave('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/cardiac_genes_one_side_p_nontf.pdf', width=10, height=10)

### CALCULATE MEANS OF EXPRESSION DATA (PAIGE)
file_ = read.table('/Users/woojunshim/Research/Data/Paige/expression/expression_data.txt', stringsAsFactors = F)
table_ = matrix(nrow=nrow(file_), ncol=ncol(file_)/2)
colnames(table_) = c('day5','day9','day14','day2')
rownames(table_) = rownames(file_)
table_[,1] = (file_[,1]+file_[,2]) /2
table_[,2] = (file_[,3]+file_[,4]) /2
table_[,3] = (file_[,5]+file_[,6]) /2
table_[,4] = (file_[,7]+file_[,8]) /2
write.table(table_, '/Users/woojunshim/Research/Data/Paige/expression/expression_data_ave.txt',quote=F, sep='\t')

### CALCULATE PEARSON'S R 

file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3/All_widths_H3K4me3.txt', stringsAsFactors = F)
width_cor = round(cor(t(file_)), 4)
write.table(width_cor, '/Users/woojunshim/Research/Data/broadPeaks/H3K4me3/width_pearson_H3K4me3.txt',sep='\t',quote=F)

file_1 = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/All_widths_H3K27me3.txt', stringsAsFactors = F)
width_cor1 = round(cor(t(file_1)), 4)
write.table(width_cor1, '/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/width_pearson_H3K27me3.txt',sep='\t',quote=F)

gene = 'GATA4'
temp = width_cor1[,gene]
idx = order(temp, decreasing = T)
temp = temp[idx]
temp[1:100]

temp = width_cor[,gene]
idx = order(temp, decreasing = T)
temp = temp[idx]
temp[1:100]


g1 = 'ID2'
g2 = 'MYH6'
width_cor[g1,g2]
width_cor1[g1,g2]
key_regulators = c('TBX3','TBX20','GATA4','WNT2','TBX2','NKX2-5','HAND2','TBX5','ISL1','HAND1','FOXC1','FOXC2','WNT11','WNT5A','BMP2','SMARCD3','MEF2C','PDGFRA','GATA6','MESP2')
key_regulators = c(key_regulators, 'MEIS2','MEIS1','SALL1','TWIST1','JARID2','MYH6','MYH7','PAX6','ZIC1')
result_k4 = data.frame(matrix(nrow=length(key_regulators), ncol=length(key_regulators)))
result_k27 = data.frame(matrix(nrow=length(key_regulators), ncol=length(key_regulators)))
rownames(result_k4) = key_regulators
colnames(result_k4) = key_regulators
rownames(result_k27) = key_regulators
colnames(result_k27) = key_regulators
for (row in key_regulators){
  for (col in key_regulators){
    result_k4[row,col] = width_cor[row,col]
    result_k27[row,col] = width_cor1[row,col]
  }
}

write.table(result_k4, '/Users/woojunshim/Research/Data/broadPeaks/H3K4me3/example_cor_H3K4me3.txt', sep='\t', quote=F)
write.table(result_k27, '/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/example_cor_H3K27me3.txt', sep='\t', quote=F)

file_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
file_ = file_[,which(colnames(file_)!='E000')]
exp_cor = round(cor(t(file_)), 4)
write.table(exp_cor, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/exp_pearson_cor.txt', sep='\t', quote=F)
colnames_ = colnames(file_)

file_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/All_widths.txt', stringsAsFactors = F)
file_ = file_[, colnames_]
width_cor = round(cor(t(file_)), 4)
write.table(width_cor, '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/width_pearson_cor_56celltypes.txt',sep='\t',quote=F)

### Mean expression by clusters
file_ = read.table('/Users/woojunshim/Research/scRNA/MeanExpression_by_Clusters_Separate5Day.txt', stringsAsFactors = T)
new_ = file_[,2:ncol(file_)]
write.table(new_, '/Users/woojunshim/Research/scRNA/MeanExpression_by_Clusters_Separate5Day_.txt', quote=F, sep='\t')

### GO ANALYSIS FOR DAY 5 
library(topGO)
library(GO.db)
back_ = read.table('/Users/woojunshim/Research/Data/scRNA/day5_background.txt', stringsAsFactors = F)
# Create GO terms to gene ID mapping table
all_genes = back_$V1
groups = c('day5_c1', 'day5_c2', 'day5_c3', 'day5_c4')
types = c('BP','CC','MF')
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (no in 1:length(groups)){
  group = groups[no]
  file_ = paste('/Users/woojunshim/Research/Data/scRNA/',groups[no],'_.txt',sep='')
  file__ = read.table(file_, stringsAsFactors = F)
  idx = which(file__$log2FoldChange<0)
  genelist = rownames(file__)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  for (type in types){
    TopGOdata<-new("topGOdata",ontology=type,allGenes=geneList,geneSel=genelist,nodeSize=5,
                   annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
    resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
    aa = score(resultFisher)
    aa = p.adjust(aa, method='fdr')  # FDR
    aa = sort(aa, decreasing=FALSE)  # sort by FDR
    bb = as.character(go_table[names(aa),1])
    top_no = length(aa)
    results = matrix(ncol=3,nrow=top_no)  # result table 
    results[,1] = names(aa)[1:top_no]
    results[,2] = bb[1:top_no]
    results[,3] = as.numeric(aa)[1:top_no]
    colnames(results) = c('#GO_ID','GO_Term','FDR')
    file_ = paste('/Users/woojunshim/Research/Data/scRNA/GO/',group,'_GO_',type,'.txt', sep='')
    write.table(results, file_, sep='\t', quote=F)
  }
}

### scRNA-seq data
tt=readRDS('/Users/woojunshim/Research/Data/Exprs_DCVLnorm_unlog_minus1_pos_Day15.RDS')


### COMPARE H3K4ME3 VS. H3K27ME3
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf_list = tf_list$V1
col_idx = which(colnames(width_cor) %in% tf_list)
row_idx = which(rownames(width_cor) %in% tf_list)
file_h3k4me3 = width_cor[row_idx, col_idx]
col_idx = which(colnames(width_cor1) %in% tf_list)
row_idx = which(rownames(width_cor1) %in% tf_list)
file_h3k27me3 = width_cor1[row_idx, col_idx]

table_h3k4me3 = matrix(nrow=nrow(width_cor), ncol=ncol(width_cor))
table_h3k27me3 = matrix(nrow=nrow(width_cor1), ncol=ncol(width_cor1))
table_h3k4me3 = apply(width_cor, 2, sort, decreasing=T)
table_h3k4me3 = data.frame(table_h3k4me3)
rownames(table_h3k4me3) = rownames(width_cor)
colnames(table_h3k4me3) = colnames(width_cor)
idx = which(rownames(table_h3k4me3) %in% tf_list)
mean_tf_h3k4me3 = mean(as.matrix(table_h3k4me3[idx,2:101]))
mean_nontf_h3k4me3 = mean(as.matrix(table_h3k4me3[-idx,2:101]))

### GO ANALYSIS FOR GENES WITHIN THE ELBOW POINT
elbows = read.table('/Users/woojunshim/Research/Data/broadPeaks/elbow_points.txt', stringsAsFactors = F)
epigenomes = read.table('/Users/woojunshim/Research/Data/Roadmap_IDs_.txt', stringsAsFactors = F)
epigenomes = epigenomes$V1
pathway = '/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/'
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (epi in epigenomes){
  cutoff = elbows[epi, 2]
  file__ = paste(pathway, epi, '_H3K27me3_genes.txt', sep='')
  file_ = read.table(file__, stringsAsFactors = F)
  genelist = file_[1:cutoff+1,2]
  all_genes = file_[1:nrow(file_),2]
  top_no = 100 # top 100 GO terms by FDR
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology='BP',allGenes=geneList,geneSel=genelist,nodeSize=5,annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  top_no = length(aa)
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste(pathway, 'GO/', epi, '_H3K27me3_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### EXTRACT EXPRESSED GENES FROM A SINGLE CELL
### FOCUSING ON DATA FROM SINGLE CELL LEVEL 
data_ = readRDS()
data_ = readRDS('/Users/woojunshim/Research/scRNA/Exprs_DCVLnorm_unlog_minus1_pos_Day30.RDS') # or
data_ = read.table('/Users/woojunshim/Research/scRNA/Exprs_Day30_RPKM_cells.txt', stringsAsFactors = F)
genes = data_[,1]
temp = unlist(strsplit(rownames(data_), '_'))
names(genes) = temp[seq(1, length(temp), 2)]
genes = genes[order(genes, decreasing = T)]
length(which(genes!=0))
genes = genes[which(genes>0.0)]

# FILTER OUT NON-INFORMATIVE GENES
threshold_ = 0.01 #less than x% of cell populations

dim(data_)
counts = apply(data_, 1, function(x) length(which(x>2)))
no_cells = ncol(data_) * threshold_
genes = names(counts)[which(counts>no_cells)]
data_ = data_[genes,]

counts = apply(data_, 1, function(x) length(which(x>0)))
no_cells = ncol(data_) * (1 - threshold_)
genes = names(counts)[which(counts<no_cells)]
data_ = data_[genes,]

genes = unlist(strsplit(genes, '_'))[seq(1,length(genes)*2, 2)]

# COEFFICIENT OF VARIANCE
sd_ = apply(data_, 1, sd)
mean_ = apply(data_, 1, mean)
cv = sd_ / abs(mean_)
names(cv) = genes
cv = cv[order(cv, decreasing=T)]

# 

### MM10 ANALYSIS
### OR MM10 GENES
file_ = read.table('/Users/woojunshim/Research/Data/mm10/h3k4me3/genes/H3K4me3_widths.txt', stringsAsFactors = F)
genes = rownames(file_)

### SINGLUALR VALUE DECOMPOSITION
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/combined_table.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/All_widths_H3K27me3.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/mm10/h3k4me3/genes/H3K4me3_widths.txt', stringsAsFactors = F)
file_ = file_ + 1
file__ = read.table('/Users/woojunshim/Research/Data/broadPeaks/Test_Day30_C2.txt', stringsAsFactors = F)
file_ = file_[rownames(file__),]
sorted_ = file_[names(genes),]
sorted_ = file_[genes, ]
sorted_ = sorted_[complete.cases(sorted_),]
data_ = scale(sorted_, center=T, scale=T)
data_ = scale(file_, center=T, scale=T)
#data_ = sorted_
aa = svd(data_)
D = diag(aa$d)
plot(1:length(aa$d), aa$d, xlab='Eigen component', ylab='Eigen value')  # Eigenvalues 
plot(1:length(aa$d), aa$d**2/sum(aa$d)**2, xlab='Number of Eigen values', main='Scree plot', ylab='Variance (i) / Total variance',type='o', col='dark blue')
cor_u = cor(t(aa$u[,1:10]))
# Generate a random data set

pk = aa$d**2/sum(aa$d)**2
entropy_ = 0
for (i in 1:length(pk)){
  entropy_ = entropy_ + pk[i]*log(pk[i])
}
entropy_ = -(entropy_)/log(length(pk))  # Entropy of the data (Alter et al. 2000)
results = vector()
sum_ = sum(aa$d)**2
current = 0
for (i in 1:length(aa$d)){
  current = current + aa$d[i]**2
  results = c(results, current/sum_)
}
results = results / max(results)
plot(1:length(aa$d),results, xlab='Number of Eigen values', ylab='', ylim=c(0,1.0), main='Cumulative relative variance', type='o', col='dark blue')
0.7/111 # Everitt and Dunn method gives only 1 singular value to use
U = aa$u
rownames(U) = rownames(data_)
colnames(U) = colnames(data_)
write.table(D, '/Users/woojunshim/Research/Data/broadPeaks/svd/D_10.txt', sep='\t', quote=F)
V = aa$v
rownames(V) = rownames(data_)
colnames(V) = colnames(data_)
write.table(V, '/Users/woojunshim/Research/Data/broadPeaks/svd/V_10.txt', sep='\t', quote=F)


plot(seq(1,nrow(aa$u)),aa$u[,20])  #Eigen-celltype
plot(seq(1,111),aa$v[5,])  #Eigen-gene


# transcriptional responses of genes
# reconstructing genes using singular values
no = ncol(aa$u) # Number of singular values to use 
us = as.matrix(aa$u[,1:no])
vs = as.matrix(aa$v[,1:no])
ds = as.matrix(D[1:no, 1:no])
result = us %*% ds %*% t(vs)

rownames(result) = rownames(data_)
colnames(result) = colnames(data_)
rownames(us) = rownames(data_)
cor_10 = cor(t(result))

cor_111 = cor(t(data_))

# Calculate correlation of genes to eigengenes 
# <Rk | Gn> / <Gn | Gn> (Alter et al. 2010)
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf_list = tf_list$V1
tf_list = c('NKX2-5','TBX5','HAND1')
tf_list = c('TNNI3','MYH6','MYH7')
eg_g = t(vs) %*% t(result[,1:no])
g_g = result[,1:no] %*% t(result[,1:no])
for (col in colnames(eg_g)){
  eg_g[,col] = eg_g[,col] / g_g[col,col]
}
write.table(eg_g, '/Users/woojunshim/Research/Data/broadPeaks/svd/relative_cor_eg_g_d30_c2_cell1.txt', sep='\t', quote=F)
idx = which(colnames(eg_g) %in% tf_list)
col_ = rep('dark blue', ncol(eg_g))
col_[idx] = 'red'
pch_ = rep(1, ncol(eg_g))
pch_[idx] = 19
plot(eg_g[1,],eg_g[2,], col=col_)
positive = colnames(eg_g)[which(eg_g[1,]>0)]
negative = colnames(eg_g)[which(eg_g[1,]<=0)]
write.table(positive, '/Users/woojunshim/Research/Data/broadPeaks/svd/positive.txt', sep='\n', quote=F)
cor_genes = round(cor(eg_g[1:5,]), 4)
write.table(cor_genes, '~/mm10_cor.txt', sep='\t', quote=F)

plot(seq(1,10), eg_g[1:10,'NKX2-5'], col='red')
points(seq(1,10), eg_g[1:10,'GATA4'], col='green')
points(seq(1,10), eg_g[1:10,'SDF4'], col='purple')
dist_ = as.matrix(dist(t(eg_g[1:10,])))
rownames(dist_) = rownames(data_)
colnames(dist_) = rownames(data_)


rownames(us) = rownames(data_)
plot(seq(1,10), us['NKX2-5',1:10], col='red')
points(seq(1,10), us['GATA4',1:10], col='green')
points(seq(1,10), us['SDF4',1:10], col='purple')
heatmap.2(as.matrix(us), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, labRow='')

dist_ = as.matrix(dist(us[,1:10]))
rownames(dist_) = rownames(data_)
colnames(dist_) = rownames(data_)

library(gplots)
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_genes_d30_c2_cell1_h3k27me3.pdf', width=10, height=10)
pdf('~/mm10_eg_g.pdf', width=10, height=10)
heatmap.2(as.matrix(t(eg_g)), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, labRow='')
dev.off()

# Same can be applied for cell types to eigen-cell types
ec_c = t(us) %*% result[,1:no]
c_c = t(result[,1:no]) %*% result[,1:no]
for (col in colnames(ec_c)){
  ec_c[,col] = ec_c[,col] / c_c[col,col]
}
write.table(ec_c, '/Users/woojunshim/Research/Data/broadPeaks/svd/relative_cor_ec_c_d30_c2_cell1.txt', sep='\t', quote=F)
plot(ec_c[1,],ec_c[2,])
cor_celltypes = round(cor(ec_c[1:5,]), 4)
write.table(cor_celltypes, '/Users/woojunshim/Research/Data/broadPeaks/svd/cor_celltypes_d30_c2_cell1_10.txt', sep='\t', quote=F)

pdf('~/mm10_ec_c.pdf', width=10, height=10)
heatmap.2(as.matrix(ec_c), col=greenred(75), scale='none', density.info ='none', trace='none', Rowv=F, labRow='')
dev.off()

# Plot for selected genes 
selected = c('PLN','MYH7','MYH6','NKX2-5','HAND1','MEF2C')
temp = eg_g[,selected]
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_selected_genes_d30_c2_cell1_.pdf', width=10, height=10)
heatmap.2(as.matrix(t(temp[1:10,])), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, cexRow = 0.8)
dev.off()

# Plot for selected cell-types 
selected = c('E095','E065','E083','E104','E105','E070')
temp = ec_c[,selected]
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_selected_celltypes_d30_c2_cell1_.pdf', width=10, height=10)
heatmap.2(as.matrix(t(temp[1:10,])), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, cexRow = 0.8)
dev.off()

# Characteristic modes (Holter et al. 2010)
r = 1  # the first r rows
mm = ds[r,] %*% t(vs[r,])
plot(seq(1,nrow(ds)), mm[r,])

mm = ds %*% t(vs)
heatmap.2(as.matrix(t(mm)), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, labRow='')

# characteristic mode coefficient (contribution)
tt = aa$u ** 2
rownames(tt) = rownames(data_)
write.table(tt, '/Users/woojunshim/Research/Data/broadPeaks/svd/aa_square_d30_c2_cell1.txt', sep='\t', quote=F)

# Visualisation
# 1. U (genes vs eigen-cell types)
library('gplots')
tt = aa$u %*% D
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/U_D30_C2_cell1_h3k27me3.pdf', width=10, height=10)
heatmap.2(as.matrix(tt), col=greenred(75), scale='row', density.info ='none', trace='none', Colv=F, labRow='')
dev.off()

# 2. V (eigengenes vs cell-types)
tt = aa$v %*% D
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/V_D30_C2_cell1_.pdf', width=10, height=10)
heatmap.2(as.matrix(t(tt)), col=greenred(75), scale='row', density.info ='none', trace='none', Rowv=F, labCol=colnames(data_), cexCol = 0.5)
dev.off()

# 3. line plots
plot(seq(1, length()))

# projection scatter plot
# first eigengene
idx = 1# highly correlated gene indexes (e.g. GATA4, NKX2-5 etc.)
tt = result[idx,]
coor1 = tt %*% aa$v[1,]
coor2 = tt %*% aa$v[2,]
plot(coor1, coor2)


# eigen cell type analysis
no = 10
us = as.matrix(aa$u[,1:no])
vs = as.matrix(aa$v[,1:no])
ds = as.matrix(D[1:no, 1:no])
result_ = vs %*% ds %*% t(us)

library(corrplot)
rownames(result_) = colnames(data_)
cor_cell_types = round(cor(t(result_)), digits=2)
corrplot(cor_cell_types, shade.col=NA, tl.col='black', tl.srt=45)

# relative regulatory contributions of eigengene values
values = apply(aa$v, 1, sum)
values = (values)**2 / (sum(aa$v))**2  # Alter et al. (2000)
plot(seq(1,nrow(aa$v)), values)

plot(seq(1,length(aa$d)), (aa$d**2)/(sum(aa$d)**2))

# exploring eigengene values
plot(seq(1,nrow(aa$v)), aa$u[56,], type='o')
plot(seq(1,nrow(aa$u)), aa$u[,10])

# Correlation between genes and eigen-genes (or cell types and eigen-cell types)
gene = 'NKX2-5'
idx = which(rownames(data_)==gene)
test_ = cor(x=data_[idx,], y=t(aa$v))
eigengene = order(test_, decreasing = T)

gene = 'HAND2'
idx = which(rownames(data_)==gene)
test_ = cor(x=data_[idx,], y=t(aa$v))
eigengene1 = order(test_, decreasing = T)

gene = 'PAX6'
idx = which(rownames(data_)==gene)
test_ = cor(x=data_[idx,], y=t(aa$v))
eigengene2 = order(test_, decreasing = T)

gene = 'MYH6'
idx = which(rownames(data_)==gene)
test_ = cor(x=data_[idx,], y=t(aa$v))
eigengene3 = order(test_, decreasing = T)

plot(seq(1,length(eigengene)), eigengene, col='red')
points(seq(1,length(eigengene1)), eigengene1, col='blue')
points(seq(1,length(eigengene2)), eigengene2, col='green')
points(seq(1,length(eigengene3)), eigengene3, col='orange')

results = data.frame(matrix(ncol=ncol(data_), nrow=nrow(data_)))
rownames(results) = rownames(data_)
row_ = rownames(data_)
for (row in row_){
  results[row,] = cor(x=data_[row,], y=t(aa$v))
}

results = results[,1:20]
final = cor(t(results))
write.table(final, '/Users/woojunshim/Research/Data/broadPeaks/svd/cor_eigengenes_20.txt',sep='\t',quote=F)

# characterising cell type profiles 
no = 10
j = 1
us = as.matrix(aa$u[,1:no])
vs = as.matrix(aa$v[j,1:no])
ds = as.matrix(D[1:no, 1:no])
result2 = t(vs) %*% ds %*% us

idx = which(rownames(file_)=='NKX2-5')
epi_ = which(colnames(file_)=='E014')
result[idx, ]
data_[idx,]
result[idx, epi_]
data_[idx, epi_]

aa$u %*% D %*% t(aa$v)  # X = UDV'

# PCA
epigenomes = colnames(file_)
tissue_ = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/epigenomes_list_.txt', stringsAsFactors = F)
rownames(tissue_) = tissue_[,1]
tissues = as.numeric(as.factor(tissue_[epigenomes, 2]))
mycol_ = rainbow(length(unique(tissues)))
mycol = mycol_[tissues]
data_ = t(file_)
data_table = prcomp(data_, center=T, scale.=T)
plot(data_table$x[,1], data_table$x[,2], col=tissues)

### CREATE non-overlapping exon length file for genes
library(GenomicFeatures)
data_ = makeTxDbFromGFF('/Users/woojunshim/Research/Data/Homo_sapiens.GRCh38.90.gtf', format='gtf')
exons.list.per.gene <- exonsBy(data_,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
results = matrix(nrow=length(exonic.gene.sizes), ncol=2)
names_ = names(exonic.gene.sizes)
values_ = as.numeric(exonic.gene.sizes)
results[,1] = names_
results[,2] = values_
write.table(results, '/Users/woojunshim/Research/Data/gene_lengths_ensembl_GRCh38.txt', sep='\t', quote=F)

### SINGLE CELL RDSDATA
### CALCULATE RPKM 
exp_ = readRDS('/Users/woojunshim/Research/scRNA/Exprs_DCVLnorm_unlog_minus1_pos_Day30.RDS')
names_ = readRDS('/Users/woojunshim/Research/scRNA/CellNames_ClusterIDs_day0.RDS')
length_ = read.table('/Users/woojunshim/Research/Data/gene_lengths_ensembl_GRCh38.txt')
library('edgeR')
genes = rownames(exp_)
yy= vector()
for (item in genes){
  temp = unlist(strsplit(item, '_'))
  yy = c(yy, temp[2])
}
names = yy
lengths = vector()
idx = which(names %in% length_$V1)
exp_ = exp_[idx,]
idx = match(names, length_$V1)
values = length_$V2[idx]
results = rpkm(exp_, values)
saveRDS(results, '/Users/woojunshim/Research/scRNA/Exprs_Day30_RPKM.RDS')
write.table(results[,1:100], '/Users/woojunshim/Research/scRNA/Exprs_Day30_RPKM_cells.txt', sep='\t', quote=F)

### Test function
data_ = read.table('/Users/woojunshim/Research/Data/combined_results_roadmap.txt', stringsAsFactors = F, row.names=NULL)
temp = t(data_[,c('E095','E104','E105','E065','E070','E071','E068','E069','E072','E013','E012','E011','E004','E005')])
temp$group_ = rownames(temp)
pca_ = prcomp(as.matrix(data_))
library(ggfortify)
autoplot(prcomp(temp), label = T)

### Correlation 
u = aa$u
rownames(u) = rownames(data_)
cor_u = cor(t(aa$u[,1:10]))
rownames(cor_u) = rownames(data_)
colnames(cor_u) = rownames(data_)

cor_eg = cor(eg_g[1:10,])
ss = apply(cor_eg, 1, sum)
names(ss) = rownames(cor_eg)
ss = ss[order(ss, decreasing = T)]
var_ = apply(cor_eg, 1, var)
names(var_) = rownames(cor_eg)

### Drawing dot plots between two eigen components
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/genes_two_ecomponents_h3k27me3_cor_selected.pdf')
for (i in 1:10){
  for (j in i:10){
    if (i!=j){
      plot(eg_g[i,],eg_g[j,], xlab=i, ylab=j, col=col_, pch=pch_)
      points(eg_g[i,idx], eg_g[j,idx], col='red', pch=19)
      points(eg_g[i,idx2], eg_g[j,idx2], col='green', pch=19)
    }
  }
}
dev.off()

### Plot selected genes 
selected1 = c('NKX2-5','HAND1','TBX5')
selected2 = c('MYH6','TNNI3','MYH7')
y_range = range(eg_g[1:10,c(selected1,selected2)])
plot(seq(1,10), eg_g[1:10,selected1[1]], col='red', type='o', ylim=y_range, xlab='Eigen component', ylab='Correlation coordinate')
lines(seq(1,10), eg_g[1:10,selected1[2]], col='red', type='o')
lines(seq(1,10), eg_g[1:10,selected1[3]], col='red', type='o')
lines(seq(1,10), eg_g[1:10,selected2[1]], col='dark blue', type='o')
lines(seq(1,10), eg_g[1:10,selected2[2]], col='dark blue', type='o')
lines(seq(1,10), eg_g[1:10,selected2[3]], col='dark blue', type='o')
lines(seq(1,10), eg_g[1:10,selected2[4]], col='dark blue', type='o')

### SATURATION PLOT FOR FUNCTION TESTING
combined_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/saturation_combined_log10.txt', stringsAsFactors = F)
combined_ = combined_[,2:ncol(combined_)]
combined_ = apply(combined_, 2, mean)
h3k27me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/saturation_H3K27me3_rank.txt', stringsAsFactors = F)
h3k27me3_ = h3k27me3_[,2:ncol(h3k27me3_)]
h3k27me3_ = apply(h3k27me3_, 2, mean)
h3k4me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/saturation_H3K4me3_rank.txt', stringsAsFactors = F)
h3k4me3_ = h3k4me3_[,2:ncol(h3k4me3_)]
h3k4me3_ = apply(h3k4me3_, 2, mean)
exp_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/saturation_expression.txt', stringsAsFactors = F)
exp_ = exp_[,2:ncol(exp_)]
exp_ = apply(exp_, 2, mean)

random_ = seq(0.01, 1.0, 0.01)
col_ = rainbow(4)
range_ = range(0,1.0)
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/comparison_functions.pdf', width=5, height=5)
plot(seq(1,100), combined_, type='l', col=col_[1], ylim=range_, xlab='Top rank (%)', ylab='Ave.proportion')
lines(seq(1,100), h3k27me3_, type='l', col=col_[2])
lines(seq(1,100), h3k4me3_, type='l', col=col_[3])
lines(seq(1,100), exp_, type='l', col=col_[4])
lines(seq(1,100), random_, type='l', lty=2)
legend(x=0, y=1.0, legend=c('Combined','H3K27me3','H3K4me3','RPKM','Random'), col=c(col_,'black'), lty=c(1,1,1,1,2), cex=0.6)
dev.off()

### PLOT JUST USING LEFT EIGEN VECTORS FOR GENES
eg_g = t(us)
idx = which(colnames(eg_g) %in% tf_list)
col_ = rep('dark blue', ncol(eg_g))
col_[idx] = 'red'
plot(eg_g[2,],eg_g[3,], col=col_)
positive = colnames(eg_g)[which(eg_g[1,]>0)]
negative = colnames(eg_g)[which(eg_g[1,]<=0)]
write.table(positive, '/Users/woojunshim/Research/Data/broadPeaks/svd/positive.txt', sep='\n', quote=F)
cor_genes = round(cor(eg_g[1:10,]), 4)

### CHECKING CLUSTERS FROM COR_GENES 
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/svd/d30_c2_cell1_10.txt', stringsAsFactors = F)
gene = 'IRX4'
idx = which(file_[gene,]==1)
colnames(file_[gene,idx])
idx1 = which(apply(file_,1,sum)!=0)
file_ = file_[idx1,]
library(gplots)
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/d30_c2_cell1_10.pdf', width=10, height=10)
heatmap.2(as.matrix(file_), col=hmcol, trace='none', labCol='', labRow='')
dev.off()

### ASSESSING VARIANCE OF GENES (EXPRESSION) AT A TIME POINT
var_ = apply(file_,1,var)

### CORRELATION OF SELETED GENES
tf_ = c('TBX5','GATA4','NKX2-5','HAND2')
structural_ = c('TNNI3','MYH6','MYH7')
house_ = c('ACTB','GAPDH','HPRT1','B2M')
exp_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
epigenomes = colnames(exp_)
house_ %in% rownames(h3k27me3_)
h3k4me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3/All_widths_H3K4me3.txt', stringsAsFactors = F)
h3k27me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/All_widths_H3K27me3.txt', stringsAsFactors = F)
k4k27_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/combined_table.txt', stringsAsFactors = F)
combined_ = c(tf_, structural_, house_)
cor_exp = cor(t(exp_[combined_,]))
cor_h3k4me3 = cor(t(h3k4me3_[combined_,epigenomes]))
cor_h3k27me3 = cor(t(h3k27me3_[combined_,epigenomes]))
cor_k4k27 = cor(t(k4k27_[combined_,epigenomes]))

### PLOT VS EXPRESSION VALUES
h3k4me1_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me1_expression_mean.txt', stringsAsFactors = F)
h3k4me1_ = h3k4me1_[,2:ncol(h3k4me1_)]
h3k4me1_ = apply(h3k4me1_, 2, mean)
h3k9me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K9me3_expression_mean.txt', stringsAsFactors = F)
h3k9me3_ = h3k9me3_[,2:ncol(h3k9me3_)]
h3k9me3_ = apply(h3k9me3_, 2, mean)
h3k36me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K36me3_expression_mean.txt', stringsAsFactors = F)
h3k36me3_ = h3k36me3_[,2:ncol(h3k36me3_)]
h3k36me3_ = apply(h3k36me3_, 2, mean)
h3k27me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_expression_mean.txt', stringsAsFactors = F)
h3k27me3_ = h3k27me3_[,2:ncol(h3k27me3_)]
h3k27me3_ = apply(h3k27me3_, 2, mean)
h3k4me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3_expression_mean.txt', stringsAsFactors = F)
h3k4me3_ = h3k4me3_[,2:ncol(h3k4me3_)]
h3k4me3_ = apply(h3k4me3_, 2, mean)
h3k27ac_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27ac_expression_mean.txt', stringsAsFactors = F)
h3k27ac_ = h3k27ac_[,2:ncol(h3k27ac_)]
h3k27ac_ = apply(h3k27ac_, 2, mean)
exp_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/saturation_expression.txt', stringsAsFactors = F)
exp_ = exp_[,2:ncol(exp_)]
exp_ = apply(exp_, 2, mean)

wilcox.test(h3k27ac_[1:10], h3k4me3_[11:100], alternative = 'greater')

random_ = seq(0.01, 1.0, 0.01)
col_ = rainbow(6)
range_ = range(h3k4me1_, h3k4me3_, h3k27me3_, h3k36me3_, h3k9me3_,h3k27ac_) 
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/comparison_hms_exp__.pdf', width=5, height=5)
plot(seq(1,100), h3k4me1_, type='l', col=col_[1], ylim=range_, xlab='Percentile rank of the peak width', ylab='RPKM',xaxt='n')
axis(1, at=seq(0,100,20), labels=seq(100,0,-20), tick=T)
lines(seq(1,100), h3k27me3_, type='l', col=col_[2])
lines(seq(1,100), h3k4me3_, type='l', col=col_[3])
lines(seq(1,100), h3k36me3_, type='l', col=col_[4])
lines(seq(1,100), h3k9me3_, type='l', col=col_[5])
lines(seq(1,100), h3k27ac_, type='l', col=col_[6])
#lines(seq(1,100), exp_, type='l', col=col_[4])
#lines(seq(1,100), random_, type='l', lty=2)
legend(x=70, y=120, legend=c('H3K4me1','H3K27me3','H3K4me3','H3K36me3','H3K9me3','H3K27ac'), col=c(col_), lty=c(1,1,1,1,1,1), cex=0.6)
dev.off()


### PLOT VS EXPRESSION VALUES
### CUMULATIVE PLOT
h3k4me1_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me1_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k4me1_ = h3k4me1_[,2:ncol(h3k4me1_)]
h3k4me1_ = apply(h3k4me1_, 2, mean)
h3k9me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K9me3_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k9me3_ = h3k9me3_[,2:ncol(h3k9me3_)]
h3k9me3_ = apply(h3k9me3_, 2, mean)
h3k36me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K36me3_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k36me3_ = h3k36me3_[,2:ncol(h3k36me3_)]
h3k36me3_ = apply(h3k36me3_, 2, mean)
h3k27me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k27me3_ = h3k27me3_[,2:ncol(h3k27me3_)]
h3k27me3_ = apply(h3k27me3_, 2, mean)
h3k4me3_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k4me3_ = h3k4me3_[,2:ncol(h3k4me3_)]
h3k4me3_ = apply(h3k4me3_, 2, mean)
h3k27ac_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27ac_expression_mean_cumulative.txt', stringsAsFactors = F)
h3k27ac_ = h3k27ac_[,2:ncol(h3k27ac_)]
h3k27ac_ = apply(h3k27ac_, 2, mean)

cor(h3k36me3_, seq(100,1,-1))
cor.test(h3k36me3_, seq(100,1,-1))

col_ = rainbow(6)
range_ = range(h3k4me1_, h3k4me3_, h3k27me3_, h3k36me3_, h3k9me3_, h3k27ac_) 
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/comparison_hms_exp_cumulative.pdf', width=10, height=10)
plot(seq(1,100), h3k4me1_, type='l', col=col_[1], ylim=range_, xlab='Rank position (top %)', ylab='Ave.RPKM')
#axis(1, at=seq(1,100,10), labels=seq(100,1,-10), tick=T)
lines(seq(1,100), h3k27me3_, type='l', col=col_[2])
lines(seq(1,100), h3k4me3_, type='l', col=col_[3])
lines(seq(1,100), h3k36me3_, type='l', col=col_[4])
lines(seq(1,100), h3k9me3_, type='l', col=col_[5])
lines(seq(1,100), h3k27ac_, type='l', col=col_[6])
#lines(seq(1,100), random_, type='l', lty=2)
legend(x=80, y=120, legend=c('H3K4me1','H3K27me3','H3K4me3','H3K36me3','H3K9me3','H3K27ac'), col=c(col_), lty=c(1,1,1,1,1,1), cex=1.0)
dev.off()

file_ = rbind(h3k27ac_, h3k4me3_, h3k4me1_, h3k36me3_, h3k27me3_, h3k9me3_)
rownames(file_) = c('H3K27ac','H3K4me3','H3K4me1','H3K36me3','H3K27me3','H3K9me3')
colnames(file_) = seq(99,0,by=-1)
file__ = t(file_)
library(ggplot2)
library(reshape)
table = melt(file__)
range_ = range(h3k27ac_, h3k4me3_, h3k4me1_, h3k36me3_, h3k27me3_, h3k9me3_)
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(mid='grey70', limits=range_) + ggtitle('Odd ratios of expressed TFs (among all genes with RPKM>1.0)') + labs(x='Percentile Rank (RPKM)', y='Cell type') + scale_x_reverse(breaks = seq(0,95,by=5)) 

### Second histone mark ###
cor_genes1 = cor_genes

### SINGLUALR VALUE DECOMPOSITION
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/combined_table.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3/All_widths_H3K4me3.txt', stringsAsFactors = F)
file_ = file_ + 1
file__ = read.table('/Users/woojunshim/Research/Data/broadPeaks/Test_Day30_C2.txt', stringsAsFactors = F)
file_ = file_[rownames(file__),]
sorted_ = file_[names(genes),]
sorted_ = file_[genes, ]
sorted_ = sorted_[complete.cases(sorted_),]
data_ = scale(sorted_, center=T, scale=T)
data_ = scale(file_, center=T, scale=T)
#data_ = sorted_
aa = svd(data_)
D = diag(aa$d)
plot(1:length(aa$d), aa$d)  # Eigenvalues 
plot(1:length(aa$d), aa$d**2/sum(aa$d)**2, xlab='Number of singular values', main='Scree plot', type='o', col='dark blue')
cor_u = cor(t(aa$u[,1:10]))
# Generate a random data set

pk = aa$d**2/sum(aa$d)**2
entropy_ = 0
for (i in 1:length(pk)){
  entropy_ = entropy_ + pk[i]*log(pk[i])
}
entropy_ = -(entropy_)/log(length(pk))  # Entropy of the data (Alter et al. 2000)
results = vector()
sum_ = sum(aa$d)**2
current = 0
for (i in 1:length(aa$d)){
  current = current + aa$d[i]**2
  results = c(results, current/sum_)
}
results = results / max(results)
plot(1:length(aa$d),results, xlab='Number of singular values', ylab='', ylim=c(0,1.0), main='Cumulative relative variance', type='o', col='dark blue')
0.7/111 # Everitt and Dunn method gives only 1 singular value to use
U = aa$u
rownames(U) = rownames(data_)
colnames(U) = colnames(data_)
write.table(D, '/Users/woojunshim/Research/Data/broadPeaks/svd/D_10.txt', sep='\t', quote=F)
V = aa$v
rownames(V) = rownames(data_)
colnames(V) = colnames(data_)
write.table(V, '/Users/woojunshim/Research/Data/broadPeaks/svd/V_10.txt', sep='\t', quote=F)


plot(seq(1,nrow(aa$u)),aa$u[,20])  #Eigen-celltype
plot(seq(1,111),aa$v[5,])  #Eigen-gene


# transcriptional responses of genes
# reconstructing genes using singular values
no = 111 # Number of singular values to use 
us = as.matrix(aa$u[,1:no])
vs = as.matrix(aa$v[,1:no])
ds = as.matrix(D[1:no, 1:no])
result = us %*% ds %*% t(vs)

rownames(result) = rownames(data_)
colnames(result) = colnames(data_)
rownames(us) = rownames(data_)
cor_10 = cor(t(result))

cor_111 = cor(t(data_))

# Calculate correlation of genes to eigengenes 
# <Rk | Gn> / <Gn | Gn> (Alter et al. 2010)
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf_list = tf_list$V1
tf_list = c('NKX2-5','TBX5','HAND1')
tf_list = c('TNNI3','MYH6','MYH7')
eg_g = t(vs) %*% t(result[,1:no])
g_g = result[,1:no] %*% t(result[,1:no])
for (col in colnames(eg_g)){
  eg_g[,col] = eg_g[,col] / g_g[col,col]
}
write.table(eg_g, '/Users/woojunshim/Research/Data/broadPeaks/svd/relative_cor_eg_g_d30_c2_cell1.txt', sep='\t', quote=F)
idx = which(colnames(eg_g) %in% tf_list)
col_ = rep('dark blue', ncol(eg_g))
col_[idx] = 'red'
pch_ = rep(1, ncol(eg_g))
pch_[idx] = 19
plot(eg_g[1,],eg_g[2,], col=col_)
positive = colnames(eg_g)[which(eg_g[1,]>0)]
negative = colnames(eg_g)[which(eg_g[1,]<=0)]
write.table(positive, '/Users/woojunshim/Research/Data/broadPeaks/svd/positive.txt', sep='\n', quote=F)
cor_genes = round(cor(eg_g[1:10,]), 4)
write.table(cor_genes, '/Users/woojunshim/Research/Data/broadPeaks/svd/cor_genes_d30_c2_cell1_10.txt', sep='\t', quote=F)

library(gplots)
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_genes_d30_c2_cell1_h3k27me3.pdf', width=10, height=10)
heatmap.2(as.matrix(t(eg_g)), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, labRow='')
dev.off()

# Same can be applied for cell types to eigen-cell types
ec_c = t(us) %*% result[,1:no]
c_c = t(result[,1:no]) %*% result[,1:no]
for (col in colnames(ec_c)){
  ec_c[,col] = ec_c[,col] / c_c[col,col]
}
write.table(ec_c, '/Users/woojunshim/Research/Data/broadPeaks/svd/relative_cor_ec_c_d30_c2_cell1.txt', sep='\t', quote=F)
plot(ec_c[1,],ec_c[2,])
cor_celltypes = round(cor(ec_c[1:no,]), 4)
write.table(cor_celltypes, '/Users/woojunshim/Research/Data/broadPeaks/svd/cor_celltypes_d30_c2_cell1_10.txt', sep='\t', quote=F)

pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_celltypes_d30_c2_cell1_h3k27me3.pdf', width=10, height=10)
heatmap.2(as.matrix(ec_c), col=greenred(75), scale='none', density.info ='none', trace='none', Rowv=F, labRow='')
dev.off()

# Plot for selected genes 
selected = c('PLN','MYH7','MYH6','NKX2-5','HAND1','MEF2C')
temp = eg_g[,selected]
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_selected_genes_d30_c2_cell1_.pdf', width=10, height=10)
heatmap.2(as.matrix(t(temp[1:10,])), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, cexRow = 0.8)
dev.off()

# 
write.table(cor_genes, '/Users/woojunshim/Research/Data/broadPeaks/svd/D30_filtered_genes_H3K4me3_cor.txt',quote=F, sep='\t')
write.table(cor_genes1, '/Users/woojunshim/Research/Data/broadPeaks/svd/D30_filtered_genes_H3K27me3_cor.txt',quote=F, sep='\t')

#
tt = (cor_genes['NKX2-5','HAND2'] + cor_genes1['NKX2-5','HAND2']) / abs(cor_genes['NKX2-5','HAND2'] - cor_genes1['NKX2-5','HAND2'])
yy = (cor_genes['NKX2-5','MYH6'] + cor_genes1['NKX2-5','MYH6']) / abs(cor_genes['NKX2-5','MYH6'] - cor_genes1['NKX2-5','MYH6'])

gg = (cor_genes + cor_genes1) / (cor_genes-cor_genes1)

### combined correlation table for D30 
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/svd/D30_filtered_genes_combined_binary.txt', stringsAsFactors = F)

### BOXPLOTS FOR TOP5% GENES ACROSS DIFFERENT HMS
type = 'All_expressed_TF'
pathway = '/Users/woojunshim/Research/Data/highly_expressed_genes/DEG/new/'
file1 = read.table(paste(pathway,'density_ratio_',type,'_H3K4me3_a.txt', sep=''), stringsAsFactors = F)
file2 = read.table(paste(pathway,'density_ratio_',type,'_H3K27me3_a.txt', sep=''), stringsAsFactors = F)
file3 = read.table(paste(pathway,'density_ratio_',type,'_H3K36me3_a.txt', sep=''), stringsAsFactors = F)
file4 = read.table(paste(pathway,'density_ratio_',type,'_H3K9me3_a.txt', sep=''), stringsAsFactors = F)
file5 = read.table(paste(pathway,'density_ratio_',type,'_H3K4me1_a.txt', sep=''), stringsAsFactors = F)
file6 = read.table(paste(pathway,'density_ratio_',type,'_H3K27ac_a.txt', sep=''), stringsAsFactors = F)
file1$group = rep('H3K4me3',nrow(file1))
file2$group = rep('H3K27me3',nrow(file2))
file3$group = rep('H3K36me3',nrow(file3))
file4$group = rep('H3K9me3',nrow(file4))
file5$group = rep('H3K4me1',nrow(file5))
file6$group = rep('H3K27ac',nrow(file6))
file1 = file1[,c(1,2,ncol(file1))]
file2 = file2[,c(1,2,ncol(file2))]
file3 = file3[,c(1,2,ncol(file3))]
file4 = file4[,c(1,2,ncol(file4))]
file5 = file5[,c(1,2,ncol(file5))]
file6 = file6[,c(1,2,ncol(file6))]
table_ = rbind(file1,file2,file3,file4,file5,file6)
table_$method = rep('broadest',nrow(table_))
file1_ = file1
file2_ = file2
file3_ = file3
file4_ = file4
file5_ = file5
file6_ = file6

file1 = read.table(paste(pathway,'density_ratio_',type,'_H3K4me3_sum_a.txt', sep=''), stringsAsFactors = F)
file2 = read.table(paste(pathway,'density_ratio_',type,'_H3K27me3_sum_a.txt', sep=''), stringsAsFactors = F)
file3 = read.table(paste(pathway,'density_ratio_',type,'_H3K36me3_sum_a.txt', sep=''), stringsAsFactors = F)
file4 = read.table(paste(pathway,'density_ratio_',type,'_H3K9me3_sum_a.txt', sep=''), stringsAsFactors = F)
file5 = read.table(paste(pathway,'density_ratio_',type,'_H3K4me1_sum_a.txt', sep=''), stringsAsFactors = F)
file6 = read.table(paste(pathway,'density_ratio_',type,'_H3K27ac_sum_a.txt', sep=''), stringsAsFactors = F)
file1$group = rep('H3K4me3',nrow(file1))
file2$group = rep('H3K27me3',nrow(file2))
file3$group = rep('H3K36me3',nrow(file3))
file4$group = rep('H3K9me3',nrow(file4))
file5$group = rep('H3K4me1',nrow(file5))
file6$group = rep('H3K27ac',nrow(file6))
file1 = file1[,c(1,2,ncol(file1))]
file2 = file2[,c(1,2,ncol(file2))]
file3 = file3[,c(1,2,ncol(file3))]
file4 = file4[,c(1,2,ncol(file4))]
file5 = file5[,c(1,2,ncol(file5))]
file6 = file6[,c(1,2,ncol(file6))]
table__ = rbind(file1,file2,file3,file4,file5,file6)
table__$method = rep('sum',nrow(table__))
table_ = rbind(table_, table__)

colnames(table_)[1:2] = c('cell type','value')
#table = melt(table_)
#table$value = -log10(table$value)
table_$lev = factor(table_$group, levels=c('H3K4me3','H3K36me3','H3K27me3','H3K9me3','H3K4me1','H3K27ac'))
#colnames(table) = c('cell type','group','method','value','lev')
p = ggplot(table_, aes(lev, value, fill=method)) 
p + geom_boxplot() + ggtitle('All expressed TFs (RPKM>1.0)') +  ylab('Enrichment score') + xlab('')
filename = paste(pathway,'hm_comparison_top5_',type,'.pdf',sep='')
ggsave(filename, width=10, height=10)

for (g in unique(table_$method)){
  temp = table_[which(table_$method==g),]
  results = data.frame(matrix(nrow=length(unique(table_$group)),ncol=length(unique(table_$group))))
  rownames(results) = unique(table_$group)
  colnames(results) = unique(table_$group)
  for (m1 in unique(table_$group)){
    for (m2 in unique(table_$group)){
      if (m1 != m2){
        yy=wilcox.test(temp[temp$group==m1,2], temp[temp$group==m2,2])
        results[m1,m2] = yy$p.value
      } else {results[m1,m2]=1.0}
    }
  }
  filename_ = paste(pathway,'stat_',type,'_',g,'_top5.txt', sep='')
  write.table(results, filename_, sep='\t', quote=F)
}

# GO analysis for top5% genes selected by HMs
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = mrna_$X.gene
type = 'H3K4me3'
types = c('H3K36me3','H3K27me3','H3K9me3','H3K4me1','H3K27ac')
types = c('H3K27ac')
groups = c('E095','E100','E070','E013','E032','E043','E003','E057')
groups = c('E095','E100','E071','E013','E032','E043','E003','E056')
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (type in types){
  file__ = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/',type,'_widths.txt', sep=''),stringsAsFactors = F)
  for (no in 1:length(groups)){
    group = groups[no]
    temp = file__[,group]
    names(temp) = rownames(file__)
    temp = temp[order(temp, decreasing=T)]  
    temp = temp[which(temp!=0)]
    all_genes = names(temp)
    cutoff = round(length(temp) * 0.05)
    genelist = names(temp[1:cutoff])
    geneList<-factor(as.integer(all_genes %in% genelist))
    names(geneList)<-all_genes
    TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                   annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
    resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
    aa = score(resultFisher)
    aa = p.adjust(aa, method='fdr')  # FDR
    aa = sort(aa, decreasing=FALSE)  # sort by FDR
    bb = as.character(go_table[names(aa),1])
    results = matrix(ncol=3,nrow=top_no)  # result table 
    results[,1] = names(aa)[1:top_no]
    results[,2] = bb[1:top_no]
    results[,3] = as.numeric(aa)[1:top_no]
    colnames(results) = c('#GO_ID','GO_Term','FDR')
    file_ = paste('/Users/woojunshim/Research/Data/broadPeaks/GO/',group,'_',type,'_broadest_top5.txt', sep='')
    write.table(results, file_, sep='\t', quote=F)
  }
}


for (no in 1:length(groups)){
  group = groups[no]
  temp = file__[,group]
  names(temp) = rownames(file__)
  temp = temp[order(temp, decreasing=T)]  
  temp = temp[which(temp!=0)]
  all_genes = names(temp)
  cutoff = round(length(temp) * 0.05)
  genelist = names(temp[1:cutoff])
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/broadPeaks/GO/',group,'_',type,'_broadest_top5.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### HEATMAP FOR SIGNIFICANT GO TERMS ACROSS DEFINED GROUPS
### RUN THIS FOR PLOTTING ANY PREDEFINED GROUPS 

groups = colnames(file__)
groups = c('E095','E100','E070','E013','E032','E043','E003','E057')
groups = c('E095','E100','E071','E013','E032','E043','E003','E056')
pathway_ = '/Users/woojunshim/Research/Data/Paige/1/analysis/clustering/hclust/7clusters/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/Tissues/'
pathway_ = '/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/entropy/GO/'
pathway_ = '/Users/woojunshim/Research/Data/broadPeaks/GO/'
top_no = 5  # Number of GO terms to be included
terms = vector()
epi_='H3K27ac'
for (group in groups){
  filename_ = paste(pathway_,group,'_',epi_,'_broadest_top5.txt', sep='')
  file_ = read.table(filename_)
  for (i in 1:top_no){
    terms = c(terms, as.character(file_$V3[i]))
  }
}
terms= unique(terms)
table_ = data.frame(matrix(data=rep(1,length(groups)*length(terms)), ncol=length(groups), nrow=length(terms)))
colnames(table_) = groups
rownames(table_) = terms
for (no in 1:ncol(table_)){
  group = colnames(table_)[no]
  filename_ = paste(pathway_,group,'_',epi_,'_broadest_top5.txt', sep='')
  file_ = read.table(filename_)
  for (m in 1:nrow(table_)){
    t = rownames(table_)[m]
    idx = which(file_$V3==t)
    if (length(idx)!=0){
      table_[t, group] = file_[idx,4] 
    }
  }
}
po = read.table('/Users/woojunshim/Research/Data/Roadmap_IDs_.txt', stringsAsFactors = F)
for (i in 1:ncol(table_)){
  epi = colnames(table_)[i]
  colnames(table_)[i] = paste(gsub('_', ' ', po[which(po$V1==epi), 2]), '(',epi,')',sep='')
}

table_ = -log10(table_)
library(ggplot2)
library(reshape)
table_$tt = gsub('_', ' ', rownames(table_))
table_$tt <- factor(table_$tt, levels = table_$tt)
table = melt(table_)
colnames(table)[3] = 'p.value'

p = ggplot(table, aes(y=tt, x=variable, fill=p.value)) 
p + geom_tile()+  labs(x='', y='', fill='-log10(FDR)') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
filename_ = paste(pathway_,epi_,'_top5_broadest_GO_BP.png',sep='')
ggsave(filename_,width=10,height=10)

### Overlap graph of top 500 genes ranked by HMs
table_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/hm_top500_overlap_average.txt', stringsAsFactors = F)
table = melt(t(table_))
ggplot(table, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value)) + geom_text(aes(label=round(value,3))) + scale_fill_gradient(low='white', high='red') + labs(x='', y='', fill='Overlap coefficient')
ggsave('/Users/woojunshim/Research/Data/broadPeaks/hm_top500_overlap_average.png', width=10, height=10)

### DRAW HIERARCHICAL CLUSTERING FOR SVD
library(gplots)
pdf('/Users/woojunshim/Research/Data/broadPeaks/svd/plots/relative_cor_genes_d30_h3k27me3.pdf', width=10, height=10)
heatmap.2(as.matrix(t(eg_g[1:10,])), col=greenred(75), scale='none', density.info ='none', trace='none', Colv=F, labRow='')
dev.off()

### PLOT DISTRIBUTION OF H3K27ME3 / H3K9ME3 
tf_list_ = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf_list = tf_list_$V1
file_ = read.table('/Users/woojunshim/Research/Data/Paige/original/H3K27me3/H3K27me3_stats.txt',stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/Paige/H3K27me3_stats_new.txt',stringsAsFactors = F)
lev = ifelse(file_[,1] %in% tf_list, 1, 0)
labels = ifelse(lev==1, 'TF', 'NonTF')
tf = file_[which(lev==1), 3]
nontf = file_[which(lev==0), 3]

library(ggplot2)

#Sample data
dat = data.frame(values=file_[,3], label=labels)
dat$values = log10(dat$values)  # log transformation

#Plot
ggplot(dat, aes(x = values, fill = label)) + geom_density(alpha = 0.5) + labs(x='log10(median width)') +ggtitle('H3K27me3 Median')
ggsave('/Users/woojunshim/Research/Data/Paige/H3K27me3_log_median_dist.png',width=10,height=10)

### CREATING CLUSTERS OF GENES (XU ET AL. 2015)
dist_ = as.matrix(dist(t(eg_g[1:10,])))
rownames(dist_) = rownames(data_)
colnames(dist_) = rownames(data_)
k = 20
snn = dist_
write.table(round(dist_,6), '/Users/woojunshim/Research/Data/broadPeaks/svd/dist_eg_g.txt', sep='\t', quote=F)
SNN<-function(data, outfile, k, distance){
  
  if(missing(data)){
    stop(paste("Input data missing.",help,sep="\n"))
  }
  if(missing(outfile)){
    stop(paste("Output file name missing.",help,sep="\n"))
  }
  if(missing(k)){
    k=3
  }
  if(missing(distance)){
    distance<-"euclidean"  # other distance options refer to dist() in R
  }
  m<-as.data.frame(data)
  numSpl<-dim(data)[1]
  m<-dist(data, distance, diag=TRUE, upper=TRUE)
  x<-as.matrix(m)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  
  edges<-list()              # SNN graph
  for (i in 1:numSpl){
    j<-i
    while (j<numSpl){
      j<-j+1
      shared<-intersect(IDX[i,], IDX[j,])
      if(length(shared)>0){			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
        strength<-max(s)
        if (strength>0)
          edges<-rbind(edges, c(i,j,strength))
      }				
    }
  }
  write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}
SNN(data=t(eg_g[1:10,]), 'test.txt', k=3, distance='euclidean')
name_table = data.frame(no=seq(1,ncol(eg_g)), name=colnames(eg_g))
write.table(name_table, '/Users/woojunshim/Research/Data/broadPeaks/svd/clusters/gene_tables_d30.txt',sep='\t',quote=F)
idx1 = which(colnames(eg_g)=='NKX2-5')
idx2 = which(colnames(eg_g)=='GATA4')

## NOW RUN Cliq.py
## THEN COMBINE WITH SINGLE CELLS 
results = read.table('/Users/woojunshim/Research/Data/broadPeaks/svd/clusters/clusters_d30.txt', stringsAsFactors = F)
rownames(results) = rownames(data_)
data__ = readRDS('/Users/woojunshim/Research/scRNA/Exprs_DCVLnorm_unlog_minus1_pos_Day30.RDS') # or
temp = unlist(strsplit(rownames(data__), '_ENSG'))
tt = temp[seq(1,length(temp),2)]
rownames(data__) = tt

output_ = data.frame(matrix(nrow=ncol(data__), ncol=length(unique(results$V1))))
colnames(output_) = unique(results$V1)
rownames(output_) = colnames(data__)
data__ = data__[rownames(sorted_),]
for (no in 1:ncol(data__)){
  
}

table_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths.txt', stringsAsFactors = F)
#heatmap(as.matrix(table_))
tf_ = ifelse(rownames(table_) %in% tf_list, 1, 0)
tf__ = which(tf_==1)
nontf__ = which(tf_==0)
tf_mean = apply(table_[tf__,], 1, mean)
nontf_mean = apply(table_[nontf__,], 1, mean)
t.test(tf_mean, nontf_mean)

### PLOT BINARY HEATMAP
file_ = read.table('/Users/woojunshim/Research/Data/Paige/Paige_day14.txt', stringsAsFactors = F)
tt = read.table('/Users/woojunshim/Research/Data/cardiac_regulators__.txt', stringsAsFactors = F)
tt = tt$V1
file_$gene = rownames(file_)
library(ggplot2)
library(reshape)
table = melt(file_)
table$gene = factor(table$gene, levels=rev(unique(tt)))
p = ggplot(table, aes(x=variable, y=gene, fill=value))
p + geom_tile(color='grey')+scale_fill_gradient2(low='white',high='red') + ggtitle('Human embryonic cardiac cells at day 14') + labs(x='', y='', subtitle='Top 200 genes') + scale_alpha(guide = 'none')
ggsave('/Users/woojunshim/Research/Data/Paige/Paige_day14_plot.png',width=10,height=10)

### PLOT PEAK WIDTHS (AS Z-SCORES) FOR SELECTED GENES
library(ggplot2)
library(reshape2)
order = c('E003','E016','E004','E005','E006','E011','E012','E013','E007','E065','E066','E061','E062','E109','E047','E100','E106','E104','E105','E113','E087','E085','E084','E038','E112','E037','E079','E058','E059','E055','E056','E071','E050','E098','E094','E095','E096','E097')
mark = 'H3K27me3'
genes = c('NKX2-5','TBX5','GATA4','GATA6','ISL1','HAND1','HAND2','MEIS1','MEIS2','MYH6','MYH7','TNNI3')
genes = c('NKX2-5','TBX5','GATA4','MYH6','MYH7','TNNI3')
genes = c('ZNF567','ZNF678','ZNF345','ZNF595','ZNF382','ZNF544')
table_ = vector()
for (gene in genes){
  file__ = paste('/Users/woojunshim/Research/Data/broadPeaks/',gene,'_features_raw.txt',sep='')
  file_ = read.table(file__, stringsAsFactors = F)
  table_ = rbind(table_, file_[mark, order])
}
rownames(table_) = genes
table_$names = genes
table = melt(table_)
table$names = reorder(table$names, as.numeric(factor(table$names,genes)))
ggplot(table, aes(y=value, x=factor(variable), group=names, colour=names)) + geom_line() + geom_point() + ggtitle(mark) + xlab('cell type') + ylab('z-score') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
output_ = paste('/Users/woojunshim/Research/Data/broadPeaks/', mark,'_width_ZNF_dist.png',sep='')
ggsave(output_, height=10, width=10)

### PLOT HEAT MAP FOR H3K4ME3, H3K27ME3, EXP 
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/E095_g.txt', stringsAsFactors = F)
heatmap(as.matrix(file_))
k = kmeans(as.matrix(file_[,c(2,3)]), centers=10)
k$cluster[which(names(k$cluster)=='NKX2-5')]
group = k$cluster[which(names(k$cluster)=='NKX2-5')]
members = names(k$cluster)[which(k$cluster==group)]
library(cluster)
library(fpc)
plotcluster(as.matrix(file_), k$cluster)
genes = c('GATA4','GATA6','HAND2','HAND1','TBX5','TBX20','TBX3','NKX2-5','MEIS1','MEIS2','MYH6','MYH7','TNNI3')
results = data.frame(matrix(nrow=length(genes), ncol=3))
colnames(results) = c('H3K4me3_exp','H3K27me3_exp','H3K4me3_H3K27me3_exp')
rownames(results) = genes

for (name in genes){
  results[name, 'H3K4me3_exp']= ifelse(name %in% members, 1, 0)
}
write.table(results, '/Users/woojunshim/Research/Data/broadPeaks/E095_kmeans_NKX2.5.txt', sep='\t',quote=F)

### COMPARE TF VS NON-TF
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_median_summary.txt', stringsAsFactors = F)
tf = file_[which(file_$V4==1), 2]
nontf = file_[which(file_$V4==0),2]
file_ = file_[order(file_$V2, decreasing=T),]
library(ggplot2)

#Sample data
dat = data.frame(values=file_[,2], label=file_[,4])
dat$values = log10(dat$values)  # log transformation

#Plot
ggplot(dat, aes(x = values, fill = label, group=label)) + geom_density(alpha = 0.5) + labs(x='log10(median width)') +ggtitle('H3K27me3 Median')


### PLOT TF PROPORTIONS 
table_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/Combined_TF_prop_200_excl.txt', stringsAsFactors = F)
table_ = t(table_)
#table_$name = seq(1,6)
table = melt(table_)
colnames(table) = c('HM','bin','value')
ggplot(table, aes(y=value, x=factor(bin), group=HM, colour=HM)) + geom_line() + geom_point() + xlab('bin number') + ylab('TF proportion') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(breaks=seq(1, 176, by=10))
output_ = paste('/Users/woojunshim/Research/Data/broadPeaks/Combined_TF_prop_200_excl.png',sep='')
ggsave(output_, height=10, width=10)

### ANALYSE 'tf_ranked_RPKM.txt'
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf = tf_list$V1
file_tf = read.table('/Users/woojunshim/Research/Data/tf_ranked_RPKM.txt', stringsAsFactors = F)
file_exp = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths.txt', stringsAsFactors = F)
jj_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_mad_table.txt', stringsAsFactors = F)
rownames(jj_)=jj_$V1

sd(file_['HAND2',])
table_ = apply(as.matrix(file_), 1, sd)
table_ = apply(as.matrix(file_), 1, median)
table_ = apply(as.matrix(file_), 1, mean)
oo = apply(as.matrix(file_), 1, mean)
cv = table_/oo
ref = rownames(file_tf)[order(file_tf$E095)][1:1569]
hist(table_[ref])
new_ = sort(table_[ref], decreasing = T)[1:500]
uu = ifelse(ref %in% names(new_), 1, 0)
pp = which(uu==1)
ref[pp[1:50]]
281 %in% pp
tt = data.frame(values=table_)
aa = ifelse(rownames(tt) %in% tf, 1, 0)
tt$TF = aa
tt = tt[order(tt$values, decreasing=T),]
write.table(tt, '/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_mean_table.txt',sep='\t',quote=F)
write.table(ref[pp], '/Users/woojunshim/Research/Data/broadPeaks/E095_high_sd_expressed_genes.txt', quote=F, sep='\t')

file_exp = file_exp[order(file_exp$E095, decreasing = T),]
exp_genes = rownames(file_exp)[which(file_exp$E095>1.0)]
mad_ = jj_[exp_genes, 2] 
names(mad_) = exp_genes
mad_ = sort(mad_, decreasing = T)
cul_mad = mad_ / sum(mad_)
selected_genes = mad_[1:500]
final_ = file_exp[names(selected_genes), 'E095']
names(final_) = names(selected_genes)
yy = log10(final_) * selected_genes
yy = sort(yy, decreasing = T)
result_E095 = data.frame(value=yy)
write.table(result_E095, '/Users/woojunshim/Research/Data/broadPeaks/E095_results_500.txt', quote=F, sep='\t')

# Palpant data or Paige data
file_exp = read.table('/Users/woojunshim/Research/Data/Palpant/TPM_ave.txt', stringsAsFactors = F)
file_exp = read.table('/Users/woojunshim/Research/Data/Paige/expression/expression_data_ave_symbol.txt', stringsAsFactors = F)
file_exp = file_exp[order(file_exp$C.EC_RNAseq, decreasing = T),]
#rownames(file_exp) = file_exp$V1
exp_genes = rownames(file_exp)[which(file_exp$CPC_RNA.seq>1.0)]
mad_ = jj_[exp_genes, 2] 
names(mad_) = exp_genes
mad_ = sort(mad_, decreasing = T)
cul_mad = mad_ / sum(mad_)
plot(x=seq(1,length(cul_mad)), y=cul_mad, main='CPC', xlab='rank', ylab='MAD/sum(mad)')
plot(x=seq(1,length(cul_mad)), y=cul_mad, main='CPC', xlab='rank', ylab='Median/sum(Median)')

selected_genes = mad_[1:150]
final_ = file_exp[names(selected_genes), 'CPC_RNA.seq']
names(final_) = names(selected_genes)
#yy = log10(final_) * selected_genes
yy = final_
yy = sort(yy, decreasing = T)
result_ = data.frame(value=yy)
aa = ifelse(rownames(result_) %in% tf, 1, 0)
result_$TF = aa
result_ = result_[order(result_$TF, decreasing = T),]
write.table(result_, '/Users/woojunshim/Research/Data/Palpant/C.EC_results.txt', sep='\t', quote=F)

aa = ifelse(exp_genes %in% tf, 1, 0)
ss = ifelse(exp_genes %in% names(selected_genes), 1, 0)
ty = data.frame(genes=exp_genes, TF=aa, exp=file_exp[exp_genes, 3], selected=ss)

### USING ALL GENES THAT HAVE AT LEAST ONE BROAD H3K27ME3 PEAKS
broad_genes = read.table('/Users/woojunshim/Research/Data/Paige/original/H3K27me3/H3K27me3_broad_genes.txt', stringsAsFactors = F)
ii = ifelse(rownames(file_exp) %in% broad_genes$V1, 1, 0)
selected_genes = file_exp[which(ii==1),3]
names(selected_genes) = rownames(file_exp)[which(ii==1)]
i = data.frame(value=selected_genes)
write.table(i,'/Users/woojunshim/Research/Data/Palpant/C.EC_all_broad_genes.txt', sep='\t',quote=F)

### 
CPC_broad = read.table('/Users/woojunshim/Research/Data/Palpant/CPC_all_broad_genes.txt', stringsAsFactors = F)
rownames(CPC_broad) = CPC_broad$V1
exp_genes = rownames(CPC_broad)[which(CPC_broad$V2>1.0)]
mad_ = jj_[exp_genes, 2] 
names(mad_) = exp_genes
mad_ = sort(mad_, decreasing = T)
cul_mad = mad_ / sum(mad_)
plot(x=seq(1,length(cul_mad)), y=cul_mad, main='CPC', xlab='rank', ylab='MAD/sum(mad)')
aa = ifelse(exp_genes %in% tf, 1, 0)
result_ = data.frame(exp=CPC_broad[exp_genes, 2], tf=aa)
rownames(result_) = exp_genes

### PLOT ENRICHMENT PLOT FOR MAD DENSITY (TF VS NON-TF)
file_ = read.table('/Users/woojunshim/Research/Data/mm10/h3k27me3/mad_density_table_mm10.txt', stringsAsFactors = F)
file_ = file_[,2:ncol(file_)]
file_ = t(file_)
colnames(file_) = c('TF','Non-TF')
rownames(file_) = seq(99,0, by=-1)
table_ = melt(file_)
colnames(table_) = c('Var1','group','value')
ggplot(table_, aes(y=value, x=factor(Var1), group=group, colour=group)) + geom_line() + geom_point() + ggtitle('Distribution of genes by MAD H3K27me3') + xlab('Percentile rank') + ylab('Enrichment score') + scale_x_discrete(breaks=seq(0, 100, by=10))
ggsave('/Users/woojunshim/Research/Data/mm10/h3k27me3/Dist_mad_mm10.png', width=10, height=10)

### GO ANALYSIS FOR MAD SELECTED GENES (ROADMAP)
# GO analysis for each groups of significant genes
library(topGO)
library(GO.db)
t_ = read.table('/Users/woojunshim/Research/Data/MAD_tfs_roadmap.txt', stringsAsFactors = F)
t_ = read.table('/Users/woojunshim/Research/Data/combined_results_skipped_seed_roadmap.txt', stringsAsFactors = F)
# Create GO terms to gene ID mapping table
all_genes = rownames(t_)
groups = colnames(t_)
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (no in 1:ncol(t_)){
  group = colnames(t_)[no]
  idx = which(t_[,no]!=0)
  genelist = rownames(t_)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  top_no = length(which(aa<0.05))
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/mad_analysis/GO/skipped_seed/',group,'_skipped_seed_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### PERFORM PCA W/O FILTERING VS FILTERING
fil_ = read.table('/Users/woojunshim/Research/Data/selected_genes_all_TF_roadmap.txt', stringsAsFactors = F)
exp_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
exp_ = read.table('/Users/woojunshim/Research/melanoma/verf_exp.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/Data/combined_results_skipped_seed_roadmap.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/melanoma/combined_genes_filtered_verf.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/melanoma/score_table_verfaillie.txt', stringsAsFactors = F)
epigenomes = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/epigenomes_list_.txt', stringsAsFactors = F)

exp_ = exp_[which(rownames(exp_) %in% tf),]
combined_ = combined_[which(rownames(exp_) %in% tf),]

exp_ = combined_
file_ = t(exp_[,intersect(epigenomes$V1,colnames(exp_))])
exp_ = t(exp_)
exp_ = exp_[,2:ncol(exp_)]
tissues = rep('p',nrow(exp_))
i1 = grep('221',rownames(exp_))
i2 = grep('225',rownames(exp_))
tissues[c(i1,i2)] = 'i'
new_ = exp_
new_ = file_


new_ = file_[c('E071','E070','E082','E055','E056','E062','E037','E047','E038','E059','E061','E057','E058','E028','E027','E104','E105','E095','E016','E003','E024'),]
y = apply(new_, 2, mean)
y_ = which(y!=0)
new__ = new_[,y_]
tissues = vector()
for (i in 1:nrow(new__)){
  idx = which(epigenomes$V1==rownames(new__)[i])
  tissues = c(tissues, epigenomes[idx, 2])
}
mycol_ = rainbow(length(unique(tissues)))
mycol = mycol_[as.numeric(as.factor(tissues))]
pca_ = prcomp(as.matrix(new__), scale.=T, center=T)
pdf('/Users/woojunshim/Research/Data/PCA_after_roadmap_skipped_seed.pdf',width=5, height=5)
plot(x=seq(1,length(pca_$sdev)), y=pca_$sdev/sum(pca_$sdev), type='o')
#plot(pca_$x[,1:2], col=mycol, main='After filtering (without seed analysis)')
legend('topright', col=unique(mycol), unique(tissues), pch=1)
dev.off()

library(ggbiplot)
g = ggbiplot(pca_, obs.scale = 1, var.scale = 1,groups = tissues, ellipse = TRUE,circle = TRUE)
g = g + scale_color_discrete(name = '')
g = g + theme(legend.direction = 'horizontal', legend.position = 'top')
print (g)

### PLOT HEAT MAP FOR JACCARD INDEX FOR TFS IN CPC (PALPANT)
file_ = read.table('/Users/woojunshim/Research/Data/Palpant/correlated_genes_CPC_RNAseq.txt', stringsAsFactors = F)
library(ggplot2)
library(reshape)
table = melt(file_, id='TF')
p = ggplot(table, aes(x=variable, y=TF, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=0.5, mid='grey70') + ggtitle('Normalised TF ChIP-seq counts per 1 kbp') + labs(x='Percentile Rank', y='', subtitle='Ranked by the FPKM value')
p + geom_tile()
new_ = ifelse(file_>0.564814814815, 1, 0)

### FIND DEGS IN PALPANT
library(DEGseq)
file_ = read.table('/Users/woojunshim/Research/Data/Palpant/TPM_ave.txt', stringsAsFactors = F)
file_ = cbind(file_, gene=rownames(file_))
DEGexp(geneExpMatrix1 = file_, geneCol1=12, expCol1=c(5,9), geneExpMatrix2 = file_, geneCol2=12, expCol2=c(1,2,3,4,6,7,8,10,11), groupLabel1 = 'Invasive', groupLabel2 = 'Proliferative', outputDir='/Users/woojunshim/Research/melanoma/DEG/')

### PLOT ROC CURVE
fil_ = read.table('/Users/woojunshim/Research/Data/Paige/ROC_seed_paige_d14_18tf.txt', stringsAsFactors = F)
deg_ = read.table('/Users/woojunshim/Research/Data/Palpant/ROC_output_score_CPC_genes_.txt', stringsAsFactors = F)
exp_ = read.table('/Users/woojunshim/Research/Data/Paige/ROC_exp_d14_18tf.txt', stringsAsFactors = F)
ski_ = read.table('/Users/woojunshim/Research/Data/Paige/ROC_skipped_paige_d14_18tf.txt', stringsAsFactors = F)
pre_ = read.table('/Users/woojunshim/Research/Data/Paige/ROC_pre-defined_paige_d14_18tf.txt', stringsAsFactors = F)
bac_ = read.table('/Users/woojunshim/Research/Data/Palpant/ROC_background_CPC_genes.txt', stringsAsFactors = F)
k20_ = read.table('/Users/woojunshim/Research/Data/Paige/ROC_20kb_paige_d14_18tf.txt', stringsAsFactors = F)
pdf('/Users/woojunshim/Research/Data/Paige/ROC_day14.pdf', width=10, height=10)
plot(1, type="n", xlab="f[False positive rate (1-specificity)]", ylab="True positive rate (sensitivity)", main = 'CROC curve (alpha=14)', xlim=c(0, 1), ylim=c(0, 1))
plot(1, type="n", xlab="False positive rate (1-specificity)", ylab="True positive rate (sensitivity)", main = 'ROC curve', xlim=c(0, 1), ylim=c(0, 1))
len = nrow(fil_)
lines(x=seq(0,1,length.out = len), y=seq(0,1,length.out = len), col='black', lwd=1, lty=3)
#lines(x=bac_$V2, y=bac_$V1, col='black', lwd=1, lty=3)
#lines(x=deg_$V2, y=deg_$V1, col='chartreuse3', lwd=2, lty=1)
lines(x=k20_$V2, y=k20_$V1, col='chartreuse3', lwd=2, lty=1)
lines(x=exp_$V2, y=exp_$V1, col='darkgoldenrod3', lwd=2, lty=1)
lines(x=fil_$V2, y=fil_$V1, col='red', lwd=2, lty=1)
lines(x=ski_$V2, y=ski_$V1, col='purple', lwd=2, lty=1)
lines(x=pre_$V2, y=pre_$V1, col='darkblue', lwd=2, lty=1)
legend(x=0.75, y=0.15, c('Seed (0.961)','Pre-defined(0.897)','Skipped (0.896)', 'DEG (0.769)','TPM (0.746)'), col=c('red','darkblue','purple','chartreuse3','darkgoldenrod3'), lwd=2)
legend(x=0.75, y=0.15, c('Seed (0.861)','Skipped (0.642)','Pre-defined (0.667)', '20kb (0.678)','TPM (0.817)'), col=c('red','darkblue','purple','chartreuse3','darkgoldenrod3'), lwd=2)
legend(x=0.75, y=0.15, c('Seed (0.869)','Skipped (0.669)','Pre-defined (0.669)', '20kb (0.674)','TPM (0.735)'), col=c('red','darkblue','purple','chartreuse3','darkgoldenrod3'), lwd=2)
dev.off()

### DISTBITUTION OF NUMBERS OF FILTERED GENES
counts = read.table('/Users/woojunshim/Research/Data/filtered_genes_count_roadmap.txt')
ind <- which.max(abs(diff(diff(x)/diff(y))/diff(y[-1])))

### PLOT RANK POSITION CHANGES FOR MELANOMA DATA
library(ggplot2)
library(reshape2)
order = c('FPKM.217','FPKM.218','FPKM.219','FPKM.220','FPKM.222','FPKM.223','FPKM.224','FPKM.226','FPKM.227','FPKM.221','FPKM.225')
genes = c('MITF','NFIB','POU3F2')
table_ = read.table('/Users/woojunshim/Research/melanoma/results_rank.txt', stringsAsFactors = F)
table_ = t(table_)
table_ = subset(table_, select=genes )
table = melt(table_)
table$Var1 = reorder(table$Var1, as.numeric(factor(table$Var1,order)))
ggplot(table, aes(y=value, x=factor(Var1), group=Var2, colour=Var2)) + geom_line() + geom_point() + ggtitle('Proliferative vs. Invasive') + xlab('sample') + ylab('rank') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(trans = "reverse")
output_ = paste('/Users/woojunshim/Research/melanoma/rank_comparison.png',sep='')
ggsave(output_, height=10, width=10)

### CREAT MAD TABLE
temp = read.table('/Users/woojunshim/Research/Data/mm10/h3k27me3/genes/H3K27me3_widths.txt', stringsAsFactors = F)
tf_list_ = read.table('/Users/woojunshim/Research/Data/TF_list_mm10_symbol.txt', stringsAsFactors = F)
tf_list = tf_list_$V1
tt = apply(temp, 1, mad)
uu = which(rownames(temp) %in% tf_list)
uu_ = rep(0, nrow(temp))
uu_[uu] = 1
results = data.frame(MAD=tt, TF=uu_)
rownames(results) = rownames(temp)
results = results[order(results$MAD, decreasing=T),]
write.table(results, '/Users/woojunshim/Research/Data/H3K27me3_mad_table_mm10.txt', quote=F, sep='\t')

### GET GENE CONVERSION TABLE USING biomRt
library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl") # for mouse
annot<-getBM(c("ensembl_gene_id", "mgi_symbol"), mart=ensembl)
write.table(annot, '/Users/woojunshim/Research/Data/mm10_ensembl_symbols.txt', quote=F, sep='\t')

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # for human
annot<-getBM(c("ensembl_gene_id", "hgnc_symbol"), mart=ensembl)
write.table(annot, '/Users/woojunshim/Research/Data/hg19_ensembl_symbols.txt', quote=F, sep='\t')

### PLOT ROC CURVE 
### FOR OUTPUT FROM 'PERFORMANCE_ANALYSIS1'
pathway_ = '/Users/woojunshim/Research/Data/Pax6/'
tpr_ = read.table(paste(pathway_, '_new_ROC_tpr.txt', sep=''), stringsAsFactors = F)
fpr_ = read.table(paste(pathway_, '_new_ROC_fpr.txt', sep=''), stringsAsFactors = F)
auc_ = read.table(paste(pathway_, '_new_ROC_auc.txt', sep=''), stringsAsFactors = F)
labels = vector()
for (i in 1:nrow(tpr_)){
  aa = paste(tpr_$V1[i],' ','(',round(auc_$V1[i],3),')', sep='')
  labels = c(labels, aa)
}
pdf(paste(pathway_, 'ROC_curve_new.pdf', sep=''), width=10,height=10)
plot(1, type="n", xlab="False positive rate (1-specificity)", ylab="True positive rate (sensitivity)", main = 'ROC curve', xlim=c(0, 1), ylim=c(0, 1))
len = ncol(tpr_)-1
col_ = rainbow(nrow(tpr_))
lines(x=seq(0,1,length.out = len), y=seq(0,1,length.out = len), col='black', lwd=1, lty=3)
for (i in 1:nrow(tpr_)){
  lines(x=fpr_[i, 2:ncol(fpr_)], y=tpr_[i, 2:ncol(tpr_)], col=col_[i], lwd=2, lty=1)
}
legend('bottomright', labels, col=col_, lwd=2)
dev.off()

### SPEARMAN'S CORRELATION BETWEEN CELL TYPES
table_ = read.table('/Users/woojunshim/Research/Data/57epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
table_ = table_[,c('E071','E070','E082','E055','E056','E062','E037','E047','E038','E059','E061','E057','E058','E028','E027','E104','E105','E095','E016','E003','E024')]
order_ = c('E071','E070','E082','E055','E056','E062','E037','E047','E038','E059','E061','E057','E058','E028','E027','E104','E105','E095','E016','E003','E024')
cor_ = cor(as.matrix(table_), method='spearman')
uu = melt(cor_)
ggplot(uu, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value)) + scale_fill_gradient2(midpoint=0.5, low='blue', mid='white', high='red') + ggtitle("Spearman's rho (before)") + xlab('Cell type') + ylab('Cell type') + scale_x_discrete(limits=order_) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('')

table1 = read.table('/Users/woojunshim/Research/Data/combined_results_roadmap_filtered_TF.txt', stringsAsFactors = F)
table1 = table1[,c('E071','E070','E082','E055','E056','E062','E037','E047','E038','E059','E061','E057','E058','E028','E027','E104','E105','E095','E016','E003','E024')]
cor1 = cor(as.matrix(table1), method='spearman')
uu1 = melt(cor1)
ggplot(uu1, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value)) + scale_fill_gradient2(midpoint=0.5, low='blue', mid='white', high='red') + ggtitle("Spearman's rho (after)") + xlab('Cell type') + ylab('Cell type') + scale_x_discrete(limits=order_) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

### INVESTIGATING H3K27ME3 AVERAGE
table_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_mean_table.txt', stringsAsFactors = F)
y= log10(table_$V2)
names(y) = table_$V1
tf = table_[which(table_$V3==1),]
nontf = table_[which(table_$V3==0),]
y_tf = log10(tf$V2)
y_nontf = log10(nontf$V2)
hist(y_tf, breaks=100)
hist(log10(table_$V2), breaks=100)

### PLOT ENRICHMENT (FET) DNA BINDING AND TFS FOR H3K27ME3 AVE PEAK GENES
library(ggplot2)
library(reshape2)
f1 = read.table('/Users/woojunshim/Research/Data/broadPeaks/fet_tf.txt')
f2 = read.table('/Users/woojunshim/Research/Data/broadPeaks/fet_GO.0003677.txt')
table_ = cbind(f1$V1, f2$V1)
table_ = -log10(table_)
rownames(table_) = seq(nrow(f1)-1,0,-1)
colnames(table_) = c('TF','DNA binding (GO:0003677)')
table = melt(table_)
colnames(table) = c('percentile','group','value')
ggplot(table, aes(x=percentile, y=value, group=group, colour=group)) + geom_line() + geom_point() + ylab('-log10(p-value)') + xlab('Percentile rank')
ggsave('/Users/woojunshim/Research/Data/broadPeaks/enrichment_fet_tf_go.0003677.png', width=10, height=5)

### FIGURE: RPKM VS H3K27ME3 WIDTH (E047)
exp_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
wid_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths.txt', stringsAsFactors = F)
genes = rownames(wid_)
x = wid_[genes, 'E047']
y = exp_[genes, 'E047']
plot(log2(x), log2(y), xlab='H3K27me3 peak width (bp)', ylab='RPKM')

### BINARY TABLE FOR GENES WITH BROAD H3K27ME3 PEAKS ACROSS CELL TYPES
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/genes_broad_h3k27me3_binary.txt', stringsAsFactors = F)
genes_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/genes_broad_h3k27me3_significant.txt', stringsAsFactors = F)
genes = genes_$V1
table_ = temp[genes,]
table__ = ifelse(table_==1, 0, 1)
table__ = melt(table__)
ggplot(table__, aes(x=Var2, y=Var1)) + geom_tile(aes(fill=value)) + scale_fill_gradient2(midpoint=0.5, low='blue', mid='white', high='red')  + xlab('Cell type') + ylab('Cell type') + theme(axis.text.x = element_text(angle = 45, hjust = 1))

### PLOT DISTRIBUTION OF H3K27ME3 COUNTS 
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/broad_gene_count_.txt', stringsAsFactors = F)
library(ggplot2)
temp$nontf = temp$V2 - temp$V3
ggplot(temp, aes(x=V1, y=V3)) + geom_bar(stat="identity") + xlab('Number of broad H3K27me3 peaks') + ylab('Number of TF genes')
ggsave('/Users/woojunshim/Research/Data/broadPeaks/dist_broad_H3K27me3_peaks_tf.png', width=10, height=10)

### PLOT CUMULATIVE CURVES 
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/broad_gene_count_.txt')
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/cumulative_prop_broad_peaks.pdf', width=5, height=5)
plot(1, type="n", xlab='Number of the broad H3K27me3 peak(s)', ylab="Cumulative proportion", xlim=c(0, 111), ylim=c(0, 1), xaxt='n', yaxt='n')
axis(side=2, at=c(0,0.5,1), las=1)
axis(side=1, at=seq(0,110,10), labels=seq(110,0,-10),las=2)
lines(x=seq(0,111), y=temp$V6, col='blue', type='o')
#lines(x=seq(0,111), y=temp$V7, col='green', type='o')
#lines(x=seq(0,111), y=temp$V8, col='red', type='o')
dev.off()

### GO ANALYSIS ON EXPRESSED TFs WITH (OR WITHOUT) BROAD PEAK(S)
tf_list = read.table('/Users/woojunshim/Research/Data/TF/TF_combined.txt', stringsAsFactors = F)
tf = tf_list$V1
library(topGO)
library(GO.db)
epis = c('E038','E095','E082')
for (epi in epis){
  t_ = read.table(paste('/Users/woojunshim/Research/Data/',epi,'_expressed_broad_vs_narrow_genes.txt', sep=''), stringsAsFactors = F)
  # Create GO terms to gene ID mapping table
  all_genes = t_[,1]
  top_no = 100 # top 100 GO terms by FDR
  go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
  idx = which((t_[,3]==0) & (t_[,2]==0))
  genelist = t_[,1][idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  top_no = length(which(aa<0.05))
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/',epi,'_expressed_narrow_nontf_GO.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### HEAT MAP FOR OVERLAP TABLES
library(ggplot2)
library(reshape)
file_ = read.table('/Users/woojunshim/Research/Data/overlap_prop_narrow_tf_.txt', stringsAsFactors = F)
table = melt(t(file_))
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=0.5, mid='white',low="blue",high="red", limits=c(0,1)) + ggtitle('TF genes (without broad H3K27me3 peak)') + labs(x='', y='', fill='Jaccard Index') 
ggsave('/Users/woojunshim/Research/Data/overlap_prop_narrow_tf_.png', width=10, height=10)

### HORIZONTAL BAR PLOTS
file_ = read.table('/Users/woojunshim/Research/Data/extracted_go_table_tf_E095.txt', stringsAsFactors = F)
colnames(file_) = c('H3K27me3 positive', 'H3K27me3 negative')
rownames(file_) = gsub('_', ' ', rownames(file_))
table_ = melt(t(file_))
group = unique(table_$X1)
ggplot(table_, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat='identity', width=.5, position='dodge') + coord_flip() + labs(ylab='-log10(FDR)', xlab='', fill='') + ylab('-log10(FDR)') + xlab('') + ggtitle('Left ventricle (E095)') + theme(axis.text=element_text(size=12, face='bold'))
output_= paste('/Users/woojunshim/Research/Data/E095_tf_FDR_.png',sep='')
ggsave(output_, width=10, height=5)

### ENRICHMENT PLOTS FOR SELECTED GO TERMS (BROAD GENES)
file_ = read.table('/Users/woojunshim/Research/Data/go_enrichment_all.txt', stringsAsFactors = F)
file_ = file_[,2:ncol(file_)]
file_ = t(file_)
colnames(file_) = c('heart development (GO:0007507)','brain development (GO:0007420)','hemopoiesis (GO:0030097)')
rownames(file_) = seq(95,0, by=-5)
table_ = melt(file_)
colnames(table_) = c('Var1','group','value')
ggplot(table_, aes(y=value, x=factor(Var1), group=group, colour=group)) + geom_line() + geom_point() + ggtitle('') + xlab('Percentile rank') + ylab('Enrichment score') + scale_x_discrete(breaks=seq(0, 100, by=10))
ggsave('/Users/woojunshim/Research/Data/go_enrichment_all.png', width=10, height=5)

### PLOT H3K27ME3 WIDTHS DIST
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths.txt', stringsAsFactors = F)
temp_ = temp$E095
idx = which(temp_!=0)
temp_ = temp_[idx]
temp_ = sort(temp_, decreasing = T)
pdf('/Users/woojunshim/Research/Data/broadPeaks/E095_H3K27me3_width.pdf', width=10, height=5)
plot(x=seq(1, length(temp_)), y=temp_, type='o', ylab='H3K27me3 peak width', xlab='Rank position')
dev.off()

### PLOT LOCATIONAL PREFERENCE TABLE 
file_ = read.table('/Users/woojunshim/Research/Data/overlap_prop_narrow_tf_.txt', stringsAsFactors = F)
table = melt(t(file_))
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=0.5, mid='white',low="blue",high="red", limits=c(0,1)) + ggtitle('TF genes (without broad H3K27me3 peak)') + labs(x='', y='', fill='Jaccard Index') 
ggsave('/Users/woojunshim/Research/Data/overlap_prop_narrow_tf_.png', width=10, height=10)

### HEAT MAP FOR OVERLAP TABLES
library(ggplot2)
library(reshape)
file_ = read.table('/Users/woojunshim/Research/Data/odds_ratio_locational_preference.txt', stringsAsFactors = F)
table = melt(t(file_))
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=1.0, mid='white',low="blue",high="red", limits=c(0,max(table$value))) + ggtitle('') + labs(x='', y='', fill='Odds ratio') + geom_text(aes(label=round(table$value, 4))) + scale_x_discrete(limits=c("upstream","TSS","downstream")) + scale_y_discrete(limits=c("non_TF","TF","not_broad","broad")) + theme(axis.text=element_text(size=12, face='bold'))
ggsave('/Users/woojunshim/Research/Data/odds_ratio_locational_preference.png', width=5, height=5)

### FROM 'COMBINED_BINARY_46_STAT.TXT'
temp = read.table('/Users/woojunshim/Research/Data/combined_binary_46_stat.txt', stringsAsFactors = F)
idx = which(temp$TF==0)
t = temp[idx,]
ref = read.table('/Users/woojunshim/Research/Data/GO.0030097_tf_9606_.txt', stringsAsFactors = F)
ref = ref$V1
idx = which(rownames(temp) %in% ref)
t = temp[idx,]
cor(t$broad, t$exp, method='spearman')

### SIDE-BY-SIDE BOX PLOTS FOR EXP VS BROAD FROM 'COMBINED_BINARY_46_STAT.TXT'
ggplot(temp, aes(x=broad, y=exp, group=broad)) + geom_boxplot() + labs(y='No.expression', x='No.broad')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/exp_vs_broad.png', width=5, height=5)

### Venn diagram for constitutively expressed genes
library(VennDiagram)
a = read.table('/Users/woojunshim/Research/Data/constitutively_expressed_genes.txt', stringsAsFactors = F)
a = a$V1
b = read.table('/Users/woojunshim/Research/Data/not_broad_genes.txt', stringsAsFactors = F)
b = b$V1
c = intersect(a,b)
grid.newpage()
png('/Users/woojunshim/Research/Data/venn_exp_vs_broad.png', width=500, height=500)
draw.pairwise.venn(length(a), length(b), length(c), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

### HEAT MAP FOR TABLES
library(ggplot2)
library(reshape)
file_ = read.table('/Users/woojunshim/Research/Data/odds_ratio_locational_preference_peaks.txt', stringsAsFactors = F)
file_ = log2(file_)
table = melt(t(file_))
p = ggplot(table, aes(x=X1, y=X2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=0.0, mid='white',low="blue",high="red") + ggtitle('') + labs(x='', y='', fill='log2(odds ratio)') + geom_text(aes(label=round(table$value, 2))) + scale_x_discrete(limits=c("upstream","TSS","downstream")) + scale_y_discrete(limits=c("narrow","broad")) + theme(axis.text=element_text(size=12, face='bold'))
ggsave('/Users/woojunshim/Research/Data/odds_ratio_locational_preference_domains.png', width=5, height=3)

### PLOT VARIOUS GRAPHS FOR 'GENE_SUMMARY.TXT'
temp = read.table('/Users/woojunshim/Research/Data/genes_summary.txt', stringsAsFactors = F)
temp$V6 = temp$V3 / temp$V4
# 1. Box plot
table_ = subset(temp, select=c(V4, V5))
ggplot(table_, aes(x=V5, y=log10(V4), group=V5)) + geom_boxplot()
# 2. Side-by-side box plots
ggplot(temp, aes(x=V2, y=V3, group=V2)) + geom_boxplot() + labs(y='Median H3K27me3 peak width (bp)', x='No.broad peaks') 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/median_vs_broad.png', width=5, height=5)
# 3. Enrichment plot 
file_ = read.table('/Users/woojunshim/Research/Data/genes_density_table.txt', stringsAsFactors = F)
file_ = file_[,2:ncol(file_)]
file_ = t(file_)
colnames(file_) = c('median width','gene size','normalised median width')
rownames(file_) = seq(99,0, by=-1)
table_ = melt(file_)
colnames(table_) = c('Var1','method','value')
ggplot(table_, aes(y=value, x=factor(Var1), group=method, colour=method)) + geom_line() + geom_point() + ggtitle('Enrichment of TFs') + xlab('Percentile rank') + ylab('Enrichment score') + scale_x_discrete(breaks=seq(0, 100, by=10))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/TF_enrichment_median_width.png', width=8, height=5)

# 4. median width distribution (histogram)
table_ = subset()
ggplot(dat, aes(x = values, fill = label)) + geom_density(alpha = 0.5) + labs(x='log10(median width)') +ggtitle('H3K27me3 Median')

### PLOT PILEUP FILES
temp = read.table('/Users/woojunshim/Research/Data/pileup_100bp_genes_E038_combined.txt', stringsAsFactors = F)
#colnames(temp) = c('TF','non-TF')
colnames(temp) = c('Broad genes','Not-broad genes')
table_ = melt(t(temp))
colnames(table_) = c('group','X2','value')
ggplot(table_, aes(x=X2/10, y=value, col='red')) + ylim(0,1) + labs(x='Distance to TSS (kb)', y='Proportion of peaks') + geom_line() + ggtitle('E038')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/pileup_genes_E038.png', width=5, height=5)

### PLOT PILEUPS (AGGREGATION PLOT)
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/1B_high.png', height=5, width=5, units='in', res=600)
plot(1, type="n", xlab="", ylab="", xlim=c(-50, 50), ylim=c(0, 0.7), xaxt='n',yaxt='n')
title(xlab='Position relative to TSS (kb)', ylab='Density')
axis(1, at=c(-50,-25,0,25,50), labels=c(-50,-25,0,25,50))
axis(2, at=seq(0,0.7, 0.1), labels=seq(0,0.7, 0.1), las=1)
file_ = read.table('/Users/woojunshim/Research/Data/pileup/pileup_100bp_broad_peaks_combined.txt', stringsAsFactors = F)
idx = order(as.numeric(rownames(file_)))
file_ = file_[idx,]
cols = c('E038','E082','E095','E104','E105','E070','E071','E037','E047')
cols = colnames(file_)
for (n in cols){
  lines(x=seq(-50,50,0.1), y=file_[,n], col='red')
}
file_ = read.table('/Users/woojunshim/Research/Data/pileup/pileup_100bp_not_broad_peaks_combined.txt', stringsAsFactors = F)
idx = order(as.numeric(rownames(file_)))
file_ = file_[idx,]
cols = c('E038','E082','E095','E104','E105','E070','E071','E037','E047')
cols = colnames(file_)
for (n in cols){
  lines(x=seq(-50,50,0.1), y=file_[,n], col='blue')
}
legend("topleft", c('Broad','Not-broad'), col=c('red','blue'), lty=c(1,1))
graphics.off()

### BOXPLOTS FOR EXP. PROPORTIONS AND LOCATIONS
t1 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_tf.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_nontf.txt', stringsAsFactors = F)
t1$group = rep('TF',nrow(t1))
t2$group = rep('non-TF',nrow(t2))
table_ = rbind(t1[,c(2,3,4,5)], t2[,c(2,3,4,5)])
colnames(table_) = c('upstream','TSS','downstream','group')
table = melt(table_)
table$group = factor(table$group, levels=c('TF','non-TF'))
ggplot(table, aes(x=variable, y=value, fill=group)) + geom_boxplot() + scale_x_discrete(limits=c("upstream","TSS","downstream")) + labs(x='',y='Proportion of expressed genes')
a = subset(table, variable=='downstream')
wilcox.test(a$value ~ a$group)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/exp_prop_vs_location.png', width=5, height=5)

### BOXPLOTS FOR EXP. PROPORTIONS AND BROADNESS
t1 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_tf.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_nontf.txt', stringsAsFactors = F)
t1$group = rep('TF',nrow(t1))
t2$group = rep('non-TF',nrow(t2))
table_ = rbind(t1[,c(2,3,4)], t2[,c(2,3,4)])
colnames(table_) = c('not-broad','broad','group')
table = melt(table_)
table$group = factor(table$group, levels=c('TF','non-TF'))
ggplot(table, aes(x=variable, y=value, fill=group)) + geom_boxplot() + scale_x_discrete(limits=c("broad","not-broad")) + labs(x='',y='Proportion of expressed genes')
a = subset(table, variable=='downstream')
wilcox.test(a$value ~ a$group)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/exp_prop_vs_broadness.png', width=5, height=5)


### COMPLIANCE TABLE (LOCATION 2 AND EXP OUTCOME)
temp = read.table('/Users/woojunshim/Research/Data/location2_not_exp_compliance_genes.txt', stringsAsFactors = F)
wilcox.test(temp$V4 ~ temp$V5, alternative='less')
a = temp[,c(4,5)]

### ENRICHMENT PLOT FOR COMPLIANCE TABLE
file_ = read.table('/Users/woojunshim/Research/Data/enrichment_tf_broadness_exp_compliance.txt', stringsAsFactors = F)
file_ = file_[,2:ncol(file_)]
file_ = t(file_)
colnames(file_) = c('Broad','Not-broad')
rownames(file_) = seq(99,0, by=-1)
table_ = melt(file_)
colnames(table_) = c('Var1','method','value')
ggplot(table_, aes(y=value, x=factor(Var1), group=method, colour=method)) + geom_line() + geom_point() + ggtitle('') + xlab('Percentile rank') + ylab('Enrichment score') + scale_x_discrete(breaks=seq(0, 100, by=10))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/enrichment_tf_compliance_by_broadness_.png', width=8, height=5)

### WILCOXON TEST COMPARING TF VS NON-TF 
t = read.table('/Users/woojunshim/Research/Data/location_3_notexp_compliance_genes_.txt', stringsAsFactors = F)
a = t[,c(2,5)]
wilcox.test(a$V2 ~ a$V5, alternative='less')

### CALCULATE COEFFICIENCE OF VARIANCE FOR GENES 
t = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_exp.txt', stringsAsFactors = F)
p = apply(t, 1, sd) / apply(t, 1, mean) 
names(p) = rownames(t)
aa = data.frame(row.names=names(p), value=p, TF=ifelse(names(p) %in% tf, 1, 0))
write.table(aa, '/Users/woojunshim/Research/Data/cv_genes.txt', sep='\t', quote=F)

### DENSITY PLOT FOR REGULATED TFS ACROSS 46 CELL TYPES 
library(ggplot2)
library(reshape2)
file_ = read.table('/Users/woojunshim/Research/Data/density_regulated_tf_exp.txt', stringsAsFactors = F)
rownames(file_) = file_[,1]
file_ = file_[,2:ncol(file_)]
colnames(file_) = seq(99, 0, -1)
table_ = melt(t(file_))
colnames(table_) = c('Percentile','Cell','value')
ggplot(table_, aes(y=value, x=factor(Percentile), group=Cell)) + geom_line(colour='coral', show.legend = F)  + xlab('Percentile rank') + ylab('Enrichment score') + scale_x_discrete(breaks=seq(0, 100, by=10))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/enrichment_regulated_tf_by_exp.png', width=8, height=5)

### CREATE COLOURED-GRID TABLE 
file_ = read.table('/Users/woojunshim/Research/Data/cv_table.txt', stringsAsFactors = F)
library(gridExtra)
library(grid)
t1 = ttheme_default(core=list(bg_params = list(fill=c('coral1','deepskyblue1','darkgoldenrod2','green3'))))
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/cv_table.pdf', width=2, height=1)
grid.table(file_, theme = t1)
dev.off()

### BOX PLOTS FOR PROPORTIONS
temp = read.table('/Users/woojunshim/Research/Data/peak_location_stat.txt', stringsAsFactors = F)
rownames(temp) = temp[,1]
temp = temp[,c(2,3,4)]
colnames(temp) = c('upstream','TSS','downstream')
sum_ = apply(temp, 1, sum)
temp$sum_ = sum_
temp[,c(1,2,3)] = temp[,c(1,2,3)] / temp$sum_
library(ggplot2)
library(reshape2)
table_ = temp[,c(1,2,3)]
table = melt(t(table_))
ggplot(table, aes(x=Var1, y=value)) + ylim(0,1) + geom_boxplot() + labs(x='', y='Proportion')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/peak_location_proportion.png', width=5, height=5)

temp = read.table('/Users/woojunshim/Research/Data/broad_TSS_counts_combined.txt', stringsAsFactors = F)
rownames(temp) = temp[,1]
temp = temp[,c(2,3)]
colnames(temp) = c('broad','narrow')
sum_ = apply(temp, 1, sum)
temp$sum_ = sum_
temp[,c(1,2)] = temp[,c(1,2)] / temp$sum_
library(ggplot2)
library(reshape2)
table_ = temp[,c(1,2)]
table = melt(t(table_))
ggplot(table, aes(x=Var1, y=value)) + ylim(0,1) + geom_boxplot() + labs(x='', y='Proportion')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/peak_broadness_proportion.png', width=5, height=5)

temp = read.table('/Users/woojunshim/Research/Data/broad_TSS_counts_combined.txt', stringsAsFactors = F)
rownames(temp) = temp[,1]
temp = temp[,c(13,15,17,19,21,23)]
temp = t(temp)
#temp = cbind(temp, c('broad','not-broad','broad','not-broad','broad','not-broad'))
#temp = cbind(temp, c('TSS','TSS','upstream','upstream','downstream','downstream'))
#colnames(temp)[112] = 'broadness'
#colnames(temp)[113] = 'location'
table = melt(temp, id=c("broadness","location"))
table = cbind(table, rep(c('TSS','TSS','upstream','upstream','downstream','downstream'),111))
table = cbind(table, rep(c('broad','narrow','broad','narrow','broad','narrow'),111))
colnames(table)[c(4,5)] = c('location','broadness')
ggplot(table, aes(x=location, y=value, colour=broadness)) + ylim(0,1) + geom_boxplot() + labs(x='', y='Proportion', colour='') + scale_x_discrete(limits=c('upstream','TSS','downstream'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/peak_location_proportion_cond_broadness.png', width=5, height=5)

temp = read.table('/Users/woojunshim/Research/Data/broad_TSS_counts_combined.txt', stringsAsFactors = F)
rownames(temp) = temp[,1]
temp = temp[,c(12,14,16,18,20,22)]
temp = t(temp)
#temp = cbind(temp, c('broad','not broad','broad','not broad','broad','not broad'))
#temp = cbind(temp, c('TSS','TSS','upstream','upstream','downstream','downstream'))
#colnames(temp)[112] = 'broadness'
#colnames(temp)[113] = 'location'
table = melt(temp, id=c("broadness","location"))
table = cbind(table, rep(c('TSS','TSS','upstream','upstream','downstream','downstream'),111))
table = cbind(table, rep(c('broad','narrow','broad','narrow','broad','narrow'),111))
colnames(table)[c(4,5)] = c('location','broadness')
table$location = factor(table$location,levels(table$location)[c(3,2,1)])
ggplot(table, aes(x=broadness, y=value, colour=location)) + ylim(0,1) + geom_boxplot() + labs(x='', y='Proportion', colour='') 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/peak_broadness_proportion_cond_location.png', width=5, height=5)

# BOXPLOTS FOR FET 
# BROADNESS
t1 = read.table('/Users/woojunshim/Research/Data/fet_not_broad_s.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/fet_broad_s.txt', stringsAsFactors = F)
t1 = cbind(t1, rep('narrow', nrow(t1)))
t2 = cbind(t2, rep('broad', nrow(t2)))
colnames(t1)[6] = 'group'
colnames(t2)[6] = 'group'
table_ = rbind(t1, t2)
colnames(table_)[c(1:5)] = c('Cell','Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
table = melt(table_[,c(2,3,4,5,6)], id='group')
table$group = factor(table$group, levels(table$group)[c(2,1)])
ggplot(table, aes(x=group, y=value, colour=variable)) + geom_boxplot() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='', y='Odds ratio', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/odds_ratio_broadness.png', width=5, height=5)

# LOCATION
t1 = read.table('/Users/woojunshim/Research/Data/fet_upstream_s.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/fet_TSS_s.txt', stringsAsFactors = F)
t3 = read.table('/Users/woojunshim/Research/Data/fet_downstream_s.txt', stringsAsFactors = F)
t1 = cbind(t1, rep('upstream', nrow(t1)))
t2 = cbind(t2, rep('TSS', nrow(t2)))
t3 = cbind(t3, rep('downstream', nrow(t2)))
colnames(t1)[6] = 'group'
colnames(t2)[6] = 'group'
colnames(t3)[6] = 'group'
table_ = rbind(t1, t2, t3)
colnames(table_)[c(1:5)] = c('Cell','Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
table = melt(table_[,c(2,3,4,5,6)], id='group')
table$group = factor(table$group, levels(table$group)[c(3,2,1)])
ggplot(table, aes(x=group, y=value, colour=variable)) + geom_boxplot() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='', y='Odds ratio', colour='') + scale_x_discrete(limits=c('upstream','TSS','downstream'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/odds_ratio_location.png', width=8, height=5)

### ENRICHMENT PLOT FOR GROUPS OF GENES
file_ = read.table('/Users/woojunshim/Research/Data/enrichment_group_genes_by_no.broad_peaks.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/enrichment_group_genes_by_no.tss_peaks.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/enrichment_group_genes_by_no.broad_tss_peaks.txt', stringsAsFactors = F)
file_ = file_[,2:ncol(file_)]
file_ = t(file_)
colnames(file_) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
rownames(file_) = seq(99,0, by=-1)
table_ = melt(file_)
colnames(table_) = c('Var1','group','value')
ggplot(table_, aes(y=value, x=factor(Var1), group=group, colour=group)) + geom_line() + geom_point() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(y='Enrichment score', x='Percentile rank', colour='') + scale_x_discrete(breaks=seq(0, 100, by=10)) + ylim(0,16)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/enrichment_group_genes_by_no.broad_tss_peaks.png', width=8, height=5)

### WILCOXON TEST COMPARING CV BETWEEN GROUPS OF GENES 
temp = read.table('/Users/woojunshim/Research/Data/cv_genes_.txt', stringsAsFactors = F)
a= temp$value[which(temp$gene_class==1)]
b= temp$value[which(temp$gene_class==2)]
c= temp$value[which(temp$gene_class==3)]
d= temp$value[which(temp$gene_class==4)]
wilcox.test(c,b)

### BOX PLOT (GROUPS OF GENES VS CV)
temp = read.table('/Users/woojunshim/Research/Data/cv_genes_.txt', stringsAsFactors = F)
table_ = subset(temp, select=c('gene_class','value'))
idx = which(table_$gene_class==1)
table_$gene_class[idx] = 'Regulated TF'
idx = which(table_$gene_class==2)
table_$gene_class[idx] = 'Regulated non-TF'
idx = which(table_$gene_class==3)
table_$gene_class[idx] = 'Stable TF'
idx = which(table_$gene_class==4)
table_$gene_class[idx] = 'Stable non-TF'
table_$gene_class = factor(table_$gene_class, levels(as.factor(table_$gene_class))[c(2,1,4,3)])
ggplot(table_, aes(x=gene_class, y=value, colour=gene_class)) + geom_boxplot() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='', y='Coefficient of variation', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/cv_group_genes.png', width=8, height=5)

### BOX PLOT (CV VS NUMBER OF EXPRESSED CELL TYPES)
temp = read.table('/Users/woojunshim/Research/Data/cv_genes_.txt', stringsAsFactors = F)
table_ = subset(temp, select=c('no.expressed','value'))
ggplot(table_, aes(x=no.expressed, y=value, group=no.expressed)) + geom_boxplot() + labs(x='No.cell types(expressed)', y='Coefficient of variation')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/cv_vs_no.expressed_cell_types.png', width=8, height=5)

### BOX PLOT (NO.BROAD PEAKS VS NUMBER OF EXPRESSED CELL TYPES)
temp = read.table('/Users/woojunshim/Research/Data/exp_broad_stat.txt', stringsAsFactors = F)
ggplot(temp, aes(x=exp, y=broad, group=exp)) + geom_boxplot() + labs(x='No.cell types(expressed)', y='No.broad peaks') 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/no.broad_peaks_vs_no.expressed_cell_types.png', width=8, height=5)

### BOX PLOT (GROUPS OF GENES VS NUMER OF EXPRESSED CELL TYPES)
temp = read.table('/Users/woojunshim/Research/Data/exp_broad_stat_.txt')
idx = which(temp$group==1)
temp$group[idx] = 'Variably expressed TF'
idx = which(temp$group==2)
temp$group[idx] = 'Variably expressed non-TF'
idx = which(temp$group==3)
temp$group[idx] = 'Stably expressed TF'
idx = which(temp$group==4)
temp$group[idx] = 'Stably expressed non-TF'
temp$group = factor(temp$group, levels(as.factor(temp$group))[c(4,3,2,1)])
ggplot(temp, aes(x=group, y=exp, fill=group)) + geom_boxplot() + labs(x='', y='Number of cell types (expressed)',fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure4A.png', width=8, height=5, dpi=300)

### BOXPLOT FOR EXP. PROPORTIONS OF GENES FOR 4 GROUPS
temp = read.table('/Users/woojunshim/Research/Data/exp_prop_group_genes.txt')
colnames(temp) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
table = melt(t(temp))
ggplot(table, aes(x=Var1, y=value, colour=Var1)) + geom_boxplot() + labs(x='',y='Proportion of expression', colour='') + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + ylim(0,1)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/group_genes_vs_exp.prop.png', width=5, height=5)

### BOXPLOTS FOR EXP. PROPORTIONS AND LOCATIONS
### FOR 4 GROUPS OF GENES (LOCATION)
t1 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_regulated_tf.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_regulated_nontf.txt', stringsAsFactors = F)
t3 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_stable_tf.txt', stringsAsFactors = F)
t4 = read.table('/Users/woojunshim/Research/Data/exp_location_prop_stable_nontf.txt', stringsAsFactors = F)
t1$group = rep('Regulated TF',nrow(t1))
t2$group = rep('Regulated non-TF',nrow(t2))
t3$group = rep('Stable TF',nrow(t3))
t4$group = rep('Stable non-TF',nrow(t4))
table_ = rbind(t1[,c(2,3,4,5)], t2[,c(2,3,4,5)], t3[,c(2,3,4,5)], t4[,c(2,3,4,5)])
colnames(table_) = c('upstream','TSS','downstream','group')
table = melt(table_)
table$group = factor(table$group, levels=c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF'))
ggplot(table, aes(x=variable, y=value, colour=group)) + geom_boxplot() + scale_x_discrete(limits=c("upstream","TSS","downstream")) + labs(x='',y='Proportion of expression', colour='') + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + ylim(0,1)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/group_genes_vs_location_exp.prop.png', width=8, height=5)

### FOR 4 GROUPS OF GENES (BROADNESS)
t1 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_regulated_tf.txt', stringsAsFactors = F)
t2 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_regulated_nontf.txt', stringsAsFactors = F)
t3 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_stable_tf.txt', stringsAsFactors = F)
t4 = read.table('/Users/woojunshim/Research/Data/exp_broad_prop_stable_nontf.txt', stringsAsFactors = F)
t1$group = rep('Regulated TF',nrow(t1))
t2$group = rep('Regulated non-TF',nrow(t2))
t3$group = rep('Stable TF',nrow(t3))
t4$group = rep('Stable non-TF',nrow(t4))
table_ = rbind(t1[,c(2,3,4)], t2[,c(2,3,4)], t3[,c(2,3,4)], t4[,c(2,3,4)])
colnames(table_) = c('narrow','broad','group')
table = melt(table_)
table$group = factor(table$group, levels=c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF'))
ggplot(table, aes(x=variable, y=value, colour=group)) + geom_boxplot() + scale_x_discrete(limits=c("broad","narrow")) + labs(x='',y='Proportion of expression', colour='') + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + ylim(0,1)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/group_genes_vs_broad_exp.prop.png', width=8, height=5)

### CALCULATING STATISTICS
colnames(table_) = c('not_broad','broad','group')
a1 = table_$broad[which(table_$group=='Regulated TF')]
b1 = table_$broad[which(table_$group=='Regulated non-TF')]
c1 = table_$broad[which(table_$group=='Stable TF')]
d1 = table_$broad[which(table_$group=='Stable non-TF')]
a2 = table_$not_broad[which(table_$group=='Regulated TF')]
b2 = table_$not_broad[which(table_$group=='Regulated non-TF')]
c2 = table_$not_broad[which(table_$group=='Stable TF')]
d2 = table_$not_broad[which(table_$group=='Stable non-TF')]

a = temp$exp[which(temp$group==1)]
b = temp$exp[which(temp$group==2)]
c = temp$exp[which(temp$group==3)]
d = temp$exp[which(temp$group==4)]

ub = table$value[which(table$location=='upstream' & table$broadness=='broad')]
un = table$value[which(table$location=='upstream' & table$broadness=='not_broad')]
tb = table$value[which(table$location=='TSS' & table$broadness=='broad')]
tn = table$value[which(table$location=='TSS' & table$broadness=='not_broad')]
db = table$value[which(table$location=='downstream' & table$broadness=='broad')]
dn = table$value[which(table$location=='downstream' & table$broadness=='not_broad')]

### BOXPLOT FOR EXP. PROPORTIONS BETWEEN BROAD AND NARROW PEAKS
temp = read.table('/Users/woojunshim/Research/Data/prop_exp_vs_broadness.txt', stringsAsFactors = F)
rownames(temp) = temp$V1
table_ = temp[,c(2,3)]
colnames(table_) = c('Broad','Narrow')
table_ = melt(table_)
ggplot(table_, aes(x=variable, y=value, colour=variable)) + geom_boxplot() + labs(y='Proportion of expression', x='', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/prop_exp_vs_broadness.png', width=5, height=5)

### HEATMAP FOR JACCARD INDEX BETWEEN BROAD & NARROW GENES
epis = c('E038','E082','E095','E104','E105','E070','E071','E037','E047')
epis = sort(epis)
file_ = read.table('/Users/woojunshim/Research/Data/jaccard_btw_narrow_genes.txt', stringsAsFactors = F)
table_ = subset(file_, select=epis)
table_ = subset(t(table_), select=epis)
table = melt(t(table_))
p = ggplot(table, aes(x=Var1, y=Var2, fill=value))
p + geom_tile()+scale_fill_gradient2(midpoint=0.5, mid='white',low="blue",high="red", limits=c(0,1)) + labs(x='', y='', fill='Jaccard Index') 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/jaccard_btw_narrow_genes.png', width=6, height=5)

### HEAT MAP FOR TF ENRICHMENT BETWEEN BROAD AND NARROW GENES
file_ = read.table('/Users/woojunshim/Research/Data/fet_tf_btw_broadness.txt', stringsAsFactors = F)
rownames(file_) = file_$V1
file_ = file_[,c(1,2)]
file_$V2 = -log10(file_$V2)
table_ = file_[which(file_$V1 %in% epis),]
ggplot(table_, aes(x=V2, y=V1)) + geom_point(colour='red')  + labs(x='-log10(p-value)', y='') + xlim(0, max(table_$V2))

### BOX PLOT 
file_ = read.table('/Users/woojunshim/Research/Data/prop_broad_over_repressed_.txt', stringsAsFactors = F)
rownames(file_) = file_$V1
file_ = file_[,c(2,3)]
colnames(file_) = c('TF','All')
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1, y=value, colour=Var1)) + geom_boxplot() + labs(x='', y='', colour='') + ggtitle('') + labs(y='proportion')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/prop_broad_over_repressed_.png', width=5, height=5)

### PROPORTION FOR 4 GENE GROUPS (MEDIAN - CURRENT) FOR ALL EXPRESSED GENES 
file_ = read.table('/Users/woojunshim/Research/Data/gene_group_prop_H3K27me3_width_percentile_all.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/gene_group_prop_expression_percentile_only_expressed.txt', stringsAsFactors = F)
file_ = file_[,c(2:ncol(file_))]
rownames(file_) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
colnames(file_) = seq(99,0,-1)
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1+0.5, y=value, fill=Var2)) + geom_bar(stat='identity') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='Percentile rank', y='Proportion', fill='') + theme_bw() 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure4B_high.png', width=8, height=5, dpi=600)


### GO BP ANALYSIS FOR REGULATED TF
file__ = read.table('/Users/woojunshim/Research/Data/binary_regulated_tf.txt', stringsAsFactors = F)
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes_ = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)
all_genes = all_genes$V1
groups = colnames(file__)
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (no in 1:ncol(file__)){
  group = colnames(file__)[no]
  idx = which(file__[,no]==1)
  genelist = rownames(file__)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/regulated_tf/go',group,'_GO_BP_regulated_tf.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### BOXPLOT FOR NUMBER OF BROAD H3K27ME3 PEAKS (4 GENE GROUPS)
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/genes_broad_h3k27me3_stats_gene_groups.txt', stringsAsFactors = F)
g1 = temp$V2[which(temp$V5==1)]
g2 = temp$V2[which(temp$V5==2)]
g3 = temp$V2[which(temp$V5==3)]
g4 = temp$V2[which(temp$V5==4)]
table_ = data.frame(value=c(g1,g2,g3,g4), group=c(rep('Regulated TF',length(g1)),rep('Regulated non-TF',length(g2)),rep('Stable TF',length(g3)),rep('Stable non-TF',length(g4))))
table_$group = factor(table_$group, levels=c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF'))
ggplot(table_, aes(x=group, y=value, colour=group)) + geom_boxplot() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='', y='Number of broad H3K27me3 peak(s)', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/no.broad_peaks_gene_groups.png', width=8, height=5)

### PROPORTION FOR 4 GENE GROUPS NUMBER OF THE BROAD PEAKS 
file_ = read.table('/Users/woojunshim/Research/Data/gene_group_prop_no.broad_peaks.txt', stringsAsFactors = F)
file_ = file_[,c(2:ncol(file_))]
rownames(file_) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
colnames(file_) = seq(110,0,-2)
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1+1, y=value, fill=Var2)) + geom_bar(stat='identity') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='Number of the broad H3k27me3 peak(s)', y='', fill='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/gene_group_prop_no.broad_peaks.png', width=8, height=5)

### PLOT HEAT MAP ('/Users/woojunshim/Research/Data/odds_ratio_broad_peaks_for_gene_groups.txt')
file_ = read.table('/Users/woojunshim/Research/Data/odds_ratio_broad_peaks_for_gene_groups.txt', stringsAsFactors = F)
file_ = read.table('/Users/woojunshim/Research/Data/odds_ratio_top10_exp_for_gene_groups.txt', stringsAsFactors = F)
rownames(file_) = file_$V1
file_ = file_[,c(2,3,4,5)]
colnames(file_) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1, y=Var2, fill=log2(value))) + geom_tile() + scale_fill_gradient2(midpoint=0, low = "blue", high = "red") + labs(y='', x='', fill='log2(odds-ratio)') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure6_high.png', width=8, height=10, units='in', dpi=600)

### PLOTS FOR 'COND_PROP_BROAD_TSS.TXT' DEPENDENCY 
file_ = read.table('/Users/woojunshim/Research/Data/cond_prop_broad_tss.txt', stringsAsFactors = F)
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/jaccard_broad_tss.png', width=500, height=500)
boxplot(file_[,'Jaccard_index'])
dev.off()

table_ = data.frame(TSS=file_$p.TSS.broad.-file_$p.TSS., broad=file_$p.broad.TSS.-file_$p.broad.)
colnames(table_) = c('p(')
table = melt(t(table_))
       
# Venn diagram
grid.newpage()
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/venn_broad_tss.png', width=250, height=250)
draw.pairwise.venn(29.2, 89.9, 19.1, lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2))
dev.off()

### VENN DIAGRAM FOR PROPORTION OF GENE GROUPS
# BROAD, TSS OR BROAD+TSS 
t1 = read.table('/Users/woojunshim/Research/Data/broad_tss_gene_group_prop.txt')
t2 = read.table('/Users/woojunshim/Research/Data/broad_gene_group_prop.txt')
t3 = read.table('/Users/woojunshim/Research/Data/tss_gene_group_prop.txt')
table_ = rbind(BOTH=t1[112, c(2,3,4,5)], BROAD=t2[112, c(2,3,4,5)], TSS=t3[112, c(2,3,4,5)])
table_ = rbind(table_, BG=c(0.0408,0.4608,0.0503,0.4479))
colnames(table_) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
table = melt(t(table_))
ggplot(table, aes(x=factor(1), y=value, fill=Var1)) + geom_bar(width=1, stat='identity') + facet_grid(facets=. ~ Var2) + labs(y='Proportion',x='', fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/4_gene_groups_composition.png', width=5, height=5)

### 
table_ = rbind(c(0.257746734764,0.67339564366,0.00593495770838,0.0629226638674), c(0.219249607245,0.701567938884,0.00710039121461,0.0720820626559), c(0.077596299956,0.690578007324,0.0243237483151,0.207501944405), c(0.0436999658283,0.528913240373,0.0448088872964,0.382577906502))
rownames(table_) = c('B+T','BROAD','TSS','ALL')
colnames(table_) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
t_ = melt(table_)
ggplot(t_, aes(x=factor(1), y=value, fill=Var2)) + geom_bar(width=1, stat='identity') + facet_grid(facets=. ~ Var1) + labs(y='Proportion',x='', fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + theme(strip.text.x = element_text(size = 5, colour = "orange", angle = 90)) + theme_bw()
ggsave('/Users/woojunshim/Desktop/new/4_gene_groups_composition_new.png', width=5, height=5, dpi=1200)

### heatmap for spearman's correlation
file_ = read.table('/Users/woojunshim/Research/Data/spearman_expressed_regulated_nontf.txt', stringsAsFactors = F)
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(midpoint=0.5, low='blue',high='red', limits=c(0,1)) + labs(y='', x='', fill="Spearman's rho") +theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/spearman_expressed_regulated_nontf.png', height=5, width=8)

### DISTRIBUTION PLOT OF BROAD H3K27ME3 PEAKS
### CUMULATIVE LINE 
temp = read.table('/Users/woojunshim/Research/Data/count_gene_broad_tss.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V1, y=V2)) + geom_bar(stat='identity') + scale_x_reverse() + labs(x='Number of the broad H3K27me3 domains over the TSS', y='Gene count')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/count_gene_broad_tss.png', width=5, height=5)

cum_ = cumsum(rev(temp$V2))
cum_prop = cum_ / 17630
table_ = data.frame(number=rev(temp$V1), cumul=cum_prop)
ggplot(table_, aes(x=number, y=cumul)) + geom_line() + scale_x_reverse() + labs(x='Number of the broad H3K27me3 peaks over the TSS', y='Cumulative proportion')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/count_gene_broad_tss_cumulative.png', width=5, height=5)

### PROPORTION PLOTS FOR BROAD + TSS
file_ = read.table('/Users/woojunshim/Research/Data/gene_group_prop_no.broad_TSS_peaks_new.txt', stringsAsFactors = F)
file_ = file_[,c(2:ncol(file_))]
rownames(file_) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
colnames(file_) = seq(100,0,-10)
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1+5, y=value, fill=Var2)) + geom_bar(stat='identity') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(x='Number of broad-H3K27me3-TSS domains', y='Proportion', fill='') + theme_bw()
ggsave('/Users/woojunshim/Desktop/new/gene_group_prop_no.broad_TSS_peaks_new_high.png', width=8, height=5, dpi=1200)

### NEGATIVE AND OVERALL POPULATIONS
t1 = read.table('/Users/woojunshim/Research/Data/gene_group_prop_no.broad_TSS_peaks_negative.txt')
t2 = read.table('/Users/woojunshim/Research/Data/gene_group_prop_no.broad_TSS_peaks_all.txt')
table_ = cbind(NEGATIVE=t1[,2], ALL=t2[,2])
rownames(table_) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
table = melt(t(table_))
ggplot(table, aes(x=factor(1), y=value, fill=Var2)) + geom_bar(width=1, stat='identity') + facet_grid(facets=. ~ Var1) + labs(y='',x='', fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/gene_group_prop_no.broad_TSS_peaks_negative_all.png', width=5, height=5)

### FET PLOT FOR GENE GROUPS, WHEN RANKED BY NUMBER OF H3K27ME3-BROAD-TSS PEAKS
file_ = read.table('/Users/woojunshim/Research/Data/fet_gene_groups_for_broad_tss.txt', stringsAsFactors = F)
rownames(file_) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
file_ = file_[, c(2:ncol(file_))]
file_ = -log10(file_)
file_ = file_[c(1,2),]
colnames(file_) = seq(99,0,-1)
table_ = melt(t(file_))
ggplot(table_, aes(x=Var1, y=value, colour=Var2)) + geom_line(stat='identity', alpha=0.4) + geom_point() + scale_colour_manual(values=c('coral1','darkgoldenrod2')) + labs(y='-log10(p-value)', x='Percentile rank', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/fet_gene_groups_for_broad_tss.png', width=5, height=5)

### BOXPLOT OF SPEARMAN'S CORRELATION BETWEEN BROAD AND NARROW PEAKS 
file_ = read.table('/Users/woojunshim/Research/Data/broadPeaks/cor_btw_broad_narrow.txt', stringsAsFactors = F)
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/cor_btw_broad_narrow.png', width=500, height=500)
boxplot(file_$V2)
dev.off()


### CLUSTERING BEFORE AND AFTER FILTERING
fil_ = read.table('/Users/woojunshim/Research/Data/selected_genes_all_TF_roadmap.txt', stringsAsFactors = F)
exp_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
exp_ = read.table('/Users/woojunshim/Research/melanoma/verf_exp.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/Data/combined_results_skipped_seed_roadmap.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/melanoma/combined_genes_filtered_verf.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/melanoma/score_table_verfaillie.txt', stringsAsFactors = F)
combined_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_regulated_genes.txt', stringsAsFactors = F)
epigenomes = read.table('/Users/woojunshim/Research/Data/Gapped_peaks/H3K4me3/Background_models/epigenomes_list_.txt', stringsAsFactors = F)

exp_ = exp_[which(rownames(exp_) %in% tf),]
combined_ = combined_[which(rownames(exp_) %in% tf),]

exp_ = combined_
file_ = t(exp_[,intersect(epigenomes$V1,colnames(exp_))])
exp_ = t(exp_)
exp_ = exp_[,2:ncol(exp_)]
tissues = rep('p',nrow(exp_))
i1 = grep('221',rownames(exp_))
i2 = grep('225',rownames(exp_))
tissues[c(i1,i2)] = 'i'
new_ = t(exp_)
new_ = t(combined_)


new_ = new_[c('E071','E070','E082','E055','E056','E062','E037','E047','E038','E059','E061','E057','E058','E028','E027','E104','E105','E095','E016','E003','E024'),]
y = apply(new_, 2, mean)
y_ = which(y!=0)
new__ = new_[,y_]
tissues = vector()
for (i in 1:nrow(new__)){
  idx = which(epigenomes$V1==rownames(new__)[i])
  tissues = c(tissues, epigenomes[idx, 2])
}
mycol_ = rainbow(length(unique(tissues)))
mycol = mycol_[as.numeric(as.factor(tissues))]
pca_ = prcomp(as.matrix(new__), scale.=T, center=T)
pdf('/Users/woojunshim/Research/Data/PCA_after_roadmap_skipped_seed.pdf',width=5, height=5)
plot(x=seq(1,length(pca_$sdev)), y=pca_$sdev/sum(pca_$sdev), type='o')
plot(pca_$x[,1:2], col=mycol, main='After filtering (without seed analysis)')
legend('topright', col=unique(mycol), unique(tissues), pch=1)
dev.off()

### GTEX DATA
### GO BP ANALYSIS
groups = read.table('/Users/woojunshim/Research/Data/GTEX/GTEX_48_tissue_groups.txt', stringsAsFactors = F)
groups = groups$V1
library(topGO)
library(GO.db)
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (group in groups){
  filename_ = paste('/Users/woojunshim/Research/Data/GTEX/results_new/',group,'_rank_100_counts_product.txt',sep='')
  file_ = read.table(filename_)
  file_ = file_[order(file_$V2, decreasing=T),]
  genelist = file_[1:100,1]
  all_genes = file_[,1]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  top_no = length(which(aa<0.05))
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/Data/GTEX/results_new/GO/',group,'_BP_product.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### HORIZONTAL BAR PLOTS
file_ = read.table('/Users/woojunshim/Research/Data/GTEX/results/GO/GO_table_Whole_blood.txt')
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt', stringsAsFactors = F)
names = vector()
file_ = -log10(file_)
for (i in 1:nrow(file_)){
  ai = go_table[which(rownames(go_table)==rownames(file_)[i]),1]
  names = c(names, ai)
}
rownames(file_) = names
rownames(file_) = gsub('_', ' ', rownames(file_))
table_ = melt(t(file_))
table_$Var1 = factor(table_$Var1, levels=(c('Weighted','TPM')))
ggplot(table_, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat='identity', width=.5, position='dodge') + coord_flip() + labs(ylab='-log10(FDR)', xlab='', fill='') + ylab('-log10(FDR)') + xlab('') + ggtitle('Whole blood') + theme(axis.text=element_text(size=12, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/whole_blood.png', width=10, height=5)

### HIERARCHICAL CLUSTERING FOR 'COMBINED_TOP_100_PRODUCT.TXT'
file_ = read.table('/Users/woojunshim/Research/Data/GTEX/results/combined_top_100_product.txt', stringsAsFactors = F)
col__ = read.table('/Users/woojunshim/Research/Data/GTEX/GTEX_tissue_groups_.txt', stringsAsFactors = F)
col_ = vector()
for (i in 1:ncol(file_)){
  aa = col__$V2[which(col__$V1==colnames(file_)[i])]
  col_ = c(col_, aa)
}
colnames(file_) = gsub('_', ' ', colnames(file_))
col_ = gsub('_', ' ', col_)
rcol = rainbow(length(unique(col_)))
mycol = rcol[as.numeric(as.factor(col_))]
hmcol = colorRampPalette(c("white","orange", "red"))(256)
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/combined_top_100_product_full.pdf', height=10, width=10)
heatmap.2(as.matrix(file_), col=hmcol, trace='none', labRow=rep('', nrow(file_)), cexCol = 0.35, srtCol=45)
dev.off()

### CREATE A HIERARCHICAL PLOT FOR SELECTED TISSUE TYPES
file_ = read.table('/Users/woojunshim/Research/Data/GTEX/results/combined_top_100_product.txt', stringsAsFactors = F)
col__ = read.table('/Users/woojunshim/Research/Data/GTEX/GTEX_tissue_groups_.txt', stringsAsFactors = F)
col_ = vector()
colp = vector()
for (i in 1:ncol(file_)){
  aa = col__$V2[which(col__$V1==colnames(file_)[i])]
  col_ = c(col_, aa)
  aa = col__$V3[which(col__$V1==colnames(file_)[i])]
  colp = c(colp, aa)
}
col_ = gsub('_', ' ', col_)
idx = which(col_ %in% c("Skin","Heart","Brain","Artery","Adipose"))
file_ = file_[, idx]
col_ = vector()
colp = vector()
for (i in 1:ncol(file_)){
  aa = col__$V2[which(col__$V1==colnames(file_)[i])]
  col_ = c(col_, aa)
  aa = col__$V3[which(col__$V1==colnames(file_)[i])]
  colp = c(colp, aa)
}
colp = gsub('_', ' ', colp)
colnames(file_) = colp
yy = apply(file_, 1, sum)
y = which(yy!=0)
file_ = file_[y,]
rcol = rainbow(length(unique(col_)))
mycol = rcol[as.numeric(as.factor(col_))]
hmcol = colorRampPalette(c("white","orange", "red"))(256)
par(mar=c(7,4,4,2)+0.1) 
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/combined_top_100_product_selected.png', height=10, width=10, units='in', res=300)
heatmap.2(as.matrix(file_), col=hmcol, trace='none', labRow=rep('', nrow(file_)), cexCol = 1, srtCol=45, ColSideColors = mycol, margins=c(12,8))
graphics.off()

### LEGEND PLOT
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/combined_top_100_product_selected_legend.png', height=3, width=3, units='in', res=300)
plot.new()
legend("center", unique(col_), pch=15, 
       col=rcol[unique(as.numeric(as.factor(col_)))], 
       bty="n")
graphics.off()

### ROC CURVE 
pathway = '/Users/woojunshim/Research/Data/NCC/'
name = 'Neural crest cells'
name_ = 'Neural crest cells'
auc = read.table(paste(pathway, 'test_auc_new.txt', sep=''), stringsAsFactors = F)
auc[,2] = round(auc[,2], 3)
tpr = read.table(paste(pathway, 'test_tpr_new.txt', sep=''), stringsAsFactors = F)
fpr = read.table(paste(pathway, 'test_fpr_new.txt', sep=''), stringsAsFactors = F)
colnames(tpr) = tpr[1,]
colnames(fpr) = fpr[1,]
tpr = tpr[2:nrow(tpr),]
fpr = fpr[2:nrow(fpr),]
png(paste('//Users/woojunshim/Desktop/new/','ROC_',name_,'_new_high.png',sep=''), width=8.6, height=8.6, units='cm', res=600)
plot(1, type="n", xlab="False positive rate", ylab="True positive rate", main = name, xlim=c(0, 1), ylim=c(0, 1), las=1)
len = nrow(tpr)
lines(x=seq(0,1,length.out = len), y=seq(0,1,length.out = len), col='black', lwd=1, lty=3)
lines(x=fpr$Exp, y=tpr$Exp, col='green', lwd=2, lty=1)
lines(x=fpr$H3K4me3, y=tpr$H3K4me3, col='blue', lwd=2, lty=1)
lines(x=fpr$Corrected, y=tpr$Corrected, col='red', lwd=2, lty=1)

#labels = c(paste('Discordant score (AUC =',auc[2,2],')', sep=''), paste('H3K4me3 breadth (AUC =',auc[1,2],')', sep=''), paste('Expresion value (AUC =',auc[3,2],')', sep=''))
labels = c(as.character(auc[2,2]), as.character(auc[1,2]), as.character(auc[3,2]))
legend('bottomright', labels, col=c('red','blue','green'), lwd=2, cex=0.8, bty = "n")
graphics.off()

### CREATE A LEGEND FOR ROC 
png('/Users/woojunshim/Desktop/new/legend.png', width=10, height=10, units='cm', res=600)
#plot.new()
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
labels = c('Discordant score', 'H3K4me3 breadth', 'Expression value')
legend('bottomright', labels, col=c('red','blue','green'), lwd=2, cex=0.8)
graphics.off()


### BOX PLOT FOR CORRELATION BETWEEN BROADEST / SUM OF DOMAINS VS. EXP. LEVEL
a = read.table('/Users/woojunshim/Research/Data/spearman_cor_btw_broadest_width_vs_exp.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/spearman_cor_btw_sum_width_vs_exp.txt', stringsAsFactors = F)
a$group = rep('Broadest',nrow(a))
b$group = rep('Sum',nrow(b))
table_ = rbind(a,b)
table = data.frame(value=table_$V2, group=table_$group)
ggplot(table, aes(x=group, y=value, colour=group)) + geom_boxplot(show.legend = F) + labs(y="Spearman's rho", x='', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/spearman_broadest_sum_vs_exp.png', width=5, height=5)

### PLOT DOT PLOT FOR BROAD-TSS + REPRESSED CELL TYPES FOR GENES 
temp = read.table('/Users/woojunshim/Research/Data/broad_tss_effect_count.txt')
rownames(temp) = temp[,1]
temp = temp[,c(2,3)]
colnames(temp) = c('count','group_')
temp$group_ = factor(temp$group_)
table_ = melt(temp)
ggplot(temp, aes(x=group_, y=count, colour=group_)) + geom_boxplot() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3'))
a = temp$count[which(temp$group_==1)]
b = temp$count[which(temp$group_==2)]
c = temp$count[which(temp$group_==3)]
d = temp$count[which(temp$group_==4)]
temp = temp[order(temp$count, decreasing = T),]
cols = rep('green3', nrow(temp))
cols[which(temp$group_==1)] = 'coral1'
cols[which(temp$group_==2)] = 'darkgoldenrod2'
cols[which(temp$group_==3)] = 'deepskyblue1'
temp$x = rev(seq(1, nrow(temp)))
temp$count = temp$count / 46
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/repressed_broad_tss_prop_.png', height=5, width=20, units='in', res=300)
plot(1, type="n", xlab="", ylab="", xlim=c(1, nrow(temp)), ylim=c(0, max(temp$count)), xaxt='n')
points(x=temp$x, y=temp$count, col=cols, pch='|')
plot(1, type="n", xlab="", ylab="", ylim=c(1, nrow(temp)), xlim=c(0, max(temp$count)), xaxt='n')
points(x=temp$count, y=temp$x, col=cols, pch='-')
graphics.off()

### PLOT CARDIAC SPECIFIC GENES 
temp = read.table('/Users/woojunshim/Research/Data/z_table_selected_cardiac_genes_exp.txt', stringsAsFactors = F)
temp = round(temp, 3)
table_ = melt(t(temp))
temp = temp[,order(colnames(temp))]
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red") + labs(y='', x='', fill='Exp. class') + scale_y_discrete(limits=c('GATA4','GATA6','NKX2-5','TBX5','TBX20','MYH6','MYH7','MYL2','MYL3','TNNI3')) +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(limits=colnames(temp))
#ggplot(table_, aes(x=Var1, y=Var2, fill=as.factor(value))) + geom_tile() +  labs(y='', x='', fill='Exp. class') + scale_y_discrete(limits=c('GATA4','GATA6','NKX2-5','TBX5','TBX20','MYH6','MYH7','MYL2','MYL3','TNNI3')) +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(limits=colnames(temp)) + scale_fill_manual(values = c('red','orange','yellow'), limits = c("4", "3","2"), labels = c("High", 'Moderate', "Low"))

ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/z_score_cardiac_genes_disc_exp.tif', width=10, height=5, dpi=600)

### PLOT ENRICHMENT SCORE FOR REPRESSED GENES BY BROAD-TSS DOMAIN
temp = read.table('/Users/woojunshim/Research/Data/broad_tss_effect_count_enrichment.txt', stringsAsFactors = F)
rownames(temp) = c('Regulated TF','Regulated non-TF','Stable TF','Stable non-TF')
temp = temp[,2:ncol(temp)]
colnames(temp) = seq(99,0,-1)
table_ = melt(t(temp))
ggplot(table_, aes(x=Var1+0.5, y=value, colour=Var2)) + geom_line() + geom_point() + scale_colour_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + labs(y='Enrichment score',x='Percentile rank', colour='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/figures/broad_tss_effect_count_enrichment.png', width=8, height=8, dpi=600)

### FIGURE 1A
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/E095_H3K27me3_genes.txt', stringsAsFactors = F)
png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/1A_high.png',width=8, height=5, units='in', res=600)
plot(x=seq(1, nrow(temp)), y=temp$V3, type='o', xlab='Rank position', ylab='H3K27me3 domain breadth')
graphics.off()

### GO BP ANALYSIS
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
file_ = read.table('/Users/woojunshim/Research/Data/H3K27me3_negative_tf_top100.txt', stringsAsFactors = F)
all_genes = rownames(read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F))
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
for (no in 1:ncol(file_)){
  group = colnames(file_)[no]
  idx = which(file_[,no]==1)
  genelist = rownames(file_)[idx]
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  top_no = length(which(aa < 0.05))
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file__ = paste('/Users/woojunshim/Research/Data/',group,'_GO_BP_H3k27me3_negative_top100_tf.txt', sep='')
  write.table(results, file__, sep='\t', quote=F)
}

### PLOT MERGED GO RESULTS HEATMAP
temp = read.table('/Users/woojunshim/Research/Data/GTEX/results_new/GO/Brain_-_Cerebellum_combined.txt', stringsAsFactors = F)
yy = gsub('_', ' ', rownames(temp))
rownames(temp) = yy
colnames(temp) = c('Discordant score','TPM')
temp = -log10(temp)
temp = temp[order(temp$"Discordant score", decreasing = T),]
table_ = melt(t(temp))

#png('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure1A_high.png', width=5, height=6, units='in', res=600)
#temp = temp[order(temp$positive, decreasing = T),]
#table_ = melt(t(temp))
#library('RColorBrewer')
#ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() 

p = ggplot(table_, aes(y=Var2, x=Var1, fill=value)) 
p + geom_tile()+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill='-log10(FDR)') + ggtitle('Cerebellum, Brain') 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/cerebellum_brain.png', height=30, width=10, units='in', dpi=300)

### FIGURE CARDIAC
temp = read.table('/Users/woojunshim/Research/Data/z_table_selected_cardiac_genes_h3k27me3.txt', stringsAsFactors = F)
temp = temp[,c('E003','E004','E005','E006','E016','E024','E037','E038','E047','E057','E058','E059','E070','E071','E082','E095','E104','E105')]
colnames(temp)[which(temp['NKX2-5',]==min(temp['NKX2-5',]))]
colnames(temp)[order(temp['NKX2-5',])]
table_ = melt(t(temp))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red") + labs(y='', x='', fill='z-score') + scale_y_discrete(limits=c('GATA4','GATA6','NKX2-5','TBX5','TBX20','MYH6','MYH7','MYL2','MYL3','TNNI3')) +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(limits=c('E105','E104','E095','E082','E071','E070','E059','E058','E057','E047','E038','E037','E024','E016','E003','E006','E005','E004'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/cardiac_high.png', height=5, width=7, units='in', dpi=600)

### SUPP. FIGURE 7A
a = read.table('/Users/woojunshim/Research/Data/exp_prop_broad_gene_groups_new.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/exp_prop_non_broad_gene_groups_new.txt', stringsAsFactors = F)
rownames(a) = a$V1
rownames(b) = b$V1
a = a[,c(2,3,4,5)]
b = b[,c(2,3,4,5)]
colnames(a) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
colnames(b) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
a_ = melt(t(a))
a_ = cbind(a_, class=rep('broad',nrow(a_)))
b_ = melt(t(b))
b_ = cbind(b_, class=rep('non-broad',nrow(b_)))
table_ = rbind(a_, b_)

table_$Var1 = factor(table_$Var1, levels=c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF'))
ggplot(table_, aes(x=class, y=value, fill=Var1)) + geom_boxplot() + scale_x_discrete(limits=c("broad","non-broad")) + labs(x='',y='Proportion of expression', fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + ylim(0,1) + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure7A.png', height=6, width=8, units='in', dpi=300)




### SUPP. FIGURE 7B
a = read.table('/Users/woojunshim/Research/Data/exp_prop_tss_gene_groups_new.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/exp_prop_non_tss_gene_groups_new.txt', stringsAsFactors = F)
rownames(a) = a$V1
rownames(b) = b$V1
a = a[,c(2,3,4,5)]
b = b[,c(2,3,4,5)]
colnames(a) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
colnames(b) = c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF')
a_ = melt(t(a))
a_ = cbind(a_, class=rep('TSS',nrow(a_)))
b_ = melt(t(b))
b_ = cbind(b_, class=rep('not-TSS',nrow(b_)))
table_ = rbind(a_, b_)

table_$Var1 = factor(table_$Var1, levels=c('Variably expressed TF','Variably expressed non-TF','Stably expressed TF','Stably expressed non-TF'))
ggplot(table_, aes(x=class, y=value, fill=Var1)) + geom_boxplot() + scale_x_discrete(limits=c("TSS","not-TSS")) + labs(x='',y='Proportion of expression', fill='') + scale_fill_manual(values=c('coral1','darkgoldenrod2','deepskyblue1','green3')) + ylim(0,1) + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure7B_high.png', height=6, width=8, units='in', dpi=600)

### ENTROPY FIGURES
a = read.table('/Users/woojunshim/Research/Data/GTEX/results_new/entropy_scores_product.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/GTEX/results_new/entropy_scores_tpm.txt', stringsAsFactors = F)
a$V1 = gsub('_', ' ', a$V1)
b$V1 = gsub('_', ' ', b$V1)
a$group = rep('Discordant score', nrow(a))
b$group = rep('TPM', nrow(b))
table_ = rbind(a, b)
table = melt(table_)
ggplot(table, aes(x=V1, y=value, fill=group)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))  + labs(x='', y='Entropy', fill='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/entropy_GTEX_high.png', width=20, height=10, units='in', dpi=600)

### WILCOXON TEST FOR ENTROPY
results = vector()
for (i in 1:nrow(a)){
  name = a$V1[i]
  yy = wilcox.test(as.numeric(a[i,2:101]), as.numeric(b[i,2:101]), alternative='less')
  results = c(results, yy$p.value)
}
results = p.adjust(results, method='fdr')
names(results) = a$V1
tt=data.frame(FDR=results)
write.table(tt, '/Users/woojunshim/Research/Data/GTEX/results_new/wilcoxon_FDR.txt', sep='\t', quote=F)

### SUPPLEMENTART FIGURE 8 (UNSTABLE GENES)
temp = read.table('/Users/woojunshim/Research/Data/signal_saturation_broad_tss_iter1000_bs3.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_mean_score_change.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15))
ggplot(table_, aes(x=no+1.5, y=mean)) + geom_point() + geom_line() + labs(x='Number of included samples', y='Mean proportion of stably ranked genes') + geom_errorbar(data=table_, mapping=aes(x=no+1.5, ymin=high, ymax=low), width=2, size=1)  + scale_y_continuous(breaks=seq(0.1, 0.9, by=0.1))+ theme_bw() 
ggplot(table_, aes(x=no+1.5, y=mean)) + geom_point() + geom_line() + labs(x='Number of included samples', y='Mean proportion of stably ranked genes') + geom_errorbar(data=table_, mapping=aes(x=no+1.5, ymin=high, ymax=low), width=2, size=1)  + theme_bw() 
ggplot(table_, aes(x=no+1.5, y=mean)) + geom_point() + geom_line() + labs(x='Number of included cell-types', y='Mean repressive tendency score change') + geom_errorbar(data=table_, mapping=aes(x=no+1.5, ymin=high, ymax=low), width=2, size=1)  + scale_y_continuous(breaks=seq(0, max(table_$mean), by=0.001))+ theme_bw() 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/updated/figures/Supp_Figure8_.png', , width=8, height=5, units='in', dpi=300)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/saturation_prop_stable_genes.pdf', device='pdf', width=8, height=5, units='in', dpi=600)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/saturation_mean_score_changes.pdf', height=5, width=8, units='in', dpi=600)

### SCATTER PLOT (H3K27ME3 WIDTH VS RNA-SEQ)
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_logz.txt', stringsAsFactors = F)
wid = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths_z.txt', stringsAsFactors = F)
cols = colnames(exp)
gene = 'MESP1'
a = exp[gene,cols]
b = wid[gene,cols]
exp_range = range(exp)
wid_range = range(wid)

plot(x=as.numeric(b), y=as.numeric(a), main=gene, xlab='H3K4me3', ylab='log10(RNA-seq)', xlim=wid_range, ylim=exp_range)

### LINEAR REGRESSION MODELS
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_logz.txt', stringsAsFactors = F)
wid = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths_z.txt', stringsAsFactors = F)
cols = colnames(exp)
genes = intersect(rownames(exp), rownames(wid))

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

results = data.frame(matrix(nrow=length(genes), ncol=3))
rownames(results) = genes
colnames(results) = c('intercept','slope','p_value')
for (gene in genes){
  table_ = data.frame(exp=as.numeric(exp[gene,cols]), wid=as.numeric(wid[gene,cols]))
  model_ = lm(wid~exp, data=table_)
  results[gene,1] = model_$coefficients[1]
  results[gene,2] = model_$coefficients[2]
  results[gene,3] = lmp(model_)
}
write.table(results, '/Users/woojunshim/Research/Data/lm_y_h3k27me3_x_exp', quote=F, sep='\t')

### ASSIGN HUMAN H3K27ME3 BED FILES TO THE NEAREST GENES USING ChIPpeakAnno package
library(ChIPpeakAnno)
#bed <- system.file("extdata", "/Users/woojunshim/Research/Data/E003-H3K27me3.broadPeak", package="ChIPpeakAnno")
gr1 = toGRanges("/Users/woojunshim/Research/Data/E003-H3K27me3.broadPeak", format="broadPeak", header=FALSE) 

library(biomaRt)
hg19 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
data(TSS.human.GRCh37)
results_ = annotatePeakInBatch(gr1, AnnotationData = TSS.human.GRCh37)
results = data.frame(chr=as.character(results_@seqnames@values), start=results_@elementMetadata@listData$start_position, end=results_@elementMetadata@listData$end_position)
rownames(results) = results_@ranges@NAMES

### DENSITY PLOT FOR BROAD DOMAINS
temp = read.table('/Users/woojunshim/Research/Data/non_TF_dist.txt', stringsAsFactors = F)
barplot(temp$V2, main='non_TF distribution', xlab='number of broad H3K27me3', ylab='probability')

### ENAKSHI 
days = c('d5')
days = c('d15','d30')
ref_ = read.table('/Users/woojunshim/Research/Data/Clayton/GO/go_list_combined_short.txt', stringsAsFactors = F)
rownames(ref_) =ref_$V1
ref_$V2 = gsub('_', ' ', ref_$V2)
for (d in days){
  dc1 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/discordance/',d,'c1_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  dc2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/discordance/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  dc3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/discordance/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  ec1 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/expression/',d,'c1_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  ec2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/expression/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  ec3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/expression/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  oc1 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c1_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  oc2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  oc3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  oc4 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c4_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  terms = rownames(dc1)
  dc1 = dc1[,order(as.numeric(colnames(dc1)))]
  dc2 = dc2[,order(as.numeric(colnames(dc2)))]
  dc3 = dc3[,order(as.numeric(colnames(dc3)))]
  ec1 = ec1[,order(as.numeric(colnames(ec1)))]
  ec2 = ec2[,order(as.numeric(colnames(ec2)))]
  ec3 = ec3[,order(as.numeric(colnames(ec3)))]
  oc1 = oc1[,order(as.numeric(colnames(oc1)))]
  oc2 = oc2[,order(as.numeric(colnames(oc2)))]
  oc3 = oc3[,order(as.numeric(colnames(oc3)))]
  oc4 = oc4[,order(as.numeric(colnames(oc4)))]
  for (t in terms){
    name_ = ref_[t,2]
    max_ = max(dc1[t,],dc2[t,],dc3[t,],ec1[t,],ec2[t,],ec3[t,],oc1[t,],oc2[t,],oc3[t,],oc4[t,])
    #max_ = max(dc1[t,],dc2[t,],dc3[t,],ec1[t,],ec2[t,],ec3[t,],oc1[t,],oc2[t,],oc3[t,])
    #max_ = max(dc1[t,],dc2[t,],ec1[t,],ec2[t,],oc1[t,],oc2[t,])
    dtable_ = data.frame(matrix(nrow=3, ncol=100))
    etable_ = data.frame(matrix(nrow=3, ncol=100))
    otable_ = data.frame(matrix(nrow=4, ncol=100))
    rownames(dtable_) = c('c1','c2','c3')
    rownames(etable_) = c('c1','c2','c3')
    rownames(otable_) = c('c1','c2','c3','c4')
    colnames(dtable_) = seq(0,99)
    colnames(etable_) = seq(0,99)
    colnames(otable_) = seq(0,99)
    dtable_['c1',] = dc1[t,]
    dtable_['c2',] = dc2[t,]
    dtable_['c3',] = dc3[t,]
    etable_['c1',] = ec1[t,]
    etable_['c2',] = ec2[t,]
    etable_['c3',] = ec3[t,]
    otable_['c1',] = oc1[t,]
    otable_['c2',] = oc2[t,]
    otable_['c3',] = oc3[t,]
    otable_['c4',] = oc4[t,]
    dtable = melt(t(dtable_))
    etable = melt(t(etable_))
    otable = melt(t(otable_))
    ggplot(dtable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
    ggsave(paste('/Users/woojunshim/Research/Data/Clayton/genes/figures/',d,'_',t,'_discordance.png',sep=''), width=8, height=5, units='in', dpi=600)
    ggplot(etable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
    ggsave(paste('/Users/woojunshim/Research/Data/Clayton/genes/figures/',d,'_',t,'_expression.png',sep=''), width=8, height=5, units='in', dpi=600)
    ggplot(otable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
    ggsave(paste('/Users/woojunshim/Research/Data/Clayton/genes/figures/',d,'_',t,'_old.png',sep=''), width=8, height=5, units='in', dpi=600)
  }
}

### PLOT H3K4ME3 AND H3K27ME3 CORRELATION BAR PLOT
temp = read.table('/Users/woojunshim/Research/Data/H3K4me3_H3K27me3_spearman.txt', stringsAsFactors = F)
temp$V4[which(temp$V4==1)] = 'variably expressed TF'
temp$V4[which(temp$V4==2)] = 'variably expressed non-TF'
temp$V4[which(temp$V4==3)] = 'stably expressed TF'
temp$V4[which(temp$V4==4)] = 'stably expressed non-TF'
ggplot(temp, aes(x=V4, y=V2, fill=V4)) + geom_boxplot() + labs(main='Correlation between H3K4me3 nad H3K27me3 widths', x='', y="Spearman's rho", fill='') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
a = temp$V2[which(temp$V4=='variably expressed TF')]
b = temp$V2[which(temp$V4=='variably expressed non-TF')]
c = temp$V2[which(temp$V4=='stably expressed TF')]
d = temp$V2[which(temp$V4=='stably expressed non-TF')]
ggsave('/Users/woojunshim/Research/BN/figures/H3K4me3_H3K27me3_spearman.png', width=8, height=5, units='in', dpi=600)

### JACCARD INDEX BOX PLOT (H3K4ME3 TOP 5% VS H3K27ME3-POSITIVE TOP 5%)
temp = read.table('/Users/woojunshim/Research/Data/H3K4me3_H3K27me3-positive_jaccard.txt', stringsAsFactors = F)
ggplot(temp, aes(x=factor(1), y=V2)) +geom_boxplot() + labs(y='Jaccard Index', x='')
ggsave('/Users/woojunshim/Research/BN/figures/H3K4me3_H3K27me3-positive_jaccard.png', width=3, height=5, units='in', dpi=600)


### COMPETITION
days = c('d5')
days = c('d15','d30')
days = ('LT-HSC_')
ref_ = read.table('/Users/woojunshim/Research/Data/competition/go_list_combined.txt', stringsAsFactors = F)
rownames(ref_) =ref_$V1
ref_$V2 = gsub('_', ' ', ref_$V2)
for (d in days){
  dc1 = -log10(read.table(paste('/Users/woojunshim/Research/Data/competition/',d,'discord_fet_fdr_by_discordance.txt', sep=''), stringsAsFactors = F, check.names = F))
#  dc2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/discordance/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  dc3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/discordance/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  dc2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/competition/',d,'exp_fet_fdr_by_expression.txt', sep=''), stringsAsFactors = F, check.names = F))
#  ec2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/expression/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  ec3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/expression/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  oc1 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c1_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  oc2 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c2_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  oc3 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c3_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
#  oc4 = -log10(read.table(paste('/Users/woojunshim/Research/Data/Clayton/genes/old/',d,'c4_fet_fdr.txt', sep=''), stringsAsFactors = F, check.names = F))
  terms = rownames(dc1)
  dc1 = dc1[,order(as.numeric(colnames(dc1)))]
#  dc2 = dc2[,order(as.numeric(colnames(dc2)))]
#  dc3 = dc3[,order(as.numeric(colnames(dc3)))]
  dc2 = dc2[,order(as.numeric(colnames(dc2)))]
#  ec2 = ec2[,order(as.numeric(colnames(ec2)))]
#  ec3 = ec3[,order(as.numeric(colnames(ec3)))]
#  oc1 = oc1[,order(as.numeric(colnames(oc1)))]
#  oc2 = oc2[,order(as.numeric(colnames(oc2)))]
#  oc3 = oc3[,order(as.numeric(colnames(oc3)))]
#  oc4 = oc4[,order(as.numeric(colnames(oc4)))]
  for (t in terms){
    name_ = ref_[t,2]
    max_ = max(dc1[t,],dc2[t,])
    #max_ = max(dc1[t,],dc2[t,],dc3[t,],ec1[t,],ec2[t,],ec3[t,],oc1[t,],oc2[t,],oc3[t,])
    #max_ = max(dc1[t,],dc2[t,],ec1[t,],ec2[t,],oc1[t,],oc2[t,])
    dtable_ = data.frame(matrix(nrow=2, ncol=100))
#    etable_ = data.frame(matrix(nrow=1, ncol=100))
#    otable_ = data.frame(matrix(nrow=4, ncol=100))
    rownames(dtable_) = c('discordance', 'expression')
#    rownames(etable_) = c('expression')
#    rownames(otable_) = c('c1','c2','c3','c4')
    colnames(dtable_) = seq(0,99)
#    colnames(etable_) = seq(0,99)
#    colnames(otable_) = seq(0,99)
    dtable_['discordance',] = dc1[t,]
#    dtable_['c2',] = dc2[t,]
#    dtable_['c3',] = dc3[t,]
    dtable_['expression',] = dc2[t,]
#    etable_['c2',] = ec2[t,]
#    etable_['c3',] = ec3[t,]
#    otable_['c1',] = oc1[t,]
#    otable_['c2',] = oc2[t,]
#    otable_['c3',] = oc3[t,]
#    otable_['c4',] = oc4[t,]
    dtable = melt(t(dtable_))
#    etable = melt(t(etable_))
#    otable = melt(t(otable_))
    ggplot(dtable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
    ggsave(paste('/Users/woojunshim/Research/Data/competition/',d,'_',t,'_discordance.png',sep=''), width=8, height=5, units='in', dpi=600)
#    ggplot(etable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
#    ggsave(paste('/Users/woojunshim/Research/Data/competition/',d,'_',t,'_expression.png',sep=''), width=8, height=5, units='in', dpi=600)
#    ggplot(otable, aes(x=as.numeric(Var1)+1, y=value, colour=Var2)) + geom_point() + geom_line() + labs(title = name_, x='Rank position', y='-log10(FDR)', colour='') + ylim(0, max_) + theme_bw()
#    ggsave(paste('/Users/woojunshim/Research/Data/Clayton/genes/figures/',d,'_',t,'_old.png',sep=''), width=8, height=5, units='in', dpi=600)
  }
}

### ENRICHMENT PLOTS OF DIFFERENT HM LATENCY SCORES
### ONE PLOT FOR EACH HM
ymax_ = 37.7855
hm = 'H3K27ac'
hms = c('H3K4me1','H3K4me3','H3K9me3','H3K27ac','H3K27me3','H3K36me3')
temp1 = read.table(paste('/Users/woojunshim/Research/Data/latency_tables/variably_expressed_tf/discordance_fet_',hm,'.txt',sep=''), stringsAsFactors = F)
rownames(temp1) = temp1$V1
temp1 = temp1[2:ncol(temp1)]
colnames(temp1) = seq(1, 100)
temp1 = -log10(temp1)
mean1_ = apply(temp1, 2, mean)
sd1_ = apply(temp1, 2, sd)
sem1_low = mean1_ - (sd1_ / sqrt(nrow(temp1)))
sem1_high = mean1_ + (sd1_ / sqrt(nrow(temp1)))
table1_ = data.frame(matrix(nrow=100, ncol=5))
colnames(table1_) = c('mean','low','high','method','rank')
table1_$mean = mean1_
table1_$low = sem1_low
table1_$high = sem1_high
table1_$method = rep('Discordance', 100)
table1_$rank = seq(1, 100)

temp2 = read.table(paste('/Users/woojunshim/Research/Data/latency_tables/variably_expressed_tf/expression_fet_',hm,'.txt',sep=''), stringsAsFactors = F)
rownames(temp2) = temp2$V1
temp2 = temp2[2:ncol(temp2)]
colnames(temp2) = seq(1, 100)
temp2 = -log10(temp2)
mean2_ = apply(temp2, 2, mean)
sd2_ = apply(temp2, 2, sd)
sem2_low = mean2_ - (sd2_ / sqrt(nrow(temp2)))
sem2_high = mean2_ + (sd2_ / sqrt(nrow(temp2)))
table2_ = data.frame(matrix(nrow=100, ncol=5))
colnames(table2_) = c('mean','low','high','method','rank')
table2_$mean = mean2_
table2_$low = sem2_low
table2_$high = sem2_high
table2_$method = rep('Expression', 100)
table2_$rank = seq(1, 100)

table_ = rbind(table1_, table2_)
ggplot(table_, aes(x=rank, y=mean, colour=method)) + ylim(0, ymax_) + geom_line() + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black') + labs(y='-log10(p-value)', x='Rank position', colour='', title=hm) + theme_bw()
ggsave(paste('/Users/woojunshim/Research/Data/latency_tables/variably_expressed_tf/', hm,'.png'), width=8, height=5, units='in', dpi=600)

### FIND YMAX
hms = c('H3K4me1','H3K4me3','H3K9me3','H3K27ac','H3K27me3','H3K36me3')
pathway = '/Users/woojunshim/Research/Data/latency_tables/variably_expressed_tf/'
ymax_ = 0
for (hm in hms){
  temp = read.table(paste(pathway,'discordance_fet_',hm,'.txt', sep=''), stringsAsFactors = F)
  pp = max(-log10(temp[,2:ncol(temp)]))
  if (pp > ymax_){
    ymax_ = pp
  }
}


### ENRICHMENT PLOTS OF DIFFERENT HM LATENCY SCORES
### ALL LINES IN ONE PLOT
hm = 'H3K27ac'
hms = c('H3K4me1','H3K4me3','H3K9me3','H3K27ac','H3K27me3','H3K36me3')
colours = rainbow(length(hms)+1)
pathway = '/Users/woojunshim/Research/Data/latency_tables/variably_expressed_tf/'

temp2 = read.table(paste(pathway,'expression_fet_H3K27me3.txt',sep=''), stringsAsFactors = F)
rownames(temp2) = temp2$V1
temp2 = temp2[2:ncol(temp2)]
colnames(temp2) = seq(1, 100)
temp2 = -log10(temp2)
mean2_ = apply(temp2, 2, mean)
sd2_ = apply(temp2, 2, sd)
sem2_low = mean2_ - (sd2_ / sqrt(nrow(temp2)))
sem2_high = mean2_ + (sd2_ / sqrt(nrow(temp2)))
table2_ = data.frame(matrix(nrow=100, ncol=5))
colnames(table2_) = c('mean','low','high','method','rank')
table2_$mean = mean2_
table2_$low = sem2_low
table2_$high = sem2_high
table2_$method = rep('Expression', 100)
table2_$rank = seq(1, 100)
table_ = table2_
for (hm in hms){
  temp1 = read.table(paste(pathway,'discordance_fet_',hm,'.txt',sep=''), stringsAsFactors = F)
  rownames(temp1) = temp1$V1
  temp1 = temp1[2:ncol(temp1)]
  colnames(temp1) = seq(1, 100)
  temp1 = -log10(temp1)
  mean1_ = apply(temp1, 2, mean)
  sd1_ = apply(temp1, 2, sd)
  sem1_low = mean1_ - (sd1_ / sqrt(nrow(temp1)))
  sem1_high = mean1_ + (sd1_ / sqrt(nrow(temp1)))
  table1_ = data.frame(matrix(nrow=100, ncol=5))
  colnames(table1_) = c('mean','low','high','method','rank')
  table1_$mean = mean1_
  table1_$low = sem1_low
  table1_$high = sem1_high
  table1_$method = rep(hm, 100)
  table1_$rank = seq(1, 100)
  table_ = rbind(table_, table1_)
}
pdf(paste(pathway,'enrichment.pdf',sep=''), width=10, height=8)
ggplot(table_, aes(x=rank, y=mean, colour=method)) + ylim(0, ymax_) + geom_line() + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black') + labs(y='-log10(p-value)', x='Rank position', colour='') + theme_bw()
dev.off()
ggsave(paste(pathway,'enrichment.png', sep=''), width=8, height=5, units='in', dpi=600)

### PLOT HISTOGRAM OF PAIR-WISE H3K27ME3 DISSIMILARITIES 
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_negative_broad_background_specificity_roadmap.txt', stringsAsFactors = F)
temp1 = temp[temp!=0]
pdf('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/figures/Dissimilarity_between_cell_types_H3k27me3.pdf', width=8, height=5)
pdf('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/figures/figure1.pdf', width=8, height=5)
hist(temp1, breaks=1000, xlab='Specificity score', main='Absence of H3K27me3', xlim=c(0.95, 1.05))

dev.off()

### FIND THE NEAREST GENE AND THE DISTANCE 
library(ChIPpeakAnno)
data(TSS.human.GRCh37)

### 
library(topGO)
library(GO.db)
# Create GO terms to gene ID mapping table
all_genes = read.table('/Users/woojunshim/Research/Data/mRNA_genes.txt', stringsAsFactors = F)
all_genes = all_genes$V1
top_no = 100 # top 100 GO terms by FDR
go_table = read.table('/Users/woojunshim/Research/Data/GO_Terms/GO_table_.txt')
groups = c('E083','E081','E070','E095')
for (no in 1:length(groups)){
  group = groups[no]
  temp = read.table(paste('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/',group,'_significant_genes.txt',sep=''), stringsAsFactors = F)
  genelist = temp$V1
  geneList<-factor(as.integer(all_genes %in% genelist))
  names(geneList)<-all_genes
  TopGOdata<-new("topGOdata",ontology="BP",allGenes=geneList,geneSel=genelist,nodeSize=5,
                 annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = "fisher")
  aa = score(resultFisher)
  aa = p.adjust(aa, method='fdr')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),1])
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','FDR')
  file_ = paste('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/',group,'_GO_BP.txt', sep='')
  write.table(results, file_, sep='\t', quote=F)
}

### STACKED BAR GRAPH 
library(ggplot2)
library(reshape2)
table_ = as.data.frame(matrix(ncol=3, nrow=3))
neg = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_broad_domains_significant_stat_negatives.txt', stringsAsFactors = F)
table_[1,] = neg$V1
pos = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_broad_domains_significant_stat_positives.txt', stringsAsFactors = F)
table_[2,] = pos$V1
tot = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_broad_domains_significant_stat_total.txt', stringsAsFactors = F)
table_[3,] = tot$V1
colnames(table_) = c('promoter','intragenic','intergenic')
rownames(table_) = c('Low','High','All')
table = melt(t(table_))
ggplot(table, aes(x=Var2,y=value, fill=Var1)) + geom_bar(stat='identity') + labs(fill='', x='', y='') + theme_bw()
ggsave('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/figures/stacked_bar_H3K27me3_absence_negatives_positives.png', width=4, height=5, units='in', dpi=600)

### BOX PLOT FOR DISTANCES TO THE NEAREST GENES (NEGATIVE VS POSITIVE)
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_broad_domains_significant_distances.txt', stringsAsFactors = F)
neg = as.numeric(temp[,1][temp$V2=='negative'])
pos = as.numeric(temp[,1][temp$V2=='positive'])
ggplot(temp, aes(x=V2, y=V1)) + geom_boxplot()

### HEATMAP FOR DISSIMILARITY BROAD H3K27ME3 DOMAINS
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_gene_dissimilarity.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/RNA-seq/spearman_dissimilarity_roadmap_rpkm.txt', stringsAsFactors = F)
table_ = melt(t(temp))
cell_order = read.table('/Users/woojunshim/Research/cell-type-specificity/data/roadmap_cell_type_order1.txt', stringsAsFactors = F)
cell_order = read.table('/Users/woojunshim/Research/cell-type-specificity/data/roadmap_cell_type_order.txt', stringsAsFactors = F)
cell_order = cell_order$V1
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(guide='colourbar',midpoint = (max(as.matrix(temp)) + min(as.matrix(temp))) /2, low = "blue", mid='white', high = "red") + labs(x='', y='', fill="Spearman's dissimilarity", title='H3K27me3', cex=0.8) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_discrete(limits=cell_order) + scale_x_discrete(limits=cell_order)
ggsave('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_gene_dissimilarity.png', width=11, height=10, units='in', dpi=600)
       
### PLOT DISTRIBUTION OF RANDOMLY SAMPLES H3K27ME3 DISSIMILARITY 
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/test/100_H3K27me3_random_samples.txt', stringsAsFactors = F)
samples = temp$V1
pdf('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/test/100_H3K27me3_random_samples.pdf', width=8, height=5)
hist(samples, breaks=100, xlab='Dissimilarity', main='Permuted samples (m=100, n=1000)')
dev.off()

### HISTOGRAM OF EMPIRICAL H3K27ME3 
temp = read.table('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/H3K27me3_gene_empirical_background.txt', stringsAsFactors = F)
values = temp$V1
pdf('/Users/woojunshim/Research/cell-type-specificity/data/H3K27me3/figures/H3K27me3_gene_empirical_background.pdf', width=5, height=5)
hist(values, breaks=100, prob=T, main='', xlab='Dissimilarity')
dev.off()

### SATURATION PLOTS
### FORFIVE DIFFERENT THRESHOLDS
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes_top1.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table1_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('< 1%', length(high_)))

temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes_top2.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table2_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('< 2%', length(high_)))

temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes_top3.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table3_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('< 3%', length(high_)))

temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes_top4.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table4_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('< 4%', length(high_)))

temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/saturation_prop_stable_genes_top5.txt', stringsAsFactors = F)
temp = t(temp)
rownames(temp) = 3*seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(999)
high_ = results + 1.96 * sd_ / sqrt(999)
table5_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('< 5%', length(high_)))

table_ = rbind(table1_, table2_, table3_, table4_, table5_)
ggplot(table_, aes(x=no+1.5, y=mean, colour=label)) + geom_point() + geom_line() + labs(x='Number of included cell-types', y='Mean proportion of stably ranked genes', colour='Rank change') + geom_errorbar(data=table_, mapping=aes(x=no+1.5, ymin=high, ymax=low), width=2, size=1)  + scale_y_continuous(breaks=seq(0.1, 1.0, by=0.1))+ theme_bw() 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/saturation_prop_stable_genes.pdf', width=10, height=8, units='in', dpi=600)

### BOOTSTRAPPING PEARSON 
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/bootstrapping_pearson.txt', stringsAsFactors = F)

### SUMA PLOT
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/suma_roadmap.txt', stringsAsFactors = F)
temp = temp[order(temp$V2, decreasing = T),]

### BOXPLOTS FOR MEAN ASSOCIATION SCORES 
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3/H3K27me3_widths_mean.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V2, y=V1, group=V2)) + geom_boxplot() +  labs(x='Rank position', y='Mean association score') + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/H3K27me3_widths_mean.pdf', width=10, height=8, units='in', dpi=600)

### SIMILARITY PLOTS (JACCARD)
epi = 'E038'
#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K4me1.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K4me1.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K4me1_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table1_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K4me1', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K4me3.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K4me3.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K4me3_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table2_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K4me3', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K9me3.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K9me3.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K9me3_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table3_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K9me3', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K27ac.txt', sep=''), stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K27ac.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27ac_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table4_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K27ac', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K27me3.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K27me3.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table5_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K27me3', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_foreground_jaccard_H3K36me3.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K36me3.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K36me3_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/H3K27me3_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table6_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('H3K36me3', length(high_)))

#temp = read.table(paste('/Users/woojunshim/Research/Data/broadPeaks/similarity_analysis/',epi,'_background_similarity.txt', sep=''), stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_random.txt', stringsAsFactors = F)
#temp = read.table('/Users/woojunshim/Desktop/enrichment/Expression_fet_house_keeping_genes.txt')
temp = read.table('/Users/woojunshim/Desktop/enrichment/Expression_fet_variably_expressed_tf.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(111)
high_ = results + 1.96 * sd_ / sqrt(111)
table7_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('Random', length(high_)))

#col_ = rainbow(7)
col_ = c('red','yellow1','green1','cyan','blue','#FF00DBFF','darkgray')
#table_ = rbind(table1_, table2_, table3_, table4_, table5_, table6_, table7_)
table_ = rbind(table1_, table2_, table3_, table4_, table5_, table6_)
#ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line(show.legend = F, size=1) + labs(x='', y='', colour='') +  scale_colour_manual(values=col_) + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black')+ scale_x_discrete(limits=seq(1,10))+theme_bw() +theme(axis.text=element_text(size=16,face="bold")) 
#ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line() + labs(x='Rank position', y='Mean Jaccard similarity index', colour='') +  scale_colour_manual(values=col_)  + theme_bw() + theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=16,face='bold'), legend.text=element_text(size=16,face="bold"), legend.position='top')
ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line(size=2) + labs(x='Percentile rank', y='-log10(p-value)', colour='')  + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black') + scale_colour_manual(values=col_[1:6]) + theme_bw()+theme(legend.text=element_text(size=16,face="bold"), legend.position='top', axis.title=element_text(size=16, face='bold'), axis.text=element_text(siz=16, face='bold'))
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/mean_jaccard_similarity_top1000_genes_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/mean_jaccard_similarity_all_genes_new_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/fet_house_keeping_genes_new_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/fet_variably_expressed_tf_new_new_.pdf',sep=''), width=8, height=5, units='in', dpi=600)

### GO ANALYSIS ON TOP 100 COMMON ELEMENTS 
temp = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K4me3_top5_prop.txt', stringsAsFactors = F)
temp = temp[order(temp$V2, decreasing = T),]
results = go_analysis(temp$V1[1:100], temp$V1)

### BOX PLOT FOR JACCARD SIMILARITY INDEX OF TOP 10000 GENES (H3K27ME3)
hm = 'H3K27ac'
temp = read.table(paste('/Users/woojunshim/Research/Data/new_assigned/jaccard_10000_genes_by_',hm,'.txt', sep=''), stringsAsFactors = F)
temp = temp[,2:ncol(temp)]
table_ = melt(t(temp))
ggplot(table_, aes(x=Var2, y=value, group=Var2)) + geom_boxplot() + labs(x='Rank position', y='Jaccard similarity index', title=hm) + theme_bw()
ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/jaccard_10000_genes_by_',hm,'.pdf', sep=''), width=10, height=8, units='in', dpi=600)



### CALCULATE CORRELATION BETWEEN EXP AND H3K27ME3 
a = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_z.txt',stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_z.txt',stringsAsFactors = F)
genes = intersect(rownames(a), rownames(b))
a = subset(a, select=colnames(b))
a = a[genes,]
b = b[genes,]
results = vector()
for (i in 1:nrow(a)){
  q1 = as.numeric(a[i,])
  q2 = as.numeric(b[i,])
  r = cor(q1,q2, method='spearman')
  results = c(results, r)
}
output_ = data.frame(matrix(nrow=nrow(a), ncol=2))
output_$X1 = genes
output_$X2 = results
write.table(output_, '/Users/woojunshim/Research/Data/new_assigned/spearman_K27_exp.txt', sep='\t', quote=F)

### FET ENRICHMENT
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/HOX_genes_rank.txt', stringsAsFactors = F)
temp = temp[,2:ncol(temp)]
rownames(temp) = c('New+MAD','New+Old','Old')
colnames(temp) = seq(1, ncol(temp))
table_ = melt(t(temp))
table_$value = -log10(table_$value)

ggplot(table_, aes(x=Var1, y=value, group=Var2, colour=Var2)) + geom_line() + geom_point() + labs(x='Rank', y='-log10(p-value)', colour='') + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/GO_0030017_genes_enrichment.pdf', height=10, width=10, dpi=600)

ggplot(table_, aes(x=Var2, y=value, group=Var2, colour=Var2)) + geom_boxplot() + labs(x='', y='Rank', colour='') + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/HOX_genes_rank.pdf', height=10, width=10, dpi=600)


### TOP 100 or 200 COMMONLY ASSOCIATED GENES 
a1 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K4me1_top5_prop_all_genes.txt', stringsAsFactors = F)
a2 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K4me3_top5_prop_all_genes.txt', stringsAsFactors = F)
a3 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K9me3_top5_prop_all_genes.txt', stringsAsFactors = F)
a4 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27ac_top5_prop_all_genes.txt', stringsAsFactors = F)
a5 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_top5_prop_all_genes.txt', stringsAsFactors = F)
a6 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K36me3_top5_prop_all_genes.txt', stringsAsFactors = F)
a1 = a1[order(a1$V2, decreasing = T),]
a2 = a2[order(a2$V2, decreasing = T),]
a3 = a3[order(a3$V2, decreasing = T),]
a4 = a4[order(a4$V2, decreasing = T),]
a5 = a5[order(a5$V2, decreasing = T),]
a6 = a6[order(a6$V2, decreasing = T),]
table_ = data.frame(H3K4me1 = a1[1:200, 2],H3K4me3 = a2[1:200, 2],H3K9me3 = a3[1:200, 2],H3K27ac = a4[1:200, 2],H3K27me3 = a5[1:200, 2],H3K36me3 = a6[1:200, 2])
rownames(table_) = seq(1,200)
table = melt(t(table_))
ggplot(table, aes(x=Var1, y=value, colour=Var1)) + geom_boxplot() + labs(x='', y='Proportion of cell-types with the broad peak', colour='', title='') + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), legend.position='top', axis.title=element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(angle = 45, hjust = 1,size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/top200_common_genes_.pdf', height=6, width=6, dpi=600)

temp = read.table('/Users/woojunshim/Research/Data/new_assigned/top200_common_genes_jaccard.txt', stringsAsFactors = F)
table = melt(t(temp))
ggplot(table, aes(x=Var1, y=Var2, fill=value, label=round(value, 3))) + geom_tile(colour='black') + scale_fill_gradient(low = "white", high = "steelblue") + labs(x='', y='', title='', fill='Jaccard similarity')  + theme_bw()+ theme(axis.text=element_text(size=16, face='bold'), axis.text.x=element_text(angle = 45, hjust = 1,size=16, face='bold'), legend.text=element_text(size=14, face='bold'), legend.title=element_text(size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/top200_common_genes_jaccard_new.pdf', height=5, width=8, dpi=600)

temp = read.table('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/manuscript/script/day15_.txt', stringsAsFactors = F)
all = temp$V1
results = go_analysis(input=a6$V1[1:200], all=all)
write.table(results, '/Users/woojunshim/Research/Data/new_assigned/H3K36me3_top200_GO.txt', sep='\t', quote=F)

temp = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_top200_selected_GO_all.txt', stringsAsFactors = F)
t = gsub('_', ' ', rownames(temp))
rownames(temp) = t
table = melt(t(temp))
ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low = "white", high = "steelblue")  + labs(x='', y='', title='', fill='-log10(P)') + theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12, face='bold'), axis.text.y=element_text(size=12, face='bold'),legend.title=element_text(size=12, face='bold'), legend.text=element_text(size=10))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/H3K27me3_top200_selected_GO_all_.pdf', height=6, width=6, dpi=600)

### PERFORMANCE PLOTS
#pathway = '/Users/woojunshim/Research/Data/NCC/new/'
pathway = '/Users/woojunshim/Research/Data/Palpant/new/'
#pathway = '/Users/woojunshim/Research/Data/Paige/new/'
x_ = read.table(paste(pathway, 'roc_x_.txt', sep=''), stringsAsFactors = F)
y_ = read.table(paste(pathway, 'roc_y_.txt', sep=''), stringsAsFactors = F)
uu_ = read.table(paste(pathway, 'auc.txt', sep=''), stringsAsFactors = F)
auc_ = as.numeric(uu_$V2)
names(auc_) = uu_$V1
plot_performance_roc1(x_, y_, 'FPR', 'TPR', 'ROC', auc_)

ggsave(paste(pathway,'roc.pdf',sep=''), width=8, height=5, dpi=600)

x_ = read.table(paste(pathway, 'prc_x_.txt', sep=''), stringsAsFactors = F)
y_ = read.table(paste(pathway, 'prc_y_.txt', sep=''), stringsAsFactors = F)
#col = colnames(x_)[5]
#x_ = subset(x_, select=col)
#y_ = subset(y_, select=col)
plot_performance_prc(x_, y_, 'Recall', 'Precision', 'PRC', 386)
ggsave(paste(pathway,'prc.pdf',sep=''), width=8, height=5, dpi=600)
ggsave(paste(pathway,'prc_',col,'.pdf',sep=''), width=8, height=5, dpi=600)

### HEATMAP FOR SELECTED CARDIAC GENES
cols = c('E095','E104','E105','E082','E071','E070','E059','E058','E057','E047','E038','E037','E024','E016','E003','E006','E005','E004')
rows = c('TNNI3','MYL3','MYL2','MYH7','MYH6','TBX20','TBX5','NKX2-5','GATA6','GATA4')
mark = 'H3K27ac'
#input_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
input_ = read.table(paste('/Users/woojunshim/Research/Data/new_assigned/',mark,'_tss_2.5kb_combined.txt', sep=''), stringsAsFactors = F)
#input_ = log10(input_+1)
tt = c('E095','E104','E105','E071','E059','E058','E047','E038','E037','E016','E003','E006','E005','E004')
x_labs = c(rep('Heart',3),rep('Brain',1),rep('Epithelial',2),rep('Blood',3),rep('ES cell',2),rep('ES-deriv.',3))
cols = tt
input_ = input_[rows, cols]
table_ = melt(t(input_))
#x_labs = c(rep('Heart',3),rep('Brain',3),rep('Epithelial',3),rep('Blood',3),rep('ES cell',3),rep('ES-deriv.',3))
#ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_x_discrete(limits=cols) + scale_y_discrete(limits=rows) + scale_fill_gradient2(low='white', high='red') + labs(x='', y='', title=mark, fill='breadth') + theme_bw()+theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y=element_text(size=16, face='bold'), axis.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=14), plot.title=element_text(size=18, face='bold'), legend.position = 'top')
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_x_discrete(limits=cols, labels=x_labs) + scale_y_discrete(limits=rows) + scale_fill_gradient2(low='white', high='red') + labs(x='', y='', title=mark, fill='breadth') + theme_bw()+theme(axis.ticks.x=element_blank(), axis.text.y=element_text(size=16, face='bold'), axis.text.x=element_text(size=16, face='bold', angle = 90, hjust = 1), axis.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=14), plot.title=element_text(size=18, face='bold'))
#ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/Exp_cardiac_new.pdf', width=8, height=5, dpi=600)
ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/',mark,'_breadth_cardiac_new.pdf',sep=''), width=8, height=5, dpi=600)

### HEATMAP FOR VARIABLY EXPRESSED GENES
library(RColorBrewer)
library(gplots)
col = colorRampPalette(c("green","black", "red"))(n = 100)
col = colorRampPalette(c("green","white","red"))(n=100)
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/Exp_gene_centred.txt', stringsAsFactors = F)
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/Exp_gene_centred_new.pdf', height=8, width=8)
heatmap.2(as.matrix(temp), col=col, trace='none', labRow='', labCol='')
dev.off()

### PLOT MAGNIFIED VIEW DOWN TO RANK POSITION 25 (I.E. 2500 GENES FROM EACH CELL-TYPE)
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K27me3.txt', stringsAsFactors = F)
temp = temp[1:25,]

### plot expressed TF with heart development term (E095)
### highly expressed TF vs discordance score
temp = read.table('/Volumes/OMICS2018-A1171/WorkingSpace/WooJun/E082_exp_discordance_positive_test_tf.txt',stringsAsFactors = F)
rownames(temp) = temp[,1]
temp = temp[,2:ncol(temp)]
wilcox.test(as.numeric(temp[1,]), as.numeric(temp[2,]))
table_ = melt(t(temp))
ggplot(table_, aes(x=Var2, y=value, colour=Var2, group=Var2)) + geom_boxplot() + labs(x='', y='Rank', title='TFs with brain development (GO:0007420)', colour='') + scale_y_reverse() + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/TF_rank_with_brain_development_E082.pdf', width=7, height=5, dpi=600)

plot(1, type="n", xlab="Rank", ylab="", xlim=c(0, 1422), ylim=c(0, 2), yaxt = 'n')
points(x=as.numeric(tt[1,]), y=rep(1,ncol(tt)), col='blue')
points(x=as.numeric(tt[1,]), y=rep(1,ncol(tt)), col='blue')

### PLOT ROC 
cols = c('Discordance','H3K4me3','DE','Expression')  # order of plotting

### GENERATE INDIVIDUAL PLOTS FOR EACH GENE AND EACH HM (PEAKS SURROUNDING THE TSS)
hm = c('H3K4me1','H3K4me3','H3K9me3','H3K27ac','H3K27me3','H3K36me3')
genes = c('NKX2-5','MYH7')
pathway = '/Users/woojunshim/Research/Data/width_plots/'
for (h in hm){
  for (g in genes){
    file_ = paste('/Users/woojunshim/Research/Data/width_plots/',g,'/',g,'_',h,'_coordinates.txt', sep='')
    temp = read.table(file_, stringsAsFactors = F)
    temp$y = seq(1, nrow(temp))
    pdf(paste('/Users/woojunshim/Research/Data/width_plots/',g,'/',g,'_',h,'_2.5kb_upstream_25kb_downstream_.pdf', sep=''), height=5, width=2)
    empty_plot(c(0,27500), c(1,111))
    horizontal_lines(temp, 2,3,5, colour='dark blue')
    dev.off()
  }
}

### Exp_gene_centred_variably_expressed_tf.pdf

### PERFORMANCE PLOTS
pathway = '/Users/woojunshim/Research/Data/roadmap_test/'
epi = c('E095','E070','E096','E098','E100')
names = c('Left ventricle','Brain germinal matrix','Lung','Pancreas','Psoas muscle')
for (i in 1:length(epi)){
  x_ = read.table(paste(pathway,epi[i],'/roc_x.txt',sep=''), stringsAsFactors = F)
  y_ = read.table(paste(pathway,epi[i],'/roc_y.txt',sep=''), stringsAsFactors = F)
  x = add_roc_random(x_)
  y = add_roc_random(y_)
  plot_performance_roc(x, y, xlab='FPR',ylab='TPR', title=names[i])
  ggsave(paste(pathway,epi[i],'/ROC','_',epi[i],'_.pdf',sep=''), width=8, height=5, dpi=600)
  
  x_ = read.table(paste(pathway,epi[i],'/prc_x.txt',sep=''), stringsAsFactors = F)
  y_ = read.table(paste(pathway,epi[i],'/prc_y.txt',sep=''), stringsAsFactors = F)
  temp = read.table(paste(pathway,epi[i],'/positives.txt',sep=''), stringsAsFactors = F)
  no_ = nrow(temp)
  plot_performance_prc(x_, y_, xlab='Recall',ylab='Precision',title=names[i], no_positives=no_)
  ggsave(paste(pathway,epi[i],'/PRC','_',epi[i],'.pdf',sep=''), width=8, height=5, dpi=600)
}

pathway = '/Users/woojunshim/Research/Data/Palpant/'
name = 'Cardiac progenitor cells'
x_ = read.table(paste(pathway,'roc_x.txt',sep=''), stringsAsFactors = F)
y_ = read.table(paste(pathway,'roc_y.txt',sep=''), stringsAsFactors = F)
x = add_roc_random(x_)
y = add_roc_random(y_)
plot_performance_roc(x, y, xlab='FPR',ylab='TPR', title=name)
ggsave(paste(pathway,'ROC_',name,'.pdf',sep=''), width=8, height=5, dpi=600)
  
x_ = read.table(paste(pathway,'prc_x.txt',sep=''), stringsAsFactors = F)
y_ = read.table(paste(pathway,'prc_y.txt',sep=''), stringsAsFactors = F)
temp = read.table(paste(pathway,'positives.txt',sep=''), stringsAsFactors = F)
no_ = nrow(temp)
plot_performance_prc(x_, y_, xlab='Recall',ylab='Precision',title=name, no_positives=no_)
ggsave(paste(pathway,'PRC_',name,'.pdf',sep=''), width=8, height=5, dpi=600)


pathway = '/Users/woojunshim/Research/Data/roadmap_test/'
epi = c('E095','E070','E096','E098','E100')
names = c('Left ventricle','Brain germinal matrix','Lung','Pancreas','Psoas muscle')
for (i in 1:length(epi)){
  x_ = read.table(paste(pathway,epi[i],'/roc_x.txt',sep=''), stringsAsFactors = F)
  y_ = read.table(paste(pathway,epi[i],'/roc_y.txt',sep=''), stringsAsFactors = F)
  auc_ = read.table(paste(pathway,epi[i],'/auc_roc.txt',sep=''), stringsAsFactors = F )
  auc = auc_$V2
  names(auc) = auc_$V1
  pdf(paste(pathway,epi[i],'/ROC','_',epi[i],'.pdf',sep=''), width=6, height=5)
  plot_performance_roc1(x_, y_, xlab='FPR',ylab='TPR', title=names[i], auc)
  dev.off()
  
  x_ = read.table(paste(pathway,epi[i],'/prc_x.txt',sep=''), stringsAsFactors = F)
  y_ = read.table(paste(pathway,epi[i],'/prc_y.txt',sep=''), stringsAsFactors = F)
  temp = read.table(paste(pathway,epi[i],'/positives.txt',sep=''), stringsAsFactors = F)
  no_ = nrow(temp)
  plot_performance_prc(x_, y_, xlab='Recall',ylab='Precision',title=names[i], no_positives=no_)
  ggsave(paste(pathway,epi[i],'/PRC','_',epi[i],'.pdf',sep=''), width=8, height=5, dpi=600)
  
}

### HEATMAP FOR FIGURE2
temp = read.table('/Users/woojunshim/Research/Data/GTEX/GTEX_tissue_exp_subset.txt', stringsAsFactors = F)
o = c("Muscle_._Skeletal","Heart_._Left_Ventricle","Lung","Brain_._Cortex","Pancreas")
temp = temp[,o]
q=gsub('_._',' ',colnames(temp))
colnames(temp) = q
rows= c('MYOD1','TNNI1','NKX2-5','TNNI3','TBX4','SFTPB','NEUROD2','GFAP','PTF1A','REG1B','GAPDH','EEF2')
heat_map1(temp, rows, q, 'ln(Expression)', fill='z-score', map_col='blue')
#heat_map1(temp, rows, cols, 'Discordance', fill='z-score', map_col='blue')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/GTEX_expression__.pdf', width=5, height=5, dpi=600)




### INTER-SPECIES FET
pathway = '/Users/woojunshim/Research/Data/interspecies/'
pathway = '/Users/woojunshim/Research/Data/ciona/'
pathway = '//Users/woojunshim/Research/Data/Clayton/'
menu = c('Sus scrofa','Mus musculus','Homo sapiens','Galus Galus','Danio rerio','Ciona intestinalis','Cavia porcellus')
menu = c('Ciona intestinalis')
menu = c('GO.0048738','GO.0007507','GO.0030017','GO.0072358')
for (n in (1:length(menu))){
  file_ = paste(pathway,'d30c2_',menu[n],'_fet.txt',sep='')
  temp = read.table(file_, stringsAsFactors = F)
  temp = -log10(temp)
  plot_fet(temp, 'Rank bin position','-log10(p-value)',title=menu[n])
  ggsave(paste(pathway,menu[n],'_fet.pdf', sep=''),width=7,height=5,dpi=600)
}

### WILCOXON RANK-SUM TEST (BOVINE DATASET)
temp = read.table('/Users/woojunshim/Research/bovine/bovine_discordance_score_updated__.txt', stringsAsFactors = F)
results = vector()
negative = c('NE001407','NE001572','NE001634','NE001735','NE001887','NE001981')
positive = c('NE001455','NE001513','NE001636','NE001783','NE001831','NE003837')
genes = rownames(temp)
for (i in (1:nrow(temp))){
  g = genes[i]
  a = temp[g, negative]
  b = temp[g, positive]
  qq = wilcox.test(as.numeric(a),as.numeric(b))
  results = c(results, qq$p.value)
}
names(results) = genes
tt = matrix(ncol=2,nrow=length(genes))
tt[,1] = genes
tt[,2] = results
tt = data.frame(tt)
write.table(tt, '/Users/woojunshim/Research/bovine/rank_comparison_.txt', sep='\t', quote=F)

### BOXPLOTS FOR BREADTHS OF H3K27ME3 
temp = read.table('/Users/woojunshim/Research/Data/mean_H3K27me3_breadths.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V2, y=log2(V1), fill=V2)) + geom_boxplot(show.legend = F) + labs(title='Mean H3K27me3 peak breadth', y='log2(H3K27me3 peak breadth)', x='') + theme_bw() 

###
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
genes = c('GATA4','GATA6','NKX2-5','TBX5','TBX20','MYH6','MYH7','MYL2','MYL3','TNNI3')
cols = c('E095','E104','E105','E082','E071','E070','E059','E058','E057','E047','E038','E037','E024','E016','E003','E006','E005','E004')
temp = temp[which(temp$V1 %in% genes),]
temp$label = rep('1',nrow(temp))
ggplot(temp, aes(x=label, y=V1, fill=V2)) + geom_tile(color='black', show.legend = F) + scale_fill_gradient2(low='white', high='blue') + scale_y_discrete(limits=rev(genes)) + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/hg19_rt_cardiac.pdf', width=1, height=3, dpi=600)

temp = read.table('/Users/woojunshim/Research/Data/discordance_table.txt', stringsAsFactors = F)

### FIGURE 4. MULTI-OMICS APPLICATIONS
order_ = read.table('/Users/woojunshim/Research/Data/proteomics/proteomics_cell_types.txt', stringsAsFactors = F)
order_ = order_$V1
temp = read.table('/Users/woojunshim/Research/Data/proteomics/HPM_protein_exp_selected.txt', stringsAsFactors = F)
temp = temp[,order_]
t1 = scale(temp, center=T, scale=T)
t2 = scale(t(t1), center=T, scale=T)
library(RColorBrewer)
library(gplots)
hmcol = colorRampPalette(c("green","white","red"))(256)
q = gsub("_"," ", rownames(t2))
rownames(t2) = q
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/HPM_protein_exp_selected_normalised.pdf', width=8, height=9)
heatmap.2(as.matrix(t2), col=hmcol, trace='none', Rowv=F, dendrogram='column', labCol=rep('',ncol(t2)))
dev.off()
tt = t(t2)
write.table(tt, '/Users/woojunshim/Research/Data/proteomics/HPM_protein_exp_selected_normalised_final.txt', sep='\t', quote=F)

### HPM DIS VS EXP PLOT
temp = read.table('/Users/woojunshim/Research/Data/Tabula_Muris/mouse_scrna_droplet_10x_fet.txt', stringsAsFactors = F, header=T)
temp$p.value = -log10(temp$p.value)
colnames(temp)[1] = 'significance'
q = gsub('_',' ',temp$colname)
temp$tissue = q
#ggplot(temp, aes(x=rank, y=significance, colour=method, shape=tissue)) + geom_point(size=3,position=position_jitter(width=1, height=.5)) + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold')) + labs(x='Best rank position', y='-log10(p-value)')
ggplot(temp, aes(x=rank, y=significance, colour=method)) + geom_point(size=3,position=position_jitter(width=1, height=.5), alpha=0.5) + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold')) + labs(x='Best rank position', y='-log10(p-value)')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/mouse_scrna_smartseq2_fet.pdf', width=6, height=4, dpi=600)

a=temp$significance[which(temp$method=='Discordance')]
b=temp$significance[which(temp$method=='Expression')]

a=temp$rank[which(temp$method=='Discordance')]
b=temp$rank[which(temp$method=='Expression')]

wilcox.test(a,b)

### ENTANGLEMENT 
library('dendextend')
library('magrittr')
dis = read.table('/Users/woojunshim/Research/Data/proteomics/HPM_discordance_protein.txt', stringsAsFactors = F)
exp = read.table('/Users/woojunshim/Research/Data/proteomics/HPM_protein_exp_table_final.txt', stringsAsFactors = F)
exp = log10(exp+1)
genes = rownames(exp)
dis = dis[genes, ]
dis = scale(dis, center=T, scale=T)
exp = scale(exp, center=T, scale=T)
d1 = as.dendrogram(hclust(dist(t(dis)), method='average'))
d2 = as.dendrogram(hclust(dist(t(exp)), method='average'))

### RANK PRODUCT FOR BOVEIN DATASET
library(RankProd)
temp = read.table('/Users/woojunshim/Research/bovine/bovine_discordance_score_updated__.txt', stringsAsFactors = F)
table_ = temp[,c('NE001407','NE001572','NE001634','NE001735','NE001887','NE001981','NE001455','NE001513','NE001636','NE001783','NE001831','NE003837')]
cl = rep(0,6)
cl = c(cl, rep(1,6))
a = RankProducts(table_, cl)
results = data.frame(genes=rownames(table_), positive_up_p=a$pval[,1], negative_up_p=a$pval[,2], positive_up_pfp=a$pfp[,1], negative_up_pfp=a$pfp[,2], positive_up_bh=p.adjust(a$pval[,1]), negative_up_bh=p.adjust(a$pval[,2]))
write.table(results, '/Users/woojunshim/Research/bovine/rank_product_table.txt', sep='\t', quote=F)

### SELECTED SAMPLES FOR CAGE TABLE
temp = read.table('/Users/woojunshim/Research/Data/FANTOM/development_sample_ids.txt', stringsAsFactors = F)
cols = temp$V1
qq = read.table('/Users/woojunshim/Research/Data/FANTOM/hg19_cage_expression_table.txt', stringsAsFactors = F)
ww = qq[,cols]
write.table(ww, '/Users/woojunshim/Research/Data/FANTOM/hg19_cage_expression_selected.txt', sep='\t', quote=F)

### CAGE DATA ANALYSIS 
samples = read.table('/Users/woojunshim/Research/Data/FANTOM/development_sample_ids.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/FANTOM/hg19_cage_expression_selected.txt', stringsAsFactors = F)
values = read.table('/Users/woojunshim/Research/Data/FANTOM/hg19_cage_dis_vs_exp_development_top5.txt', stringsAsFactors = F, header = T)
values1 = values[which(values$method=='Expression'),]
table_ = temp[,samples$V1]
table_ = scale(table_, center=T, scale=T)
pca = prcomp(t(table_))
cf = rainbow(max(samples$V2))
colours = cf[samples$V2]
qq = pca$x[,c('PC1', 'PC2')]
qq = as.data.frame(qq)
qq$cell_type = samples$V3
qq$value = -log10(values1$p.value)
qq$cell_type = gsub('_', ' ', qq$cell_type)
ggplot(qq, aes(x=PC1, y=PC2)) + geom_point(aes(colour=cell_type, size=value, alpha=0.5)) + scale_size(limits=c(0,85))+theme_bw()+theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/hg19_cage_expression_development_fet_top5.pdf', width=12, height=7, dpi=600)

### DENDROGRAM FOR INTERSPECIES ANALYSIS
library(lsa)
temp = read.table('/Users/woojunshim/Research/Data/interspecies/heart/interspecies_discordance_table.txt', stringsAsFactors = F)
temp = normalise1(temp)
cs = cosine(as.matrix(temp))
w = as.dist(1-cs)
ww = hclust(w)
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/Dendrogram_discordance.pdf', width=5, height=5)
plot(ww, horiz=TRUE)
dev.off()

### SCATTER PLOT FOR CORRELATION BETWEEN EXPRESSION AND DISCORDANCE ACROSS CELL-TYPES
temp = read.table('/Users/woojunshim/Research/Data/pearson_dis_vs_exp.txt', stringsAsFactors = F)
table_ = melt(t(temp))
ggplot(table_, aes(x=Var2, y=value)) + geom_point(size=3, color='red') + geom_hline(yintercept = 0, linetype='dashed', color='black') + scale_y_continuous(limits=c(0,1)) + theme_bw()+theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=10, face='bold',angle = 45, hjust = 1), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold')) + labs(y="Pearson correlation coefficient", x='')
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/pearson_dis_vs_exp_.pdf', width=8, height=5, dpi=600)

### GO ANALYSIS FOR VARIABLY EXPRESSED GENES 
temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
ref = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)
ref = ref$V1
cols = colnames(temp)
for (c in cols){
  temp = temp[order(temp[,c], decreasing=T),]
  all = rownames(temp)[which(temp[,c]>1)]
  genes = all[which(all %in% ref)][1:50]
  results = go_analysis(genes, all)
  write.table(results, paste('/Users/woojunshim/Research/Data/variably_expressed_tf/',c,'_go.txt', sep=''), quote=F, sep='\t')
}

### HEATMAP FOR TOP 50 VARIABLY EXPRESSED TF FOR SELECTED CELL-TYPES
temp = read.table('/Users/woojunshim/Research/Data/variably_expressed_tf/variably_expressed_tf_top50_selected_go_terms.txt', stringsAsFactors = F)
colnames(temp) = c('Brain Germinal Matrix','Pancreatic Islets','Left Ventricle','Primary T helper naive cells','H1 BMP4 Derived Mesendoderm')
t = gsub('_', ' ', rownames(temp))
rownames(temp) = t
table = melt(t(temp))
ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low = "white", high = "steelblue")  + labs(x='', y='', title='', fill='-log10(P)') + theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12, face='bold'), axis.text.y=element_text(size=12, face='bold'),legend.title=element_text(size=12, face='bold'), legend.text=element_text(size=10)) + scale_y_discrete(limits=c('forebrain development', 'generation of neurons','central nervous system neuron differentiation', 'pancreas development','endocrine pancreas development','pancreatic A cell differentiation','heart morphogenesis','embryonic heart tube morphogenesis','heart development','alpha-beta T cell differentiation','leukocyte differentiation','lymphocyte differentiation','endodermal cell fate specification','somatic stem cell population maintenance','gastrulation'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/variably_expressed_tf_top50_selected_go_terms.pdf', height=6, width=7, dpi=600)

### CALCULATE AVERAGE FOR EACH TISSUE (MOUSE SCRNA-SEQ)
library(Matrix)
ref = read.table('/Users/woojunshim/Research/Data/Tabula_Muris/SmartSeq2_samples.txt', stringsAsFactors = F)
cols = unique(ref$V2)
dis = readRDS('/Users/woojunshim/Research/Data/Tabula_Muris/TabulaMuris_SmartSeq2_discordance.RDS')
genes = rownames(dis)
results = data.frame(matrix(nrow=length(genes), ncol=length(cols)))
rownames(results) = genes
colnames(results) = cols
for (c in cols){
  cat (c)
  ids = ref[,1][which(ref$V2==c)]
  tt = dis[,ids]
  qq = apply(tt, 1, mean)
  results[,c] = qq
}
write.table(results, '/Users/woojunshim/Research/Data/Tabula_Muris/droplet_10x_discordance_for_tissues.txt')
  
### SPATIALLY RESOLVED RNA-SEQ 
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_fet_for_COSMIC_tf_rank_ratable.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_fet_for_kegg_top300_genes_tf_table.txt', stringsAsFactors = F)
colnames(temp) = c('sig','method','x','y')
temp$sig = -log10(temp$sig)
qq = temp[which(temp$method=='Discordance'),]
qq = temp[which(temp$method=='Expression'),]
ggplot(qq, aes(x=x, y=y, fill=sig)) + geom_tile(color='black') + scale_fill_gradient2(low='white', high='red',limits=c(0, max(temp$sig))) + labs(x='', y='', title='Enrichment of KEGG prostate cancer genes', fill='-log10(p-value)') + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'))
#ggplot(qq, aes(x=x, y=y, fill=sig)) + geom_tile(color='black') + scale_fill_continuous(low='red', high='white',limits=c(1, 100)) + labs(x='', y='', title='Enrichment of COSMIC cancer genes', fill='Rank position') + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/prostate_fet_for_kegg_top300_genes_tf_table_exp.pdf', width=10, height=8, dpi=600)

# CORRELATION BETWEEN SIGNIFICANCE AND CELL MARKER EXPRESSION
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_marker_genes.txt', stringsAsFactors = F)
ref = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_fet_for_COSMIC_at_top5_tf.txt', stringsAsFactors = F, header=T)
cols = ref$colname[which(ref$method=='Expression')]
gene = 'NPY'
a = temp[gene, cols]
b = ref$p.value[which(ref$method=='Expression')]
b = -log10(b)
cor(as.numeric(a),as.numeric(b))

### CORRELATION BETWEEN EXPRESSION AND H3K27ME3 BREADTHS
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
k27 = read.table('/Users/woojunshim/Research/Data/H3K27me3_width_table.txt', stringsAsFactors = F)
#k27 = read.table('/Users/woojunshim/Research/Data/broadPeaks/H3K27me3_widths_sum.txt', stringsAsFactors = F)
r1 = rownames(exp)
r2 = rownames(k27)
genes = intersect(r1, r2)
cols = colnames(exp)
exp = exp[genes, ]
k27 = k27[genes, cols]
k27 = scale(k27, center=T, scale=T)
exp = scale(exp, center=T, scale=T)
results = vector()
for (g in genes){
  a = exp[g,cols]
  b = k27[g,cols]
  results = c(results, cor(a,b))
}
names(results) = genes

var_tf = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)
house_ = read.table('/Users/woojunshim/Research/Data/house_keeping_genes.txt', stringsAsFactors = F)
all = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
var_tf = var_tf$V1
house_ = house_$V1
all_ = rownames(all)
stable_tf = read.table('/Users/woojunshim/Research/Data/stable_tf.txt', stringsAsFactors = F)
stable_tf = stable_tf$V1
stable_str = read.table('/Users/woojunshim/Research/Data/regulated_nontf.txt', stringsAsFactors = F)
stable_str = stable_str$V1

a = results[which(names(results) %in% var_tf)]
b = results[which(names(results) %in% house_)]
c = results[which(names(results) %in% all_)]
d = results[which(names(results) %in% stable_tf)]
e = results[which(names(results) %in% stable_str)]

table_ = data.frame(value=c(a,b,c,d,e), label=c(rep('Variably expressed TFs', length(a)), rep('House-keeping genes', length(b)), rep('Protein-coding genes', length(c)), rep('Non-variably expressed TFs', length(d)), rep('Variably expressed non-TFs', length(e))))
ggplot(table_, aes(x=label, y=value, fill=label)) + geom_boxplot(alpha=0.5) + scale_fill_manual(values=c('blue','green','yellow','purple', 'red')) + theme_bw() + labs(x='', y="Pearson's correlation coefficient", fill='') + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_blank(), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/correlation_btw_exp_and_h3k27me3_breadth_new.pdf', width=7, height=6, dpi=600)

# CORRELATION WITHIN A CELL-TYPE
results1 = vector()
for (c in cols){
  a = exp[genes,c]
  b = k27[genes,c]
  results1 = c(results1, cor(a,b, method='spearman'))
}
names(results1) = genes

### DEMARCATING PROSTATE CANCER CELLS USING THE MARKER GENE (SPINK1)
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_normalized_spatial_transcriptomics_.txt', stringsAsFactors = F)
temp = scale(temp, center=T, scale=T)
write.table(temp, '/Users/woojunshim/Research/Data/prostate_cancer/prostate_normalized_spatial_transcriptomics_z.txt', sep='\t', quote=F)


### FET OF HOUSE-KEEPING GENES (EXP VS DIS)
temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.discordance.variably_expressed_tf.distribution.txt')
temp = read.table('/Users/woojunshim/Research/Data/rt_vs_exp_new.txt', stringsAsFactors = F)
temp = temp[,2:ncol(temp)]
#temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(46)
high_ = results + 1.96 * sd_ / sqrt(46)
table1_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('Discordance', length(high_)))

temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.expression.variably_expressed_tf.distribution.txt')
temp = temp[,2:ncol(temp)]
temp[temp==0] = min(temp[temp!=0])
#temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(46)
high_ = results + 1.96 * sd_ / sqrt(46)
table2_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15), label=rep('Expression', length(high_)))

#col_ = rainbow(7)
#table_ = rbind(table1_, table2_, table3_, table4_, table5_, table6_, table7_)
table_ = rbind(table1_, table2_)
#ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line(show.legend = F, size=1) + labs(x='', y='', colour='') +  scale_colour_manual(values=col_) + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black')+ scale_x_discrete(limits=seq(1,10))+theme_bw() +theme(axis.text=element_text(size=16,face="bold")) 
#ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line() + labs(x='Rank position', y='Mean Jaccard similarity index', colour='') +  scale_colour_manual(values=col_)  + theme_bw() + theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=16,face='bold'), legend.text=element_text(size=16,face="bold"), legend.position='top')
ggplot(table_, aes(x=no, y=mean, colour=label)) +  geom_line(size=2) + labs(x='Rank position', y='Proportion', colour='')  + geom_errorbar(aes(ymin=low, ymax=high), width=.5, colour='black') + theme_bw()+theme(legend.text=element_text(size=16,face="bold"), legend.position='top', axis.title=element_text(size=16, face='bold'), axis.text=element_text(siz=16, face='bold'))
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/mean_jaccard_similarity_top1000_genes_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/mean_jaccard_similarity_all_genes_new_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
#ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/fet_house_keeping_genes_new_new.pdf',sep=''), width=8, height=5, units='in', dpi=600)
ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/proportion_variably_expressed_tf_roadmap.pdf',sep=''), width=8, height=5, units='in', dpi=600)

### DISTRIBUTION OF VARIABLY EXPRESSED TF BY RT SCORE
temp = read.table('/Volumes/backup/Research/Data/new_assigned/prop_variably_expressed_tf_by_rt.txt', stringsAsFactors = F)
temp$x = seq(1, 100)
ggplot(temp, aes(x=x, y=V1)) + geom_line(size=1) + geom_point() + labs(x='Rank position by RT', y='Proportion of variably expressed TFs') + theme_bw() + theme(legend.text=element_text(size=10,face="bold"), legend.position='top', axis.title=element_text(size=12, face='bold'), axis.text=element_text(siz=16, face='bold')) + geom_hline(yintercept = 0.01, linetype='dashed', colour='red', size=1)
ggsave('/Volumes/backup/Research/Documents/manuscripts/H3K27me3/new_figures/prop_variably_expressed_tf_by_rt.pdf', width=8, height=5, units='in', dpi=600)

### RT VS EXPRESSION (LOG2) WITH CONFIDENCE INTERVAL (95%)
temp = read.table('/Users/woojunshim/Research/Data/rt_vs_exp_ori.txt', stringsAsFactors = F)
aa = temp[order(temp$V2, decreasing = T),]
aa$x = seq(1, nrow(aa))
ggplot(aa, aes(x=x, y=V3)) + geom_line() + geom_ribbon(aes(ymin=aa$V4, ymax=aa$V5), fill='grey70') 
ggplot(temp, aes(x=V7, y=log2(V3), group=V7)) + geom_boxplot()

temp = read.table('/Volumes/backup/Research/Data/rt_vs_prop_expression_15000_new.txt', stringsAsFactors = F)
temp = temp[,2:ncol(temp)]
#temp = -log10(temp)
temp = t(temp)
#temp = temp[1:10,2:ncol(temp)]
rownames(temp) = seq(1, nrow(temp)) 
results = apply(temp, 1, mean)
sd_ = apply(temp, 1, sd)
low_ = results - 1.96 * sd_ / sqrt(46)
high_ = results + 1.96 * sd_ / sqrt(46)
table1_ = data.frame(no=as.numeric(names(results)), mean=round(results, 15), low=round(low_,15), high=round(high_,15))
ggplot(table1_, aes(x=no, y=mean)) + geom_line(size=1) + geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5) + labs(x='Rank bin (RT)', y='Log2(RPKM)') + annotate('text', x=30, y=4.3, label="Spearman's rho = 0.928", colour='red', size=5)+ theme_bw() +theme(legend.text=element_text(size=6,face="bold"), legend.position='top', axis.title=element_text(size=16, face='bold'), axis.text=element_text(siz=16, face='bold'))
ggsave('/Volumes/backup/Research/Documents/manuscripts/H3K27me3/new_figures/rt_vs_exp_new_15000.pdf', width=8, height=5, dpi=600)

### PLOT RT VS PROP. OF EXPRESSION 
temp = read.table('/Volumes/backup/Research/Data/rt_vs_prop_expression_15000.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V1, y=V2)) + geom_line(size=1) + geom_ribbon(aes(ymin=V4, ymax=V3), alpha=0.5) + labs(x='Rank bin (RT)', y='Expressed / all cell-types') + annotate('text', x=35, y=0.95, label="Spearman's rho = 0.974", colour='red', size=5)+ theme_bw() +theme(legend.text=element_text(size=6,face="bold"), legend.position='top', axis.title=element_text(size=16, face='bold'), axis.text=element_text(siz=16, face='bold'))
ggsave('/Volumes/backup/Research/Documents/manuscripts/H3K27me3/new_figures/rt_vs_prop_expression_15000.pdf', width=8, height=5, dpi=600)

### PLOT KEGG PATHWAY ENRICHMENT FOR PROSTATE CANCER DATA
temp = read.table('//Users/woojunshim/Research/Data/prostate_cancer/SPINK1_KEGG_full_table.txt', stringsAsFactors = F)
colnames(temp) = c('High SPINK1, Dis', 'High SPINK1, Exp','Low SPINK1, Dis','Low SPINK1, Exp')
order__ = read.table('/Users/woojunshim/Research/Data/prostate_cancer/kegg_order.txt', stringsAsFactors = F)
order_ = gsub('_', ' ', order__$V1)
temp = temp[,c(2,4)]
t = gsub('_', ' ', rownames(temp))
rownames(temp) = t
table = melt(t(temp))
ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low = "white", high = "steelblue")  + labs(x='', y='', title='', fill='-log10(P)') + theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) + scale_y_discrete(limits=order_)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/SPINK1_KEGG_exp_table.pdf', height=10, width=10, dpi=600)

# SPINK1 VS KEGG GENE ENRICHMENT CORRELATION
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/SPINK1.txt', stringsAsFactors = F, header=T)
input_ = read.table('/Users/woojunshim/Research/Data/KEGG/kegg_fet_table_dis.txt', stringsAsFactors = F)
input_ = read.table('/Users/woojunshim/Research/Data/KEGG/kegg_fet_table_exp.txt', stringsAsFactors = F)
#aa = input_['hsa05200-Pathways_in_cancer',temp$names]
table_ = melt(t(input_))
ggplot(table_, aes(x=Var1, y=Var2, fill=-log10(value))) + geom_tile() + scale_x_discrete(limits=temp$names) + scale_fill_gradient2(low='white',high='red') 
ggsave('/Users/woojunshim/Research/Data/KEGG/kegg_fet_table_exp.pdf', height=10, width=10, dpi=600)

### MAPPING GENE SYMBOLS TO ENTREZ IDS 
temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
genes = rownames(temp)
a = unlist(mget(genes, revmap(org.Hs.egSYMBOL),ifnotfound=NA))
a_ = data.frame(entrez=a)
write.table(a_, '/Users/woojunshim/Research/Data/KEGG/symbol_to_entrez.txt', sep='\t', quote=F)

d1 = data.frame(KEGGPATHID2NAME)
d2 =  data.frame(KEGGPATHID2EXTID)
write.table(d1, '/Users/woojunshim/Research/Data/KEGG/kegg_description.txt', sep='\t', quote=F)
write.table(d2, '/Users/woojunshim/Research/Data/KEGG/kegg_to_entrez.txt', sep='\t', quote=F)

### FET FOR KEGG TERMS (OR BP TERMS)
temp = read.table('/Users/woojunshim/Research/Data/KEGG/prop_bp_pluripotency_.txt', stringsAsFactors = F)
rr = gsub('_', ' ', temp$V1)
rownames(temp) = rr
temp = temp[, 2:ncol(temp)]
colnames(temp) = seq(10,100,10)
table_ = melt(t(temp))
no_terms = nrow(temp)
ggplot(table_, aes(x=Var1, y=value, colour=Var2)) + geom_line(size=1) + geom_point() + labs(x='Rank bin (by RT)', y='Proportion', color='') + scale_x_continuous(breaks=seq(10,100,10)) + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) + geom_hline(yintercept=0.1, linetype='dashed', colour='dark gray', size=1)
ggplot(table_, aes(x=Var1, y=value, colour=Var2)) + geom_line(size=1) + geom_point() + labs(x='Rank bin (by RT)', y='Proportion', color='') + scale_x_continuous(breaks=seq(10,100,10)) + scale_colour_discrete(label=rep('',no_terms)) +theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) + geom_hline(yintercept=0.1, linetype='dashed', colour='dark gray', size=1)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/prop_bp_pluripotency_.pdf', width=10, height=8, dpi=600)

### GET GO TERMS FOR DEVELOPMENT 
a=c('GO:0019827','GO:0007492','GO:0007498','GO:0007398','GO:0072358','GO:0007417','GO:0055123','GO:0035270','GO:0035272','GO:0002520','GO:0007399','GO:0060541','GO:0061458','GO:0001501')
a=c('GO:0072089')
for (i in a){
  result = get_go_terms(i)
  write.table(result, paste('/Users/woojunshim/Research/Data/genes_with_go_term/',i,'.txt',sep=''), sep='\t', quote=F)
}

### DRAW SPINK1 HIGH AND LOW MAP
temp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/SPINK1_high_low_map.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V1, y=V2, fill=V3)) + geom_tile(color='black', show.legend = F) + scale_fill_continuous(low='white', high='red') + labs(title='SPINK1',x='',y='') + theme_bw()
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/SPINK1_map.pdf',width=7, height=6, dpi=600)

### RANK PRODUCT TEST 
library(RankProd)
temp = read.table('/Users/woojunshim/Research/bovine/new/discordance_by_TPM.txt', stringsAsFactors = F)
low = c('NE001887','NE001634','NE001735','NE001572','NE001981','NE001407')
high = c('NE003837','NE001513','NE001831','NE001783','NE001636','NE001455')
table_ = temp[,c(low,high)]
cl = c(rep(0,length(low)),rep(1,length(high)))
a = RankProducts(table_, cl)
results = data.frame(genes=rownames(table_), high_up_p=a$pval[,1], low_up_p=a$pval[,2], high_up_pfp=a$pfp[,1], low_up_pfp=a$pfp[,2], high_up_bh=p.adjust(a$pval[,1]), low_up_bh=p.adjust(a$pval[,2]))
write.table(results, '/Users/woojunshim/Research/bovine/new/results/rank_product_result_for_TPM.txt', sep='\t', quote=F)

### PLOT RT DOT PLOT
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
a = sort(as.numeric(temp$V2), decreasing=T)
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/RT_distribution.pdf', width=5, height=5)
plot(x=seq(1,length(a)), y=a, ylab='RT', xlab='')
dev.off()

### BOVINE 
#temp = read.table('/Users/woojunshim/Research/bovine/new/NEV_homologues.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/bovine/new/expression_data_.txt', stringsAsFactors = F)
genes = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
genes = genes$V1
#new_ = temp[which(temp$human_gene_name %in% genes),]
new_ = temp[which(rownames(temp) %in% genes),]
rownames(new_) = new_$human_gene_name
write.table(new_, '/Users/woojunshim/Research/bovine/new/expression_data_by_read_count.txt', sep='\t', quote=F)

### CALCULATE EUCLIDEAN DISTANCE FOR CARDIAC SC DATA
pathway = '/Users/woojunshim/Research/cardiac_sc/ds/'
prefix = 'ds'
menu = c('d0','d2','d5','d15','d30')
for (i in 1:length(menu)){
  file = paste(pathway,prefix,'_',menu[i],'.txt',sep='')
  cat(file)
  temp = read.table(file, stringsAsFactors = F)
  a= dist(t(temp), method='euclidean',diag=F)
  b= as.matrix(a)
  write.table(b, paste(pathway,prefix,'_',menu[i],'_dist.txt',sep=''), sep='\t', quote=F)
}

### HEATMAP 
day = 'd30'
ds = read.table(paste('/Users/woojunshim/Research/cardiac_sc/ds/ds_',day,'_dist.txt',sep=''),stringsAsFactors = F)
exp = read.table(paste('/Users/woojunshim/Research/cardiac_sc/exp/exp_',day,'_dist.txt',sep=''),stringsAsFactors = F)
plot_heatmap(ds, 'DS (day30)', 'Euclidean distance',range(ds))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/ds_d30_dist.pdf', width=6, height=5, dpi=600)
plot_heatmap(exp, 'Expression (day30)', 'Euclidean distance',range(ds))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/exp_d30_dist.pdf', width=6, height=5, dpi=600)

### CORRELATION BETWEEN SPP (or HOMER) AND MACS
#a = read.table('/Users/woojunshim/Research/Data/spp_peaks/H3K27me3/processed/spp_rt.txt', stringsAsFactors = F)
a = read.table('/Users/woojunshim/Research/Data/homer_peaks/H3K27me3/processed/homer_rt.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/macs2/macs2_rt.txt', stringsAsFactors = F)
#b = read.table('/Users/woojunshim/Research/Data/spp_peaks/H3K27me3/processed/spp_peaks_rt.txt', stringsAsFactors = F)
#a = read.table('/Users/woojunshim/Research/Data/macs2/macs2_peaks_rt.txt', stringsAsFactors = F)
genes = intersect(a$V1, b$V1)
rownames(a) = a$V1
rownames(b) = b$V1
v1 = a[genes,2]
v2 = b[genes,2]
value=round(cor(v1,v2,method='pearson'),3)
table_ = data.frame(macs2=v1, homer=v2)
ggplot(table_, aes(x=macs2, y=homer)) + geom_point(aes(colour='red', alpha=0.5), show.legend = F) + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) + annotate('text', x=0.75, y=0.90, label=paste("Pearson's"," rho = ",value,sep=''), colour='black', size=5)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/homer_vs_macs2.pdf',width=5, height=5, dpi=600)

### GO ANALYSIS FOR TOP 1359 VS ...
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
temp = temp[order(temp$V2, decreasing = T),]
input_ = temp$V1[(nrow(temp)-1359):nrow(temp)]
all_ = temp$V1
bp = go_analysis(input_, all_)
cc = go_analysis(input_, all_, onto='CC')
write.table(bp, '/Users/woojunshim/Research/Data/1359/bottom_1359_bp.txt',sep='\t',quote=F)
write.table(cc, '/Users/woojunshim/Research/Data/1359/bottom_1359_cc.txt',sep='\t',quote=F)

### ENRICHMENT PLOT FOR 1359 GENES
temp = read.table('/Users/woojunshim/Research/Data/1359/go_enrichment_table_selected.txt', stringsAsFactors = F)
terms = read.table('/Users/woojunshim/Research/Data/1359/selected_terms.txt', stringsAsFactors = F)
terms = terms$V1
temp = temp[terms, c('Top_1359','Next_1359','Bottom_1359')]
pdf('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/go_enrichment_table_selected.pdf', width=8, height=5)
sig_plot(temp)
dev.off()

ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_y_discrete(limits=terms) + scale_fill_gradient2(low='white',high='red') + labs(x='', y='') + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/go_enrichment_table_selected.pdf', width=5, height=5) 

### HIERARCHICAL CLUSTERING OF VARIABLY EXPRESSED TFS
qq = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
cols = colnames(qq)
regulated = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)
regulated = regulated$V1
temp = read.table('/Users/woojunshim/Research/Data/H3K27me3_width_table.txt', stringsAsFactors = F)
temp = temp[,cols]
temp = scale(temp, center=T, scale=T)
table_ = temp[which(rownames(temp) %in% regulated),]
table2_ = scale(t(table_), center=T, scale=T)
table2_ = t(table2_)
write.table(table2_, '/Users/woojunshim/Research/Data/H3K27me3_relative_absence.txt', sep='\t', quote=F)

### LOGIT REGRESSION
temp = read.table('/Volumes/backup/Research/bigdata/roadmap/broad_peaks/assigned/hm_data_integrated.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/logit_regression/inputs.txt', stringsAsFactors = F)
train = temp[which(temp$V8!='E095'),c('V2','V3','V4','V5','V6','V7','V9')]
test = temp[which(temp$V8=='E095'),c('V2','V3','V4','V5','V6','V7','V9')]
rownames(test) = temp$V1[as.numeric(rownames(test))]
temp = temp[,2:ncol(temp)]
model = glm(V9~., family=binomial(link='logit'), data=train)
summary(model)
anova(model, test='Chisq')
results = predict(model, newdata=subset(tt, select=c('V2','V3','V4','V5','V6','V7')), type='response')

### CONVERT WIDTHS TO Z-SCORE 
file_ = 'H3K36me3'
temp = read.table(paste('/Users/woojunshim/Research/Data/new_assigned/', file_, '_tss_2.5kb_combined.txt', sep=''), stringsAsFactors = F)
temp = scale(temp, center=T, scale=T)
write.table(temp, paste('/Users/woojunshim/Research/Data/new_assigned/', file_, '_tss_2.5kb_combined_z.txt', sep=''), sep='\t', quote=F)

### FET FOR DE FOR SC CARDIAC DATASETS
temp = read.table('/Users/woojunshim/Research/cardiac_sc/comparison/d2c1_endoderm_fet.txt', stringsAsFactors = F)
title_ = 'Day2 Cluster 1, Endoderm development (GO:0007492)'
temp$V1 = gsub('_', ' ', temp$V1)
rownames(temp) = temp$V1
temp = temp[,2:ncol(temp)]
colnames(temp) = seq(1,100)
table_ = melt(t(temp))
ggplot(table_, aes(x=Var1, y=-log10(value), color=Var2)) + geom_line(size=1) + geom_point(size=1) + labs(title=title_, x='Rank position', y='-log10(p-value)', color='') + scale_colour_manual(values=c("#F8766D","#619CFF","#7CAE00","#C77CFF"), limits=c('Discordance','Expression','DE C2','DE C3'))  + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/d2c1_endoderm_fet.pdf', width=8, height=6, dpi=600)

### COMPARING JACCARD SIMILARITY 
k27me3 = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K27me3.txt', stringsAsFactors = F)
k9me3 = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K9me3.txt', stringsAsFactors = F)
k4me3 = read.table('/Users/woojunshim/Research/Data/new_assigned/jaccard_15000_all_genes_by_H3K4me3.txt', stringsAsFactors = F)

### QUICK CHECK FOR PROSTATE CANCER DATASETS
### NKX2-2, IRF2
map = read.table('/Users/woojunshim/Research/Data/prostate_cancer/SPINK1.txt', stringsAsFactors = F)
dis = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_discordance.txt', stringsAsFactors = F)
exp = read.table('/Users/woojunshim/Research/Data/prostate_cancer/prostate_normalized_spatial_transcriptomics_.txt', stringsAsFactors = F)
coor = map$names[which(map$value>4)]
rest = map$names[which(map$value<=4)]
a1=as.numeric(dis['TBX3', coor])
a2=as.numeric(dis['TBX3', rest])
b1=as.numeric(exp['TBX3', coor])
b2=as.numeric(exp['TBX3', rest])
wilcox.test(b1,b2)

### PLOT DISTRIBUTION PLOTS 
library(ggplot2)
library(reshape2)
a = read.table('/Users/woojunshim/Research/cardiac_sc/fet/ds_d30_heart_development_tf_fet.txt', stringsAsFactors = F)
a = as.numeric(a[1,2:ncol(a)])
b = read.table('/Users/woojunshim/Research/cardiac_sc/fet/ds_d30_sarcomere_fet.txt', stringsAsFactors = F)
b = as.numeric(b[1,2:ncol(b)])
c = read.table('/Users/woojunshim/Research/cardiac_sc/fet/ds_d30_heart_signaling_fet.txt', stringsAsFactors = F)
c = as.numeric(c[1,2:ncol(c)])
d = read.table('/Users/woojunshim/Research/cardiac_sc/fet/ds_d30_house_keeping_fet.txt', stringsAsFactors = F)
d = as.numeric(d[1,2:ncol(d)])
table_ = data.frame(heart_development_tf=a, sarcomere=b, heart_signaling=c, house_keeping=d)
table = melt(t(table_))
ggplot(table, aes(x=Var2, y=-log10(value), colour=Var1)) + geom_line(size=1) + labs(x='Rank position', y='Proportion', colour='') + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/d30_ds_.pdf', width=8, height=5, dpi=600)

library(RColorBrewer)
cols <- rev(brewer.pal(11, 'RdYlBu'))
a = read.table('/Users/woojunshim/Research/cardiac_sc/fet/exp_d30_heart_development_tf_fet.txt', stringsAsFactors = F)
a = as.numeric(a[1,2:ncol(a)])
temp = order(a)
a[temp] = seq(100,1,-1)
b = read.table('/Users/woojunshim/Research/cardiac_sc/fet/exp_d30_sarcomere_fet.txt', stringsAsFactors = F)
b = as.numeric(b[1,2:ncol(b)])
temp = order(b)
b[temp] = seq(100,1,-1)
c = read.table('/Users/woojunshim/Research/cardiac_sc/fet/exp_d30_heart_signaling_fet.txt', stringsAsFactors = F)
c = as.numeric(c[1,2:ncol(c)])
temp = order(c)
c[temp] = seq(100,1,-1)
d = read.table('/Users/woojunshim/Research/cardiac_sc/fet/exp_d30_house_keeping_fet.txt', stringsAsFactors = F)
d = as.numeric(d[1,2:ncol(d)])
temp = order(d)
d[temp] = seq(100,1,-1)
table_ = data.frame(heart_development_tf=a, sarcomere=b, heart_signaling=c, house_keeping=d)
colnames(table_) = c('Heart development TF','Sarcomere','Heart signaling','House keeping')
table = melt(t(table_))
ggplot(table, aes(x=Var2, y=Var1, fill=value)) + geom_tile(size=0.25) + scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=50) + labs(x='Rank position (gene)', y='Rank position (significance)', fill='Rank percentile')  + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
#ggplot(table, aes(x=Var2, y=value, colour=Var1)) + geom_line(size=1) + labs(x='Rank position (gene)', y='Rank position (significance)', colour='') + theme_bw() + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/d30_exp_rank.pdf', width=7, height=3, dpi=600)

### HEATMAP FOR H3K27ME3 BREADTH 
library(ggplot2)
library(reshape2)
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_broad_table.txt', stringsAsFactors = F)
a=sort(rownames(temp))
temp=temp[a,]
#temp = scale(temp, scale=T, center=T)
temp = log2(temp+1)
table_ = melt(t(temp))
#hmcol = colorRampPalette(c("white","red"))(256)
#heatmap.2(as.matrix(temp), col=hmcol, dendrogram='none',trace='none', labCol='none', labRow='none')

ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(show.legend = F) + scale_fill_gradient2(low='white', high='black') + labs(x='', y='', fill='H3K27me3') + scale_x_discrete(labels=rep('', ncol(temp))) + scale_y_discrete(labels=rep('', nrow(temp)), limits=rownames(temp)) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/H3K27me3_breadth_heatmap_broad.pdf', width=5, height=5)

### SCALE BAR FOR RTS
a = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
aa = sort(a$V1)
rownames(a) = a$V1
values = a[rownames(temp),]$V2
table_ = data.frame(value=values, class=rep(1,length(values)), name=rownames(temp))
ggplot(table_, aes(x=class, y=name, fill=value)) + geom_tile() + scale_fill_gradient2(low='white', high='black') + labs(x='', y='', fill='RTS') + scale_x_discrete(labels=rep('',)) + scale_y_discrete(labels=rep('', nrow(temp)),limits=rownames(temp)) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/RTS_heatmap.pdf', width=3, height=5)

### CREATE RTS HEAT MAPS FOR SELECTED ROADMAP SAMPLES 
samples = c('E003','E013','E037','E055','E112','E054','E070','E100','E095','E109','E066','E097')
ref = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
ds_ = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_ds.txt', stringsAsFactors = F)
for (s in samples){
  a=temp[,s]
  names(a) = rownames(temp)
  c=sort(a[a>0],decreasing = T)
  values=ref$V2[match(names(c),ref$V1)]
  names(values) = names(c)
  values = values[!is.na(values)]
  #order_=rev(names(values))
  table_ = data.frame(value=values, class=rep(1,length(values)), name=names(values))
  order_ = sort(names(values))
  #ggplot(table_, aes(x=class, y=name, fill=value)) + geom_tile() + scale_fill_gradient2(low='white', high='black',limits=c(0,1)) + labs(x='', y='', fill='RTS') + scale_x_discrete(labels=rep('',)) + scale_y_discrete(labels=rep('', length(values)),limits=order_) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
  #ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/RTS/RTS_heatmap_',s,'_RTS_black.pdf',sep=''), width=3, height=5*length(c)/26833)

  #ggplot(table_, aes(x=class, y=name, fill=value)) + geom_tile() + scale_fill_gradient(low='black', high='yellow',limits=c(0,1)) + labs(x='', y='', fill='RTS') + scale_x_discrete(labels=rep('',)) + scale_y_discrete(labels=rep('', length(values)),limits=order_) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
  #ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/abstract_plot/black_colour/RTS/RTS_heatmap_',s,'_RTS_yellow.pdf',sep=''), width=3, height=5*length(c)/26833)
  
  exp=temp[order_,s]
  names(exp)=order_
  #aa=exp
  #aa[order(exp, decreasing=T)]=seq(length(exp),1,-1)
  table_ = data.frame(value=exp, class=rep(1,length(exp)), name=names(exp))
  ggplot(table_, aes(x=class, y=name, fill=log(value+1))) + geom_tile() + scale_fill_gradient(low='black', high='yellow') + labs(x='', y='', fill='ln(exp)') + scale_x_discrete(labels=rep('',)) + scale_y_discrete(labels=rep('', length(exp)),limits=order_) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
  ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/abstract_plot/black_colour/EXP/heatmap_',s,'_exp_ln.pdf',sep=''), width=3, height=5*length(c)/26833)

  ds=ds_[order_,s]
  names(ds)=order_
  table_ = data.frame(value=ds, class=rep(1,length(ds)), name=names(ds))
  ggplot(table_, aes(x=class, y=name, fill=value)) + geom_tile() + scale_fill_gradient(low='black', high='yellow') + labs(x='', y='', fill='DS') + scale_x_discrete(labels=rep('',)) + scale_y_discrete(labels=rep('', length(ds)),limits=order_) + theme_bw() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) 
  ggsave(paste('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/abstract_plot/black_colour/DS/heatmap_',s,'_ds_.pdf',sep=''), width=3, height=5*length(c)/26833)

}

### CHECK RELATIONSHIP BETWEEN EXPRESSION H3K27ME3 BREADTH
a=read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
b=read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_width_table.txt', stringsAsFactors = F)
ta=scale(a, center=T, scale=T)
tb=scale(b, center=T, scale=T)
cols = intersect(colnames(ta), colnames(tb))
x=tb['NKX2-5',cols]
y=ta['NKX2-5',cols]

### EM ALGORITHM
library(EMCluster, quietly = TRUE)
set.seed(1234)
x2 <- da2$da
ret <- init.EM(x2, nclass = 2)
ret.new <- assign.class(x2, ret, return.all = FALSE)
cols = ifelse(ret.new$class==1, 'red', 'blue')
x2$colour=cols
plot(x=x2$x, y=x2$y, col=x2$colour, pch=20)

x1 = read.table('/Users/woojunshim/Research/BN/data/h3k27me3_signal_nontf.tsv', stringsAsFactors = F)
require("ggrepel")

p1=rnorm(46,100,15)
p2=rnorm(46,200,15)

x1=data.frame(V1=p1, V2=p2)

gene = 'GTF2B'
x1 = read.table(paste('/Users/woojunshim/Research/BN/data/',gene,'.tsv',sep=''), stringsAsFactors = F)
#x1 <- da1$da

x1=data.frame(exp_=as.numeric(exp[gene, cols]), k27_=as.numeric(k27[gene, cols]))
emobj <- simple.init(x1[,c(1,2)], nclass = 3)
emobj <- shortemcluster(x1[,c(1,2)], emobj)
ret <- emcluster(x1[,c(1,2)], emobj, assign.class = TRUE)
ret$llhdval

x1$colour = ifelse(ret$class==1, 'LC 1','LC 2') 
ggplot(x1, aes(x=exp_, y=k27_, colour=colour, label=colour)) + geom_point(size=1) + stat_ellipse(aes(x=exp_, y=k27_, colour=colour),type = "norm") + geom_text_repel(size=2.5) + labs(x='log2(RPKM+1)', y='H3K27me3(bp)', title=gene, colour='') + theme_bw()
ggsave(paste('/Users/woojunshim/Research/BN/figures/',gene,'.pdf',sep=''), width=5, height=5, dpi=600)

### SLOPE OF THE LINEAR REGRESSION LINE
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
exp = log2(exp+1)
k27 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_width_table.txt', stringsAsFactors = F)
genes = intersect(rownames(exp), rownames(k27))
cols = colnames(exp)
results = vector()

for (g in genes){
  table_ = data.frame(exp_=as.numeric(exp[g,cols]), k27_=as.numeric(k27[g,cols]))
  rr = lm(k27_ ~ exp_, data=table_)
  if (nrow(summary(rr)$coefficients)==2){
    results = rbind(results, c(g, summary(rr)$coefficients[2,1], summary(rr)$coefficients[2,4]))
  }
}
results = data.frame(results)
colnames(results) = c('#Gene','Slope','p-value')
write.table(results, '/Users/woojunshim/Research/BN/lm_results.txt', sep='\t', quote=F)

# PEARSON'S CORRELATION
results1 = vector()
for (g in genes){
  aa = cor.test(as.numeric(exp[g,cols]), as.numeric(k27[g,cols]), method='pearson')
  results = rbind(results, c(g, aa$estimate, aa$p.value))
}
results1 = data.frame(results)
colnames(results1) = c('#Gene','Pearson','p-value')
write.table(results1, '/Users/woojunshim/Research/BN/pearsons_results.txt', sep='\t', quote=F)

# COMBINE PEARSON'S CORRELATION, INTERCEPT AND SLOPE FROM LINEAR REGRESSION 
results=vector()
for (g in genes){
  aa = cor(as.numeric(exp[g,cols]), as.numeric(k27[g,cols]), method='pearson')
  table_ = data.frame(exp_=as.numeric(exp[g,cols]), k27_=as.numeric(k27[g,cols]))
  rr = lm(k27_ ~ exp_, data=table_)  
  if (nrow(summary(rr)$coefficients)==2){
    results = rbind(results, c(g, aa, summary(rr)$coefficients[1,1], summary(rr)$coefficients[2,1]))
  }
}
results = data.frame(results)
colnames(results)=c('gene','pearson','intercept','slope')
results=results[complete.cases(results), ]

results=data.frame()
for (g in genes){
  aa = cor(as.numeric(exp[g,cols]), as.numeric(k27[g,cols]), method='spearman')
  mean_ = mean(as.numeric(k27[g,cols]))
  sd_ = sd(as.numeric(exp[g,cols]))
  results = rbind(results, c(aa, mean_, sd_))
}
rownames(results) = genes
colnames(results)=c('spearman','K27_mean','Exp_sd')
results=results[complete.cases(results), ]
write.table(results, '/Users/woojunshim/Research/BN/features_new.txt',sep='\t', quote=F)

d=data.frame(pearson=as.numeric(results$pearson), intercept=as.numeric(results$intercept), slope=as.numeric(results$slope))
rownames(d)=results$gene
#rownames(results)=results$gene
#results=results[,c(2:ncol(results))]

d=read.table('/Users/woojunshim/Research/BN/features.txt', stringsAsFactors = F)
a=scale(d, center=T, scale=T)
#write.table(a, '/Users/woojunshim/Research/BN/features.txt', sep='\t', quote=F)

b=a[a[,1]<10 & a[,2]<10 & a[,3]<10, ]
b=b[b[,1]>-10 & b[,2]>-10 & b[,3]>-10, ]

# CLUSTERING
library(gplots)
hmcol = colorRampPalette(c("dark blue","white","red"))(256)

qq=hclust(dist(b), method='ward.D')
order1=qq$labels[qq$order]
ww = cutree(qq, 4)
rainbow_ = c('blue','green','yellow','red')
#rainbow_ = rainbow(4)
cols = rainbow_[as.numeric(ww)]

labels = rainbow_[ww[order1]]
names(labels) = order1
ll = data.frame(label=labels)
write.table(ll, '/Users/woojunshim/Research/BN/clustering/labels.txt', sep='\t', quote=F)


pdf('/Users/woojunshim/Research/BN/figures/features_hc.pdf', width=5, height=5)
heatmap.2(as.matrix(b), main='', trace = 'none', col=hmcol, na.color='black',labRow = F, cexCol = 0.8, labCol=colnames(b), distfun=dist, hclustfun=function(x) hclust(x, method = "ward.D"), RowSideColors = cols)
dev.off()

# PERFORM GO ANALYSIS FOR EACH CLUSTER
results = go_analysis(names(labels)[which(labels=='green')], names(labels))
write.table(results, '/Users/woojunshim/Research/BN/clustering/green_go.txt', sep='\t', quote=F)

# PLOT LEAST SQUARE LINE FOR EACH CLUSTER
# FIRST CALCULATE COORDINATES FOR INTERCEPT (I.E. WHEN LOG2(EXP+1)==0)
# AND WHEN K27=0
# GET AN AVERAGE LINE FOR EACH CLUSTER
a=read.table('/Users/woojunshim/Research/BN/features.txt', stringsAsFactors = F)
results = data.frame(x1=rep(0,nrow(a)), y1=rep(0,nrow(a)), x2=rep(0,nrow(a)), y2=rep(0,nrow(a)))
rownames(results) = rownames(a)
results[,2] = a[,2]
results[,3] = -a[,2]/a[,3]
write.table(results, '/Users/woojunshim/Research/BN/coordinates.txt', sep='\t', quote=F)

### SPEARMAN'S RHO
results = vector()
for (g in genes){
  a=cor(as.numeric(log2(exp[g,cols]+1)), as.numeric(k27[g,cols]), method='spearman')
  results = rbind(results, c(g,a))
}
write.table(results, '/Users/woojunshim/Research/BN/data/spearman.txt', sep='\t', quote=F)

### RANDOM DATA 
a=rnorm(24,100,15)
b=rnorm(24,500,15)
c=rnorm(24,900,15)
d1=data.frame(value=c(a,b,c))

emobj <- simple.init(d1, nclass = 1)
emobj <- shortemcluster(d1, emobj)
ret1 <- emcluster(d1, emobj, assign.class = TRUE)
ret1$llhdval

emobj <- simple.init(d1, nclass = 2)
emobj <- shortemcluster(d1, emobj)
ret2 <- emcluster(d1, emobj, assign.class = TRUE)
ret2$llhdval

emobj <- simple.init(d1, nclass = 3)
emobj <- shortemcluster(d1, emobj)
ret3 <- emcluster(d1, emobj, assign.class = TRUE)
ret3$llhdval

d=rnorm(46,100,15)
e=rnorm(46,200,15)
d1=data.frame(value=c(d,e))
emobj <- simple.init(d1, nclass = 4)
emobj <- shortemcluster(d1, emobj)
ret <- emcluster(d1, emobj, assign.class = TRUE)
ret$llhdval

### CALCULATE SILOUETTE FOR ALL GENES
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
k27 = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_width_table.txt', stringsAsFactors = F)
exp = log2(exp+1)
k27 = log2(k27+1)
cols = colnames(exp)
library(cluster)

non_zero = function(list_){
  return (length(which(list_>1)))
}

g1_=apply(exp, 1, non_zero)
g1 = names(g1_)[which(g1_>1)]
g2_=apply(k27, 1, non_zero)
g2 = names(g2_)[which(g2_>1)]
genes = intersect(g1, g2)

results = vector()
for (gene in genes){
  x1=data.frame(exp_=as.numeric(exp[gene, cols]), k27_=as.numeric(k27[gene, cols]))
  
  emobj <- simple.init(x1, nclass = 2)
  emobj <- shortemcluster(x1, emobj)
  ret <- emcluster(x1, emobj, assign.class = TRUE)
  if (length(unique(ret$class))==2){
    a=dist(x1,"euclidean")
    clu=ret$class
    si2 = silhouette(clu, a)
    result = sum(si2[,3])/nrow(x1)
    results = rbind(results, c(gene, result))
  }
}

### APPROACH EM ALGORITHM WITH K=1
### WE ASSUME THAT CELL-TYPE SPECIFIC GENES HAVE A POOR FIT

### OK
### FILTER OUT GENE WITH LIKELIHOOD > 0
### DO MULTIPLE TIMES FOR EACH GENE AND AVERAGE
gene='GATA4'
x1=data.frame(exp_=as.numeric(exp[gene, cols]), k27_=as.numeric(k27[gene, cols]))
emobj <- simple.init(x1, nclass = )
emobj <- shortemcluster(x1, emobj)
ret <- emcluster(x1, emobj, assign.class = TRUE)
#length(unique(ret$class))
ret$llhdval
###

### EXRACT AVERAGE LOG LIKELIHOOD FOR ALL GENES 
### 10 ITERATIONS
results = vector()
for (gene in genes){
  x1=data.frame(exp_=as.numeric(exp[gene, cols]), k27_=as.numeric(k27[gene, cols]))
  emobj <- simple.init(x1, nclass = 2)
  emobj <- shortemcluster(x1, emobj)
  ret <- emcluster(x1, emobj, assign.class = TRUE)
  if ((length(unique(ret$class)==2)) && (!is.na(ret$llhdval))) {
    if (ret$llhdval<0){
      results = rbind(results, c(gene, ret$llhdval))
    }
  }
}

### ASSUME K=1 
results = vector()
for (gene in genes){
  x1=data.frame(exp_=as.numeric(exp[gene, cols]), k27_=as.numeric(k27[gene, cols]))
  emobj <- simple.init(x1[,c(1,2)], nclass = 1)
  emobj <- shortemcluster(x1[,c(1,2)], emobj)
  ret <- emcluster(x1[,c(1,2)], emobj, assign.class = TRUE)
  results = rbind(results, c(gene, ret$llhdval))
}

write.table(results, '/Users/woojunshim/Research/BN/data/ll_one_state.txt', sep='\t', quote=F)


### GO ANALYSIS FOR TOP 1% DS GENES ACROSS ROADMAP SAMPLES
temp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols_discordance.txt', stringsAsFactors = F)
ref_ = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
ref_ = ref_[order(ref_$V2, decreasing=T),]
cols = colnames(temp)
for (c in cols){
  temp = temp[order(temp[,c], decreasing=T),]
  all_ = temp[,c]
  names(all_) = rownames(temp)
  a = names(all_[which(all_>0)])
  input_ = a[1:(length(a)*0.01)]
  result = go_analysis(input_, a)
  write.table(result, paste('/Users/woojunshim/Research/Data/roadmap_ds/',c,'_top1.txt', sep=''), sep='\t', quote=F)
}

temp = read.table('/Users/woojunshim/Research/Data/roadmap_ds/top1_ds_table.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/Data/roadmap_ds/top1_ds_table_all.txt', stringsAsFactors = F)
colnames(temp) = gsub('_', ' ', colnames(temp))
rownames(temp) = gsub('_', '', rownames(temp))
table = melt(t(temp))
qq = read.table('/Users/woojunshim/Research/Data/roadmap_ds/selected_terms_all.txt', stringsAsFactors = F)
q = gsub('_', ' ', qq$V1)
ww = read.table('/Users/woojunshim/Research/Data/roadmap_ds/desc_cols.txt', stringsAsFactors = F)
w = gsub('_', ' ', ww$V1)

ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low = "white", high = "red")  + labs(x='', y='', title='', fill='-log10(P)') + theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12, face='bold'), axis.text.y=element_text(size=12, face='bold'),legend.title=element_text(size=12, face='bold'), legend.text=element_text(size=10)) + scale_y_discrete(limits=q) + scale_x_discrete(limits=w)
ggsave('/Users/woojunshim/Research/Documents/manuscripts/H3K27me3/new_figures/test.pdf', height=20, width=15, dpi=600)

### CHECK RANDOM DISTRIBUTION LOG LIKELIHOOD
results1=vector()
no=100
for (i in 1:no){
  #x1=data.frame(exp_=rnorm(100), k27_=rnorm(100)) # Random
  exp1=rnorm(50)
  k271=rnorm(50)
  exp2=rnorm(50,1)
  k272=rnorm(50,1)
  x1=data.frame(exp_=c(exp1,exp2), k27_=c(k271,k272))
  emobj <- simple.init(x1[,c(1,2)], nclass = 1)
  emobj <- shortemcluster(x1[,c(1,2)], emobj)
  ret <- emcluster(x1[,c(1,2)], emobj, assign.class = TRUE)
  results1=c(results1, ret$llhdval)
}

### TEST LL FOR HOUSE-KEEPING AND VARIABLY EXPRESSED TFS
genes = intersect(rownames(exp), rownames(k27))
house = read.table('/Users/woojunshim/Research/Data/house_keeping_genes.txt', stringsAsFactors = F)
house = house$V1
house = intersect(house, genes)
tf = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)
tf = tf$V1
tf = intersect(tf, genes)
r1 = vector()
r2 = vector()
cols = colnames(exp)
for (g in house){
  x1=data.frame(exp_=as.numeric(exp[g, cols]), k27_=as.numeric(k27[g,cols]))
  emobj <- simple.init(x1[,c(1,2)], nclass = 1)
  emobj <- shortemcluster(x1[,c(1,2)], emobj)
  ret <- emcluster(x1[,c(1,2)], emobj, assign.class = TRUE)
  r1=c(r1, ret$llhdval)
}
names(r1)=house

for (g in tf){
  x1=data.frame(exp_=as.numeric(exp[g, cols]), k27_=as.numeric(k27[g,cols]))
  emobj <- simple.init(x1[,c(1,2)], nclass = 1)
  emobj <- shortemcluster(x1[,c(1,2)], emobj)
  ret <- emcluster(x1[,c(1,2)], emobj, assign.class = TRUE)
  r2=c(r2, ret$llhdval)
}
names(r2)=tf

table_ = data.frame(value=c(r1, r2), label=c(rep('House-keeping genes', length(house)), rep('Variably expressed TFs', length(tf))))
ggplot(table_, aes(x=value, y=..density.., fill=label, colour=label)) + geom_histogram(position='identity',alpha=0.5) + labs(x='Log likelihood', y='Density', colour='', fill='') + theme_bw() + theme(axis.text.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=12, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=12)) 
ggsave('/Users/woojunshim/Research/BN/figures/ll_comparison.pdf', width=7, height=5, dpi=600)

### BOXPLOT FOR KNN 
temp = read.table('/Users/woojunshim/Research/BN/knn_table.txt', stringsAsFactors = F)
temp$V4 = gsub('_', ' ', temp$V4)
ggplot(temp, aes(x=factor(V3), y=V2, colour=factor(V4))) + geom_boxplot(outlier.shape = 20) 
ggsave('/Users/woojunshim/Research/BN/figures/knn_table.pdf', height=10, width=10, dpi=600)

### FIND WHICH K HAS THE BETS P-VALUE
results=vector()
for (i in seq(1,46)){
  a=as.numeric(temp$V2[which(temp$V4=='House keeping genes' & temp$V3==i)])
  b=as.numeric(temp$V2[which(temp$V4=='Variably expressed TF' & temp$V3==i)])
  q=t.test(a,b)
  results=c(results, q$p.value)
}

### k27_sorted_by_exp.txt
temp = read.table('/Users/woojunshim/Research/BN/k27_sorted_by_exp.txt', stringsAsFactors = F)
#temp$V4 = gsub('_', ' ', temp$V4)
results = vector()
for (i in seq(1,46)){
  a=as.numeric(temp$V2[which(temp$V4=='House_keeping_genes' & temp$V3==i)])
  b=as.numeric(temp$V2[which(temp$V4=='Variably_expressed_TFs' & temp$V3==i)])
  q1=ci(a)
  q2=ci(b)
  results=rbind(results, c(i,as.numeric(q1[1]),as.numeric(q1[2]),mean(a),'House_keeping_gene'))
  results=rbind(results, c(i,as.numeric(q2[1]),as.numeric(q2[2]),mean(b),'Variably_expressed_TF'))
}
results = data.frame(results)
write.table(results, '/Users/woojunshim/Research/BN/ci_table.txt', sep='\t', quote=F)

table_ = read.table('/Users/woojunshim/Research/BN/ci_table.txt', stringsAsFactors = F)
ggplot(table_, aes(x=X1, y=X4, group=as.factor(X5)), colour=as.factor(X5)) + geom_line(colour='black') + geom_ribbon(aes(ymin=table_$X2, ymax=table_$X3), fill='grey70', alpha=0.5) + ylim(-0.6,0.6) + theme_bw() + labs(x='Cell types sorted by expression value', y='z-score(H3K27me3 breadth)')
ggsave('/Users/woojunshim/Research/BN/figures/k27_sorted_by_exp.pdf', width=5, height=5, dpi=600)

### CLUSTERING 
rho = read.table('/Users/woojunshim/Research/BN/data/spearman.txt', stringsAsFactors = F)
rownames(rho) = rho$V1
exp_ = read.table('/Users/woojunshim/Research/BN/data/ave_exp_table.txt', stringsAsFactors = F)
exp_ = log2(exp_+1)
k27_ = read.table('/Users/woojunshim/Research/BN/data/ave_k27_table.txt', stringsAsFactors = F)
genes = intersect(rho$V1, rownames(exp_))
genes = intersect(genes, rownames(k27_))
table_ = data.frame(k27_high=k27_[genes, 1], exp_high=exp_[genes, 1], k27_low=k27_[genes, 9], exp_low=exp_[genes, 9])
#table_ = scale(table_, center=T, scale=T)
rownames(table_) = genes

qq = hclust(dist(table_), method='ward.D')
plot(qq)
ee = cutree(qq, 5)
rb =  

library(gplots)
hmcol = colorRampPalette(c("blue","white","red"))(256)
colnames(table_) = c("K27 (High)","Exp (High)","K27 (Low)","Exp (Low)")
pdf('/Users/woojunshim/Research/BN/figures/K27_exp_hc.pdf', height=5, width=5)
heatmap.2(as.matrix(table_), distfun=dist, hclustfun=function(x) hclust(x, method = "ward.D"), trace='none', col=hmcol, labRow = F, cexCol = 0.8)
dev.off()

###CALCULATE AVERAGE LINES FOR EACH CLUSTER
temp = read.table('/Users/woojunshim/Research/BN/coordinates.txt', stringsAsFactors = F)
ref = read.table('/Users/woojunshim/Research/BN/clustering/labels.txt', stringsAsFactors = F)
results = data.frame()
groups = c('red','blue','green','yellow')
for (g in groups){
  genes = ref$V1[which(ref$V2==g)]
  idx = which(rownames(temp) %in% genes)
  x = mean(temp$x2[idx])
  y = mean(temp$y1[idx])
  results = rbind(results, c(x,y))
}

### AVERAGE AND 95% CI FOR COORDINATES
temp = read.table('/Users/woojunshim/Research/BN/coordinates_full.txt', stringsAsFactors = F)
colnames(temp) = seq(0,10)
ref = read.table('/Users/woojunshim/Research/BN/clustering/labels.txt', stringsAsFactors = F)
groups = unique(ref$V2)
means = data.frame()
mins = data.frame()
maxs = data.frame()
for (group in groups){
  rr = ref$V1[which(ref$V2 %in% group)]
  idx = which(rownames(temp) %in% rr)
  trun = temp[idx,]
  mean_ = apply(trun, 2, mean)
  sd_ = apply(trun, 2, sd)
  min_ = mean_ - 1.96*sd_/sqrt(nrow(trun))
  max_ = mean_ + 1.96*sd_/sqrt(nrow(trun))
  means = rbind(means, mean_)
  mins = rbind(mins, min_)
  maxs = rbind(maxs, max_)
}
rownames(means) = groups
rownames(mins) = groups
rownames(maxs) = groups
colnames(means) = seq(0,10)
colnames(mins) = seq(0,10)
colnames(maxs) = seq(0,10)

#table_ = melt(t(means['red',]))
#table_ = melt(t(means))
#h = ggplot(table_, aes(x=Var1, y=value, colour=Var2)) + geom_line(size=1, col='red') + scale_colour_manual(values=c('red'))
#h = h + geom_ribbon(data=table_, aes(ymin=as.numeric(mins['red',]), ymax=as.numeric(maxs['red',])), fill='grey70', alpha=0.5)
#print (h)

t1 = melt(t(means))

t2 = melt(t(mins))

t3 = melt(t(maxs))

table_ = cbind(t1,t2$value,t3$value)
colnames(table_) = c('x','cluster','mean','min','max')

ggplot(table_, aes(x=x, y=mean, colour=cluster)) + geom_line(size=1)+ geom_errorbar(aes(ymin=min, ymax=max), size=0.5, width=0.3, colour='black', alpha=0.5) + theme_bw() + scale_colour_manual(values=c('blue','green','red','yellow')) + labs(x='log2(RPKM+1)', y='H3K27me3 breadth (bp)', colour='') + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
ggsave('/Users/woojunshim/Research/BN/figures/coordinates_full.pdf', width=6, height=4, dpi=600)

ggplot(table_, aes(x=x, y=mean, ymin=min, ymax=max, group=cluster)) + geom_line() + geom_ribbon(aes(fill=cluster), alpha=0.5, show.legend = F) + theme_bw() + scale_fill_manual(values=c('blue','green','red','yellow')) + labs(x='log2(RPKM+1)', y='H3K27me3 breadth (bp)', colour='') + theme(axis.text.x=element_text(size=16, face='bold'), axis.text.y=element_text(size=16, face='bold'),legend.title=element_text(size=16, face='bold'), legend.text=element_text(size=16)) 
ggsave('/Users/woojunshim/Research/BN/figures/coordinates_full.pdf', width=6, height=4, dpi=600)

### PLOT DENSITY HEATMAP PLOTS FOR EACH GENE CLUSTER
ref = read.table('/Users/woojunshim/Research/BN/clustering/labels.txt', stringsAsFactors = F)
temp = read.table('/Users/woojunshim/Research/BN/data/all_data.tsv', stringsAsFactors = F)
groups = unique(ref$V2)

library(RColorBrewer)
cols = rev(brewer.pal(9, 'YlOrRd'))
cols = rev(brewer.pal(11, 'Spectral'))

aa = ref$V1[which(ref$V2=='red')]
qq = temp[which(temp$V3 %in% aa),]

group = 'red'
id_ = ref$V1[which(ref$V2==group)]
temp$V3[which(temp$V3 %in% id_)] = rep('category 1', length(which(temp$V3 %in% id_)))
group = 'yellow'
id_ = ref$V1[which(ref$V2==group)]
temp$V3[which(temp$V3 %in% id_)] = rep('category 2', length(which(temp$V3 %in% id_)))
group = 'green'
id_ = ref$V1[which(ref$V2==group)]
temp$V3[which(temp$V3 %in% id_)] = rep('category 3', length(which(temp$V3 %in% id_)))
group = 'blue'
id_ = ref$V1[which(ref$V2==group)]
temp$V3[which(temp$V3 %in% id_)] = rep('category 4', length(which(temp$V3 %in% id_)))
temp = temp[which(temp$V3 %in% c('category 1','category 2','category 3','category 4')),]
temp = temp[temp$V2!=0,]

ggplot(temp) + geom_hex(aes(x=V1, y=log2(V2), fill=..density..)) + scale_fill_gradientn(colours = cols)+ labs(xlab='log2(RPKM)', ylab='log2(H3K27me3)')+ theme_bw(base_size = 16) + facet_grid(facets=. ~ V3) + stat_smooth(data = temp, method = "lm", aes(x=V1, y=log2(V2)))
ggsave('/Users/woojunshim/Research/BN/figures/density_map.pdf', width=10, height=6, dpi=600)

### HISTOGRAM FOR LOG2(H3K27ME3) GIVEN RPKM < 1
ww=temp[which(temp$V1<1), c(2,3)]
ggplot(ww, aes(x=log2(V2), y=..density.., fill=V3, colour=V3)) + geom_histogram(position='identity',alpha=0.5,bins=50) + scale_fill_manual(values=c('red','yellow','green','blue')) + theme_bw(base_size = 16) + labs(x='log2(H3K27me3)', fill='', colour='')
ggsave('/Users/woojunshim/Research/BN/figures/H3K27me3_histogram_given_exp_less_than_1.pdf', width=8, height=5, dpi=600)

### HEATMAP FOR TOP 10 GO TERMS FOR 4 GROUPS
temp = read.table('/Users/woojunshim/Research/BN/clustering/top10_terms.txt', stringsAsFactors = F)
terms = read.table('/Users/woojunshim/Research/BN/clustering/selected_terms.txt', stringsAsFactors = F)
terms = gsub('_', ' ', terms$V1)

rownames(temp) = gsub('_', ' ', rownames(temp))
colnames(temp) = gsub('_', ' ', colnames(temp))
table_ = melt(t(temp))
ggplot(table_, aes(y=Var2, x=Var1, fill=value)) + geom_tile(colour='black')+  labs(x='', y='') +scale_fill_gradient(low = "white", high = "steelblue") + theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill='-log10(FDR)') 
ggsave('/Users/woojunshim/Research/BN/figures/top10_terms.pdf', width=9, height=10, dpi=600)

### mclust
library(mclust)
exp=read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
k27=read.table('/Users/woojunshim/Research/intergenic/data/peak_assigned/H3K27me3_broadest_table.txt', stringsAsFactors = F)
k4=read.table('/Users/woojunshim/Research/intergenic/data/peak_assigned/H3K4me3_broadest_table.txt', stringsAsFactors = F)
cols = colnames(exp)
genes = intersect(rownames(exp), rownames(k27))
genes = intersect(genes, rownames(k4))

c='E070'
x1=data.frame(exp_=log2(as.numeric(exp[genes, c]+1)), k27_=log2(as.numeric(k27[genes,c]+1)), k4_=log2(as.numeric(k4[genes,c]+1)))
rownames(x1) = genes
fit = Mclust(x1, G=1:20)
summary(fit)
plot(fit, what='classification')
fit$classification[which(rownames(x1)=='NKX2-5')]

#

### SIMILARIY ANALYSIS 
labels = read.table('/Users/woojunshim/Research/modelling/clustering/labels.txt', stringsAsFactors = F)
genes = labels$V1[which(labels$V2=='red')]
temp = exp[which(rownames(exp) %in% genes),]
pca = prcomp(as.matrix(temp))

### mclust
gene = 'GAPDH'
x1=data.frame(exp_=as.numeric(exp[gene,cols]), k27_=as.numeric(k27[gene,cols]))
rownames(x1) = colnames(exp)
fit_house = Mclust(x1, G=1:2)
fit_regulatory = Mclust(x1, G=1:2)
fit_structural = Mclust(x1, G=1:2)

### GO ANALYSIS ON TOP 50 HIGHLY EXPRESSED VEG
exp = read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F)
veg = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)$V1
pathway = '/Users/woojunshim/Research/modelling/data/top50_veg/'
for (c in cols){
  all_ = rownames(exp)[which(exp[,c] > 0)]
  input_ = sort(exp[all_,c], decreasing=T)
  names(input_) = all_
}

### HISTOGRAM 
top5 = read.table('/Users/woojunshim/Research/intergenic/data/location/top5_counts.txt', stringsAsFactors = F)
top1 = read.table('/Users/woojunshim/Research/intergenic/data/location/top1_counts.txt', stringsAsFactors = F)
pdf('/Users/woojunshim/Research/intergenic/figures/broadest.pdf', width=5, height=5)
plot(1, type="n", xlab="", ylab="", xlim=range(top5$V1), ylim=range(top5$V2, top1$V2))
points(x=top5$V1, y=top5$V2, col='blue')
points(x=top1$V1, y=top1$V2, col='red')
dev.off()

### CUMULATIVE ANALYSIS
ref = read.table('/Users/woojunshim/Research/intergenic/data/result/H3K27me3_all_bonferroni_final.bed', stringsAsFactors = F)  
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/cumulative_table.txt', stringsAsFactors = F)
colnames(temp) = gsub('X', '', colnames(temp))

# HISTOGRAM
ggplot(ref, aes(x=V5)) + geom_histogram(bins=20, colour='black', fill='blue', aes(y=..density..)) + theme_bw(base_size=16) + labs(x='RTS') + annotate('text', x=80, y=0.07, label='Total number of enriched regions = 338,354')
ggsave('/Users/woojunshim/Research/intergenic/figures/RTS_hist.pdf', width=6, height=5, dpi=600)

ggplot(ref, aes(x=V7)) + geom_histogram(bins=100, colour='black', fill='blue', aes(y=..density..)) + theme_bw(base_size=16) + labs(x='Region breadth (bp)') + annotate('text', x=200000, y=0.00017, label='Total number of enriched regions = 338,354')
ggsave('/Users/woojunshim/Research/intergenic/figures/region_breadth_hist.pdf', width=6, height=5, dpi=600)

# HEATMAP
a = data.frame(RTS=ref$V6, Breadth=ref$V8)
b = data.frame(table(a) / nrow(a))
ggplot(b, aes(x=RTS, y=Breadth, fill=Freq)) + geom_tile() + scale_fill_gradient(low='white', high='red') + labs(x='RTS rank', y='Breadth rank', fill='density') + scale_x_discrete(breaks=seq(0,100,25)) + scale_y_discrete(breaks=seq(0,100,25))
ggsave('/Users/woojunshim/Research/intergenic/figures/RTS_vs_breadth_without_grid.pdf', width=8, height=6)

# CUMULATIVE GRAPHS
# UP TO TOP 10 % REGIONS 
no = 10
result = data.frame(matrix(nrow=no, ncol=100))
colnames(result) = seq(1,100)

for (i in (1:no)){
    ids = ref$V4[which(ref$V8==i)]
    a = temp[which(rownames(temp) %in% ids),]
    b = apply(a, 2, mean)
    result[i,names(b)] = b
}
table_ = melt(t(result))
ggplot(table_, aes(x=Var1, y=value, colour=as.factor(Var2))) + geom_line(size=1) + geom_abline(intercept=0, slope=0.01, colour='black', linetype='dashed', size=1) + theme_bw(base_size = 10) + labs(x='H3K27me3 peaks included (rank threshold)', y='RTS recovered (proportion)', colour='Region breadth (rank)')
ggsave('/Users/woojunshim/Research/intergenic/figures/RTS_recovered_by_region_breadth.pdf', width=6, height=4, dpi=600)

result = data.frame(matrix(nrow=10, ncol=100))
colnames(result) = seq(1,100)
mins = seq(1,91,10)
maxs = seq(10,100,10)
for (i in (1:10)){
  ids = ref$V4[which(ref$V8>=mins[i] & ref$V8<=maxs[i])]
  a = temp[which(rownames(temp) %in% ids),]
  b = apply(a, 2, mean)
  result[i,names(b)] = b
}
table_ = melt(t(result))
ggplot(table_, aes(x=Var1, y=value, colour=as.factor(Var2))) + geom_line(size=1) + geom_abline(intercept=0, slope=0.01, colour='black', linetype='dashed', size=1) + theme_bw(base_size = 10) + labs(x='H3K27me3 peaks included (rank threshold)', y='RTS recovered (proportion)', colour='Region breadth (rank)')
ggsave('/Users/woojunshim/Research/intergenic/figures/RTS_recovered_by_region_breadth_overall.pdf', width=6, height=4, dpi=600)

# OR HEATMAP
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low='white', high='blue') + labs(x='H3K27me3 peaks included (rank threshold)', y='RTS rank', fill='density') + scale_x_discrete(breaks=seq(0,100,25)) + scale_y_discrete(breaks=seq(0,100,25))

#mins = seq(1,91,10)
#maxs = seq(10,100,10)
#for (i in (1:10)){
#  ids = ref$V4[which(ref$V8>=mins[i] & ref$V8<=maxs[i])]
#  a = temp[which(rownames(temp) %in% ids),]
#  b = apply(a, 2, mean)
#  result[i,names(b)] = b
#}


           
library(RColorBrewer)
cols = rev(brewer.pal(11, 'Spectral'))
table_ = data.frame(V1=ref$V5, V2=ref$V6)
ggplot(table_) + geom_hex(aes(x=V1, y=V2, fill=..density..)) + scale_fill_gradientn(colours = cols)+ theme_bw(base_size = 16) + stat_smooth(data = table_, method = "lm", aes(x=V1, y=V2)) + labs(x='RTS', y='Breadth') + annotate('text', x=50, y=350000, label="Spearman's rho = 0.934", colour='black', size=3)
# Pearson's r = 0.6812486
ggsave('/Users/woojunshim/Research/intergenic/figures/RTS_vs_breadth.pdf', width=6, height=5, dpi=600)

### CORRELATION BETWEEN RTS
old = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
new = read.table('/Users/woojunshim/Research/intergenic/data/tss_2.5kb_score.txt', stringsAsFactors = F)
rownames(old) = old$V1
rownames(new) = new$V1
genes = intersect(old$V1, new$V1)
old = old[genes,]
new = new[genes,]
cor(old$V2, new$V2, method='spearman')

### FISHER'S EXACT TEST
x=matrix(c(3510,30330,9710,294804), nrow=2)  # top 10 % scores, proximal region
fisher.test(x, alternative='greater')

### TEST FOR MATRIX FACTORISATION
a=matrix(c(1,1,0,1,1,1,0,1,1), nrow=3,ncol=3)

### DENSITY PLOTS 
temp = read.table('/Users/woojunshim/Research/modelling/data/predictive/vet_5kb_table_cage.txt', stringsAsFactors = F)
colnames(temp) = seq(-5000, 4999)

table_ = melt(t(temp))
ggplot(table_) + geom_hex(aes(x=Var1, y=value, fill=..density..)) + scale_fill_gradientn(colours = cols)+ theme_bw(base_size = 16) + labs(x='', y='H3K27me3 count', title='Heart development (GO:007507)') + scale_x_continuous(breaks=c(-5000,-4000,-3000,-2000,-1000,0,1000,2000,3000,4000,4999), labels=c('-5kb','-4kb','-3kb','-2kb','-1kb','TSS','1kb','2kb','3kb','4kb','5kb')) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/woojunshim/Research/modelling/figures/GO_007507_density.pdf', width=10, height=8, dpi=600)

### AVERAGE + 95% CI 

temp = read.table('/Users/woojunshim/Research/modelling/data/predictive/GO.0032502_5kb_table.txt', stringsAsFactors = F)
colnames(temp) = seq(-5000, 4999)
mean_=apply(temp, 2, mean)
sd_ = apply(temp, 2, sd)
low_ = mean_ - 1.96 * sd_ / sqrt(nrow(temp))
high_ = mean_ + 1.96 * sd_ / sqrt(nrow(temp))
a = data.frame(mean_=mean_, high_=high_, low_=low_, label=rep('Developmental process (GO:0032502)'), x=seq(-5000, 4999))

temp = read.table('/Users/woojunshim/Research/modelling/data/predictive/GO.0030198_5kb_table.txt', stringsAsFactors = F)
colnames(temp) = seq(-5000, 4999)
mean_=apply(temp, 2, mean)
sd_ = apply(temp, 2, sd)
low_ = mean_ - 1.96 * sd_ / sqrt(nrow(temp))
high_ = mean_ + 1.96 * sd_ / sqrt(nrow(temp))
b = data.frame(mean_=mean_, high_=high_, low_=low_, label=rep('Extracellular matrix organisation (GO:0030198)'), x=seq(-5000, 4999))

temp = read.table('/Users/woojunshim/Research/modelling/data/predictive/GO.0016071_5kb_table.txt', stringsAsFactors = F)
colnames(temp) = seq(-5000, 4999)
mean_=apply(temp, 2, mean)
sd_ = apply(temp, 2, sd)
low_ = mean_ - 1.96 * sd_ / sqrt(nrow(temp))
high_ = mean_ + 1.96 * sd_ / sqrt(nrow(temp))
c = data.frame(mean_=mean_, high_=high_, low_=low_, label=rep('mRNA metabolic process (GO:0016071)'), x=seq(-5000, 4999))

temp = read.table('/Users/woojunshim/Research/modelling/data/predictive/GO.0007507_5kb_table.txt', stringsAsFactors = F)
colnames(temp) = seq(-5000, 4999)
mean_=apply(temp, 2, mean)
sd_ = apply(temp, 2, sd)
low_ = mean_ - 1.96 * sd_ / sqrt(nrow(temp))
high_ = mean_ + 1.96 * sd_ / sqrt(nrow(temp))
d = data.frame(mean_=mean_, high_=high_, low_=low_, label=rep('Heart development (GO:0007507)'), x=seq(-5000, 4999))

table_ = rbind(a,b,c,d)

ggplot(table_, aes(x=x, y=mean_, ymin=low_, ymax=high_, group=label)) + geom_line() + geom_ribbon(aes(fill=label), alpha=0.5) + theme_bw(base_size=16) + labs(x='', y='H3K27me3 count' ) + scale_x_continuous(breaks=c(-5000,-4000,-3000,-2000,-1000,0,1000,2000,3000,4000,4999), labels=c('-5kb','-4kb','-3kb','-2kb','-1kb','TSS','1kb','2kb','3kb','4kb','5kb')) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/woojunshim/Research/modelling/figures/aa.pdf', width=10, height=8, dpi=600)

### TSS QUANTIFICATION
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K27me3_quantification_table_refseq.txt', stringsAsFactors = F)
ref = read.table('/Users/woojunshim/Research/intergenic/data/go/GO.0016071_cage.txt', stringsAsFactors = F)
ref = ref$V1
table_ = temp[ref,]
order_ = read.table('/Users/woojunshim/Research/intergenic/data/epigenome_list.txt', stringsAsFactors = F)$V1
#ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient(low='white', high='red') + scale_x_discrete(limits=order_$V1)
table_ = table_[,order_$V1]
ww = apply(table_, 2, mean)
ee = apply(table_, 2, sd)
plot(x=seq(1,length(ww)), y=ww)

### TABLE TEST
tf = read.table('/Users/woojunshim/Research/Scripts/TF_combined.txt', stringsAsFactors = F)$V1
ref = read.table('/Users/woojunshim/Research/intergenic/data/tss_2.5kb_score_top10.txt', stringsAsFactors = F)$V1
k27_ = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K27me3_quantification_table_refseq.txt', stringsAsFactors = F)
k4_ = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K4me3_quantification_table_refseq.txt', stringsAsFactors = F)
k9_ = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K9me3_quantification_table_refseq.txt', stringsAsFactors = F)
order_ = read.table('/Users/woojunshim/Research/intergenic/data/epigenome_list.txt', stringsAsFactors = F)$V1
labels = read.table('/Users/woojunshim/Research/intergenic/data/tissue_list.txt', stringsAsFactors = F)$V1
term = 'GO.0007507'
go = read.table(paste('/Users/woojunshim/Research/Data/genes_with_go_term/',term,'.txt',sep=''), stringsAsFactors = F)$x
genes = intersect(ref,go)
genes = intersect(genes,tf)
k27 = k27_[genes,]

table_ = melt(t(k27))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_x_discrete(limits=order_, labels=labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + labs(x='', y='', fill='H3K27me3') 
ggsave(paste('/Users/woojunshim/Research/intergenic/figures/',term,'_top10_tf_H3K27me3.pdf',sep=''), width=6, height=6, dpi=600)

### MLE, with assumption of Gaussian distribution
gos = c('0007507', '0007420', '0030097')

### HEATMAP
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/genes/GO.0007420_tf_combined.txt', stringsAsFactors = F)
temp = temp / max(temp)
rownames(temp) = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/genes/epi_list_GO.0007420_tf.txt', stringsAsFactors = F)$V1
order_ = read.table('/Users/woojunshim/Research/intergenic/data/epigenome_list.txt', stringsAsFactors = F)$V1
labels = read.table('/Users/woojunshim/Research/intergenic/data/tissue_list.txt', stringsAsFactors = F)$V1
colnames(temp) = seq(-2500,2500)
table_ = melt(t(temp))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_y_discrete(limits=order_, labels=labels) + labs(x='Distance from TSS', y='', fill='', title='Brain development TFs')
ggsave('/Users/woojunshim/Research/intergenic/data/result/predictive/genes/GO.0007420_tf_combined.txt.pdf', width=10, height=8, dpi=600)

### SVD VISUALISATION
name = 'eigen_cell_'
for (i in seq(1,111)){
  temp = read.table(paste('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/patterns_activated/',name,i,'_activated.txt', sep=''), stringsAsFactors = F)
  temp = temp / max(temp)
  rownames(temp) = read.table(paste('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/patterns_activated/epi_list_',name,i,'_activated.txt', sep=''), stringsAsFactors = F)$V1
  order_ = read.table('/Users/woojunshim/Research/intergenic/data/epigenome_list.txt', stringsAsFactors = F)$V1
  labels = read.table('/Users/woojunshim/Research/intergenic/data/tissue_list.txt', stringsAsFactors = F)$V1
  colnames(temp) = seq(-2500,2500)
  table_ = melt(t(temp))
  ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_y_discrete(limits=order_, labels=labels) + labs(x='Distance from TSS', y='', fill='', title=paste(name,i,' top 5% activated genes', sep=''))
  ggsave(paste('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/figures/activated/final_',name,i,'_activated.pdf', sep=''), width=10, height=8, dpi=600)
}

### GO FUNCTIONAL ANALYSIS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K27me3_tissue_specific.txt', stringsAsFactors = F)
cols = colnames(temp)
background = rownames(temp)
for (c in cols){
  input_=background[which(temp[,c]==1)]
  result = go_analysis(input_, background)
  write.table(result, paste('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/go/',c,'_GO.txt',sep=''), sep='\t', quote=F)
}

### HEAT FOR 5 TISSUE GROUPS
table_ = temp[,1:5]
idx = which(apply(table_, 1, sum)>0)
table_ = table_[idx,]
idx = vector()
for (i in 1:ncol(table_)){
  temp = which(table_[,i]==1)
  idx = c(idx, temp)
}
idx = unique(idx)
table_ = table_[idx,]
table = melt(t(table_))
ggplot(table, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +scale_y_discrete(limits=rownames(table_))
ggsave('/Users/woojunshim/Research/intergenic/data/result/predictive/class/significance_heatmap.pdf', width=8, height=8, dpi=600)

#### 
k27_ = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/tss_H3K27me3_quantification_table_refseq.txt', stringsAsFactors = F)
ref = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/class/tss_H3K27me3_tissue_specific.txt', stringsAsFactors = F)
genes = rownames(ref)[which(ref$Brain==1)]
genes = intersect(genes, rownames(k27_))
table_ = melt(t(k27_[genes,]))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_x_discrete(limits=order_)

q = svd(k27_)
result = q$u
rownames(result) = rownames(k27_)
colnames(result) = paste('eigen_cell_',seq(1,111),sep='')
write.table(result, '/Users/woojunshim/Research/intergenic/data/result/predictive/svd/H3K27me3_svd_u_pre_final.txt', sep='\t', quote=F)
table_ = melt(t(result[genes,]))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_x_discrete(limits=order_)

### GO ANALYSIS FOR EACH EIGEN CELL TYPE
result = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/H3K27me3_svd_u_final.txt', stringsAsFactors = F)
all_ = rownames(result)
for (c in colnames(result[,40:ncol(result)])){
  input_ = rownames(result)[order(result[,c], decreasing=T)[1:921]] # top 5% genes in an ascending order
  zz = go_analysis(input_, all_)
  write.table(zz, paste('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/go_activated/go_',c,'_final_activated.txt', sep=''), sep='\t', quote=F)
}

### SVD VISUALISATION 
ref = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/top5_repressed_binary_svd_u.txt', stringsAsFactors = F)
col = 'eigen_cell_1'
genes = rownames(ref)[which(ref[,col]==1)]

### HEATMAP FOR TOP 20 TERMS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/go_activated/go_top_20_merged_short_activated.txt', stringsAsFactors = F)
temp = -log10(temp)
table_ = melt(t(temp))
aa=table_[order(table_$Var1),]
idx = vector()
for (i in (1:20)){
  aa = aa[order(aa$value, decreasing=T),]
  qq = which(aa$Var1==paste('eigen_feature_',i,sep=''))
  idx = c(idx, qq)
}
terms = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/go_activated/go_top_20_merged_short_terms_activated.txt',stringsAsFactors = F)$V1
terms = unique(terms)
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient(low='white', high='red') + scale_x_discrete(limits=paste('eigen_feature_',seq(1,20),sep='')) + scale_y_discrete(limits=terms) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='', y='', fill='-log10(FDR)')
ggsave('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/figures/activated/go_top_20_merged_short_activated.pdf', width=10, height=30, dpi=600)

### TEST SVD FOR H3K27ME3 WIDTH TABLE
### LINE PLOT FOR EXAMPLE GENES
temp = read.table('/Users/woojunshim/Research/Data/new_assigned/H3K27me3_width_table.txt', stringsAsFactors = F)
q = svd(temp)
nkx2.5 = q$u[which(rownames(temp)=='NKX2-5'),]
gata4 = q$u[which(rownames(temp)=='GATA4'),]
tnni3 = q$u[which(rownames(temp)=='TNNI3'),]
myh7 = q$u[which(rownames(temp)=='MYH7'),]
gapdh = q$u[which(rownames(temp)=='GAPDH'),]
gtf2b = q$u[which(rownames(temp)=='GTF2B'),]
table_ = data.frame(NKX2.5=nkx2.5, GATA4=gata4, TNNI3=tnni3, MYH7=myh7, GAPDH=gapdh, GTF2B=gtf2b)
table = melt(t(table_))
ggplot(table, aes(x=Var2, y=value, colour=Var1)) + geom_line() + geom_point() + labs(x='Eigen cell-type', y='H3K27me3', colour='')
ggsave('/Users/woojunshim/Research/intergenic/figures/H3K27me3_width_svd_selected_cardiac.pdf', width=7, height=5, dpi=600)

w=data.frame(q$u)
rownames(w) = rownames(temp)
write.table(w, '/Users/woojunshim/Research/intergenic/data/H3K27me3_width_u.txt', sep='\t', quote=F)

# CORRELATION BETWEEN SVD VS ORIGINAL RTS
# SPEARMAN'S RHO = 0.906
a = read.table('/Users/woojunshim/Research/intergenic/data/H3K27me3_width_svd_top_bottom_5_percent_counts.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/Data/new_assigned/tables/repressive_hg19_cont_prop.txt', stringsAsFactors = F)
genes = intersect(a$V1, b$V1)
rownames(a) = a$V1
rownames(b) = b$V1
a = a[genes, ]
b = b[genes, ]
cor_ = cor(a$V2, b$V2, method='spearman')

### GO BP ENRICHMENT FOR LOW-H3K27ME3-GENES USING ORIGIANL ROADMAP CELL-TYPES
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/tss_H3K27me3_quantification_table_refseq.txt', stringsAsFactors = F)
genes = rownames(read.table('/Users/woojunshim/Research/Data/46epigenomes.RPKM.symbols.txt', stringsAsFactors = F))
genes = intersect(genes, rownames(temp))
temp = temp[genes,]
cols = colnames(temp)
cols = c('E095','E083','E070','E100','E038')
for (c in cols){
  input_ = rownames(temp)[order(temp[,c], decreasing=F)[1:921]] # top 5% genes in an ascending order
  zz = go_analysis(input_, genes)
  write.table(zz, paste('/Users/woojunshim/Research/intergenic/data/result/tss_quantification/go/repressed/go_',c,'_repressed_original.txt', sep=''), sep='\t', quote=F)
}

### PLOTS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/top5_activated_svd_u_final.txt', stringsAsFactors = F)
qq = apply(temp, 1, sum)
pdf('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/figures/repressed/hist_activated_H3K27me3_state.pdf', width=5, height=4)
hist(qq, breaks=25, xlab='Number of eigen cell-type', ylab='Gene count', main='Association with activated H3K27me3 state')
dev.off()

### FISHER'S EXACT TEST TO SEE IF ANY GROUPS ARE ENRICHED WITH REF GENES
table_ = table(qq)
result = vector()
ref = read.table('/Users/woojunshim/Research/Data/house_keeping_genes.txt', stringsAsFactors = F)$V1
for (i in (1:25)){
  input_ = names(qq)[which(qq==i)]
  a = length(which(input_ %in% ref))
  b = length(input_) - a
  input_ = names(qq)[which(qq!=i)]
  c = length(which(input_ %in% ref))
  d = length(input_) - c
  result = c(result, fisher.test(matrix(c(a,c,b,d), nrow=2), alternative = 'greater')$p.value)
}

### LINE PLOT FOR H3K27ME3 ACROSS EIGEN CELL-TYPES
k27 = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/tss_H3K27me3_quantification_table_refseq_coding_only.txt', stringsAsFactors = F)
table_ = data.frame(matrix(nrow=111,ncol=3))
colnames(table_)=paste('eigen_cell_',seq(1,3), sep='')
order_ = read.table('/Users/woojunshim/Research/intergenic/data/epigenome_list.txt', stringsAsFactors = F)$V1
labels = read.table('/Users/woojunshim/Research/intergenic/data/tissue_list.txt', stringsAsFactors = F)$V1
rownames(table_) = order_
ref = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/top5_repressed_svd_u_final.txt', stringsAsFactors = F)
for (c in colnames(table_)){
  pos = rownames(ref)[which(ref[,c]==1)]
  temp = k27[pos, order_]
  ww = apply(temp, 2, mean)
  table_[names(ww),c] = ww
}
table = melt(t(table_))
ggplot(table, aes(x=Var2, y=value, colour=Var1, group=Var1)) + geom_line() + scale_x_discrete(limits=order_, labels=labels) + labs(x='', y='H3K27me3 deposition', colour='', title='Low-H3K27me3 genes') + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 
ggsave('/Users/woojunshim/Research/intergenic/data/result/predictive/svd/figures/repressed/low_H3K27me3_first_3_eigen_cells.pdf', width=7, height=5, dpi=600)

### HEATMAP (K27) FOR TOP 20% GENES FOR SELECTED GO TERMS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/predictive/tss_H3K27me3_quantification_table_refseq_top_20.txt', stringsAsFactors = F)
refs = c('GO.0007507','GO.0007420','GO.0030097','GO.0007398','GO.0007498','GO.0007492','GO.0043588','GO.0001822')
names = c('Heart_development','Brain_development','Hemopoiesis','Ectoderm_development','Mesoderm_development','Endoderm_development','Skin_development','Kidney_development')
for (i in 1:length(refs)){
  ref = read.table(paste('/Users/woojunshim/Research/Data/genes_with_go_term/',refs[i],'.txt',sep=''),stringsAsFactors = F)$x
  ww = temp[which(rownames(temp) %in% ref),]
  table_ = melt(t(ww))
  ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + labs(x='', y='', value='H3K27me3', main=names[i]) + scale_x_discrete(limits=order_, labels=labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 
  ggsave(paste('/Users/woojunshim/Research/intergenic/data/result/predictive/top_20_percent/figure/',refs[i],'_top_20_percent.pdf',sep=''), width=7,height=10, dpi=600)
}

### BOXPLOT FOR FANTOM5 RTS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/H3K27me3_rts_FANTOM5_selected_genes.txt', stringsAsFactors = F)
ggplot(temp, aes(x=V2, y=V1, colour=V2, group=V2)) + geom_boxplot() + labs(x='', y='RTS', colour='')
ggsave('/Users/woojunshim/Research/intergenic/figures/FANTOM5_rts_selected_genes.pdf', width=5, height=3, dpi=600)

### CORRELATION BETWEEN ENCODE AND ROADMAP
vet = read.table('/Users/woojunshim/Research/Data/regulated_tf.txt', stringsAsFactors = F)$V1
house = read.table('/Users/woojunshim/Research/Data/house_keeping_genes.txt', stringsAsFactors = F)$V1
a= read.table('/Users/woojunshim/Research/intergenic/data/result/H3K27me3_rts_Roadmap_hg19.txt', stringsAsFactors = F)
b=read.table('/Users/woojunshim/Research/intergenic/data/result/H3K27me3_rts_ENCODE_hg19.txt', stringsAsFactors = F)
rownames(a) = a$V1
rownames(b) = b$V1
genes = intersect(a$V1, b$V1)
q=a[genes,2]
w=b[genes,2]
q1=(q-min(q))/(max(q)-min(q))
w1=(w-min(w))/(max(w)-min(w))
names(q1) = genes
names(w1) = genes
idx = which(names(q1) %in% intersect(genes, house))
cor(q1[idx],w1[idx],method='spearman')
table_ = data.frame(roadmap=q1[idx], encode=w1[idx])
ggplot(table_, aes(x=roadmap, y=encode)) + geom_point(colour='blue',alpha=0.3) + annotate(geom='text',x=0.25, y=0.9, label="Spearman's rho = 0.32 (n=3,500)")
ggsave('/Users/woojunshim/Research/intergenic/figures/scatter_plot_rts_roadmap_vs_encode_house.pdf', width=5, height=5, dpi=600)

### COMPARE RTS FOR FANTOM5 
a = read.table('/Users/woojunshim/Research/intergenic/data/result/differentially_expressed_enhancers_rts.txt', stringsAsFactors = F)
b = read.table('/Users/woojunshim/Research/intergenic/data/result/ubiquitous_enhancers_rts.txt', stringsAsFactors = F)
table_ = data.frame(value=c(a$V2, b$V2), label=c(rep('DE', nrow(a)), rep('Ubiquitous', nrow(b))))
ggplot(table_, aes(x=label, y=value, colour=label, group=label)) + geom_boxplot() + labs(x='', y='RTS', colour='', title='FANTOM5 enhancers (p=8.95e-16, Wilcoxon ranksum)') + scale_x_discrete(labels=c('DE (n=14,058)','Ubiquitous (n=446)'))
ggsave('/Users/woojunshim/Research/intergenic/figures/FANTOM5_rts_enhancers.pdf', width=5, height=4, dpi=600)

### HEATMAP FOR FANTOM5 DIFFERENTIALLY EXPRESSED ENHANCERS
temp = read.table('/Users/woojunshim/Research/intergenic/data/result/blood_specific_enhancers_H3K27me3_table.txt', stringsAsFactors = F)
table_ = melt(t(temp))
ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_y_discrete(labels=rep('', nrow(temp))) + scale_x_discrete(limits=order_, labels=labels) + labs(x='', y='') + scale_fill_gradient2(low='white', high='red') + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 
ggsave('/Users/woojunshim/Research/intergenic/figures/FANTOM5_blood_specific_enhancers.pdf', width=10, height=8, dpi=600)
