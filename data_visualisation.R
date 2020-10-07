### GO enrichment analysis
library(topGO)
library(GO.db)
library(org.Hs.eg.db)
go_table = read.table('GO_table.txt')  # requires GO to description table
rownames(go_table) = go_table$V1

go_analysis <- function(input, all, onto='BP', stat='fisher'){
  geneList <- factor(as.integer(all %in% input))
  names(geneList)<-all
  TopGOdata<-new("topGOdata",ontology=onto,allGenes=geneList,geneSel=genelist,nodeSize=5,annot=annFUN.org,mapping="org.Hs.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = "classic", statistic = stat)
  aa = score(resultFisher)
  #aa = p.adjust(aa, method='BH')  # FDR
  aa = sort(aa, decreasing=FALSE)  # sort by FDR
  bb = as.character(go_table[names(aa),2])
  top_no = length(which(aa<=0.05))
  results = matrix(ncol=3,nrow=top_no)  # result table 
  results[,1] = names(aa)[1:top_no]
  results[,2] = bb[1:top_no]
  results[,3] = as.numeric(aa)[1:top_no]
  colnames(results) = c('#GO_ID','GO_Term','p_value')
  return(results)
}

### ROC plot
plot_performance_roc <- function(x_matrix, y_matrix, xlab, ylab, title, auc){
  x_ = melt(t(x_matrix))
  y_ = melt(t(y_matrix))
  x_$value_y = y_$value
  ggplot(x_, aes(x=value, y=value_y, colour=Var1, group=Var1)) + geom_line(size=1) + labs(x=xlab, y=ylab, title=title, colour='', group='')  + theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold')) 
}

### PRC plot
plot_performance_prc <- function(x_matrix, y_matrix, xlab, ylab, title, no_positives, baseline=T){
  if (baseline==T){
    interval = nrow(x_matrix)
    x_matrix$Random = seq(0,1,length.out=interval)
    y_matrix$Random = rep(no_positives / nrow(x_matrix), interval)
  }
  x_ = melt(t(x_matrix))
  y_ = melt(t(y_matrix))
  x_$value_y = y_$value
  ggplot(x_, aes(x=value, y=value_y, colour=Var1)) + geom_line(size=1, alpha=0.7) + labs(x=xlab, y=ylab, title=title, colour='', group='') + ylim(0,1)+ theme_bw() + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'))
}

### Extract genes with a given GO term
library(GO.db)
go_id = GOID( GOTERM[ Term(GOTERM) == "sarcomere"])

get_go_terms <- function(go_id){
  allegs = get(go_id, org.Hs.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Hs.egSYMBOL))
  genes = unique(genes)
  return (genes)
}

### Figures 2B, 3F, 5E, Enrichment plot
library(corrplot)
sig_plot <- function(input){
  corrplot(as.matrix(input), is.corr=F)
}

### CALCULATE CORRELATION BETWEEN DIS AND EXP FOR EACH CELL-TYPE
cor_btw_tables <- function(t1,t2,method='pearson'){
  genes = rownames(t1)
  results = vector()
  cols = colnames(t1)
  for (c in cols){
    l1 = t1[genes, c]
    l2 = t2[genes, c]
    results = c(results, cor(l1,l2,method=method))
  }
  names(results) = cols
  return (results)
}

### PLOT HEATMAP 
library(ggplot2)
library(reshape2)
plot_heatmap <- function(input, title, legend, value_range){
  table_ = melt(t(input))
  ggplot(table_, aes(x=Var1, y=Var2, fill=value)) + geom_tile(colour='black') + scale_fill_gradient(low='white',high='red',limits=value_range) + geom_text(aes(label=round(value,3))) + labs(x='', y='', title=title, fill=legend) + theme_bw() +theme(legend.text=element_text(size=16,face="bold"), axis.title=element_text(size=16, face='bold'), axis.text=element_text(siz=16, face='bold'))
}

### STACKED BAR GRAPH Example
# STACKED H3K27ME3 DOMAINS (LIKE FIGURE 1D)
order_ = read.table('heart_sample_prioritised_list.txt', stringsAsFactors = F)$V1
heart_id = order_[107:111]
hm = c('H3K27me3','H3K4me3','H3K27ac')
genes = c('NKX2-5','IRX4','GATA5')
tsss = c(172662315,1887293,61051026)
for (h in hm){
  temp = read.table(paste(h,'_selected_genes_overlaps.txt', sep=''), stringsAsFactors = F)
  for (i in 1:length(genes)){
    gene = genes[i]
    tss = tsss[i]
    table_ = temp[temp$V4==gene,]
    table_$start = tss - table_$V3
    table_$end = tss - table_$V2
    sample = intersect(order_, table_$V5)
    cols = rep('dark blue', length(sample))
    cols[sample %in% heart_id] = 'red'
    rownames(table_) = table_$V5
    table_ = table_[sample,]
    pdf(paste(h,'_',gene,'_stacked_bar_sroted_by_tissue.pdf',sep=''), width=2.5, height=5)
    plot(1, type="n", xlab="", ylab="", main=gene, xaxt = 'n', yaxt = 'n', xlim=c(-2500, 25000), ylim=c(0, 111))
    for (i in 1:length(sample)){
      x_line = c(table_[i,6], table_[i,7])
      lines(x=x_line, y=c(i,i), col=cols[i])
    }
    axis(1, at=c(-2500,0,25000), labels=c('-2.5kb','TSS','+25kb'))
    dev.off()
  }
}

### Fisher's exact test plot
plot_fet <- function(matrix, xlab, ylab, title, legend=TRUE){
  table_ = melt(t(matrix))
  ggplot(table_ , aes(x=Var2, y=value, colour=Var1)) + geom_line(show.legend = legend, size=1) + geom_point(show.legend = legend) + theme_bw() + labs(x=xlab, y=ylab, colour='', title=title) + theme(legend.text=element_text(size=16,face="bold"), axis.title = element_text(size=16, face='bold'), axis.text.y=element_text(siz=16, face='bold'), axis.text.x=element_text(size=16, face='bold'), plot.title=element_text(size=16, face='bold'), legend.title=element_text(size=16, face='bold'))
}