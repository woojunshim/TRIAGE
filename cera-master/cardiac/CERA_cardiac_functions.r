# Jun Xu 16/04/2018
# functions used in CERA_main
# 1: matrix.to.results - generate result table from matrix
# 2: prepare.results - used in ad-hoc app, prepare output tables for shiny tabs
# 3: cardiac.results - used by app_cardiac*.R, calling matrix.to.results
# 4: filtered.results - filter results based on input ranges
# 5: update.sliders - update the sliders min, max, value range based on the matrix dynamically
# 6: prepare.files - prepare the download files based on the input and output
# 7: summary.page - output function for the summary page
# 8: tsne2d.clustcolor - prepare for tSNE plot
# 9: elbow.finder - find the elbow point

# function 1: generate result table from matrix
matrix.to.results <- function(m, species, gene, sym, reference, session, inputday) {
  if (species == 1) {  # human
    # Repression score of MAD proportion + TSS information
    w <- read.table('repressive_hg19_cont_prop.txt')  #  hybrid model (continuous variable * proportion)
    h27bp <- read.table('H3K27me3_mean_table.txt')  # means of h3k27me3 width
    cv <- read.table('cv_genes.txt')  # read coefficient of variation of the genes
  } else if (species == 2){  # mouse
    w <- read.table('genes_broad_h3k27me3_stats_mm10_.txt')
    h27bp <- read.table('H3K27me3_mean_mm10.txt')  # means of h3k27me3 width
  }
  r <- as.data.frame(matrix(0,nrow=nrow(m),ncol=14))  # build initial result table
  colnames(r) <- c('Gene','TF','k27bp','CV','Repression','Proportion',
                   'mean_exp_w_0','mean_exp_wo_0','discordance_w_0',
                   'discordance_wo_0','SD','CC','DE','CO')
  rownames(r) <- rownames(m)

  h27bps <- h27bp[rownames(r),]
  h27bps[is.na(h27bps)] <- 0  # clear NA
  r$Gene <- rownames(r)
  r$TF <- h27bps[,2]  # TF flag
  r$k27bp <- as.integer(h27bps[,1])  # h3k27me3 width in bp
  cv <- cv[rownames(r),]
  cv[is.na(cv)] <- 0  # clear NA
  r$CV <- round(cv$value,3)
  rownames(w) <- w[,1]  # set the row names for Repression scores
  Repression <- w[rownames(r),2]  # RT score for the genes
  Repression[is.na(Repression)] <- 0  # clear NA
  r$Repression <- round(Repression,3)  # calculated chromatin Repression score with pseudo count
  
  if (!is.null(gene)) {
    if (sym == 1) {
      # subset with all selected genes expressed
      m <- m[,which(is.finite(colSums(log(m)[gene,,drop=FALSE]))),drop=FALSE]
    } else if (sym == 2) {
      # subset with any selected genes expressed
      m <- m[,which(colSums(m[gene,,drop=FALSE])!=0),drop=FALSE]   
    }
  }
  session$userData$g <- c('', sort(rownames(m)))   # build the list of available genes from the data for selection
  
  if (ncol(m)>0) {   # in case filtering gene is expressed in no cells, do not calculate metrics
    r$Proportion <- round(1-apply(m==0,1,sum)/apply(m>=0,1,sum),3)  # percent of expression
    r$mean_exp_w_0 <- round(apply(m,1,mean),3)  # mean of original expression with zero expressed cells
    r$mean_exp_wo_0 <- round(apply(m,1,function(x){mean(x[x>0])}),3)  # mean of original expression w/ot zeros
    r$SD <- round(apply(m,1,sd),3)  # standard deviation of original expression with zero expressed cells
    if (!is.null(reference)) {
      if (reference != '') {
        r$CC <- t(cor(m[reference,],t(m)))
        
        n <- m[,which(m[reference,]>0)]    # select the cells expressing the reference gene
        r$CO <- 1-apply(n==0,1,sum)/apply(n>=0,1,sum)

        # # re-generate all cluster matrix for differnetial analysis
        # exp <- get(paste0('m',inputday,1))
        # for (j in 2:length(table(get(paste0('dat3d_cluster',inputday,'$data.cluster'))))) {
        #   exp <- cbind(exp, get(paste0('m',inputday,j)))
        # }
        
        # Abbi Helfer (adjusted)
        # count positive and negative cells 
        gene.pos <- m[, which(m[reference,] > 0)]
        gene.neg <- m[, which(m[reference,] == 0)]
        # count positive and negative cells
        pos <- ncol(gene.pos)
        neg <- ncol(gene.neg)
        tot <- ncol(m)
        # create gene +/- condition
        cond <- rep('neg', ncol(m))
        cond[which(m[reference, ] >  0)] <- 'pos'  
        cond <- as.factor(cond)
        # pseudocount of 1 is added
        # DESeq
        cds <- newCountDataSet(round(m + 1), cond)
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds, fitType = 'local')
        res <- nbinomTest(cds, "neg", "pos" )
        rm(cds)
        bonf <- 0.05/nrow(res)
        res.new <- res[order(rownames(r)), ]
        res.new$adjFoldChange <- (res.new$baseMeanB - 1)/(res.new$baseMeanA - 1)
        r$DE <- log2(res.new$adjFoldChange)
        # r$DE[which(res.new$pval > bonf)] <- NA
      }
    }
    
    # calculate the discordance score for each gene/cell
    # pseudo 1 added to avoid 0 expressions and negative log-transformed values
    log.weighted <- rep(Repression, ncol(m)) * log(m + 1)

    # mean of the log-weighted NOT counting 0s
    r$discordance_wo_0 <- round(apply(log.weighted, 1, 
                                        function(x){mean(x[x>0])}),3)

    # mean of the log-weighted counting 0s
    r$discordance_w_0 <- round(apply(log.weighted, 1, mean),3)

  }
  r
}


# function 2: used in ad-hoc app, prepare output tables for shiny tabs (calling matrix.to.results)
prepare.results <- function(session, input, n) {
  req(input[[paste0('file',n)]])
  session$userData$g0 <- input$filter
  if (is.null(session$userData$g0)) {session$userData$g0 <- ''}  # transform NULL input
  
  # reload and recalculate if the file name or file size of gene selection has been changed
  if ((get(paste0('fn',n),session$userData) != 
       input[[paste0('file',n)]]$name) |
      (get(paste0('z',n),session$userData) != 
       input[[paste0('file',n)]]$size) |
      !all(get(paste0('g',n),session$userData) == session$userData$g0) |
      (get(paste0('s',n),session$userData) != input$sym)) {
    # not necessary to reload if gene selection is not changed
    # but considering the memory limit on shinyapps.io, do not keep m  
    m <- readRDS(input[[paste0('file',n)]]$datapath)
    eval(parse(text=paste0('session$userData$r',n, '<- matrix.to.results(m, input$species, 
                           input$filter, input$sym, input$reference, session, 0)')))
    # if the refreshing is NOT triggered by new gene selection
    # summary table might have different gene list from different matrices
    # here we assume the matrices to be compared have similar gene list
    if (all(get(paste0('g',n),session$userData) == session$userData$g0)) {
      updateSelectizeInput(session, 'filter', choices = session$userData$g,
                             selected =input$filter, options = list(maxItems=3))
    # g was updated inside matrix.to.results
      updateSelectizeInput(session, 'reference', choices = session$userData$g,  
                           selected =input$reference, options = list(maxItems=1))
      updatesliders(session, input)
    }
  }
  
  # record the current selection
  eval(parse(text=paste0('session$userData$fn',n,'<- input[[paste0(\'file\',n)]]$name')))
  eval(parse(text=paste0('session$userData$z',n,'<- input[[paste0(\'file\',n)]]$size')))
  eval(parse(text=paste0('session$userData$g',n,'<- session$userData$g0'))) # record day
  eval(parse(text=paste0('session$userData$s',n,'<- input$sym'))) # record symbol
  }


# function 3: used by app_cardiac*.R, calling matrix.to.results
cardiac.results <- function(session, inputday, species, input, m, n) {
  # m: matrix to process
  # n: cluster number
  session$userData$g0 <- input$filter     # record the input to <variable>0
  if (is.null(session$userData$g0)) {session$userData$g0 <- ''}  # transform NULL input
  session$userData$rg0 <- input$reference
  if (is.null(session$userData$rg0)) {session$userData$rg0 <- ''}
  
  # calculate the result if the gene or day selection has been changed
  if ((!all(get(paste0('g',n),session$userData) == session$userData$g0) | # filtering gene changed
       get(paste0('rg',n),session$userData) != input$reference |      # or reference changed
       get(paste0('d',n),session$userData) != inputday |              # or day changed
       get(paste0('s',n),session$userData) != input$sym) &            # or sym changed
       exists(paste0('m',inputday,n)) &                          # and source matrix does exist
    # do not re-calculate if only reference gene is changed for summary tab, as no output will need to change
      !(input$tab == 'Summary' &
         all(get(paste0('g',n),session$userData) == session$userData$g0) & 
         get(paste0('d',n),session$userData) == inputday &
         get(paste0('s',n),session$userData) == input$sym)) {
      
      eval(parse(text=paste0('session$userData$r',n,' <- matrix.to.results(m, species, 
                             input$filter, input$sym, input$reference, session, inputday)')))
      
      if (input$tab == 'Summary') {session$userData$rs0 <- 0}
      if (input$tab == 'Cluster 1') {session$userData$rs1 <- 0}
      if (input$tab == 'Cluster 2') {session$userData$rs2 <- 0}
      if (input$tab == 'Cluster 3') {session$userData$rs3 <- 0}
      if (input$tab == 'Cluster 4') {session$userData$rs4 <- 0}
      if (input$tab == 'Multiple') {session$userData$rsm <- 0}
      
      # summary table might have different gene list from different clusters/matrices
      # here we have assumed the matrices to be compared have similar gene list
      # g: from matrix.to.results
      updateSelectizeInput(session, 'filter', choices = session$userData$g,
                           selected =input$filter, options = list(maxItems=3))
      updateSelectizeInput(session, 'reference', choices = session$userData$g,
                           selected =input$reference, options = list(maxItems=1))
      updatesliders(session, input)
      
      # record the current selection
      eval(parse(text=paste0('session$userData$g',n,'<- session$userData$g0'))) # record gene
      eval(parse(text=paste0('session$userData$rg',n,'<- session$userData$rg0'))) # record reference
      eval(parse(text=paste0('session$userData$d',n,'<- inputday'))) # record day
      eval(parse(text=paste0('session$userData$s',n,'<- input$sym'))) # record symbol
    }
  }


# function 4: filter results based on input ranges
filtered.results <- function(session, r, inputday, input, matrix.name) {
  update <- 0
  if (input$zeroexpression == 1) {
    r <- r[,-c(8,10)]   # including 0 by dropping 'without zero expression' columns
  } else if (input$zeroexpression == 2) {
    r <- r[,-c(7,9)]   # excluding 0 by dropping 'with zero expression' columns
  }
  colnames(r)[7] <- 'mean_exp'
  colnames(r)[8] <- 'discordance'
  
  # set 0 for the NA values in the result table to ease filtering
  if (nrow(r[is.na(r$mean_exp),]) > 0) {
    r[is.na(r$mean_exp),]$mean_exp <- 0
  }
  if (nrow(r[is.na(r$discordance),]) > 0) {
    r[is.na(r$discordance),]$discordance <- 0
  }
  # if resetting, different initialisation for including/excluding null expressions
  if ((input$tab == 'Summary' & session$userData$rs0 == 0)|
      (input$tab == 'Cluster 1' & session$userData$rs1 == 0) |
      (input$tab == 'Cluster 2' & session$userData$rs2 == 0) |
      (input$tab == 'Cluster 3' & session$userData$rs3 == 0) |
      (input$tab == 'Cluster 4' & session$userData$rs4 == 0) |
      (input$tab == 'Multiple' & session$userData$rsm == 0)) {
    k <- kmeans(r$discordance,2)
    r <- r[which(k$cluster == which(k$centers==max(k$centers))),]
    r <- r[order(r$discordance, decreasing=TRUE),]
    
      session$userData$ranges <- c(min(session$userData$ranges[1],min(r$k27bp)),
                                   max(session$userData$ranges[2],max(r$k27bp)),
                                   min(session$userData$ranges[3],min(r$CV)),
                                   max(session$userData$ranges[4],max(r$CV)),
                                   min(session$userData$ranges[5],min(r$Repression)),
                                   max(session$userData$ranges[6],max(r$Repression)),
                                   min(session$userData$ranges[7],min(r$Proportion*100)),
                                   max(session$userData$ranges[8],max(r$Proportion*100)),
                                   min(session$userData$ranges[9],min(r$mean_exp)),
                                   max(session$userData$ranges[10],max(r$mean_exp)),
                                   min(session$userData$ranges[11],min(r$discordance)),
                                   max(session$userData$ranges[12],max(r$discordance)),
                                   min(session$userData$ranges[13],min(r$CC)),
                                   max(session$userData$ranges[14],max(r$CC)),
                                   min(session$userData$ranges[15],min(r$CO)),
                                   max(session$userData$ranges[16],max(r$CO)))

    # mark reset complete and update min/max
    if (input$tab == 'Summary') {
      if ((session$userData$count == length(table(get(paste0('dat3d_cluster',
                                    inputday), envir=.GlobalEnv)$data.cluster)))) {
        session$userData$rs0 <- 1
        session$userData$count <- 1
        update <- 1
      } else {
        # count for all clusters
        session$userData$count <- session$userData$count + 1
      }
    } else {
      if (input$tab == 'Cluster 1') {session$userData$rs1 <- 1}
      if (input$tab == 'Cluster 2') {session$userData$rs2 <- 1}
      if (input$tab == 'Cluster 3') {session$userData$rs3 <- 1}
      if (input$tab == 'Cluster 4') {session$userData$rs4 <- 1}
      if (input$tab == 'Multiple') {session$userData$rsm <- 1}
      update <- 1
    }

    # update the slider bars
    if (update == 1) {
      updateSliderInput(session, 'range_width', 
                      value = c(session$userData$ranges[1],session$userData$ranges[2]))
      updateSliderInput(session, 'range_CV', 
                      value = c(session$userData$ranges[3],session$userData$ranges[4]))
      updateSliderInput(session, 'range_Repression', 
                      value = c(session$userData$ranges[5],session$userData$ranges[6]))
      updateSliderInput(session, 'range_percent', 
                      value = c(session$userData$ranges[7],session$userData$ranges[8]))
      updateSliderInput(session, 'range_mean', 
                      value = c(session$userData$ranges[9],session$userData$ranges[10]))
      updateSliderInput(session, 'range_discordance', 
                      value = c(session$userData$ranges[11],session$userData$ranges[12]))
      updateSliderInput(session, 'range_cc', 
                      value = c(session$userData$ranges[13],session$userData$ranges[14]))
      updateSliderInput(session, 'range_co', 
                      value = c(session$userData$ranges[15],session$userData$ranges[16]))
      update <- 0
      session$userData$ranges <- c(1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,
                                   1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10)
    }
      
  } else {
      r <- r[!is.na(r$k27bp),]
      r <- r[!is.na(r$CV),]
      r <- r[!is.na(r$Repression),]
      r <- r[r$k27bp>=input$range_width[1] & r$k27bp<=input$range_width[2],]
      r <- r[r$CV>=input$range_CV[1] & r$CV<=input$range_CV[2],]
      r <- r[r$Repression>=input$range_Repression[1] & r$Repression<=input$range_Repression[2],]
      r <- r[r$Proportion*100>=input$range_percent[1] & r$Proportion*100<=input$range_percent[2],]
      r <- r[r$mean_exp>=input$range_mean[1] & r$mean_exp<=input$range_mean[2],]
      r <- r[r$discordance>=input$range_discordance[1] & r$discordance<=input$range_discordance[2],]
      r <- r[r$CC>=input$range_cc[1] & r$CC<=input$range_cc[2],]
      r <- r[r$CO>=input$range_co[1] & r$CO<=input$range_co[2],]
  }
  
  r$Proportion <- percent(r$Proportion)   # percentage format for the expressed percentage
  r$CO <- percent(r$CO)
  
  # sort the filtered result and generate summary (rr) and detail (rd) results
  if (input$sort == 1) {
    rd <- r[order(r$Gene),]
    rr <- as.matrix(r[order(r$Gene),1])
  } else if (input$sort == 2) {
    rd <- r[order(r$TF,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$TF,decreasing=TRUE),1])
  } else if (input$sort == 3) {
    rd <- r[order(r$k27bp,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$k27bp,decreasing=TRUE),1])
  } else if (input$sort == 4) {
    rd <- r[order(r$CV,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$CV,decreasing=TRUE),1])
  } else if (input$sort == 5) {
    rd <- r[order(r$Repression,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$Laltency,decreasing=TRUE),1])
  } else if (input$sort == 6) {
    rd <- r[order(r$Proportion,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$Proportion,decreasing=TRUE),1])
  } else if (input$sort == 7) {
    rd <- r[order(r$mean_exp,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$mean_exp,decreasing=TRUE),1])
  } else if (input$sort == 8) {
    rd <- r[order(r$discordance,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$discordance,decreasing=TRUE),1])
  } else if (input$sort == 9) {
    rd <- r[order(r$CC,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$CC,decreasing=TRUE),1])
  } else if (input$sort == 10) {
    rd <- r[order(r$DE,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$DE,decreasing=TRUE),1])
  } else if (input$sort == 11) {
    rd <- r[order(r$CO,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$CO,decreasing=TRUE),1])
  }
  colnames(rr) <- matrix.name  # add matrix id to the result column
  common <- unique(Reduce(union, list(session$userData$f1[,1], session$userData$f2[,1], 
                                   session$userData$f3[,1],session$userData$f4[,1])))
  updateSelectizeInput(session, 'common', choices = common,
                       selected =input$common, options = list(maxItems=1))
  filtered.list <- list('sum' = rr, 'detail' = rd)   # returned values
}


# function 5: update the sliders min, max, value range based on the matrix dynamically
updatesliders <- function(session, input){
  
  # compare the current min/max with the pre-filtered results, for all clusters sequentially
  if (input$zeroexpression == 1) {   # including null expressions
    if (exists('r1',session$userData)) {
      min_mean <- min(session$userData$r1$mean_exp_w_0[
        is.finite(session$userData$r1$mean_exp_w_0)])
      max_mean <- max(session$userData$r1$mean_exp_w_0[
        is.finite(session$userData$r1$mean_exp_w_0)])
      min_discordance <- min(session$userData$r1$discordance_w_0[
        is.finite(session$userData$r1$discordance_w_0)])
      max_discordance <- max(session$userData$r1$discordance_w_0[
        is.finite(session$userData$r1$discordance_w_0)])
    }
    if (exists('r2',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r2$mean_exp_w_0[
        is.finite(session$userData$r2$mean_exp_w_0)]))
      max_mean <- max(max_mean, max(session$userData$r2$mean_exp_w_0[
        is.finite(session$userData$r2$mean_exp_w_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r2$discordance_w_0[
        is.finite(session$userData$r2$discordance_w_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r2$discordance_w_0[
        is.finite(session$userData$r2$discordance_w_0)]))
    }
    if (exists('r3',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r3$mean_exp_w_0[
        is.finite(session$userData$r3$mean_exp_w_0)]))
      max_mean <- max(max_mean, max(session$userData$r3$mean_exp_w_0[
        is.finite(session$userData$r3$mean_exp_w_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r3$discordance_w_0[
        is.finite(session$userData$r3$discordance_w_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r3$discordance_w_0[
        is.finite(session$userData$r3$discordance_w_0)]))
    }
    if (exists('r4',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r4$mean_exp_w_0[
        is.finite(session$userData$r4$mean_exp_w_0)]))
      max_mean <- max(max_mean, max(session$userData$r4$mean_exp_w_0[
        is.finite(session$userData$r4$mean_exp_w_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r4$discordance_w_0[
        is.finite(session$userData$r4$discordance_w_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r4$discordance_w_0[
        is.finite(session$userData$r4$discordance_w_0)]))
    }
  } else if (input$zeroexpression == 2) {   # excluding null expressions
    if (exists('r1',session$userData)) {
      min_mean <- min(session$userData$r1$mean_exp_wo_0[
        is.finite(session$userData$r1$mean_exp_wo_0)])
      max_mean <- max(session$userData$r1$mean_exp_wo_0[
        is.finite(session$userData$r1$mean_exp_wo_0)])
      min_discordance <- min(session$userData$r1$discordance_wo_0[
        is.finite(session$userData$r1$discordance_wo_0)])
      max_discordance <- max(session$userData$r1$discordance_wo_0[
        is.finite(session$userData$r1$discordance_wo_0)])
    }
    if (exists('r2',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r2$mean_exp_wo_0[
        is.finite(session$userData$r2$mean_exp_wo_0)]))
      max_mean <- max(max_mean, max(session$userData$r2$mean_exp_wo_0[
        is.finite(session$userData$r2$mean_exp_wo_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r2$discordance_wo_0[
        is.finite(session$userData$r2$discordance_wo_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r2$discordance_wo_0[
        is.finite(session$userData$r2$discordance_wo_0)]))
    }
    if (exists('r3',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r3$mean_exp_wo_0[
        is.finite(session$userData$r3$mean_exp_wo_0)]))
      max_mean <- max(max_mean, max(session$userData$r3$mean_exp_wo_0[
        is.finite(session$userData$r3$mean_exp_wo_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r3$discordance_wo_0[
        is.finite(session$userData$r3$discordance_wo_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r3$discordance_wo_0[
        is.finite(session$userData$r3$discordance_wo_0)]))
    }
    if (exists('r4',session$userData)) {
      min_mean <- min(min_mean, min(session$userData$r4$mean_exp_wo_0[
        is.finite(session$userData$r4$mean_exp_wo_0)]))
      max_mean <- max(max_mean, max(session$userData$r4$mean_exp_wo_0[
        is.finite(session$userData$r4$mean_exp_wo_0)]))
      min_discordance <- min(min_discordance, min(session$userData$r4$discordance_wo_0[
        is.finite(session$userData$r4$discordance_wo_0)]))
      max_discordance <- max(max_discordance, max(session$userData$r4$discordance_wo_0[
        is.finite(session$userData$r4$discordance_wo_0)]))
    } 
  }
  
  # update the input selections
  updateSliderInput(session, 'range_mean', value = c(min_mean, max_mean), min = min_mean, max = max_mean)
  updateSliderInput(session, 'range_discordance', value = c(min_discordance, max_discordance), 
                    min = min_discordance, max = max_discordance)
}

# function 6: prepare the download files based on the input and output
prepare.files <- function(session, input) {
  write(paste('Day:',as.character(input$day),sep='\t'),'user.inputs.txt')
  write(paste('Multiple clusters:',input$clusters,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('reference:',input$reference,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('Expressed genes:',input$filter,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('Symbol of multiple genes (1:AND, 2: OR):',
              input$sym,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('Cell type specific (1:Yes, 2: No):',
              input$specific,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('sorted by (1: Gene, 2: TF, 3: K27bp, 4: coefficient of variant, 5: Repression,
              6: percent of expression, 7: mean_exp, 8: discordance, 9: correlation):',
              input$sort,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('h3k27me3 width:',input$range_width[1],input$range_width[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('Coefficient of variation:',input$range_CV[1],input$range_CV[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('Repression:',input$range_Repression[1],input$range_Repression[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('percent of expression:',input$range_percent[1],input$range_percent[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('mean expression:',input$range_mean[1],input$range_mean[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('discordance:',input$range_discordance[1],input$range_discordance[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('Correlation:',input$range_cc[1],input$range_cc[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  files <- 'user.inputs.txt'
  if (exists('o1',session$userData)) {
    write.csv(session$userData$o1, 'cluster1.csv')
    files <- c('cluster1.csv', files)
  }
  if (exists('o2',session$userData)) {
    write.csv(session$userData$o2, 'cluster2.csv')
    files <- c('cluster2.csv', files)
  }
  if (exists('o3',session$userData)) {
    write.csv(session$userData$o3, 'cluster3.csv')
    files <- c('cluster3.csv', files)
  }
  if (exists('o4',session$userData)) {
    write.csv(session$userData$o4, 'cluster4.csv')
    files <- c('cluster4.csv', files)
  }
  if (exists('om',session$userData)) {
    write.csv(session$userData$om, 'mixed.clusters.csv')
    files <- c('mixed.clusters.csv', files)
  }
  if (exists('op',session$userData)) {
    write.csv(session$userData$op, 'selected_gene_ratio.csv')
    files <- c('selected_gene_ratio.csv', files)
  }
  if (exists('oq',session$userData)) {
    write.csv(session$userData$oq, 'selected_gene_expression.csv')
    files <- c('selected_gene_expression.csv', files)
  }
  files
}

# function 7: output function for the summary page
summary.page <- function(session, input, inputday, species, cluster) {
  
  # calculate the results table for all tabs
  cardiac.results(session, inputday, species, input, get(paste0('m',inputday,1)), 1)
  cardiac.results(session, inputday, species, input, get(paste0('m',inputday,2)), 2)
  cardiac.results(session, inputday, species, input, get(paste0('m',inputday,3)), 3)
  cardiac.results(session, inputday, species, input, get(paste0('m',inputday,4)), 4)
  
  # generate filtered results based on the input if the result table and the matrix both exist
  if (cluster == 1 & exists('r1',session$userData) & exists(paste0('m',inputday,1))) {  
    session$userData$f1 <- filtered.results(session, session$userData$r1, inputday, input, 'cluster 1')$sum
    rownames(session$userData$f1) <- session$userData$f1[,1]
  }
  if (cluster == 2 & exists('r2',session$userData) & exists(paste0('m',inputday,2))) {
    session$userData$f2 <- filtered.results(session, session$userData$r2, inputday, input, 'cluster 2')$sum
    rownames(session$userData$f2) <- session$userData$f2[,1]
  }    
  if (cluster == 3 & exists('r3',session$userData) & exists(paste0('m',inputday,3))) {
    session$userData$f3 <- filtered.results(session, session$userData$r3, inputday, input, 'cluster 3')$sum
    rownames(session$userData$f3) <- session$userData$f3[,1]
  }
  if (cluster == 4 & exists('r4',session$userData) & exists(paste0('m',inputday,4))) {
    session$userData$f4 <- filtered.results(session, session$userData$r4, inputday, input, 'cluster 4')$sum
    rownames(session$userData$f4) <- session$userData$f4[,1]
  }
  
  # generate dummy filter table for intersection, won't be output if matrix not existing    
  if (!exists(paste0('m',inputday,4))) {session$userData$f4 <- session$userData$f1}
  if (!exists(paste0('m',inputday,3))) {session$userData$f3 <- session$userData$f1}
  if (!exists(paste0('m',inputday,2))) {session$userData$f2 <- session$userData$f1}
  
  # select the unique values if "cluster specific" sign is on
  if (input$specific==1) {
    if (exists(paste0('r',cluster),session$userData) & exists(paste0('m',inputday,cluster))) {
      eval(parse(text=paste0('session$userData$f',cluster,'[!session$userData$f',cluster,
              '[,1]%in%Reduce(intersect, list(session$userData$f1[,1], session$userData$f2[,1], 
              session$userData$f3[,1],session$userData$f4[,1])),,drop=FALSE]')))
    }
  }
  
  # if "cluster specific" sign is off
  else if (exists(paste0('r',cluster),session$userData) & exists(paste0('m',inputday,cluster))) {
    # check whether result table exists
    eval(parse(text=paste0('session$userData$f',cluster)))
    # filter the results based on the input
  }
}

# function 8: prepare for tSNE plot
# Abbi Helfer 21-06-2017
# Tailored by Jon Xu 04-04-2018
# Can generate a tSNE plot each colored by cluster

tsne2d.clustcolor <- function(gene, dat3d_all, day1=day){
  max_cluster <- max(dat3d_all$cluster) # Total number of clusters
  color_codes <- c("#D3D3D3", "#0E3D59", "#88A61B", "#F29F05", "#D92525")
  dat3d_all$color <- rep('Not Expressed', length(dat3d_all[,4]))
  legend_names <- rep('Not Expressed', max_cluster + 1)
  # Loop through expression data and find cells that express the gene and assign them to a group based on cluster
  for (i in 1:max_cluster) {
    clust <- which(dat3d_all$cluster == i & dat3d_all$gene_exp > 0)
    dat3d_all$color[+clust] <- paste(gene, paste(' in Cluster ', i, sep=""), sep="")
    legend_names[i+1] <-  paste(gene, paste(' in Cluster ', i, sep=""), sep="")
    # If there are no cells expressing gene in a cluster, remove that color from the color list
    if (length(dat3d_all$cluster[dat3d_all$cluster == i & dat3d_all$gene_exp > 0]) == 0) {
      remove <- which(dat3d_all$cluster == i)
      color_codes <- color_codes[-(i+1)]
    }
  }
  colnames(dat3d_all) <- c('tSNE1', 'tSNE2', 'tSNE3', 'cluster', 'expression', 'color')
  dat3d_all$color <- factor(dat3d_all$color, levels=legend_names)
  values <- list(dat3d_all, max_cluster, color_codes)
}

# function 9: find the elbow point

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)
  
  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}