# Jun Xu 16/04/2018
# functions used in CERA_adhoc
# 1: matrix.to.results - generate result table from matrix
# 2: prepare.results - used by app_cardiac*.R, calling matrix.to.results
# 3: filtered.results - filter results based on input ranges
# 4: update.sliders - update the sliders min, max, value range based on the matrix dynamically
# 5: prepare.files - prepare the download files based on the input and output
# 6: summary.page - output function for the summary page

# function 1: generate result table from matrix
matrix.to.results <- function(m, species, gene, sym, reference, session, inputday) {
  w <- read.table('repressive_hg19_cont_prop.txt')  # read human repressive tendency scores
  colnames(w) <- c('human', 'rts')
  h27bp <- read.table('H3K27me3_mean_table.txt')  # means of h3k27me3 width
  if (species == 2) {t <- unique(read.table('Human_C.intestinalis.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 3) {t <- unique(read.table('Human_Chicken.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 4) {t <- unique(read.table('Human_Cow.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 5) {t <- unique(read.table('Human_Fruitfly.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 6) {t <- unique(read.table('Human_Mouse.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 7) {t <- unique(read.table('Human_Pig.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 8) {t <- unique(read.table('Human_Rat.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species == 9) {t <- unique(read.table('Human_Zebrafish.txt', TRUE, sep=',')[, c(2,10)])}   # C.intestinalis
  if (species != 1) {
    t <- t[order(t[,2]),]
    t <- t[!duplicated(t[,2]),]
    rownames(t) <- t[,2]
    colnames(t) <- c('human', 'animal')
    w <- merge(t,w)[-1]
    h27bp[,3] <- rownames(h27bp)
    colnames(h27bp)[3] <- 'human'
    h27bp <- merge(t,h27bp)
    rownames(h27bp) <- h27bp[,2]
    h27bp <- h27bp[,c(3,4)]
  }
  rownames(w) <- w[,1]  # set the row names for rts scores
  
  r <- as.data.frame(matrix(0,nrow=nrow(m),ncol=12))  # build initial result table
  colnames(r) <- c('Gene','TF','k27bp','RTS','Proportion', 'mean_exp_w_0','mean_exp_wo_0',
                   'discordance_w_0', 'discordance_wo_0','CC','DE','CO')
  rownames(r) <- rownames(m)
  h27bps <- h27bp[rownames(r),]
  h27bps[is.na(h27bps)] <- 0  # clear NA
  r$Gene <- rownames(r)
  r$TF <- h27bps[,2]  # TF flag
  r$k27bp <- as.integer(h27bps[,1])  # h3k27me3 width in bp
  rts <- w[rownames(r),2]  # rts for the genes in result table
  rts[is.na(rts)] <- 0  # clear NA
  r$RTS <- round(rts,3)  # calculated chromatin rts score with pseudo count
  
  if (!is.null(gene)) {
    if (sym == 1) {
      # subset with ALL selected genes expressed
      m <- m[,which(is.finite(colSums(log(m)[gene,,drop=FALSE]))),drop=FALSE]
    } else if (sym == 2) {
      # subset with ANY selected genes expressed
      m <- m[,which(colSums(m[gene,,drop=FALSE])!=0),drop=FALSE]   
    }
  }
  session$userData$g <- c('', sort(rownames(m)))   # build the list of available genes from the data for selection
  
  if (ncol(m)>0) {   # in case filtering gene is expressed in no cells, do not calculate metrics
    r$Proportion <- round(1-apply(m==0,1,sum)/apply(m>=0,1,sum),3)  # percent of expression
    r$mean_exp_w_0 <- round(apply(m,1,mean),3)  # mean of original expression with zero expressed cells
    r$mean_exp_wo_0 <- round(apply(m,1,function(x){mean(x[x>0])}),3)  # mean of original expression
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
    log.weighted <- rep(rts, ncol(m)) * log(m + 1)
    
    # mean of the log-weighted NOT counting 0s
    r$discordance_wo_0 <- round(apply(log.weighted, 1, 
                                      function(x){mean(x[x>0])}),3)
    
    # mean of the log-weighted counting 0s
    r$discordance_w_0 <- round(apply(log.weighted, 1, mean),3)
    
  }
  rm(m)
  gc()
  r
}


# function 2: used by CERA_adhoc, calling matrix.to.results to get result tables
prepare.results <- function(session, input, n) {
  # n: matrix number
  session$userData$sp0 <- input$species     # record the input to <var0>
  if (is.null(session$userData$sp0)) {session$userData$sp0 <- ''}  # transform NULL input
  session$userData$g0 <- input$filter
  if (is.null(session$userData$g0)) {session$userData$g0 <- ''}
  session$userData$rg0 <- input$reference
  if (is.null(session$userData$rg0)) {session$userData$rg0 <- ''}
  session$userData$fn0 <- input[[paste0('file',n)]]$name
  if (is.null(session$userData$fn0)) {session$userData$fn0 <- ''}
  session$userData$z0 <- input[[paste0('file',n)]]$size
  if (is.null(session$userData$z0)) {session$userData$z0 <- 0}
  # re-calculate result if any data-changing selection has been updated - <var0> refers to the last updated selections
  if (get(paste0('sp',n),session$userData) != session$userData$sp0 | # species changed
      all(get(paste0('g',n),session$userData) != session$userData$g0) | # filtering gene changed, "all" for multiple filtering genes
      get(paste0('rg',n),session$userData) != input$reference |      # or reference changed
      get(paste0('fn',n),session$userData) != session$userData$fn0 |   # or file name changed
      get(paste0('z',n),session$userData) != session$userData$z0  |   # or file size changed
      get(paste0('s',n),session$userData) != input$sym &            # or sym changed
      !is.null(input[[paste0('file',n)]]) &                          # and file selected
      # do not re-calculate if only reference gene is changed for summary tab, as no output will need to change
      !(input$tab == 'Summary' &
        all(get(paste0('g',n),session$userData) == session$userData$g0) & # "all" for multiple filtering genes
        get(paste0('s',n),session$userData) == input$sym)) {
    # end of if conditions
    # start of if execution    
    req(input[[paste0('file',n)]])
    m <- read.csv(input[[paste0('file',n)]]$datapath, TRUE, ',', row.names=1)
    assign(paste0('m0',n), m, envir = session$userData)
    eval(parse(text=paste0('session$userData$r',n,' <- matrix.to.results(m, input$species, input$filter, input$sym, input$reference, session, 0)')))
    
    if (input$tab == 'Summary') {session$userData$rs0 <- 0}
    if (input$tab == 'Matrix 1') {session$userData$rs1 <- 0}
    if (input$tab == 'Matrix 2') {session$userData$rs2 <- 0}
    if (input$tab == 'Matrix 3') {session$userData$rs3 <- 0}
    if (input$tab == 'Matrix 4') {session$userData$rs4 <- 0}
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
    eval(parse(text=paste0('session$userData$s',n,'<- input$sym'))) # record symbol
    eval(parse(text=paste0('session$userData$sp',n,'<- input$species'))) # record species
    eval(parse(text=paste0('session$userData$fn',n,'<- input[[paste0(\'file\',n)]]$name'))) # record file names
    eval(parse(text=paste0('session$userData$z',n,'<- input[[paste0(\'file\',n)]]$size'))) # record file sizes
  }
  if (is.null(input[[paste0('file',n)]])) {assign(paste0('m0',n), NULL, envir = session$userData)}
}


# function 3: filter results based on input ranges
filtered.results <- function(session, r, inputday, input, matrix.name) {
  update <- 0
  if (input$zeroexpression == 1) {
    r <- r[,-c(7,9)]   # including 0 by dropping 'without zero expression' columns
  } else if (input$zeroexpression == 2) {
    r <- r[,-c(6,8)]   # excluding 0 by dropping 'with zero expression' columns
  }
  colnames(r)[6] <- 'mean_exp'
  colnames(r)[7] <- 'discordance'
  
  # set 0 for the NA values in the result table to ease filtering
  if (nrow(r[is.na(r$mean_exp),]) > 0) {
    r[is.na(r$mean_exp),]$mean_exp <- 0
  }
  if (nrow(r[is.na(r$discordance),]) > 0) {
    r[is.na(r$discordance),]$discordance <- 0
  }
  # if resetting, different initialisation for including/excluding null expressions
  if ((input$tab == 'Summary' & session$userData$rs0 == 0)|
      (input$tab == 'Matrix 1' & session$userData$rs1 == 0) |
      (input$tab == 'Matrix 2' & session$userData$rs2 == 0) |
      (input$tab == 'Matrix 3' & session$userData$rs3 == 0) |
      (input$tab == 'Matrix 4' & session$userData$rs4 == 0) |
      (input$tab == 'Multiple' & session$userData$rsm == 0)) {
    k <- kmeans(r$discordance,2)
    r <- r[which(k$cluster == which(k$centers==max(k$centers))),]
    r <- r[order(r$discordance, decreasing=TRUE),]
    
    session$userData$ranges <- c(min(session$userData$ranges[1],min(r$k27bp)),
                                 max(session$userData$ranges[2],max(r$k27bp)),
                                 min(session$userData$ranges[3],min(r$RTS)),
                                 max(session$userData$ranges[4],max(r$RTS)),
                                 min(session$userData$ranges[5],min(r$Proportion*100)),
                                 max(session$userData$ranges[6],max(r$Proportion*100)),
                                 min(session$userData$ranges[7],min(r$mean_exp)),
                                 max(session$userData$ranges[8],max(r$mean_exp)),
                                 min(session$userData$ranges[9],min(r$discordance)),
                                 max(session$userData$ranges[10],max(r$discordance)),
                                 min(session$userData$ranges[11],min(r$CC)),
                                 max(session$userData$ranges[12],max(r$CC)),
                                 min(session$userData$ranges[13],min(r$CO)),
                                 max(session$userData$ranges[14],max(r$CO)))
    
    # mark reset complete and update min/max
    if (input$tab == 'Summary') {
      n = substr_right(matrix.name,1)
      if (is.null(input[[paste0('file',as.integer(n)+1)]])) {
        session$userData$rs0 <- 1
        update <- 1
      }
    } else {
      if (input$tab == 'Matrix 1') {session$userData$rs1 <- 1}
      if (input$tab == 'Matrix 2') {session$userData$rs2 <- 1}
      if (input$tab == 'Matrix 3') {session$userData$rs3 <- 1}
      if (input$tab == 'Matrix 4') {session$userData$rs4 <- 1}
      if (input$tab == 'Multiple') {session$userData$rsm <- 1}
      update <- 1
    }
    
    # update the slider bars
    if (update == 1) {
      updateSliderInput(session, 'range_width',
                        value = c(session$userData$ranges[1],session$userData$ranges[2]))
      updateSliderInput(session, 'range_rts',
                        value = c(session$userData$ranges[3],session$userData$ranges[4]))
      updateSliderInput(session, 'range_percent',
                        value = c(session$userData$ranges[5],session$userData$ranges[6]))
      updateSliderInput(session, 'range_mean', 
                        value = c(session$userData$ranges[7],session$userData$ranges[8]))
      updateSliderInput(session, 'range_discordance', 
                        value = c(session$userData$ranges[9],session$userData$ranges[10]))
      updateSliderInput(session, 'range_cc', 
                        value = c(session$userData$ranges[11],session$userData$ranges[12]))
      updateSliderInput(session, 'range_co', 
                        value = c(session$userData$ranges[13],session$userData$ranges[14]))
      # reset for next round
      update <- 0
      session$userData$ranges <- c(1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10)
    }
    
  } else {
    r <- r[!is.na(r$k27bp),]
    r <- r[!is.na(r$RTS),]
    r <- r[r$k27bp>=input$range_width[1] & r$k27bp<=input$range_width[2],]
    r <- r[r$RTS>=input$range_rts[1] & r$RTS<=input$range_rts[2],]
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
    rd <- r[order(r$RTS,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$RTS,decreasing=TRUE),1])
  } else if (input$sort == 5) {
    rd <- r[order(r$Proportion,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$Proportion,decreasing=TRUE),1])
  } else if (input$sort == 6) {
    rd <- r[order(r$mean_exp,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$mean_exp,decreasing=TRUE),1])
  } else if (input$sort == 7) {
    rd <- r[order(r$discordance,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$discordance,decreasing=TRUE),1])
  } else if (input$sort == 8) {
    rd <- r[order(r$CC,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$CC,decreasing=TRUE),1])
  } else if (input$sort == 9) {
    rd <- r[order(r$DE,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$DE,decreasing=TRUE),1])
  } else if (input$sort == 10) {
    rd <- r[order(r$CO,decreasing=TRUE),]
    rr <- as.matrix(r[order(r$CO,decreasing=TRUE),1])
  }
  colnames(rr) <- matrix.name  # add matrix id to the result column
  common <- unique(Reduce(union, list(session$userData$f1[,1], session$userData$f2[,1], 
                                      session$userData$f3[,1], session$userData$f4[,1])))
  updateSelectizeInput(session, 'common', choices = common,
                       selected =input$common, options = list(maxItems=1))
  filtered.list <- list('sum' = rr, 'detail' = rd)   # returned values
}


# function 4: update the sliders min, max, value range based on the matrix dynamically
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

# function 5: prepare the download files based on the input and output
prepare.files <- function(session, input) {
  tmpdir <- tempdir()
  setwd(tempdir())
  write('User inputs:\r\n','user.inputs.txt')
  write(paste('\r\nSelected species:',input$species,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nMultiple clusters:',input$clusters,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nreference:',input$reference,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nExpressed genes:',input$filter,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nSymbol of multiple genes (1:AND, 2: OR):',
              input$sym,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nCell type specific (1:Yes, 2: No):',
              input$specific,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nsorted by (1: Gene, 2: TF, 3: K27bp, 4: coefficient of variant, 5: rts,
              6: percent of expression, 7: mean_exp, 8: discordance, 9: correlation):',
              input$sort,sep='\t'),'user.inputs.txt',append=TRUE)
  write(paste('\r\nh3k27me3 width:',input$range_width[1],input$range_width[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('\r\nrts:',input$range_rts[1],input$range_rts[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('\r\npercent of expression:',input$range_percent[1],input$range_percent[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('\r\nmean expression:',input$range_mean[1],input$range_mean[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('\r\ndiscordance:',input$range_discordance[1],input$range_discordance[2],sep='\t'),
        'user.inputs.txt',append=TRUE)
  write(paste('\r\nCorrelation:',input$range_cc[1],input$range_cc[2],sep='\t'),
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

# function 6: output function for the summary page
summary.page <- function(session, input, inputday, cluster) {
  
  # calculate the results table for all tabs
  prepare.results(session, input, cluster)
  
  # generate filtered results based on the input if result table and matrix both exist
  if (cluster == 1 & exists('r1',session$userData) & exists('m01', session$userData)) {  
    session$userData$f1 <- filtered.results(session, session$userData$r1, inputday, input, 'Matrix 1')$sum
    rownames(session$userData$f1) <- session$userData$f1[,1]
  }
  if (cluster == 2 & exists('r2',session$userData) & exists('m02', session$userData)) {
    session$userData$f2 <- filtered.results(session, session$userData$r2, inputday, input, 'Matrix 2')$sum
    rownames(session$userData$f2) <- session$userData$f2[,1]
  }    
  if (cluster == 3 & exists('r3',session$userData) & exists('m03', session$userData)) {
    session$userData$f3 <- filtered.results(session, session$userData$r3, inputday, input, 'Matrix 3')$sum
    rownames(session$userData$f3) <- session$userData$f3[,1]
  }
  if (cluster == 4 & exists('r4',session$userData) & exists('m04', session$userData)) {
    session$userData$f4 <- filtered.results(session, session$userData$r4, inputday, input, 'Matrix 4')$sum
    rownames(session$userData$f4) <- session$userData$f4[,1]
  }
  
  # generate dummy filter table for intersection, won't be output if matrix not existing    
  if (is.null(session$userData$m04)) {session$userData$f4 <- session$userData$f1}
  if (is.null(session$userData$m03)) {session$userData$f3 <- session$userData$f1}
  if (is.null(session$userData$m02)) {session$userData$f2 <- session$userData$f1}
  
  # select the unique values if "cluster specific" sign is on
  if (input$specific==1) {
    if (exists(paste0('r',cluster),session$userData) & !is.null(get(paste0('m',inputday,cluster), session$userData))) {
      eval(parse(text=paste0('session$userData$f',cluster,'[!session$userData$f',cluster,
                             '[,1]%in%Reduce(intersect, list(session$userData$f1[,1], session$userData$f2[,1], 
              session$userData$f3[,1],session$userData$f4[,1])),,drop=FALSE]')))
    }
  } else if (exists(paste0('r',cluster),session$userData) & exists(paste0('m',inputday,cluster), session$userData)) {
    # if "cluster specific" sign is off
    # check whether result table exists
    eval(parse(text=paste0('session$userData$f',cluster)))
    # filter the results based on the input
  }
}

