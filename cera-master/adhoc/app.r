# Jun Xu 16/04/2018
# shiny application to enable Chromatin Enhanced scRNA Analysis
# users can upload up to 4 matrices to analyse

library(shiny)
library(scales)
library(zip)
library(dplyr)
library(ggplot2)
require(reshape2)
library(DESeq)
library(FedData)
source('CERA_adhoc_functions.r')

cleanMem <- function(n=10) { for (i in 1:n) gc() }

# Define UI
ui <- fluidPage(
  
  # App title ----
  titlePanel('ad-hoc TRIAGE analysis platform'),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar selections
    sidebarPanel(
      h6('1) Please upload normalised expression matrix in .csv file with comma as separator'),
      h6('2) CSV content: Row Names - gene names, Column Names - barcodes/samples'),
      h6('3) Gene names should NOT be duplicated and should align with nomenclature (https://en.wikipedia.org/wiki/Gene_nomenclature)'),
      h6('4) After changing criteria, please press "Refresh" to get results updated'),
      h6('5) Use "Reset" + "Refresh" buttons to clear all criteria changes'),
      h6('6) You need to navigate through tabs before downloading results'),
      h6('7) Use "Max Range" button to clear all filtering criteria'),
      h6(),
      
      fluidRow(
        downloadButton('downloadData', 'Download',
                       style='color: #fff; background-color: #DDA0DD'),
        actionButton('reset', 'Reset', icon('paper-plane'), width = 105,
                     style='color: #fff; background-color: #87CEFA'),
        actionButton('maxrange', 'Max Range', icon('resize-horizontal'), width = 105,
                     style='color: #000; background-color: #ADFF2F'),
        actionButton('refresh', 'Refresh', icon('refresh'), width = 105,
                     style='color: #000; background-color: #FFD700'),
        HTML("<br><br><br>")
      ),
      
      fluidRow(
        column(6,
               fileInput('file1', 'matrix 1', accept = c('.csv'), width = 250)),
        column(6,
               fileInput('file2', 'matrix 2', accept = c('.csv'), width = 250)),
        column(6,
               fileInput('file3', 'matrix 3', accept = c('.csv'), width = 250)),
        column(6,
               fileInput('file4', 'matrix 4', accept = c('.csv'), width = 250)),
        
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('species', label = 'Species:', 
                           choices = c('Human' = 1, 'C.interstinalis' = 2,
                                       'Chicken' = 3, 'Cow' = 4, 'Fruitfly' = 5,
                                       'Mouse' = 6, 'Pig' = 7, 'Rat' = 8, 'Zebrafish' = 9),
                           selected = 1, options = list(maxItems=1))),
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('clusters', label = 'Multiple Matrices:', 
                           choices = c('Matrix 1', 'Matrix 2', 'Matrix 3', 'Matrix 4'),
                           selected = NULL, options = list(maxItems=4))),
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectInput('zeroexpression', 'Cells on Boxplot',
                        choices = list('including null expressions' = 1, 
                                       'excluding null expressions' = 2), 
                        selected = 1)),  # zero expressions
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('common', label = 'Gene for tSNE plot:', 
                           choices = NULL, options = list(maxItems=1))),
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('reference', label = 'Reference gene for CC/DE/CO:', 
                           choices = NULL, options = list(maxItems=1))),
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('filter', label = 'Filtering Genes:', 
                           choices = NULL, options = list(maxItems=10))),
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectInput('sym', label = 'ALL/ANY for multiple filtering genes',
                        choices = list('ALL' = 1, 'ANY' = 2), selected = 1)),  # Input: ALL/ANY
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectInput('specific', label = 'Cluster specific', 
                        choices = list('Yes' = 1, 'No' = 2), selected = 2)),  # Cell type specific
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectInput('sort', 'Sorted by',
                        choices = list('Gene' = 1, 'Transcription Factor' = 2, 'H3K27me3 width' = 3, 
                                       'RTS score' = 4, 'Cell Proportion' = 5, 'Mean Expression' = 6, 
                                       'discordance score' = 7, 'Correlation' = 8, 'Differential Expression' = 9, 
                                       'Co-expression' = 10), selected = 7))),
      fluidRow(
        sliderInput('range_width', 'H3K27me3 width:', 
                    min = 0, max = 70000, value = c(20000,70000), step = 100),
        sliderInput('range_rts', 'RTS score:',
                    min = 0, max = 1, value = c(0.25,1), step = 0.001),
        sliderInput('range_percent', 'Proportion of expression:',
                    min = 0, max = 100, value = c(30,100), step = 0.1),
        sliderInput('range_mean', 'Mean expression:',
                    min = 0, max = 750, value = c(0,750), step = 0.001),
        sliderInput('range_discordance', 'discordance score:',
                    min = -50, max = 10, value = c(-50,10), step = 0.001),
        sliderInput('range_cc', 'Correlation with the reference:',
                    min = -1, max = 1, value = c(-1,1), step = 0.01),
        sliderInput('range_co', 'Co-expression ratio with the reference:',
                    min = 0, max = 1, value = c(0,1), step = 0.01)
      )
      , width = 4),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(
        tabPanel('Summary',
                 fluidRow(
                   column(2,tableOutput('sum_m1')),
                   column(2,tableOutput('sum_m2')),
                   column(2,tableOutput('sum_m3')),
                   column(2,tableOutput('sum_m4')))),
        tabPanel('Filtering genes', 
                 fluidRow(
                   tableOutput('filter_ratio'),
                   tableOutput('filter_exp'))),
        tabPanel('Matrix 1', column(12,'', tableOutput('cluster1_table')),
                 column(12,'', plotOutput('cluster1_plot')),
                 h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Matrix 2', column(12,'', tableOutput('cluster2_table')),
                 column(12,'', plotOutput('cluster2_plot')),
                 h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Matrix 3', column(12,'', tableOutput('cluster3_table')),
                 column(12,'', plotOutput('cluster3_plot')),
                 h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Matrix 4', column(12,'', tableOutput('cluster4_table')),
                 column(12,'', plotOutput('cluster4_plot')),
                 h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Multiple',  column(12,'', tableOutput('multiple_table')),
                 column(12,'', plotOutput('multiple_plot')),
                 h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Glossary',
                 br(),
                 h5(strong('Gene:'),'Name of the genes filtered by the weighted expression 
                    with H3K27me3 breadth pattern'),
                 h5(strong('TF:'), 'Whether the gene is a known Transcription Factor (1=Yes, 0=No)'),     
                 h5(strong('k27bp:'), 'Mean width of H3K27me3 breadth near the gene, based 
                    on the 111 reference cell types in ROADMAP'),                 
                 h5(strong('RTS scores:'), 'The weight calculated from a gene\'s consistent 
                    association with broad H3K27me3 domains'),                 
                 h5(strong('Proportion:'), 'the proportion of cells with the gene expressed 
                    in each cluster'),
                 h5(strong('mean_exp:'), 'mean of the original gene expression without applying the RTS score'),                 
                 h5(strong('mean_discordance:'), 'mean of the product of the RTS score and the log-transformed gene 
                    expression'),
                 h5(strong('CC:'), 'Pearson correlation coefficient of the expression 
                    between a gene and the reference gene'),  
                 h5(strong('DE:'), 'differential expression ratio between reference gene is 
                    or not expressed'),  
                 h5(strong('CO:'), 'Co-expression ratio for the selected and the reference gene'),
                 h5(strong('Expression:'), 'Gene expression in each cell')),
        id = 'tab'
                 )
                 )
  )
  )

#------------------------------------------------------------------------

# Define server logic for slider
server <- function(input, output, session) {
  
  observe({
    # periodically collect
    invalidateLater(1000,session)
    cleanMem()
  })
  
  # selections for each matrix 
  # (g: gene, rg: regulator, s: sym, sp: species, c: matrix, rs: reset, fn: file names, z: file size)
  session$userData$g0 <- ''; session$userData$g1 <- ''; session$userData$g2 <- ''; 
  session$userData$g3 <- ''; session$userData$g4 <- ''; session$userData$rg0 <- ''; 
  session$userData$rg1 <- ''; session$userData$rg2 <- ''; session$userData$rg3 <- '';
  session$userData$rg4 <- ''; session$userData$s0 <- 0; session$userData$s1 <- 0; 
  session$userData$s2 <- 0; session$userData$s3 <- 0; session$userData$s4 <- 0;  
  session$userData$d0 <- -1; session$userData$d1 <- -1; session$userData$d2 <- -1; 
  session$userData$d3 <- -1; session$userData$d4 <- -1; session$userData$sp0 <- 0; 
  session$userData$sp1 <- 0; session$userData$sp2 <- 0; session$userData$sp3 <- 0; 
  session$userData$sp4 <- 0; session$userData$c0 <- ''; session$userData$c1 <- ''; 
  session$userData$c0p <- ''; session$userData$rs0 <- 0; session$userData$rs1 <- 0; 
  session$userData$rs2 <- 0; session$userData$rs3 <- 0; session$userData$rs4 <- 0; 
  session$userData$rsm <- 0; session$userData$fn1 <- ''; session$userData$fn2 <- ''; 
  session$userData$fn3 <- ''; session$userData$fn4 <- ''; session$userData$z1 <- 0; 
  session$userData$z2 <- 0; session$userData$z3 <- 0; session$userData$z4 <- 0;
  session$userData$ranges <- c(1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,
                               1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10)
  
  # set size limit to 500MB
  options(shiny.maxRequestSize=500*1024^2)
  
  # actions for reset button
  observeEvent(input$reset, {
    if (input$tab == 'Summary') {session$userData$rs0 <- 0}
    if (input$tab == 'Matrix 1') {session$userData$rs1 <- 0}
    if (input$tab == 'Matrix 2') {session$userData$rs2 <- 0}
    if (input$tab == 'Matrix 3') {session$userData$rs3 <- 0}
    if (input$tab == 'Matrix 4') {session$userData$rs4 <- 0}
    if (input$tab == 'Multiple') {session$userData$rsm <- 0}
  }) 
  
  observeEvent(input$maxrange, {
    updateSliderInput(session, 'range_width', 
                      min = 0, max = 70000, value = c(0,70000), step = 100)
    updateSliderInput(session, 'range_rts', 
                      min = 0, max = 1, value = c(0,1), step = 0.001)
    updateSliderInput(session, 'range_percent', 
                      min = 0, max = 100, value = c(0,100), step = 0.1)
    updateSliderInput(session, 'range_mean', 
                      min = 0, max = 750, value = c(0,750), step = 0.001)
    updateSliderInput(session, 'range_discordance', 
                      min = -50, max = 10, value = c(-50,10), step = 0.001)
    updateSliderInput(session, 'range_cc', 
                      min = -1, max = 1, value = c(-1,1), step = 0.01)
    updateSliderInput(session, 'range_co', 
                      min = 0, max = 1, value = c(0,1), step = 0.01)
  })
  
  # summary of all clusters
  sum1 <- eventReactive(input$refresh, {
    summary.page(session, input, 0, 1)}, ignoreNULL = FALSE)
  output$sum_m1 <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {sum1()})
  })
  
  sum2 <- eventReactive(input$refresh, {
    summary.page(session, input, 0, 2)}, ignoreNULL = FALSE)
  output$sum_m2 <- renderTable({sum2()})
  
  sum3 <- eventReactive(input$refresh, {
    summary.page(session, input, 0, 3)}, ignoreNULL = FALSE)
  output$sum_m3 <- renderTable({sum3()})
  
  sum4 <- eventReactive(input$refresh, {
    summary.page(session, input, 0, 4)}, ignoreNULL = FALSE)
  output$sum_m4 <- renderTable({sum4()})
  
  # expression ratio for the filtering genes
  filterratio <- eventReactive(input$refresh, {
    p <- matrix(0,length(input$filter),5)   # create the result table
    colnames(p) <- c('Gene','Matrix 1','Matrix 2','Matrix 3','Matrix 4')
    p[,1] <- input$filter
    rownames(p) <- p[,1]
    for (c in 1:4) {
      if (exists(paste0('m0',c),session$userData)) {
        cell_tot <- ncol(get(paste0('m0',c),session$userData))   # total number of cells in the matrix
        i <- 1;
        if (!is.null(input$filter)) {
          for (fg in input$filter) {
            # subset of data which has selected genes expressed
            x <- get(paste0('m0',c),session$userData)[,which(
              get(paste0('m0',c),session$userData)[fg,,drop=FALSE]!=0),drop=FALSE]
            if (is.null(x)) {
              if (is.null(cell_tot)) {
                p[i,c+1] <- '-'
              } else {
                p[i,c+1] <- 0
              }
            } else {
              xc <- ncol(x)
              p[i,c+1] <- percent(xc/cell_tot)
            }
            i <- i + 1
          }
        }
      }
    }
    session$userData$op <- p
    p
  }, ignoreNULL = FALSE)
  output$filter_ratio <- renderTable({filterratio()})
  
  # mean of expression for the filtering genes
  filterexp <- eventReactive(input$refresh, {
    q <- matrix(0,length(input$filter),5)
    colnames(q) <- c('Gene','Matrix 1','Matrix 2','Matrix 3','Matrix 4')
    q[,1] <- input$filter
    rownames(q) <- q[,1]
    for (c in 1:4) {
      if (exists(paste0('m0',c),session$userData)) {
        i <- 1;
        if (!is.null(input$filter)) {
          for (fg in input$filter) {
            # subset of data which has selected genes expressed
            m1 <- get(paste0('m0',c),session$userData)[fg,which(
              get(paste0('m0',c),session$userData)[fg,,drop=FALSE]!=0),drop=FALSE]
            if (!is.null(m1)) {
              if (input$zeroexpression == 1) {
                q[i,c+1] <- round(apply(m1,1,mean),3)  # including zero expressed cells
              } else {
                q[i,c+1] <- round(apply(m1,1,function(x){mean(x[x>0])}),3)  # excluding
              }
            } else {
              q[i,c+1] <- '-'
            }
            i <- i + 1
          }
        }
      }
    }
    session$userData$oq <- q
    q
  }, ignoreNULL = FALSE)
  output$filter_exp <- renderTable({filterexp()})
  
  # details for matrix 1
  cluster1t <- eventReactive(input$refresh, {
    prepare.results(session, input, 1)
    if (exists('r1',session$userData) & exists(paste0('m0',1),session$userData)) {
      session$userData$o1 <- filtered.results(session, session$userData$r1, 0, 
                                              input, 'Matrix 1')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster1_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster1t()
    })
  })
  
  cluster1p <- eventReactive(input$refresh, {
    prepare.results(session, input, 1)
    # need to rebuild result table as the _table and _plot sections run independently
    if (exists('r1',session$userData) & exists(paste0('m0',1),session$userData)) {
      session$userData$o1 <- filtered.results(session, session$userData$r1, 0, 
                                              input, 'Matrix 1')$detail
      # change matrix format for ggplot
      mm <- melt(as.matrix(get(paste0('m0',1),session$userData)[rownames(session$userData$o1),,drop=FALSE]))
      colnames(mm) <- c('Gene','Cell','Expression')
      if (input$zeroexpression == 2) {
        mm <- mm[which(mm$Expression!=0),]
      }
      ggplot(data=mm, aes(x=Gene, y=Expression)) + geom_jitter(alpha = 0.3, color = 'grey') +
        geom_boxplot(aes(fill=Gene), outlier.shape=NA) + 
        scale_y_continuous(limits = quantile(mm$Expression, c(0.1, 0.9)))
    }
  }, ignoreNULL = FALSE)
  output$cluster1_plot <- renderPlot({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster1p()
    })
  })
  
  # details for matrix 2
  cluster2t <- eventReactive(input$refresh, {
    prepare.results(session, input, 2)
    if (exists('r2',session$userData) & exists(paste0('m0',2),session$userData)) {
      session$userData$o2 <- filtered.results(session, session$userData$r2, 0, 
                                              input, 'Matrix 2')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster2_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster2t()
    })
  })
  
  cluster2p <- eventReactive(input$refresh, {
    prepare.results(session, input, 2)
    if (exists('r2',session$userData) & exists(paste0('m0',2),session$userData)) {
      session$userData$o2 <- filtered.results(session, session$userData$r2, 0, 
                                              input, 'Matrix 2')$detail
      
      mm <- melt(as.matrix(get(paste0('m0',2),session$userData)[rownames(session$userData$o2),,drop=FALSE]))
      colnames(mm) <- c('Gene','Cell','Expression')
      if (input$zeroexpression == 2) {
        mm <- mm[which(mm$Expression!=0),]
      }
      ggplot(data=mm, aes(x=Gene, y=Expression)) + geom_jitter(alpha = 0.3, color = 'grey') +
        geom_boxplot(aes(fill=Gene), outlier.shape=NA, width = 0.6) + 
        scale_y_continuous(limits = quantile(mm$Expression, c(0.1, 0.9)))
    }
  }, ignoreNULL = FALSE)
  output$cluster2_plot <- renderPlot({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster2p()
    })
  })
  
  # details for matrix 3
  cluster3t <- eventReactive(input$refresh, {
    prepare.results(session, input, 3)
    if (exists('r3',session$userData) & exists(paste0('m0',3),session$userData)) {
      session$userData$o3 <- filtered.results(session, session$userData$r3, 0, 
                                              input, 'Matrix 3')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster3_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster3t()
    })
  })
  
  cluster3p <- eventReactive(input$refresh, {
    prepare.results(session, input, 3)
    if (exists('r3',session$userData) & exists(paste0('m0',3),session$userData)) {
      session$userData$o3 <- filtered.results(session, session$userData$r3, 0, 
                                              input, 'Matrix 3')$detail
      mm <- melt(as.matrix(get(paste0('m0',3),session$userData)[rownames(session$userData$o3),,drop=FALSE]))
      colnames(mm) <- c('Gene','Cell','Expression')
      if (input$zeroexpression == 2) {
        mm <- mm[which(mm$Expression!=0),]
      }
      ggplot(data=mm, aes(x=Gene, y=Expression)) + geom_jitter(alpha = 0.3, color = 'grey') +
        geom_boxplot(aes(fill=Gene), outlier.shape=NA) + 
        scale_y_continuous(limits = quantile(mm$Expression, c(0.1, 0.9)))
    }
  }, ignoreNULL = FALSE)
  output$cluster3_plot <- renderPlot({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster3p()
    })
  })
  
  # details for matrix 4
  cluster4t <- eventReactive(input$refresh, {
    prepare.results(session, input, 4)
    if (exists('r4',session$userData) & exists(paste0('m0',4),session$userData)) {
      session$userData$o4 <- filtered.results(session, session$userData$r4, 0, 
                                              input, 'Matrix 4')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster4_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster4t()
    })
  })
  
  cluster4p <- eventReactive(input$refresh, {
    prepare.results(session, input, 4)
    if (exists('r4',session$userData) & exists(paste0('m0',4),session$userData)) {
      session$userData$o4 <- filtered.results(session, session$userData$r4, 0, 
                                              input, 'Matrix 4')$detail
      mm <- melt(as.matrix(get(paste0('m0',4),session$userData)[rownames(session$userData$o4),,drop=FALSE]))
      colnames(mm) <- c('Gene','Cell','Expression')
      if (input$zeroexpression == 2) {
        mm <- mm[which(mm$Expression!=0),]
      }
      ggplot(data=mm, aes(x=Gene, y=Expression)) + geom_jitter(alpha = 0.3, color = 'grey') +
        geom_boxplot(aes(fill=Gene), outlier.shape=NA) + 
        scale_y_continuous(limits = quantile(mm$Expression, c(0.1, 0.9)))
    }
  }, ignoreNULL = FALSE)
  output$cluster4_plot <- renderPlot({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster4p()
    })
  })
  
  # details for multiple clusters
  clustermt <- eventReactive(input$refresh, {
    merge <- 0; reload <- 0  # intial flags
    session$userData$g1 <- input$filter; 
    session$userData$c1 <- input$clusters
    if (is.null(session$userData$g1)) {session$userData$g1 <- ''}  # transform NULL input
    if (is.null(session$userData$c1)) {session$userData$c1 <- ''}
    
    # calculate the result if any gene, species, symbol or cluster selection has been changed
    if (session$userData$d0 != 0 |
        session$userData$sp0 != input$species |
        !all(session$userData$g0 == session$userData$g1) |
        session$userData$s0 != input$sym | 
        !all(session$userData$c0 == session$userData$c1)) {
      for (item in input$clusters) {
        i <- substr(item,8,8)   # extract the cluster numbers
        if (merge == 0) {       # first matrix
          m <- get(paste0('m0',i),session$userData)
          merge <- 1  # set flag for second matrix
        } else {
          m <- cbind(m, get(paste0('m0',i),session$userData))  # merge matrices
        }
      }
      
      session$userData$rmix <- 
        matrix.to.results(m, input$species, input$filter, input$sym, input$reference, session) 
      
      # update the selections if the refreshing is NOT triggered by new gene selection
      if (all(session$userData$g0 == session$userData$g1)) {
        # update the gene selection for new files, g is updated inside matrix.to.results
        updateSelectizeInput(session, 'filter', choices = session$userData$g,
                             selected =input$filter, options = list(maxItems=3))
        updatesliders(session, input)
      }
      session$userData$d0 <- 0
      session$userData$sp0 <- input$species
      session$userData$s0 <- input$sym
      session$userData$c0 <- session$userData$c1
      session$userData$g0 <- session$userData$g1  # record the selections
    }
    if (exists('rmix',session$userData)) {
      session$userData$om <- filtered.results(session, session$userData$rmix, 0, input,'')$detail
    }
  }, ignoreNULL = FALSE)
  output$multiple_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      clustermt()
    })
  })
  
  clustermp <- eventReactive(input$refresh, {
    merge <- 0; reload <- 0  # intial flags
    session$userData$g1 <- input$filter; 
    session$userData$c1 <- input$clusters
    if (is.null(session$userData$g1)) {session$userData$g1 <- ''}  # transform NULL input
    if (is.null(session$userData$c1)) {session$userData$c1 <- ''}
    
    # recalculate the result if any gene, day, symbol or cluster selection has been changed
    if (session$userData$d0 != 0 |
        session$userData$sp0 != input$species |
        !all(session$userData$g0 ==session$userData$g1) |
        session$userData$s0 != input$sym |
        !all(session$userData$c0p ==session$userData$c1)) {
      for (item in input$clusters) {
        i <- substr(item,8,8)   # extract the cluster numbers
        if (merge == 0) {       # first matrix
          m <- get(paste0('m0',i),session$userData)
          merge <- 1  # set flag for second matrix
        } else {
          m <- cbind(m, get(paste0('m0',i),session$userData))  # merge matrices
        }
      }
      
      session$userData$mp <- m        # keep for filtering usage
      
      session$userData$rmix <- 
        matrix.to.results(m, input$species, input$filter, input$sym, input$reference, session) 
      
      # update the selections if the refreshing is NOT triggered by new gene selection
      if (all(session$userData$g0 == session$userData$g1)) {
        # update the gene selection for new files, g is updated inside matrix.to.results
        updateSelectizeInput(session, 'filter', choices = session$userData$g,
                             selected =input$filter, options = list(maxItems=3))
        updatesliders(session, input)
      }
      session$userData$d0 <- 0
      session$userData$sp0 <- input$species
      session$userData$s0 <- input$sym
      session$userData$c0p <- session$userData$c1
      session$userData$g0 <- session$userData$g1  # record the selections
    }
    if (exists('rmix',session$userData)) {
      session$userData$om <- filtered.results(session, session$userData$rmix, 0, input,'mixed')$detail
      mm <- melt(as.matrix(session$userData$mp[rownames(session$userData$om),,drop=FALSE]))
      colnames(mm) <- c('Gene','Cell','Expression')
      if (input$zeroexpression == 2) {
        mm <- mm[which(mm$Expression!=0),]
      }
      ggplot(data=mm, aes(x=Gene, y=Expression)) + geom_jitter(alpha = 0.3, color = 'grey') +
        geom_boxplot(aes(fill=Gene), outlier.shape=NA) + 
        scale_y_continuous(limits = quantile(mm$Expression, c(0.1, 0.9)))
    }
  }, ignoreNULL = FALSE)
  output$multiple_plot <- renderPlot({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      clustermp()
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('results','zip',sep='.')
    },
    content = function(fname) {
      curdir <- getwd()
      fs = prepare.files(session, input)
      zip(zipfile=fname, files=fs)
      setwd(curdir)
    },
    contentType = 'application/zip'
  )
  
  onStart = function() {
    cat('Doing application setup\n')
    
    onStop(function() {
      cat('Doing application cleanup\n')
      rm(list = setdiff(ls(), lsf.str()))
    })
  }
}

# Create Shiny app ----
shinyApp(ui, server)
