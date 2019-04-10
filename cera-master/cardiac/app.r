# Jun Xu 16/04/2018
# shiny application to enable Chromatin Enhanced RNA Analysis
# upload cardiac scRNA data per user-selected date
# calculate different metrics based on known cell clusters
# table and plot view with slice bars to filter key regulatory genes
# two global variables: m<day><cluster> (e.g. m301), and dat3d_cluster<day>

library(shiny)
library(scales)
library(zip)
library(dplyr)
library(ggplot2)
require(reshape2)
library(DESeq)
library(FedData)
source('CERA_cardiac_functions.r')

# need to change to user session specific
inputday <- readline(prompt='Please select Day 0, 2, 5, 15 or 30: ')
species <- 1   # human

if (!exists(paste0('m',inputday,1))) {
  
  # read in clustered scRNA data
  data <- readRDS(paste0('day',inputday,'_separate_final.RDS'))
  data.expr <- data$Exprs
  dat3d <- data$tSNE # x,y,z coordinates for tSNE plot
  data.cluster <- data$Clusters
  assign(paste0('dat3d_cluster',inputday),
         as.data.frame(cbind(dat3d, data.cluster)), envir=.GlobalEnv) # tSNE and cluster combined
  rownames(data.expr) <- unlist(as.character(
                         lapply(strsplit(as.character(rownames(data.expr)), '_'), '[', 1)))
  
  rm(data)
  gc()
  
  # correct the duplicated gene names in the data
  if (inputday==0) {
    rownames(data.expr)[1201] <- 'RGS5.1'
    rownames(data.expr)[3012] <- 'CYB561D2.1'
    rownames(data.expr)[4722] <- 'MATR3.1'
    rownames(data.expr)[9683] <- 'DNAJC9-AS1.1'
    rownames(data.expr)[10123] <- 'EMG1.1'
    rownames(data.expr)[10551] <- 'LINC01481.1'
    rownames(data.expr)[12896] <- 'COG8.1'
  } else if (inputday==2) {
    rownames(data.expr)[1228] <- 'RGS5.1'
    rownames(data.expr)[3132] <- 'CYB561D2.1'
    rownames(data.expr)[4866] <- 'MATR3.1'
    rownames(data.expr)[9908] <- 'DNAJC9-AS1.1'
    rownames(data.expr)[10363] <- 'EMG1.1'
    rownames(data.expr)[10799] <- 'LINC01481.1'
    rownames(data.expr)[13216] <- 'COG8.1'
    rownames(data.expr)[13534] <- 'TMEM256-PLSCR3.1'
  } else if (inputday==5) {
    rownames(data.expr)[1265] <- 'RGS5.1'
    rownames(data.expr)[3209] <- 'CYB561D2.1'
    rownames(data.expr)[5018] <- 'MATR3.1'
    rownames(data.expr)[10166] <- 'DNAJC9-AS1.1'
    rownames(data.expr)[10630] <- 'EMG1.1'
    rownames(data.expr)[11085] <- 'LINC01481.1'
    rownames(data.expr)[13556] <- 'COG8.1'
  } else if (inputday==15) {
    rownames(data.expr)[1236] <- 'RP11-466F5.8.1'
    rownames(data.expr)[3145] <- 'RRP9.1'
    rownames(data.expr)[4883] <- 'ANKHD1.1'
    rownames(data.expr)[9831] <- 'CISD1.1'
    rownames(data.expr)[10284] <- 'FBXL14.1'
    rownames(data.expr)[10739] <- 'OS9.1'
    rownames(data.expr)[11121] <- 'ATP6V0A2.1'
    rownames(data.expr)[13110] <- 'CTD-2600O9.2.1'
  } else if (inputday==30) {
    rownames(data.expr)[1253] <- 'RGS5.1'
    rownames(data.expr)[3190] <- 'CYB561D2.1'
    rownames(data.expr)[3641] <- 'SCHIP1.1'
    rownames(data.expr)[4986] <- 'MATR3.1'
    rownames(data.expr)[10086] <- 'DNAJC9-AS1.1'
    rownames(data.expr)[10537] <- 'EMG1.1'
    rownames(data.expr)[13415] <- 'COG8.1'
    rownames(data.expr)[13737] <- 'TMEM256-PLSCR3.1'
  }
  
  gc()
  
  # separate into available clusters
  for (j in 1:length(table(data.cluster))) {
    temp <- data.expr[,data.cluster==j]
    assign(paste0('m',inputday,j), temp, envir=.GlobalEnv)
  }

  rm(data.expr, dat3d, data.cluster, temp)
  gc()
}

#------------------------------------------------------------------------

# Define UI
ui <- fluidPage(

  # App title ----
  titlePanel(paste('Chromatin Enhanced RNA Analysis - cardiac scRNA Day', inputday)),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar selections
    sidebarPanel(
      fluidRow(
        h5('Please navigate the tabs to get results before downloading'),
        downloadButton('downloadData', 'Download',
                       style='color: #fff; background-color: #DDA0DD'),
        actionButton('reset', 'Reset', icon('paper-plane'), width = 105,
                     style='color: #fff; background-color: #87CEFA'),
        actionButton('maxrange', 'Max Range', width = 105,
                     style='color: #000; background-color: #ADFF2F'),
        actionButton('refresh', 'Refresh', icon('refresh'), width = 105,
                     style='color: #000; background-color: #FFD700'),
        HTML("<br><br><br>")
        ),
      
      fluidRow(
        div(style='display: inline-block;vertical-align:top; width: 250px;',
            selectizeInput('clusters', label = 'Multiple clusters:', 
                   choices = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'),
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
                                  'CV of Chromatin Pattern' = 4, 'Repression score' = 5,
                                  'Cell Proportion' = 6, 'Mean Expression' = 7, 'discordance score' = 8, 
                                  'Correlation' = 9, 'Differential Expression' = 10, 
                                  'Co-expression' = 11), selected = 8))),
      fluidRow(
        sliderInput('range_width', 'H3K27me3 width:', 
                   min = 0, max = 70000, value = c(20000,70000), step = 100),
        sliderInput('range_CV', 'Coefficient of Variation:', 
                   min = 0, max = 7, value = c(0,7), step = 0.001),
        sliderInput('range_Repression', 'Repression score:',
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
              column(2,tableOutput('sum_m4'))),
            fluidRow(
              column(12,plotOutput('sum_plot')))),
        tabPanel('Filtering genes', 
            fluidRow(
              tableOutput('filter_ratio'),
              tableOutput('filter_exp'))),
        tabPanel('Cluster 1', column(12,'', tableOutput('cluster1_table')),
                              column(12,'', plotOutput('cluster1_plot')),
                              h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Cluster 2', column(12,'', tableOutput('cluster2_table')),
                              column(12,'', plotOutput('cluster2_plot')),
                              h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Cluster 3', column(12,'', tableOutput('cluster3_table')),
                              column(12,'', plotOutput('cluster3_plot')),
                              h5('     Note: a definition of the metrics is available on the Glossary tab')),
        tabPanel('Cluster 4', column(12,'', tableOutput('cluster4_table')),
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
                 h5(strong('CV:'), 'Coefficient of Variation for the genes\' H3K27me3 breadth 
                    among the 111 reference cell types'),
                 h5(strong('Repression scores:'), 'The weight calculated from a gene\'s consistent 
                    association with broad H3K27me3 domains'),                 
                 h5(strong('Proportion:'), 'the proportion of cells with the gene expressed 
                    in each cluster'),
                 h5(strong('mean_exp:'), 'mean of the original gene expression without applying the Repression score'),                 
                 h5(strong('mean_discordance:'), 'mean of the product of the Repression score and the log-transformed gene 
                    expression'),
                 h5(strong('CC:'), 'Pearson correlation coefficient of the expression 
                    between a gene and the reference gene'),  
                 h5(strong('DE:'), 'differential expression ratio between reference gene is 
                    or not expressed'),  
                 h5(strong('CO:'), 'Co-expression ratio for the selected and the reference gene'),
                 h5(strong('GO:'), 'related to a known Gene Ontology'),
                 h5(strong('Expression:'), 'Gene expression in each cell')),
        id = 'tab'
        )
    )
  )
)

#------------------------------------------------------------------------

# Define server logic for slider
server <- function(input, output, session) {

  # selections for each cluster (g: gene, rg: regulator, s: sym, d: day, c: cluster, rs: reset)
  session$userData$g0 <- ''; session$userData$g1 <- ''; session$userData$g2 <- ''; 
  session$userData$g3 <- ''; session$userData$g4 <- ''; session$userData$rg0 <- ''; 
  session$userData$rg1 <- ''; session$userData$rg2 <- ''; session$userData$rg3 <- '';
  session$userData$rg4 <- ''; session$userData$s0 <- 0; session$userData$s1 <- 0; 
  session$userData$s2 <- 0; session$userData$s3 <- 0; session$userData$s4 <- 0;  
  session$userData$d0 <- -1; session$userData$d1 <- -1; session$userData$d2 <- -1; 
  session$userData$d3 <- -1; session$userData$d4 <- -1; session$userData$w0 <- 0; 
  session$userData$w1 <- 0; session$userData$w2 <- 0; session$userData$w3 <- 0; 
  session$userData$w4 <- 0; session$userData$c0 <- ''; session$userData$c1 <- '';
  session$userData$c0p <- ''; session$userData$rs0 <- 0; session$userData$rs1 <- 0; 
  session$userData$rs2 <- 0; session$userData$rs3 <- 0; session$userData$rs4 <- 0; 
  session$userData$rsm <- 0; session$userData$count <- 1;
  session$userData$ranges <- c(1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10,
                               1e10,-1e10,1e10,-1e10,1e10,-1e10,1e10,-1e10)
  
  # set size limit to 100MB
  options(shiny.maxRequestSize=100*1024^2)
  
  # actions for reset button
  observeEvent(input$reset, {
    if (input$tab == 'Summary') {session$userData$rs0 <- 0}
    if (input$tab == 'Cluster 1') {session$userData$rs1 <- 0}
    if (input$tab == 'Cluster 2') {session$userData$rs2 <- 0}
    if (input$tab == 'Cluster 3') {session$userData$rs3 <- 0}
    if (input$tab == 'Cluster 4') {session$userData$rs4 <- 0}
    if (input$tab == 'Multiple') {session$userData$rsm <- 0}
  }) 
  
  observeEvent(input$maxrange, {
    updateSliderInput(session, 'range_width', 
                      min = 0, max = 70000, value = c(0,70000), step = 100)
    updateSliderInput(session, 'range_CV', 
                      min = 0, max = 7, value = c(0,7), step = 0.001)
    updateSliderInput(session, 'range_Repression', 
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
    summary.page(session, input, inputday, species, 1)}, ignoreNULL = FALSE)
  output$sum_m1 <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {sum1()})
  })
  
  sum2 <- eventReactive(input$refresh, {
    summary.page(session, input, inputday, species, 2)}, ignoreNULL = FALSE)
  output$sum_m2 <- renderTable({sum2()})
  
  sum3 <- eventReactive(input$refresh, {
    summary.page(session, input, inputday, species, 3)}, ignoreNULL = FALSE)
  output$sum_m3 <- renderTable({sum3()})
  
  sum4 <- eventReactive(input$refresh, {
    summary.page(session, input, inputday, species, 4)}, ignoreNULL = FALSE)
  output$sum_m4 <- renderTable({sum4()})
  
  sumplot <- eventReactive(input$refresh, {
    if (input$common!='') {
      data.expr <- get(paste0('m',inputday,1))
      for (j in 2:length(table(get(paste0('dat3d_cluster',inputday))$data.cluster))) {
        data.expr <- cbind(data.expr, get(paste0('m',inputday,j)))
      }
      dat3d_all <- cbind(get(paste0('dat3d_cluster',inputday), envir=.GlobalEnv),
                         data.expr[input$common,])
      rm(data.expr)
      gc()
      colnames(dat3d_all) <- c('V1','V2','V3','cluster','gene_exp')
      tSNE <- tsne2d.clustcolor(input$common,dat3d_all,inputday)
      names(tSNE) <- c('dat3d_all','max_cluster','color_codes')
      rm(dat3d_all)
      gc()
      # Create the tSNE plot
      ggplot(tSNE$dat3d_all, aes(tSNE1,  tSNE2)) + 
        geom_point(aes(color = color), alpha = 0.5, size = 1.5) + 
        coord_fixed() +
        scale_color_manual(values = tSNE$color_codes[1:(tSNE$max_cluster + 1)]) +  
        theme_bw() + theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) + 
        theme(legend.title = element_blank()) + 
        labs(title = input$common) + 
        xlab('tSNE 1') +
        ylab('tSNE 2') +
        theme(plot.title = element_text(face = 'bold', colour = '#000000', size = 16)) +
        theme(axis.text = element_text(size = 12)) +
        theme(axis.title = element_text(size = 12))
    }
  }, ignoreNULL = FALSE)
  output$sum_plot <- renderPlot({sumplot()})

  # expression ratio for the filtering genes
  filterratio <- eventReactive(input$refresh, {
    p <- matrix(0,length(input$filter),5)   # create the result table
    colnames(p) <- c('Gene','cluster1','cluster2','cluster3','cluster4')
    p[,1] <- input$filter
    rownames(p) <- p[,1]
    for (c in 1:4) {
      if (exists(paste0('m',inputday,c))) {
        cell_tot <- ncol(get(paste0('m',inputday,c)))   # total number of cells in the matrix
        i <- 1;
        if (!is.null(input$filter)) {
          for (fg in input$filter) {
            # subset of data which has selected genes expressed
            x <- get(paste0('m',inputday,c))[,which(
              get(paste0('m',inputday,c))[fg,,drop=FALSE]!=0),drop=FALSE]
            if (is.null(x)) {xc <- 0} else {xc <- ncol(x)}
            p[i,c+1] <- percent(xc/cell_tot)
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
    colnames(q) <- c('Gene','cluster1','cluster2','cluster3','cluster4')
    q[,1] <- input$filter
    rownames(q) <- q[,1]
    for (c in 1:4) {
      if (exists(paste0('m',inputday,c))) {
        i <- 1;
        if (!is.null(input$filter)) {
          for (fg in input$filter) {
            # subset of data which has selected genes expressed
            m1 <- get(paste0('m',inputday,c))[fg,which(
                get(paste0('m',inputday,c))[fg,,drop=FALSE]!=0),drop=FALSE]
            if (!is.null(m1)) {
              if (input$zeroexpression == 1) {
                q[i,c+1] <- round(apply(m1,1,mean),3)  # including zero expressed cells
              } else {
                q[i,c+1] <- round(apply(m1,1,function(x){mean(x[x>0])}),3)  # excluding
              }
              i <- i + 1
            }
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
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,1)), 1)
    if (exists('r1',session$userData) & exists(paste0('m',inputday,1))) {
      session$userData$o1 <- filtered.results(session, session$userData$r1, inputday, 
                                              input, 'cluster 1')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster1_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster1t()
    })
  })
  
  cluster1p <- eventReactive(input$refresh, {
    # need to rebuild result table as the _table and _plot sections run independently
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,1)), 1)
    if (exists('r1',session$userData) & exists(paste0('m',inputday,1))) {
      session$userData$o1 <- filtered.results(session, session$userData$r1, inputday, 
                                              input, 'cluster 1')$detail
      # change matrix format for ggplot
      mm <- melt(get(paste0('m',inputday,1))[rownames(session$userData$o1),,drop=FALSE])
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
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,2)), 2)
    if (exists('r2',session$userData) & exists(paste0('m',inputday,2))) {
      session$userData$o2 <- filtered.results(session, session$userData$r2, inputday, 
                                              input, 'cluster 2')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster2_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster2t()
    })
  })
  
  cluster2p <- eventReactive(input$refresh, {
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,2)), 2)
    if (exists('r2',session$userData) & exists(paste0('m',inputday,2))) {
      session$userData$o2 <- filtered.results(session, session$userData$r2, inputday, 
                                              input, 'cluster 2')$detail
      
      mm <- melt(get(paste0('m',inputday,2))[rownames(session$userData$o2),,drop=FALSE])
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
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,3)), 3)
    if (exists('r3',session$userData) & exists(paste0('m',inputday,3))) {
      session$userData$o3 <- filtered.results(session, session$userData$r3, inputday, 
                                              input, 'cluster 3')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster3_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster3t()
    })
  })
  
  cluster3p <- eventReactive(input$refresh, {
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,3)), 3)
    if (exists('r3',session$userData) & exists(paste0('m',inputday,3))) {
      session$userData$o3 <- filtered.results(session, session$userData$r3, inputday, 
                                              input, 'cluster 3')$detail
      mm <- melt(get(paste0('m',inputday,3))[rownames(session$userData$o3),,drop=FALSE])
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
    cardiac.results(session, inputday,species, input, get(paste0('m',inputday,4)), 4)
    if (exists('r4',session$userData) & exists(paste0('m',inputday,4))) {
      session$userData$o4 <- filtered.results(session, session$userData$r4, inputday, 
                                              input, 'cluster 4')$detail
    }
  }, ignoreNULL = FALSE)
  output$cluster4_table <- renderTable({
    withProgress(message = 'Refreshing, please wait...', value = 0, {
      cluster4t()
    })
  })
  
  cluster4p <- eventReactive(input$refresh, {
    cardiac.results(session, inputday, species, input, get(paste0('m',inputday,4)), 4)
    if (exists('r4',session$userData) & exists(paste0('m',inputday,4))) {
      session$userData$o4 <- filtered.results(session, session$userData$r4, inputday, 
                                              input, 'cluster 4')$detail
      mm <- melt(get(paste0('m',inputday,4))[rownames(session$userData$o4),,drop=FALSE])
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
    
    # calculate the result if any gene, day, symbol or cluster selection has been changed
    if (session$userData$d0 != inputday | 
        !all(session$userData$g0==session$userData$g1) |
        session$userData$s0 != input$sym | 
        !all(session$userData$c0==session$userData$c1)) {
      for (item in input$clusters) {
        i <- substr(item,9,9)   # extract the cluster numbers
        if (merge == 0) {       # first matrix
          m <- get(paste0('m',inputday,i))
          merge <- 1  # set flag for second matrix
        } else {
          m <- cbind(m, get(paste0('m',inputday,i)))  # merge matrices
        }
      }
      
      session$userData$rmix <- 
        matrix.to.results(m, species, input$filter, input$sym, input$reference, session) 
      
      # update the selections if the refreshing is NOT triggered by new gene selection
      if (all(session$userData$g0 == session$userData$g1)) {
        # update the gene selection for new files, g is updated inside matrix.to.results
        updateSelectizeInput(session, 'filter', choices = session$userData$g,
                             selected =input$filter, options = list(maxItems=3))
        updatesliders(session, input)
      }
      session$userData$d0 <- inputday
      session$userData$s0 <- input$sym
      session$userData$c0 <- session$userData$c1
      session$userData$g0 <- session$userData$g1  # record the selections
    }
    if (exists('rmix',session$userData)) {
      session$userData$om <- filtered.results(session, session$userData$rmix, inputday, input,'')$detail
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
    
    # calculate the result if any gene, day, symbol or cluster selection has been changed
    if (session$userData$d0 != inputday |
        !all(session$userData$g0==session$userData$g1) |
        session$userData$s0 != input$sym |
        !all(session$userData$c0p==session$userData$c1)) {
      for (item in input$clusters) {
        i <- substr(item,9,9)   # extract the cluster numbers
        if (merge == 0) {       # first matrix
          m <- get(paste0('m',inputday,i))
          merge <- 1  # set flag for second matrix
        } else {
          m <- cbind(m, get(paste0('m',inputday,i)))  # merge matrices
        }
      }
      
      session$userData$mp <- m        # keep for filtering usage
      
      session$userData$rmix <- 
        matrix.to.results(m, species, input$filter, input$sym, input$reference, session) 
      
      # update the selections if the refreshing is NOT triggered by new gene selection
      if (all(session$userData$g0 == session$userData$g1)) {
        # update the gene selection for new files, g is updated inside matrix.to.results
        updateSelectizeInput(session, 'filter', choices = session$userData$g,
                             selected =input$filter, options = list(maxItems=3))
        updatesliders(session, input)
      }
      session$userData$d0 <- inputday
      session$userData$s0 <- input$sym
      session$userData$c0p <- session$userData$c1
      session$userData$g0 <- session$userData$g1  # record the selections
    }
    if (exists('rmix',session$userData)) {
      session$userData$om <- filtered.results(session, session$userData$rmix, inputday, input,'mixed')$detail
      mm <- melt(session$userData$mp[rownames(session$userData$om),,drop=FALSE])
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
      paste('output','zip',sep='.')
    },
    content = function(fname) {
      fs = prepare.files(input, session)
      zip(zipfile=fname, files=fs)
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
