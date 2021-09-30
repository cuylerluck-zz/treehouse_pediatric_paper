library(shiny) #version 1.6.0
library(ggplot2) #version 3.3.5
library(pheatmap) #version 1.0.12
library(viridis) #version 0.6.1
library(plotly) #version 4.9.4.1
library(RColorBrewer) #version 1.1-2
library(shinycssloaders) #version 1.0.0
library(shinythemes) #version 1.2.0
library(data.table) #version 1.14.0
library(ggbiplot) #version 0.55
library(gridExtra) #version 2.3
library(tidyr) #version 1.1.3
library(heatmaply) #version 1.2.1
library(DT) #version 0.19
library(EnhancedVolcano) #version 1.6.0
library(ggpubr) #version 0.4.0



# Define UI 
ui <- fluidPage(
   
  navbarPage("Pediatric Cancer Cell Line Comparison",
             
             
             tabPanel("About",
                      mainPanel(
                      tags$h2("Transcriptional fidelity enhances cancer cell line selection in pediatric cancers"),
                      div(HTML("Cuyler Luck<sup>1,2,3</sup>, Katharine Yu<sup>2,3</sup>, Ross A. Okimoto<sup>1,4</sup>, and Marina Sirota<sup>2,3*</sup>")),
                      div(HTML('<font size="2">1. Department of Medicine, University of California, San Francisco, San Francisco, California, USA</font>')),
                      div(HTML('<font size="2">2. Bakar Computational Health Sciences Institute, University of California, San Francisco, San Francisco California, USA</font>')),
                      div(HTML('<font size="2">3. Department of Pediatrics, University of California, San Francisco, San Francisco California, USA</font>')),                
                      div(HTML('<font size="2">4. Helen Diller Family Comprehensive Cancer Center, University of California, San Francisco, San Francisco, California, USA</font>'))               
                      ),
                      mainPanel(tags$h5("Multi-omic technologies have allowed for comprehensive profiling of patient-derived tumor samples and the cell lines that are intended to model them. Yet, our understanding of how cancer cell lines reflect native pediatric cancers in the age of molecular subclassification remains unclear and represents a clinical unmet need. Here we use Treehouse public data to provide an RNA-seq driven analysis of 799 cancer cell lines, focusing on how well they correlate to 1,655 pediatric tumor samples spanning 12 tumor types. For each tumor type of interest, we present a ranked list of the most representative cell lines based on correlation of their transcriptomic profiles to those of the tumor. We found that most (8/12) tumor types best correlated to a cell line of the closest matched disease type. We furthermore showed that inferred molecular subtype differences in medulloblastoma significantly impacted correlation between medulloblastoma tumor samples and cell lines. Our results are available as an interactive web application to help researchers select cancer cell lines that more faithfully recapitulate pediatric cancer.", style = "position:relative;")),
                      
                      mainPanel(imageOutput("graphical_abstract", height = "100px", width = "auto"), style = "position:relative;")

             ),
             
             tabPanel("Explore the samples", 
                      mainPanel(
                        tags$h2("UMAP of all pediatric tumor samples and cell lines"),
                        tags$h4("uses the top 5000 most variable genes across all samples"),
                        withSpinner(plotlyOutput("UMAP_all", height = "600px", width = "auto"), type=5)
                      ),
                      
                      mainPanel(
                        tags$h2("UMAP of pediatric tumor samples only"),
                        tags$h4("uses the top 5000 most variable genes across all tumor samples"),
                        withSpinner(plotlyOutput("UMAP_tumor", height = "600px", width = "auto"), type=5)
                      ),
                      
                      mainPanel(
                        tags$h2("UMAP of cell lines only"),
                        tags$h4("uses the top 5000 most variable genes across all cell lines"),
                        withSpinner(plotlyOutput("UMAP_cell", height = "600px", width = "auto"), type=5)
                      )
             ),
             
             tabPanel("Summarized data",
               
                      
                      
                      mainPanel(
                        tags$h2("Correlations between matched tumor samples and cell line diseases"),
                        tags$h6("See publication for tumor-cell disease pairings."),
                        withSpinner(plotOutput("violinMaster", height = "500px"), type = 5)
                      ),
                      
                      
                      mainPanel(
                        tags$h2("Mean correlations summarised by tumor disease and cell line TCGA code"),
                        tags$h6("Correlations shown are summarized to mean values."),
                        withSpinner(plotOutput("master_heatmap"), type = 5)
                      )
               
             ),
             
             tabPanel("Query individual tumor diseases",
                      sidebarPanel(
                        selectInput("cancer", "Choose a pediatric tumor type", 
                                    choices = list("acute myeloid leukemia", "acute lymphoblastic leukemia", "alveolar rhabdomyosarcoma", "embryonal rhabdomyosarcoma", "ependymoma", "ewing sarcoma", "glioma", "medulloblastoma", "medulloblastoma_SHH_GRP3_like", "medulloblastoma_WNT_like", "neuroblastoma", "osteosarcoma", "rhabdomyosarcoma", "wilms tumor"))
                      ),
                      mainPanel(
                        downloadButton("downloadData", "Download all correlations"),
                        tags$h2(),
                        downloadButton("downloadData_cancer", "Download correlations of selected tumor type with all cell lines"),
                        tags$h2(),
                        downloadButton("downloadData_cancer_matched", "Download correlations of selected tumor type with matched cell lines"),
                        tags$h2()
                      ),
                      mainPanel(
                        tags$h2("PCAs of uncorrected (left) and ComBat corrected (right) tumor samples"),
                        withSpinner(plotlyOutput("PCAs"), type = 5)
                      ),
                      
                      mainPanel(
                        tags$h2("Violin plots ranked by median correlation with cell line TCGA Code (left) or with the top 10 individual cell lines (right)"),
                        withSpinner(plotOutput("violins", height = "500px", width = "1500px"), type = 5)
                        
                      ),
                      
                      mainPanel(
                        tags$h2("Heatmap of correlations between disease-specific tumor samples and matched cell lines"),
                        withSpinner(plotOutput("matched_heatmap_plot", height = "700px"), type = 5)
                        
                      ),
                      
                      mainPanel(
                        
                        tags$h2("Median correlations between disease-specific tumor samples and all cell lines"),
                        downloadButton("downloadData_cancer_corrs", "Download median correlations of selected tumor type with all cell lines"),
                        
                        withSpinner(DTOutput("table_all_corrs"), type = 5)
                      )
                      
                    

                      
              ),
             
              tabPanel("Investigate medulloblastoma inferred subtyping",
                       mainPanel(
                         tags$h2("PCA of medulloblastoma expected counts"),
                         tags$h4("uses the top 5000 most variable genes across all samples"),
                         withSpinner(plotlyOutput("medullo_counts", width = "700px", height = "700px"),type=5)
                        
                       ),
                       
                       mainPanel(
                         tags$h2("EnhancedVolcano of WNT-like vs. SHH/Group 3-like medulloblastoma samples"),
                         withSpinner(plotOutput("evolcano", width = "1200px", height = "800px"), type = 5)
                       ),
                       
                       
                       
                       
                       mainPanel(
                         
                         tags$h2("Look up fold change and p-value for your favorite gene"),
                         tags$h2(),
                         downloadButton("downloadData_medullo_volcano", "Download the data in this table"),
                         tags$h2(),
                         withSpinner(DTOutput("evolcano_table"), type = 5)
                       ),
                       
                       mainPanel(
                         tags$h2("Correlations between medulloblastoma tumor samples and cell lines are impacted by inferred subtype"),
                         tags$h4("Bulk comparison"),
                         withSpinner(plotOutput("bulk_medullo"), type = 5),
                         tags$h4("Broken down by cell line"),
                         withSpinner(plotOutput("by_line_medullo"), type = 5),
                         tags$h5("Formatted p-values pre-multiple comparisons correction are shown. FDR-corrected p-values are 0.024 (D283 Med), 0.024 (D341 Med), 0.018 (Daoy), and 0.024 (ONS-76). Crossbars indicate mean correlations, which are compared in statistical tests (Wilcoxon test used in all comparisons).")
                       )
                       
                       
                       
                       
              )
             
             
             
             

  )
)

# Define server logic 
server <- function(input, output) {
  
  manualcolors<-c('forestgreen', 'red2', 'orange', 'cornflowerblue', 
                  'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                  'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                  "yellow1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                  'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                  "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro", 'black')
  
  tumorMeta = fread(file = "Intermediate_Data/tCompendiumv11_pediatric_selectDiseases_forComparison_meta.csv", header = TRUE)
  tumorMeta = dplyr::select(tumorMeta, th_sampleid, comparison_disease, Database)
  
  cellMetaV2 = fread(file = "Intermediate_Data/cellMetaV2.csv", header = TRUE)
  cellMetaV2 = dplyr::select(cellMetaV2, th_sampleid, source_sample_ID, tcga_code)
  
  
  tumorColors = manualcolors[1:14]
  names(tumorColors) = c(unique(tumorMeta$comparison_disease), "medulloblastoma_SHH_GRP3_like", "medulloblastoma_WNT_like")
  tumorColors = tumorColors[order(names(tumorColors))]
  
  cellColors = manualcolors[1:29]
  names(cellColors) = unique(cellMetaV2$tcga_code)
  cellColors = cellColors[order(names(cellColors))]
  

  output$UMAP_all = renderPlotly({
    umap_all_plotting = (readRDS("results/umap_plotting_all_samples.RDS"))
    
    for(ID in umap_all_plotting[umap_all_plotting$type == "cell_line",]$th_sampleid){
      umap_all_plotting[umap_all_plotting$th_sampleid == ID,]$th_sampleid = cellMetaV2[cellMetaV2$th_sampleid == ID,]$source_sample_ID
    }
    
    
    
    umap_all = ggplot(umap_all_plotting, aes(x = x, y = y, shape = type, fill = disease, label = th_sampleid)) + 
      geom_point(data = umap_all_plotting[umap_all_plotting$type == "cell_line",], size = 3, alpha = 0.5) +
      geom_point(data = umap_all_plotting[umap_all_plotting$type == "tumor_sample",],size = 3) +
      scale_fill_manual(values = tumorColors) +
      scale_shape_manual(values = c(21,24)) +
      xlab("UMAP1") + ylab("UMAP2")
    
    umap_all_ggplotly = ggplotly(umap_all)
    umap_all_ggplotly = style(umap_all_ggplotly, showlegend = FALSE, traces = 1:29)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "acute lymphoblastic leukemia", traces = 30)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "acute myeloid leukemia", traces = 31)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "alveolar rhabdomyosarcoma", traces = 32)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "embryonal rhabdomyosarcoma", traces = 33)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "ependymoma", traces = 34)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "ewing sarcoma", traces = 35)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "glioma", traces = 36)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "medulloblastoma", traces = 37)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "neuroblastoma", traces = 38)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "osteosarcoma", traces = 39)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "rhabdomyosarcoma", traces = 40)
    umap_all_ggplotly = style(umap_all_ggplotly, name = "wilms tumor", traces = 41)
    umap_all_ggplotly$x$layout$annotations[1][[1]]$text = ""
    umap_all_ggplotly = umap_all_ggplotly %>% layout(legend=list(title=list(text='<b> Tumor Disease </b>')))
    umap_all_ggplotly
    
    })
  
  output$UMAP_tumor = renderPlotly({
    umap_tumor_plotting = readRDS(file = "results/umap_plotting_tumor.RDS")
    
    umap_tumor = ggplot(umap_tumor_plotting, aes(x = x, y = y, fill = comparison_disease, label = th_sampleid)) +   
      scale_fill_manual(values=tumorColors) +
      geom_point(size = 3, pch=21) +
      guides(fill = guide_legend(override.aes = list(shape = 21) ),
             shape = guide_legend(override.aes = list(fill = "black"))) +
      xlab("UMAP1") + ylab("UMAP2") +
      labs(fill = "Disease")
    
    umap_tumor_ggplotly = ggplotly(umap_tumor)
    umap_tumor_ggplotly$x$layout$annotations[1][[1]]$text = ""
    umap_tumor_ggplotly = umap_tumor_ggplotly %>% layout(legend=list(title=list(text='<b> Tumor Disease </b>')))
    umap_tumor_ggplotly
    
  })
  
  output$UMAP_cell = renderPlotly({
    umap_cell_plotting = readRDS(file = "results/umap_plotting_cells.RDS")
    
    umap_cell = ggplot(umap_cell_plotting, aes(x = x, y = y, fill = tcga_code, label = source_sample_ID)) + 
      geom_point(size = 3) +
      scale_fill_manual(values=cellColors) +
      scale_shape_manual(values = c(21,24)) +
      xlab("UMAP1") + ylab("UMAP2") +
      labs(shape = "Blood or Solid Cancer", fill = "TCGA Code")
    
    umap_cell_ggplotly = ggplotly(umap_cell)
    umap_cell_ggplotly$x$layout$annotations[1][[1]]$text = ""
    umap_cell_ggplotly = umap_cell_ggplotly %>% layout(legend=list(title=list(text='<b> Cell Line TCGA Code </b>')))
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "ALL", traces = 1)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "BLCA", traces = 2)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "BRCA", traces = 3)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "CLL", traces = 4)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "COAD/READ", traces = 5)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "DLBC", traces = 6)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "ESCA", traces = 7)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "GBM", traces = 8)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "HNSC", traces = 9)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "KIRC", traces = 10)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LAML", traces = 11)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LCML", traces = 12)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LGG", traces = 13)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LIHC", traces = 14)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LUAD", traces = 15)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "LUSC", traces = 16)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "MB", traces = 17)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "MESO", traces = 18)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "MM", traces = 19)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "NB", traces = 20)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "OV", traces = 21)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "PAAD", traces = 22)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "PRAD", traces = 23)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "SARC", traces = 24)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "SCLC", traces = 25)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "SKCM", traces = 26)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "STAD", traces = 27)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "THCA", traces = 28)
    umap_cell_ggplotly = style(umap_cell_ggplotly, name = "UCEC", traces = 29)
    
    
    
    umap_cell_ggplotly

  })
  
  output$medullo_counts = renderPlotly({
    medullo_counts_pca = readRDS("results/medullo_counts_pca.RDS")
    medullo_sites = readRDS("results/medullo_sites.RDS")
    
    medullo_sites = medullo_sites %>% tibble::rownames_to_column("th_sampleid")
    
    medullo_counts_pca_object = ggbiplot(medullo_counts_pca, var.axes = FALSE, groups = medullo_sites$site) + 
      scale_color_manual(values=manualcolors[10:20])

    medullo_counts_pca_plotly = ggplotly(medullo_counts_pca_object)
    medullo_counts_pca_plotly$x$layout$annotations[1][[1]]$text = ""
    medullo_counts_pca_plotly = medullo_counts_pca_plotly %>% layout(legend=list(title=list(text='<b> Sample Study of Origin </b>')))
    
    
    medullo_IDs_and_PCA_coords = as.data.frame(medullo_counts_pca_object$data)
    medullo_IDs_and_PCA_coords = medullo_IDs_and_PCA_coords %>% tibble::rownames_to_column("th_sampleid")
    medullo_IDs_and_PCA_coords$xvar = as.numeric(format(round(medullo_IDs_and_PCA_coords$xvar, 7), nsmall = 7))
    medullo_IDs_and_PCA_coords$yvar = as.numeric(format(round(medullo_IDs_and_PCA_coords$yvar, 7), nsmall = 7))
    

    for(i in 1:length(medullo_counts_pca_plotly$x$data)){
      for(j in 1:length(medullo_counts_pca_plotly$x$data[[i]]$x)){
        x_coord = as.numeric(format(round(medullo_counts_pca_plotly$x$data[[i]]$x[j], 7), nsmall = 7))
        y_coord = as.numeric(format(round(medullo_counts_pca_plotly$x$data[[i]]$y[j], 7), nsmall = 7))
        
        current_sample = medullo_IDs_and_PCA_coords[medullo_IDs_and_PCA_coords$xvar == x_coord & medullo_IDs_and_PCA_coords$yvar == y_coord,]
        
        medullo_counts_pca_plotly$x$data[[i]]$text[j] = paste(medullo_counts_pca_plotly$x$data[[i]]$text[j], "<br />ID: ", current_sample$th_sampleid)
        
      }
    }
    
    medullo_counts_pca_plotly
    
  })
  
  
  output$PCAs <- renderPlotly({
    uncorrected_PCA_object <- readRDS(paste("results/", input$cancer, "/", input$cancer, "_unCorrPCA.RDS", sep = ""))
    uncorr_withDatabase = readRDS(paste("results/", input$cancer, "/", input$cancer, "_unCorr_withDatabase.RDS", sep = ""))
    uncorr_plot = ggbiplot(uncorrected_PCA_object, groups = uncorr_withDatabase$site, ellipse = TRUE, var.axes = FALSE) +
      labs(color='Study')
    uncorr_for_annotation = uncorr_plot$data
    uncorr_for_annotation = uncorr_for_annotation %>% tibble::rownames_to_column("th_sampleid")
    uncorr_for_annotation$xvar = as.numeric(format(round(uncorr_for_annotation$xvar, 7), nsmall = 7))
    uncorr_for_annotation$yvar = as.numeric(format(round(uncorr_for_annotation$yvar, 7), nsmall = 7))
    uncorr_plotly = ggplotly(uncorr_plot)
    uncorr_plotly$x$layout$annotations[1][[1]]$text = ""
    uncorr_plotly = uncorr_plotly %>% layout(legend=list(title=list(text='<b> Sample Study of Origin </b>')))
    
    
    for(i in 1:length(uncorr_plotly$x$data)){
      for(j in 1:length(uncorr_plotly$x$data[[i]]$x)){
        
        x_coord = as.numeric(format(round(uncorr_plotly$x$data[[i]]$x[j], 7), nsmall = 7))
        y_coord = as.numeric(format(round(uncorr_plotly$x$data[[i]]$y[j], 7), nsmall = 7))
        
        current_sample = uncorr_for_annotation[uncorr_for_annotation$xvar == x_coord & uncorr_for_annotation$yvar == y_coord,]
        
        uncorr_plotly$x$data[[i]]$text[j] = paste(uncorr_plotly$x$data[[i]]$text[j], "<br />ID: ", current_sample$th_sampleid)
        
      }
    }
    
    bcorrected_PCA_object <- readRDS(paste("results/", input$cancer, "/", input$cancer, "_batchCorrPCA.RDS", sep = ""))
    batchcorr_withDatabase = readRDS(paste("results/", input$cancer, "/", input$cancer, "_batchCorr_withDatabase.RDS", sep = ""))
    batchcorr_plot = ggbiplot(bcorrected_PCA_object, groups = batchcorr_withDatabase$site, ellipse = TRUE, var.axes = FALSE) +
      labs(color='Study')
    batchcorr_for_annotation = batchcorr_plot$data
    batchcorr_for_annotation = batchcorr_for_annotation %>% tibble::rownames_to_column("th_sampleid")
    batchcorr_for_annotation$xvar = as.numeric(format(round(batchcorr_for_annotation$xvar, 7), nsmall = 7))
    batchcorr_for_annotation$yvar = as.numeric(format(round(batchcorr_for_annotation$yvar, 7), nsmall = 7))
    batchcorr_plotly = ggplotly(batchcorr_plot)
    batchcorr_plotly$x$layout$annotations[1][[1]]$text = ""

    
    for(i in 1:length(batchcorr_plotly$x$data)){
      for(j in 1:length(batchcorr_plotly$x$data[[i]]$x)){
        
        x_coord = as.numeric(format(round(batchcorr_plotly$x$data[[i]]$x[j], 7), nsmall = 7))
        y_coord = as.numeric(format(round(batchcorr_plotly$x$data[[i]]$y[j], 7), nsmall = 7))
        
        current_sample = batchcorr_for_annotation[batchcorr_for_annotation$xvar == x_coord & batchcorr_for_annotation$yvar == y_coord,]
        
        batchcorr_plotly$x$data[[i]]$text[j] = paste(batchcorr_plotly$x$data[[i]]$text[j], "<br />ID: ", current_sample$th_sampleid)
        
      }
    }
    
    subplot(uncorr_plotly, batchcorr_plotly)

  })

  
  output$violins = renderPlot({
    
    batchCorrelations = readRDS(paste("results/", input$cancer, "/", input$cancer, "_correlations.RDS", sep = ""))
    
    batchCorrelations = tibble::rownames_to_column(batchCorrelations, "th_sampleid")
    
    batchCorrelations = merge(batchCorrelations, cellMetaV2, by = "th_sampleid")
    
    batchCorrelations = tibble::column_to_rownames(batchCorrelations, "th_sampleid")
    
    #pivoting correlation matrices to make them amenable for plotting
    batchPivoted = pivot_longer(batchCorrelations, !(source_sample_ID | tcga_code), names_to = "tumor_id", values_to = "cor")
    
    #finding median correlations per cell line tcga code
    batchMedians = batchPivoted %>% dplyr::group_by(tcga_code) %>% dplyr::summarise(med = median(cor)) %>% arrange(desc(med))
    
    violins_by_tcga_code = ggplot(batchPivoted, aes(x = factor(tcga_code, levels = rev(batchMedians$tcga_code)), y = cor), fill = tcga_code) + 
      geom_violin() + 
      geom_boxplot(width = 0.1) + 
      coord_flip() + 
      xlab("Cell Line TCGA Code") +
      ylab("Correlation") + 
      theme(legend.text = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 16), axis.title.x = element_text(color = "black", size = 18), axis.title.y = element_text(color = "black", size = 18)) +
      scale_fill_manual(values = cellColors)
    
    
    #finding median correlations per cell line ID
    batchMedians_lines = batchPivoted %>% dplyr::group_by(source_sample_ID) %>% dplyr::summarise(med = median(cor)) %>% arrange(desc(med))
    
    #making subset of pivot table that just includes top 10 cell lines by median correlation
    batchCellSubset = batchPivoted[batchPivoted$source_sample_ID %in% batchMedians_lines$source_sample_ID[1:10],]
    
    violin_top10 = ggplot(batchCellSubset, aes(x = factor(source_sample_ID, levels = rev(batchMedians_lines$source_sample_ID)), y = cor, fill = tcga_code)) + 
      geom_violin() + 
      geom_boxplot(width = 0.1) + 
      coord_flip() + 
      xlab("Cell Line") +
      ylab("Correlation") +
      labs(fill='Cell Line TCGA Code') + 
      theme(legend.text = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 16), axis.title.x = element_text(color = "black", size = 18), axis.title.y = element_text(color = "black", size = 18)) +
      scale_fill_manual(values = cellColors)
    
    
    grid.arrange(violins_by_tcga_code,violin_top10, ncol = 2)   
    
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "all_correlations.csv"
    },
    content = function(file) {
      write.csv(readRDS("results/masterCorrs.rds"), file, quote = F)
    }
  )
  
  output$downloadData_cancer <- downloadHandler(
    filename = function() {
      paste(input$cancer, "_correlations.csv", sep = "")
    },
    content = function(file) {
      write.csv(readRDS(paste("results/", input$cancer,"/",input$cancer, "_correlations.RDS", sep = "")), file, quote = F)
    }
    
  )
  
  output$downloadData_cancer_matched <- downloadHandler(
    filename = function() {
      paste(input$cancer, "_correlations_matched.csv", sep = "")
    },
    content = function(file) {
      write.csv(readRDS(paste("results/", input$cancer,"/",input$cancer, "_correlations_matched.RDS", sep = "")), file, quote = F)
    }
    
  )
  
  output$downloadData_medullo_volcano <- downloadHandler(
    filename = "medullo_volcano_data.csv",
    content = function(file) {
      write.csv(readRDS(file = "results/medullo_evolcano.RDS"), file, quote = F)
    }
    
  )
  
  output$downloadData_cancer_corrs <- downloadHandler(
    filename = paste(input$cancer,"_correlations_by_cell_line.csv", sep = ""),
    content = function(file) {
      write.csv(readRDS(paste("results/", input$cancer, "/", input$cancer, "_all_cell_lines_corrs_table.RDS", sep = "")), file, quote = F)
    }
    
  )
  
  
  
  output$graphical_abstract <- renderImage({
    list(src = "Sirota Visual Methods.png",
         alt = "Graphical abstract",
         width = "auto",
         height = "650%"
    )
  }, deleteFile = FALSE)
  
  output$matched_heatmap = renderImage({
    
    list(src = paste("results/",input$cancer,"/",input$cancer,"_heatmap_matched_cell_lines_only.png", sep = ""),
         alt = "matched heatmap",
         height = 2000,
         width = 2000,
         contentType = "image/png")
    
    
  }, deleteFile = FALSE)
  
  output$matched_heatmap_plot = renderPlot({

    matchedHeatmapDF = readRDS(paste("results/", input$cancer, "/", input$cancer, "_correlations_matched.RDS", sep = ""))
    
    #replace cell line sample IDs with actual names
    matchedHeatmapDF = tibble::rownames_to_column(matchedHeatmapDF, "th_sampleid")
    matchedHeatmapDF = tibble::column_to_rownames(matchedHeatmapDF, "source_sample_ID")
    matchedHeatmapDF = dplyr::select(matchedHeatmapDF, -th_sampleid, -tcga_code)
    
    breaksList = seq(0, 1, by = 0.01) 
    plot = pheatmap(matchedHeatmapDF, show_colnames = F, color = colorRampPalette(c("white", "purple"))(length(breaksList)), breaks = breaksList)
    plot
    dev.off()
    plot
  })
  
  output$evolcano = renderPlot({
    volcano_DF = readRDS(file = "results/medullo_evolcano.RDS")
    
    evolcano = EnhancedVolcano(volcano_DF, x = "logFC", y = "p.adj", lab = volcano_DF$hgnc_symbol, labSize = 5, title = "WNT-like vs. SHH/Group 3-like MB Samples")
    
    evolcano
    
  })
  
  output$evolcano_table = renderDT({
    volcano_table = readRDS(file = "results/medullo_evolcano.RDS")
    
    datatable(volcano_table) %>% formatSignif(columns = c(3,4,5,6),3)
    
    
  })
  
  
  output$table_all_corrs <- renderDT({
    
    table_all_DF = readRDS(paste("results/", input$cancer, "/", input$cancer, "_all_cell_lines_corrs_table.RDS", sep = ""))
    
    datatable(table_all_DF, colnames = c("Cell Line", "Median Correlation", "TCGA Code")) %>% formatRound(columns = "med", 3)

  })
  
  output$by_line_medullo = renderPlot({
    medullo_pivot = readRDS("results/medullo_pivot.RDS")
    
    bulk_plot = ggplot(medullo_pivot, aes(x = source_sample_ID, y = cor, fill = PC)) + 
      geom_violin(position = position_dodge(width = 0.9)) + 
      geom_boxplot(width = 0.25, position = position_dodge(0.9)) + 
      geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9, jitter.width = 0.25)) + 
      stat_compare_means(method = "wilcox.test", size = 8, label = "p.format", label.y = 0.68) + 
      stat_summary(fun=mean, geom = "crossbar", position = position_dodge(0.9)) +
      theme(legend.text = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 16), axis.text.x = element_text(color = "black", size = 16), axis.title.x = element_text(color = "black", size = 18), axis.title.y = element_text(color = "black", size = 18))
      
    
    bulk_plot
    
  })
  
  output$bulk_medullo = renderPlot({
    medullo_pivot = readRDS("results/medullo_pivot.RDS")
    
    line_plot = ggplot(medullo_pivot, aes(x = PC, y = cor, fill = PC)) + 
      geom_violin() + 
      geom_boxplot(width = 0.25) + 
      geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9, jitter.width = 0.25)) +
      stat_compare_means(method = "wilcox.test", label.x.npc = "center", size = 8) + 
      stat_summary(fun=mean, geom = "crossbar") +      
      theme(legend.text = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16), axis.title.x = element_text(color = "black", size = 18), axis.title.y = element_text(color = "black", size = 18))

    
    line_plot
    
  })
  
  output$violinMaster = renderPlot({
    violinMaster = readRDS(file = "results/violinMaster_fig1.RDS")
    violinMedians = readRDS(file = "results/violinMedians_fig1.RDS")
    
    plot = ggplot(violinMaster, aes(x = factor(comparison_disease, levels = violinMedians$comparison_disease), y = cor, fill = comparison_disease)) + 
      geom_violin() + geom_boxplot(width = 0.1) + 
      ggtitle(paste("Matched Cell Line - Tumor Samples")) + 
      xlab("Tumor Type") +
      ylab("Correlation") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 14), legend.text = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 16), axis.title.x = element_text(color = "black", size = 18), axis.title.y = element_text(color = "black", size = 18)) + 
      scale_fill_manual(values = tumorColors)
    
    plot
    
    
  })
  
  output$master_heatmap = renderPlot({
    
    masterCorrs = readRDS(file = "results/masterCorrs_fig1_heatmap.RDS")
    
    breaksList = seq(min(masterCorrs), max(masterCorrs), by = 0.01) 
    
    plot = pheatmap(masterCorrs, fontsize = 12, color = colorRampPalette(c("white", "purple"))(length(breaksList)), breaks = breaksList) 
    
    plot
    dev.off()
    plot
    
  })
  

  
}

# Run the application 
shinyApp(ui = ui, server = server)

