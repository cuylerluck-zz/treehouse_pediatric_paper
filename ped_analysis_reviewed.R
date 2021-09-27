library(data.table)
library(pheatmap)
library(ggplot2)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggbiplot)
library(gmodels)
library(umap)
library(corrplot)
library(grDevices)
library(sva)
library(ggpubr)
library(gridExtra)
library(umap)
library(ggrepel)
library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler)
library(biomaRt)
library(stringr)
library(org.Hs.eg.db)

#ensure graphics device is cleared or else plots may fail
dev.off()

#Set working directory as needed
setwd("/Volumes/Cuyler_T5/UCSF Sirota Lab Rotation")

#Set seed for reproducibility, I arbitrarily pulled one from .Random.seed
set.seed(-1535835272)

#Tumor Compendium Metadata -- select diseases (>=30 samples) and age at dx <= 18 with comparison_disease that I defined
tumorMeta = fread(file = "Intermediate_Data/tCompendiumv11_pediatric_selectDiseases_forComparison_meta.csv", header = TRUE)
tumorMeta = dplyr::select(tumorMeta, th_sampleid, comparison_disease, Database)

#Cell Compendium Metadata v2 -- cell lines from Treehouse compendium with TCGA acronyms taken from CCLE annotation file
#Excluding those with no acronym or unclassified acronym
cellMetaV2 = fread(file = "Intermediate_Data/cellMetaV2.csv", header = TRUE)
cellMetaV2 = dplyr::select(cellMetaV2, th_sampleid, source_sample_ID, tcga_code)

#Select Pediatric Tumor TPM Data
pedData = fread("Intermediate_Data/tCompendium_v11_selectPedSamples_TPM.csv")

#Select Cell Line TPM Data
cellData = fread("Intermediate_Data/cellData_TPM.csv")



###analysis
#requires reading in pedData, tumorMeta, cellData, and cellMetaV2







###Generating UMAPs for uncorrected samples in Figure 1
#color palette derived from Ferroao on stack overflow: https://stackoverflow.com/a/46810812
manualcolors<-c('forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "yellow1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro", 'black')

tumorColors = manualcolors[1:14]
names(tumorColors) = c(unique(tumorMeta$comparison_disease), "medulloblastoma_SHH_GRP3_like", "medulloblastoma_WNT_like")
tumorColors = tumorColors[order(names(tumorColors))]

cellColors = manualcolors[1:29]
names(cellColors) = unique(cellMetaV2$tcga_code)
cellColors = cellColors[order(names(cellColors))]


options(ggrepel.max.overlaps = Inf)


#merging all samples together by gene
umap_merge = inner_join(pedData, cellData, by = "Gene")
umap_merge = tibble::column_to_rownames(umap_merge, "Gene")
umap_merge = as.data.frame(t(umap_merge))

#collapsing to just the top 5000 most variable genes across all samples (by IQR) for computational speed
umapIQRs = apply(umap_merge, 2, IQR)
ordered_umapIQRs = umapIQRs[order(umapIQRs, decreasing = T)]
top_umapGenes = names(ordered_umapIQRs)[1:5000]
umap_top = umap_merge[,top_umapGenes]

#generating UMAP object and plotting dataframe
umap_object = umap(umap_top)
umap_plotting = data.frame(x = umap_object$layout[,1], y = umap_object$layout[,2])
umap_plotting = tibble::rownames_to_column(umap_plotting, "th_sampleid")

#adding information on whether samples are cell lines or tumor samples
umap_plotting$type = NA
umap_plotting[umap_plotting$th_sampleid %in% cellMetaV2$th_sampleid, ]$type = "cell_line"
umap_plotting[umap_plotting$th_sampleid %in% tumorMeta$th_sampleid, ]$type = "tumor_sample"

#adding disease data if I want to use it to plot with later
umap_plotting$disease = NA
for(ID in umap_plotting$th_sampleid){
  if(ID %in% cellMetaV2$th_sampleid){
    umap_plotting[umap_plotting$th_sampleid == ID,]$disease = cellMetaV2[cellMetaV2$th_sampleid == ID,]$tcga_code
  }
  else{
    umap_plotting[umap_plotting$th_sampleid == ID,]$disease = tumorMeta[tumorMeta$th_sampleid == ID,]$comparison_disease
  }
}

#picking tumor samples to label
picked_samples = c()
for(value in unique(umap_plotting[umap_plotting$th_sampleid %in% tumorMeta$th_sampleid,]$disease)){
  IDs = c(umap_plotting[umap_plotting$disease == value & umap_plotting$th_sampleid %in% tumorMeta$th_sampleid,]$th_sampleid)
  picked_samples = c(picked_samples, IDs[20])
}

umap_plotting = umap_plotting %>% dplyr::mutate(label = ifelse(th_sampleid %in% picked_samples, "yes","no"))

#saving RDS object and plotting combined UMAP
saveRDS(umap_plotting, file = "results/umap_plotting_all_samples.RDS")


pdf(paste("results/UMAP_all_samples_precorrection_tumor_sample_colors.pdf", sep = ""))
print(
  ggplot(umap_plotting, aes(x = x, y = y, shape = type, label = disease, fill = disease)) + 
    geom_point(data = umap_plotting[umap_plotting$type == "cell_line",], size = 3, alpha = 0.5) +
    geom_point(data = umap_plotting[umap_plotting$type == "tumor_sample",],size = 3) +
    scale_fill_manual(values = tumorColors) +
    scale_shape_manual(values = c(21,24)) +
    geom_label_repel(data = umap_plotting[umap_plotting$label == "yes",], show.legend = F, nudge_y = 2, fontface = "bold") +
    guides(fill = guide_legend(override.aes = list(shape = 21) ),
           shape = guide_legend(override.aes = list(fill = "gray"))) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(shape = "Tumor Sample or Cell Line", fill = "Disease")
)
dev.off()


blood_cancers = c("ALL", "MM", "DLBC", "LAML", "LCML", "CLL")


#UMAP cell lines only
cellData_forUMAP = cellData
cellData_forUMAP = tibble::column_to_rownames(cellData_forUMAP, "Gene")
cellData_forUMAP = as.data.frame(t(cellData_forUMAP))

#collapsing to just the top 5000 most variable genes across all cell lines (by IQR) for computational speed
cell_umapIQRs = apply(cellData_forUMAP, 2, IQR)
cell_ordered_umapIQRs = cell_umapIQRs[order(cell_umapIQRs, decreasing = T)]
cell_top_umapGenes = names(cell_ordered_umapIQRs)[1:5000]
cell_umap_top = cellData_forUMAP[,cell_top_umapGenes]

#generating UMAP object and plotting dataframe
cell_umap_object = umap(cell_umap_top)
cell_umap_plotting = data.frame(x = cell_umap_object$layout[,1], y = cell_umap_object$layout[,2])
cell_umap_plotting = tibble::rownames_to_column(cell_umap_plotting, "th_sampleid")

#attaching disease info
cell_umap_plotting = left_join(cell_umap_plotting, cellMetaV2, by = "th_sampleid")

cell_umap_plotting = dplyr::mutate(cell_umap_plotting, blood_solid = ifelse(cell_umap_plotting$tcga_code %in% blood_cancers, "blood", "solid"))


#pick out one cell line each for different disease types to label (arbitrarily chooses second line in each type)
picked_lines = c()
for(value in unique(cell_umap_plotting$tcga_code)){
  IDs = c(cell_umap_plotting[cell_umap_plotting$tcga_code == value,]$th_sampleid)
  picked_lines = c(picked_lines, IDs[2])
}

#annotate which lines are the winners to be labeled
cell_umap_plotting = dplyr::mutate(cell_umap_plotting, label = ifelse(cell_umap_plotting$th_sampleid %in% picked_lines, "yes", "no"))

#saving RDS object for web app
saveRDS(cell_umap_plotting, file = "results/umap_plotting_cells.RDS")

#plot cell line only PCA
pdf("results/UMAP_cell_lines_only.pdf", width = 8, height = 6)
print(
ggplot(cell_umap_plotting, aes(x = x, y = y, fill = tcga_code, shape = blood_solid, label = tcga_code)) + 
  geom_point(size = 3) +
  scale_fill_manual(values=cellColors) +
  scale_shape_manual(values = c(21,24)) +
  geom_label_repel(data = cell_umap_plotting[cell_umap_plotting$label == "yes",], show.legend = F, nudge_x = 2, nudge_y = 1, fontface = "bold") +
  guides(fill = guide_legend(override.aes = list(shape = 21) ),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("UMAP1") + ylab("UMAP2") +
  labs(shape = "Blood or Solid Cancer", fill = "TCGA Code")
)
dev.off()


#UMAP Tumor Samples only

pedData_forUMAP = pedData
pedData_forUMAP = tibble::column_to_rownames(pedData_forUMAP, "Gene")
pedData_forUMAP = as.data.frame(t(pedData_forUMAP))

#collapsing to just the top 5000 most variable genes across all tumor samples (by IQR) for computational speed
ped_umapIQRs = apply(pedData_forUMAP, 2, IQR)
ped_ordered_umapIQRs = ped_umapIQRs[order(ped_umapIQRs, decreasing = T)]
ped_top_umapGenes = names(ped_ordered_umapIQRs)[1:5000]
ped_umap_top = pedData_forUMAP[,ped_top_umapGenes]


#generating UMAP object and plotting dataframe
ped_umap_object = umap(ped_umap_top)
ped_umap_plotting = data.frame(x = ped_umap_object$layout[,1], y = ped_umap_object$layout[,2])
ped_umap_plotting = tibble::rownames_to_column(ped_umap_plotting, "th_sampleid")

#attaching disease info
ped_umap_plotting = left_join(ped_umap_plotting, tumorMeta, by = "th_sampleid")

#pick out one tumor sample each for different disease types to label (chooses 20th sample in each type because it labels
#groups best
picked_samples = c()
for(value in unique(ped_umap_plotting$comparison_disease)){
  IDs = c(ped_umap_plotting[ped_umap_plotting$comparison_disease == value,]$th_sampleid)
  picked_samples = c(picked_samples, IDs[20])
}

#annotate which samples are the winners to be labeled
ped_umap_plotting = dplyr::mutate(ped_umap_plotting, label = ifelse(ped_umap_plotting$th_sampleid %in% picked_samples, "yes", "no"))

#saving RDS object for web app
saveRDS(ped_umap_plotting, file = "results/umap_plotting_tumor.RDS")

#plot tumor sample only PCA
pdf("results/UMAP_tumor_samples_only.pdf", width = 8, height = 6)
print(
  ggplot(ped_umap_plotting, aes(x = x, y = y, fill = comparison_disease, label = comparison_disease)) +   
    scale_fill_manual(values=tumorColors) +
    geom_point(size = 3, pch=21) +
    geom_label_repel(data = ped_umap_plotting[ped_umap_plotting$label == "yes",], show.legend = F, nudge_x = 2, nudge_y = 1, fontface = "bold") +
    guides(fill = guide_legend(override.aes = list(shape = 21) ),
           shape = guide_legend(override.aes = list(fill = "black"))) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(fill = "Disease")
)
dev.off()


#pulling out characters before _ or - which means pulling out site ID from th_sampleid
tumorSites = mutate(tumorMeta, site = stringr::str_match(tumorMeta$th_sampleid, ".[^_-]+"))
#View(unique(sort(tumorSites$site))) #checking to make sure the syntax is OK and turned up what I want

###
#medulloblastoma expected counts investigation to help separate medulloblastoma samples before correlation analysis

###do not need to run if the csv at the end exists already -- just pulling medullo expected counts data out of all data
#load tumor sample expected counts data and extract only these medulloblastoma samples
#tumor_expected_counts = fread("Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/TumorCompendium_v11_PolyA_ensembl_expected_count_58581genes_2020-04-09.tsv")
#medullo_samples = tumorMeta$th_sampleid[tumorMeta$comparison_disease == "medulloblastoma"]
#medullo_expected_counts = dplyr::select(tumor_expected_counts, c("Gene",all_of(medullo_samples)))
#write.csv(medullo_expected_counts,file = "Intermediate_Data/medullo_expected_counts.csv")
###


#data cleaning to read in medullo counts and format properly
medullo_expected_counts = fread("Intermediate_Data/medullo_expected_counts.csv")
medullo_expected_counts = dplyr::select(medullo_expected_counts, -V1)

medullo_expected_counts = tibble::column_to_rownames(medullo_expected_counts, "Gene")
medullo_expected_counts = as.data.frame(t(medullo_expected_counts))

#rounding counts in each column (for each gene)
medullo_expected_counts = medullo_expected_counts %>% dplyr::mutate(across(everything(),round))

#remove genes with zero variance across all samples
medullo_vars = apply(medullo_expected_counts, 2, var)
ordered_medullo_vars = medullo_vars[order(medullo_vars, decreasing = T)]
non_zero_var_medullo_genes = names(ordered_medullo_vars[ordered_medullo_vars>0])
medullo_expected_counts = dplyr::select(medullo_expected_counts, all_of(non_zero_var_medullo_genes))

#log transforming
medullo_expected_counts_log2 = log2(medullo_expected_counts + 1)

#select top 5000 most variable genes to use for a PCA
medullo_IQRs = apply(medullo_expected_counts_log2, 2, IQR)
ordered_medullo_IQRs = medullo_IQRs[order(medullo_IQRs, decreasing = T)]
top_medullo_genes = names(ordered_medullo_IQRs)[1:5000]
medullo_top = medullo_expected_counts_log2[,top_medullo_genes]

#make pca object
medullo_counts_pca = fast.prcomp(medullo_top, center = TRUE)

#get and order study site information for labelling pca plot
medullo_sites = tumorSites[tumorSites$th_sampleid %in% rownames(medullo_expected_counts_log2),]
medullo_sites = tibble::column_to_rownames(medullo_sites, "th_sampleid")
medullo_sites = medullo_sites[rownames(medullo_top),]

saveRDS(medullo_counts_pca, file = "results/medullo_counts_pca.RDS")
saveRDS(medullo_sites, file = "results/medullo_sites.RDS")

#plot pca
pdf(paste("results/medullo_pca.pdf"), width = 4, height = 4)
print(
  ggbiplot(medullo_counts_pca, var.axes = FALSE, groups = medullo_sites$site) + scale_color_manual(values=manualcolors[10:20])
)
dev.off()

#k-means analysis of PCA object to objectively prove the two groupings are a good idea
pca_df = data.frame(medullo_counts_pca$x[,1:2])
k = kmeans(pca_df, centers = 2)
pdf(paste("results/medullo_kmeans.pdf"), width = 4, height = 4)
print(
ggplot(pca_df, aes(x = PC1, y = PC2, color = k$cluster)) + geom_point(show.legend = F)
)
dev.off()
#divide MB samples into two groups based on the sign of standardized PC1 (i.e. left or right clusters in the PCA)
#the signs of original PC components should correspond to scaled PCs on the PCA, so can group this way
#even though I don't know the exact scaling method

#use this to demonstrate that the scaled PCs are just a narrower version of the original PCs, and the signs are the same
#ggplot(data = as.data.frame(medullo_counts_pca$x), aes(x=PC1, y=PC2)) + geom_point()
unscaled_PCs = as.data.frame(medullo_counts_pca$x)
neg_PC1 = rownames(unscaled_PCs[unscaled_PCs$PC1 < 0,])
pos_PC1 = rownames(unscaled_PCs[unscaled_PCs$PC1 > 0,])

#use edgeR to do differential expression analysis between samples with positive vs. negative PC1s
#using original counts (not log transformed) for this. log transformation was just for PCA visualization.

#transpose expected counts matrix to match edgeR syntax
medullo_expected_counts = as.data.frame(t(medullo_expected_counts))

group = colnames(medullo_expected_counts) %in% neg_PC1 #in the same order as samples in medullo_expected_counts
#gives TRUE if it's a sample with NEG PC1, FALSE otherwise
group = ifelse(group == "TRUE", "NEG", "POS") #converting T/F to named NEG and POS groupings
master = DGEList(counts = medullo_expected_counts, group = group)
keep = filterByExpr(master)
master = master[keep, ,keep.lib.sizes = F]
master = calcNormFactors(master)
design = model.matrix(~group)
master = estimateDisp(master, design)

et = exactTest(master)
exactResults = et$table
exactResults_FDR = dplyr::mutate(exactResults, p.adj = p.adjust(PValue, method = "BH"))

#converting gene names to hugo symbols for better EnhancedVolcano and gseGO
gene_names = rownames(exactResults_FDR)
trimmed_gene_names = stringr::str_extract(gene_names, "[^\\.]*")
mart = biomaRt::useMart('ensembl', "hsapiens_gene_ensembl")
medullo_key = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values = trimmed_gene_names, mart = mart)
medullo_key = medullo_key[!medullo_key$hgnc_symbol == "",]

medullo_key$hgnc_symbol[duplicated(medullo_key$hgnc_symbol)]
#19893 genes have matches but four symbols are duplicated -- RMRP, SNORA17B, LINC01238, POLR2J4
#removing these
medullo_key = medullo_key[!medullo_key$hgnc_symbol %in% medullo_key$hgnc_symbol[duplicated(medullo_key$hgnc_symbol)],]

medullo_key$ensembl_gene_id[duplicated(medullo_key$ensembl_gene_id)]
#three ensembl IDs duplicated -- "ENSG00000254876" "ENSG00000276085" "ENSG00000230417"
#removing these
medullo_key = medullo_key[!medullo_key$ensembl_gene_id %in% medullo_key$ensembl_gene_id[duplicated(medullo_key$ensembl_gene_id)],]

#check for duplicates of either symbol type again
medullo_key$hgnc_symbol[duplicated(medullo_key$hgnc_symbol)]
medullo_key$ensembl_gene_id[duplicated(medullo_key$ensembl_gene_id)]

#seeing no more duplicates, trimming gene symbols in exactResults_FDR and adding name information for genes that have it
#available, which is 19879 genes

exactResults_FDR = tibble::rownames_to_column(exactResults_FDR, "ensembl_gene_id")
exactResults_FDR = dplyr::mutate(exactResults_FDR, ensembl_gene_id = stringr::str_extract(ensembl_gene_id, "[^\\.]*"))
exactResults_FDR = left_join(medullo_key, exactResults_FDR, by = "ensembl_gene_id")

#saving exactResults_FDR to RDS so I can use it in the web app
saveRDS(exactResults_FDR, file = "results/medullo_evolcano.RDS")

pdf(paste("results/medullo_enhancedvolcano.pdf"), width = 6.75)
print(
  EnhancedVolcano(exactResults_FDR, x = "logFC", y = "p.adj", lab = exactResults_FDR$hgnc_symbol)
)
dev.off()

FCs = dplyr::select(exactResults_FDR, hgnc_symbol, logFC)
FCvector = FCs$logFC
names(FCvector) = FCs$hgnc_symbol
FCvector = sort(FCvector, decreasing = T)

FCs = FCs[order(-FCs$logFC),]
write.table(FCs, file = "Intermediate_Data/medullo_DE_fold_changes.txt", quote = F, col.names = F, row.names = F, sep = "\t")

#to use in GSEA preranked, change the file extension to .rnk in your file browser. Then it can be loaded into the application.
#non-default settings for GSEA preranked:
#Collapse/Remap to gene symols: No_Collapse
#Min Size = 10

#the paper I got my custom gene lists from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112909

#Using gseGO to do an unbiased GSEA 
gse = gseGO(geneList=FCvector, keyType = "SYMBOL", pAdjustMethod = "BH", OrgDb = org.Hs.eg.db, seed = T)

gseUp = gse
gseUp@result = gseUp@result[gseUp@result$NES>0,] #subsetting results with positive normalized enrichment score
#this makes emapplots nicer (split by direction of fold change)

gseDown = gse
gseDown@result = gseDown@result[gseDown@result$NES<0,] #same with negative NES


e1 = emapplot(gseUp, showCategory = 20)
e2 = emapplot(gseDown, showCategory = 20)

pdf(paste("results/medullo_emapplot.pdf"), width = 12, height = 5)
print(
  cowplot::plot_grid(e2, e1, labels = c("Negative NES", "Positive NES"))
)
dev.off()

#using WNT16 (highly negative logFC of -4) to determine which group is positive vs negative FC
WNT16_values = medullo_expected_counts["ENSG00000002745",]
rownames(WNT16_values) = "WNT16"
WNT16_values = as.data.frame(t(WNT16_values))
WNT16_values = tibble::rownames_to_column(WNT16_values, "sample")
WNT16_values = mutate(WNT16_values, PC1_Sign = ifelse(sample %in% pos_PC1, "POS", "NEG"))
WNT16_values = mutate(WNT16_values, WNT16 = log2(WNT16+1)) #TRANSFORMED WNT16 COUNTS TO LOG2(X+1) FOR VISUALIZATION

pdf(paste("results/medullo_WNT16.pdf"))
print(
  ggplot(data = WNT16_values, aes(x = PC1_Sign, y = WNT16)) + 
    geom_violin(alpha = 0.3, fill = "blue") + 
    geom_point() + 
    ylab("WNT16 log2(count+1)")
)
dev.off()

#double checking using PTCH2 (positive logFC of 2.6)
PTCH2_values = medullo_expected_counts["ENSG00000117425",]
rownames(PTCH2_values) = "PTCH2"
PTCH2_values = as.data.frame(t(PTCH2_values))
PTCH2_values = tibble::rownames_to_column(PTCH2_values, "sample")
PTCH2_values = mutate(PTCH2_values, PC1_Sign = ifelse(sample %in% pos_PC1, "POS", "NEG"))
PTCH2_values = mutate(PTCH2_values, PTCH2 = log2(PTCH2+1)) #TRANSFORMED PTCH2 COUNTS TO LOG2(X+1) FOR VISUALIZATION

pdf(paste("results/medullo_PTCH2.pdf"))
print(
  ggplot(data = PTCH2_values, aes(x = PC1_Sign, y = PTCH2)) + 
    geom_violin(alpha = 0.3, fill = "red") + 
    geom_point() + 
    ylab("PTCH2 log2(count+1)")
)
dev.off()

#these plots indicate that positive logFC = higher in samples with a positive PC1 = gene further left on GSEA barcode plot


#separate out medullo metadata into wnt-like and shh/group3-like in the tumor meta data

tumorMeta$comparison_disease[tumorMeta$th_sampleid %in% pos_PC1] = "medulloblastoma_SHH_GRP3_like"
tumorMeta$comparison_disease[tumorMeta$th_sampleid %in% neg_PC1] = "medulloblastoma_WNT_like"

###


pedData = tibble::column_to_rownames(pedData, "Gene")

#transposing expression matrix
pedData = as.data.frame(t(pedData))

#removing genes from the cell data with zero variance 
cellData = tibble::column_to_rownames(cellData, "Gene")
cellData = cellData[apply(cellData, 1, var, na.rm=TRUE) !=0,]
cellData = tibble::rownames_to_column(cellData, "Gene")


diseases = unique(tumorMeta$comparison_disease)
diseases = c(diseases, "medulloblastoma")

masterCorrs = data.frame(cellIDs = cellMetaV2$th_sampleid)

#defining function to manually choose cell line TCGA code to pair with each tumor type
getMatch = function(comparison_disease){
  if(comparison_disease %in% c("medulloblastoma_SHH_GRP3_like", "medulloblastoma_WNT_like")){return("MB")}
  else if(comparison_disease == "medulloblastoma"){return("MB")}
  else if(comparison_disease == "ewing sarcoma"){return("SARC")}
  else if(comparison_disease == "acute lymphoblastic leukemia"){return("ALL")}
  else if(comparison_disease == "alveolar rhabdomyosarcoma"){return("SARC")}
  else if(comparison_disease == "osteosarcoma"){return("SARC")}
  else if(comparison_disease == "embryonal rhabdomyosarcoma"){return("SARC")}
  else if(comparison_disease == "glioma"){return(c("LGG","GBM"))} #Glioma uses LGG and GBM cell lines
  else if(comparison_disease == "acute myeloid leukemia"){return("LAML")}
  else if(comparison_disease == "wilms tumor"){return("KIRC")} #WT uses KIRC -- no great match
  else if(comparison_disease == "neuroblastoma"){return("NB")}
  else if(comparison_disease == "ependymoma"){return(c("LGG"))}
  else if(comparison_disease == "rhabdomyosarcoma"){return("SARC")}
}

for(disease in diseases){

  
  
  
  if(disease == "medulloblastoma"){
    disease_samples = tumorMeta$th_sampleid[tumorMeta$comparison_disease == "medulloblastoma_SHH_GRP3_like" | tumorMeta$comparison_disease == "medulloblastoma_WNT_like"]
  } else{
    disease_samples = tumorMeta$th_sampleid[tumorMeta$comparison_disease == disease]
  }
  
  subset = pedData[disease_samples, ]
  
  #transposing subset for ComBat syntax
  subset = as.data.frame(t(subset))
  
  subset_meta = tumorSites[tumorSites$th_sampleid %in% disease_samples,]
  
  #ordering subset_meta samples the same as subset's samples
  subset_meta = tibble::column_to_rownames(subset_meta, "th_sampleid")
  subset_meta = subset_meta[colnames(subset),]
  
  #batch correction step
  subset_batchCorr = ComBat(subset, subset_meta$site)
  subset_batchCorr = as.data.frame(subset_batchCorr)
  
  #making a copy of un-corrected data for PCA generation
  subset_unCorr = subset
  
  #transposing corrected and uncorrected subsets for next steps
  subset_batchCorr = as.data.frame(t(subset_batchCorr))
  subset_unCorr = as.data.frame(t(subset_unCorr))
  
  #removing genes with zero variance from both datasets
  subset_batchCorr_rmNoVar = subset_batchCorr[,apply(subset_batchCorr, 2, var, na.rm=TRUE) !=0]
  
  subset_unCorr_rmNoVar = subset_unCorr[,apply(subset_unCorr, 2, var, na.rm=TRUE) !=0]
  
  #choosing the top 5000 most variable genes in each dataset
  
  geneIQRs = apply(subset_batchCorr_rmNoVar, 2, IQR)
  orderedIQRs = geneIQRs[order(geneIQRs, decreasing = T)]
  topGenes = names(orderedIQRs)[1:5000]
  
  topBatchCorr = subset_batchCorr_rmNoVar[, topGenes]
  
  
  geneIQRs_unCorr = apply(subset_unCorr_rmNoVar, 2, IQR)
  orderedIQRs_unCorr = geneIQRs_unCorr[order(geneIQRs_unCorr, decreasing = T)]
  topGenes_unCorr = names(orderedIQRs_unCorr)[1:5000]
  
  topUnCorr = subset_unCorr_rmNoVar[, topGenes_unCorr]
  
  
  #making PCA objects
  batchCorr_pca = fast.prcomp(topBatchCorr, center = TRUE)
  unCorr_pca = fast.prcomp(topUnCorr, center = TRUE)
  
  #adding metadata for grouping in PCA plots
  batchCorr_withDatabase = tibble::rownames_to_column(topBatchCorr, "th_sampleid")
  batchCorr_withDatabase = left_join(batchCorr_withDatabase, tumorSites, by = "th_sampleid")
  
  unCorr_withDatabase = tibble::rownames_to_column(topUnCorr, "th_sampleid")
  unCorr_withDatabase = left_join(unCorr_withDatabase, tumorSites, by = "th_sampleid")
  
  
  #PCA plot grouping by database (to check for institutional batch effects)
  ifelse(!dir.exists(file.path(paste("results/", disease, sep = ""))),
         dir.create(file.path(paste("results/", disease, sep = ""))), FALSE)
  
  pdf(paste("results/", disease, "/", disease, "_pca_top5000Genes_batchCorrected.pdf", sep=""), width = 4, height = 4)
  print(ggbiplot(batchCorr_pca, groups = batchCorr_withDatabase$site, ellipse = TRUE, var.axes = FALSE) +
          ggtitle(paste(disease,"b corr")) +
          labs(color='Study')) 
  dev.off()
  
  pdf(paste("results/", disease, "/", disease, "_pca_top5000Genes_unCorr.pdf", sep=""), width = 4, height = 4)
  print(ggbiplot(unCorr_pca, groups = unCorr_withDatabase$site, ellipse = TRUE, var.axes = FALSE) +
          ggtitle(paste(disease,"uncorr")) +
          labs(color='Study'))
          
  dev.off()
  
  #also saving PCA objects and metadata objects for web app
  saveRDS(batchCorr_pca, file = paste("results/", disease, "/", disease, "_batchCorrPCA.RDS", sep=""))
  saveRDS(unCorr_pca, file = paste("results/", disease, "/", disease, "_unCorrPCA.RDS", sep=""))
  saveRDS(batchCorr_withDatabase, file = paste("results/", disease, "/", disease, "_batchCorr_withDatabase.RDS", sep=""))
  saveRDS(unCorr_withDatabase, file = paste("results/", disease, "/", disease, "_unCorr_withDatabase.RDS", sep=""))
  
  
  
  #generating correlations between disease-subset tumor samples & all cell lines
  #using the batch corrected tumor data from above, before I subsetted genes -- need to repick them here
  
  #transposing tumor data to get it in the proper orientation to join with cellData
  
  subset_batchCorr_rmNoVar = as.data.frame(t(subset_batchCorr_rmNoVar))

  #putting row names into gene column for merging
  subset_batchCorr_rmNoVar = tibble::rownames_to_column(subset_batchCorr_rmNoVar, "Gene")

  #inner joining the cell data onto the tumor data so that the only genes are those with variance in both datasets
  batchMerged = inner_join(subset_batchCorr_rmNoVar, cellData, by = "Gene")

  #putting Gene back on rownames
  batchMerged = tibble::column_to_rownames(batchMerged, "Gene")

  #subsetting to the top 5000 most variable genes across disease-specific tumor samples
  #potentially different from similar calculation above because the genes remaining here is the intersect of those
  #that have non-zero variance in the individual tumor sample & cell line datasets
  batchMergedIQRs = apply(batchMerged[, disease_samples], 1, IQR)
  orderedBatchMergedIQRs = batchMergedIQRs[order(batchMergedIQRs, decreasing = T)]
  topBatchMergedGenes = names(orderedBatchMergedIQRs)[1:5000]
  batchMerged = batchMerged[topBatchMergedGenes,]
  
  #generating correlation matrix
  batchCorrelations = as.data.frame(cor(batchMerged, method = "s"))

  #subsetting rows to cell lines and columns to tumor samples
  batchCorrelations = batchCorrelations[cellMetaV2$th_sampleid, disease_samples]

  #add to master correlation dataframe
  if(!disease == "medulloblastoma"){
  toAdd = batchCorrelations
  toAdd = tibble::rownames_to_column(toAdd, "cellIDs")
  masterCorrs = inner_join(masterCorrs, toAdd, by = "cellIDs")
  }
  
  #save as RDS for web tool
  saveRDS(batchCorrelations, file = paste("results/",disease,"/",disease,"_correlations.RDS", sep = ""))
  
  #attaching cell meta data
  batchCorrelations = tibble::rownames_to_column(batchCorrelations, "th_sampleid")

  batchCorrelations = merge(batchCorrelations, cellMetaV2, by = "th_sampleid")

  batchCorrelations = tibble::column_to_rownames(batchCorrelations, "th_sampleid")

  #pivoting correlation matrices to make them amenable for plotting
  batchPivoted = pivot_longer(batchCorrelations, !(source_sample_ID | tcga_code), names_to = "tumor_id", values_to = "cor")

  #finding median correlations per cell line tcga code
  batchMedians = batchPivoted %>% dplyr::group_by(tcga_code) %>% dplyr::summarise(med = median(cor)) %>% arrange(desc(med))

  
  #plotting violin plot with x axis as cell line tcga code
  pdf(paste("results/", disease, "/", disease, "_violin_by_cell_line_tcgacode_batchCorrected.pdf", sep=""))
  print(
    ggplot(batchPivoted, aes(x = factor(tcga_code, levels = rev(batchMedians$tcga_code)), y = cor), fill = tcga_code) + 
      geom_violin() + 
      geom_boxplot(width = 0.1) + 
      coord_flip() + 
      ggtitle(paste("Batch Corrected, Tumor disease:",disease)) + 
      xlab("Cell Line TCGA Code") +
      ylab("Correlation") + 
      theme(axis.text.x = element_text(face = "bold", color = "black")) +
      scale_fill_manual(values = cellColors)
  )
  dev.off()
  
  
  #finding median correlations per cell line ID
  batchMedians_lines = batchPivoted %>% dplyr::group_by(source_sample_ID) %>% dplyr::summarise(med = median(cor)) %>% arrange(desc(med))

  #adding TCGA Code information back onto batchMedians_lines and then saving as RDS for web app
  batchMedians_lines_forWeb = batchMedians_lines
  batchMedians_lines_forWeb = inner_join(batchMedians_lines_forWeb, cellMetaV2, by = "source_sample_ID")
  batchMedians_lines_forWeb = batchMedians_lines_forWeb %>% dplyr::select(-"th_sampleid")
  saveRDS(batchMedians_lines_forWeb, file = paste("results/", disease, "/", disease, "_all_cell_lines_corrs_table.RDS", sep = ""))
  
  #also subsetting to just matched cell lines (based on TCGA Code) and saving as RDS
  batchMedians_lines_matched_forWeb = batchMedians_lines_forWeb[batchMedians_lines_forWeb$tcga_code %in% getMatch(disease),]
  saveRDS(batchMedians_lines_matched_forWeb, file = paste("results/", disease, "/", disease, "_matched_cell_lines_corrs_table.RDS", sep = ""))
  
  #making subset of pivot table that just includes top 10 cell lines by median correlation
  batchCellSubset = batchPivoted[batchPivoted$source_sample_ID %in% batchMedians_lines$source_sample_ID[1:10],]

  #plotting violin plot with x axis as top 10 cell lines, colored by cell line tcga code
  pdf(paste("results/", disease, "/", disease, "_violin_by_top10_cell_line_batchCorrected.pdf", sep=""))
  print(
    ggplot(batchCellSubset, aes(x = factor(source_sample_ID, levels = rev(batchMedians_lines$source_sample_ID)), y = cor, fill = tcga_code)) + 
      geom_violin() + 
      geom_boxplot(width = 0.1) + 
      coord_flip() + 
      ggtitle(paste("Batch corrected, Tumor disease:",disease)) + 
      xlab("Cell Line") +
      ylab("Correlation") +
      labs(fill='Cell Line TCGA Code') + 
      theme(axis.text.x = element_text(face = "bold", color = "black")) +
      scale_fill_manual(values = cellColors)
  )
  dev.off()
  
  
  #making heatmap of disease-specific tumor samples with matched disease cell lines only
  matchedHeatmapDF = batchCorrelations[batchCorrelations$tcga_code %in% getMatch(disease),]
  
  #saving matched DF as RDS for web app
  saveRDS(matchedHeatmapDF, file = paste("results/",disease,"/",disease,"_correlations_matched.RDS", sep = ""))
  
  #replace cell line sample IDs with actual names
  matchedHeatmapDF = tibble::rownames_to_column(matchedHeatmapDF, "th_sampleid")
  matchedHeatmapDF = tibble::column_to_rownames(matchedHeatmapDF, "source_sample_ID")
  matchedHeatmapDF = dplyr::select(matchedHeatmapDF, -th_sampleid, -tcga_code)
  
  
  breaksList = seq(0, 1, by = 0.01) 
  
  pdf(paste("results/", disease, "/", disease, "_heatmap_matched_cell_lines_only.pdf", sep = ""))
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  
  print(
    pheatmap(matchedHeatmapDF, show_colnames = F, color = colorRampPalette(c("white", "purple"))(length(breaksList)), breaks = breaksList)
  )
  
  setHook("grid.newpage", NULL, "replace")
  grid.text("Tumor Samples", y=-0.02, x = 0.4, gp=gpar(fontsize=16)) #y axis
  grid.text("Cell Lines", x = -0.06, rot=90, gp=gpar(fontsize=16)) #x axis
  
  dev.off()
  
  jpeg(paste("results/", disease, "/", disease, "_heatmap_matched_cell_lines_only.jpg", sep = ""), pointsize = 2, res = 300, quality = 100, height = 2000, width = 2000)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  
  print(
    pheatmap(matchedHeatmapDF, show_colnames = F, color = colorRampPalette(c("white", "purple"))(length(breaksList)), breaks = breaksList)
  )
  
  setHook("grid.newpage", NULL, "replace")
  grid.text("Tumor Samples", y=-0.02, x = 0.4, gp=gpar(fontsize=16)) #y axis
  grid.text("Cell Lines", x = -0.06, rot=90, gp=gpar(fontsize=16)) #x axis
  
  dev.off()

  
  
  rm(annotation_row)
  rm(matchedHeatmapDF)
  rm(disease_samples)
  rm(topBatchCorr)
  rm(batchCorr_pca)
  rm(batchCorr_withDatabase)
  rm(topUnCorr)
  rm(unCorr_pca)
  rm(unCorr_withDatabase)
  rm(geneIQRs)
  rm(geneIQRs_unCorr)
  rm(orderedIQRs)
  rm(orderedIQRs_unCorr)
  rm(topGenes)
  rm(topGenes_unCorr)
  rm(subset)
  rm(subset_meta)
  rm(subset_batchCorr)
  rm(subset_unCorr)
  rm(subset_batchCorr_rmNoVar)
  rm(subset_unCorr_rmNoVar)
  rm(batchCellSubset)
  rm(batchCorrelations)
  rm(batchMedians)
  rm(batchMedians_lines)
  rm(batchMerged)
  rm(batchPivoted)
  rm(batchMergedIQRs)
  rm(orderedBatchMergedIQRs)
  rm(topBatchMergedGenes)
  rm(breaksList)
  if(disease == "medulloblastoma"){
    rm(toAdd)
  }
}

diseases = diseases[!diseases == "medulloblastoma"]

masterCorrs = tibble::column_to_rownames(masterCorrs, "cellIDs")
masterCorrs_forViolin = masterCorrs
masterCorrs_forMedullo = masterCorrs

#saving masterCorrs as an RDS so people can download it from the web app
saveRDS(masterCorrs, file = "results/masterCorrs.rds")

#masterCorrs is now a 799 (cell lines) x 1655 (tumor samples) matrix of correlations between cell lines and tumor samples
#where there was batch correction done for tumor samples within each tumor type

#need to collapse each dimension into disease type. let's start with cell lines.
#attaching cell line TCGA code

masterCorrs = tibble::rownames_to_column(masterCorrs, "th_sampleid")
masterCorrs = inner_join(masterCorrs, cellMetaV2, by = "th_sampleid")
masterCorrs = dplyr::select(masterCorrs, -th_sampleid, -source_sample_ID)

#collapsing cell lines based on TCGA code, using mean
master_tumor_IDs = colnames(masterCorrs[colnames(masterCorrs) != "tcga_code"]) #for use in next line's dplyr code
masterCorrs = masterCorrs %>% group_by(tcga_code) %>% dplyr::summarise(across(all_of(master_tumor_IDs), mean))
masterCorrs = tibble::column_to_rownames(masterCorrs, "tcga_code")

#transposing and doing the same with tumor samples
masterCorrs = as.data.frame(t(masterCorrs))
masterCorrs = tibble::rownames_to_column(masterCorrs, "th_sampleid")
masterCorrs = inner_join(masterCorrs, tumorMeta, by = "th_sampleid")
masterCorrs = dplyr::select(masterCorrs, -th_sampleid, -Database)

#collapsing tumor samples based on comparison disease, using mean
master_tcga_codes = colnames(masterCorrs[colnames(masterCorrs) != "comparison_disease"])#or use in next line's dplyr code
masterCorrs = masterCorrs %>% group_by(comparison_disease) %>% dplyr::summarise(across(all_of(master_tcga_codes), mean))
masterCorrs = tibble::column_to_rownames(masterCorrs, "comparison_disease")


#saving input to heatmap as RDS for web app
saveRDS(masterCorrs, file = "results/masterCorrs_fig1_heatmap.RDS")


##Generating heatmap in figure 1
breaksList = seq(min(masterCorrs), max(masterCorrs), by = 0.01) 

pdf(paste("results/Fig1d_heatmap.pdf"), height = 8, width = 8)

  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  print(pheatmap(masterCorrs, color = colorRampPalette(c("white", "purple"))(length(breaksList)), breaks = breaksList))
  setHook("grid.newpage", NULL, "replace")
  grid.text("Cell Lines", y=0, x = 0.35, gp=gpar(fontsize=16)) #y axis
  grid.text("Pediatric Cancers", x=-0.07, rot=90, gp=gpar(fontsize=16)) #x axis

dev.off()



###Generating violin plot in figure 1
masterCorrs = masterCorrs_forViolin
masterCorrs = as.data.frame(t(masterCorrs))
masterCorrs = tibble::rownames_to_column(masterCorrs, "th_sampleid")
masterCorrs = inner_join(masterCorrs, tumorMeta, by = "th_sampleid")

#pivoting so that each row is one tumor - cell line pair with correlation and tumor disease information
pivotedMaster = pivot_longer(masterCorrs, !(comparison_disease | Database | th_sampleid), names_to = "cell_id", values_to = "cor")

#adding cell line metadata to this
colnames(cellMetaV2) = c("cell_id", "source_sample_ID","tcga_code")
pivotedMaster = inner_join(pivotedMaster, cellMetaV2, by = "cell_id")

#initialize data frame that will hold only rows from pivotedMaster that have correct tumor-cell line matches
violinMaster = pivotedMaster[0,]

#loop through all tumor types to populate data frame properly
for(disease in diseases){
  violinMaster = rbind(violinMaster, 
                       pivotedMaster[pivotedMaster$comparison_disease == disease & pivotedMaster$tcga_code %in% getMatch(disease),])
}

#determining median correlations for ordering purposes in figure
violinMedians = violinMaster %>% dplyr::group_by(comparison_disease) %>% dplyr::summarise(med = median(cor)) %>% arrange(desc(med))

#saving violinMaster and violinMedians for web app
saveRDS(violinMaster, file = "results/violinMaster_fig1.RDS")
saveRDS(violinMedians, file = "results/violinMedians_fig1.RDS")

#plot violin plot in figure 1
pdf(paste("results/Fig1c_violin.pdf"), width = 10, height = 5)
print(
  ggplot(violinMaster, aes(x = factor(comparison_disease, levels = violinMedians$comparison_disease), y = cor, fill = comparison_disease)) + 
    geom_violin() + geom_boxplot(width = 0.1) + 
    ggtitle(paste("Matched Cell Line - Tumor Samples")) + 
    xlab("Tumor Type") +
    ylab("Correlation") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"))
  + scale_fill_manual(values = tumorColors)
)
dev.off()



###Medullo plots
#subsetting cell lines down to just MB lines
masterCorrs_forMedullo = masterCorrs_forMedullo[cellMetaV2$cell_id[cellMetaV2$tcga_code == "MB"],]

#subsetting tumor samples down to just MB samples
masterCorrs_forMedullo = masterCorrs_forMedullo[,tumorMeta$th_sampleid[grepl("medulloblastoma",tumorMeta$comparison_disease)]]

masterCorrs_forMedullo = tibble::rownames_to_column(masterCorrs_forMedullo, "cell_id")
masterCorrs_forMedullo = inner_join(masterCorrs_forMedullo, cellMetaV2, by = "cell_id")
masterCorrs_forMedullo = dplyr::select(masterCorrs_forMedullo, -cell_id, -tcga_code)

#pivoting individual MB cell line - MB tumor sample correlations
medullo_pivot = pivot_longer(masterCorrs_forMedullo, cols = colnames(masterCorrs_forMedullo)[!colnames(masterCorrs_forMedullo) == "source_sample_ID"], names_to = "tumor_sample", values_to = "cor")

#add labels for whether tumor sample had positive or negative PC1
medullo_pivot = mutate(medullo_pivot, PC = ifelse(tumor_sample %in% neg_PC1, "WNT-like", "SHH/Group3-like"))

#saving medullo_pivot for web app
saveRDS(medullo_pivot, file = "results/medullo_pivot.RDS")

#subtype comparison by cell line
pdf(paste("results/medullo_PC_by_line.pdf"), width = 7, height = 5)
print(
ggplot(medullo_pivot, aes(x = source_sample_ID, y = cor, fill = PC)) + 
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(width = 0.25, position = position_dodge(0.9)) + 
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9, jitter.width = 0.25)) + 
  #stat_compare_means(method = "wilcox.test", size = 10, label = "p.adj", label.y = 0.65) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(0.9))
)
dev.off()
compare_means(data = medullo_pivot, formula = cor ~ PC, group.by = "source_sample_ID", p.adjust.method = "fdr")

#bulk subtype comparison
pdf(paste("results/medullo_PC.pdf"), width = 5, height = 5)
print(
ggplot(medullo_pivot, aes(x = PC, y = cor, fill = PC)) + 
  geom_violin() + 
  geom_boxplot(width = 0.25) + 
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9, jitter.width = 0.25)) +
  stat_compare_means(method = "wilcox.test") + 
  stat_summary(fun=mean, geom = "crossbar")
)
dev.off()






###GBM subtyping?
#Original paper, list can be found in supplemental (Gene List A_500 Gene Classifier)
#https://www.nature.com/articles/s41598-019-43173-y#MOESM1 
#file used here was made by copying gene list (without header) into a new file in Excel and saving as CSV

GBM_gene_list = fread("Raw_Data/GBM/GBM_gene_list.csv", header = F)

#let's pull out pediatric samples I classified as glioma from the expected counts data

#run this only if you haven't already written glioma counts to a CSV. faster to only read in glioma counts
all_counts = fread("Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/TumorCompendium_v11_PolyA_ensembl_expected_count_58581genes_2020-04-09.tsv")
glioma_samples = tumorMeta$th_sampleid[tumorMeta$comparison_disease == "glioma"]
glioma = dplyr::select(all_counts, Gene, all_of(glioma_samples))
write.csv(glioma, file = "Intermediate_Data/glioma_expected_counts.csv", row.names = F)

glioma = fread("Intermediate_Data/glioma_expected_counts.csv")
glioma = tibble::column_to_rownames(glioma, "Gene")

glioma = as.data.frame(t(glioma))

#round counts in each column
glioma = glioma %>% dplyr::mutate(across(everything(),round))

#remove genes with no variance across all samples
glioma_vars = apply(glioma, 2, var)
ordered_glioma_vars = glioma_vars[order(glioma_vars, decreasing = T)]
non_zero_var_glioma_genes = names(ordered_glioma_vars[ordered_glioma_vars>0])
glioma = dplyr::select(glioma, all_of(non_zero_var_glioma_genes))

#select top 5000 most var genes by IQR for a PCA
glioma_IQRs = apply(glioma, 2, IQR)
ordered_glioma_IQRs = glioma_IQRs[order(glioma_IQRs, decreasing = T)]
top_glioma_genes = names(ordered_glioma_IQRs)[1:5000]
glioma_top = glioma[,top_glioma_genes]

#log transform for PCA
glioma_expected_counts_log2 = log2(glioma_top + 1)

#make PCA object
glioma_pca = fast.prcomp(glioma_expected_counts_log2, center = TRUE)

#plot PCA
pdf(paste("results/glioma_counts_top5k_pca.pdf"), width = 4, height = 4)
print(
ggbiplot(glioma_pca, var.axes = F)
)
dev.off()


#now what if we use gene list from paper?
#first need to convert gene names to Hugo symbols

glioma_names = colnames(glioma)
trimmed_glioma_names = stringr::str_extract(glioma_names, "[^\\.]*")
mart = biomaRt::useMart('ensembl', "hsapiens_gene_ensembl")
glioma_key = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values = trimmed_glioma_names, mart = mart)
glioma_key = glioma_key[!glioma_key$hgnc_symbol == "",]

glioma_key$hgnc_symbol[duplicated(glioma_key$hgnc_symbol)]
#theres some duplicated HGNC symbols, removing these
#RMRP, POLR2J4, SNORA17B, PINX1, LINC01238, SIGLEC5, DUXAP8, SNORD3D, GOLGA8M, ITFG2-AS1, SNORA16A, SNORA50A
glioma_key = glioma_key[!glioma_key$hgnc_symbol %in% glioma_key$hgnc_symbol[duplicated(glioma_key$hgnc_symbol)],]

glioma_key$ensembl_gene_id[duplicated(glioma_key$ensembl_gene_id)]
#three duplicate ensembl IDs, remove these
glioma_key = glioma_key[!glioma_key$ensembl_gene_id %in% glioma_key$ensembl_gene_id[duplicated(glioma_key$ensembl_gene_id)],]

#check for duplicates of either symbol type again
glioma_key$hgnc_symbol[duplicated(glioma_key$hgnc_symbol)]
glioma_key$ensembl_gene_id[duplicated(glioma_key$ensembl_gene_id)]

#with no duplicates, trimming glioma expected counts to just these genes
glioma = as.data.frame(t(glioma))
glioma = tibble::rownames_to_column(glioma, "ensembl_gene_id")
glioma$ensembl_gene_id = stringr::str_extract(glioma$ensembl_gene_id, "[^\\.]*")
glioma = left_join(glioma_key, glioma, by = 'ensembl_gene_id')
glioma = dplyr::select(glioma, -ensembl_gene_id)

#trimming further to just gene set genes
colnames(GBM_gene_list) = "hgnc_symbol"
glioma_GBM_gene_list = glioma[glioma$hgnc_symbol %in% GBM_gene_list$hgnc_symbol,]
rownames(glioma_GBM_gene_list) = c()
glioma_GBM_gene_list = tibble::column_to_rownames(glioma_GBM_gene_list, "hgnc_symbol")
glioma_GBM_gene_list = as.data.frame(t(glioma_GBM_gene_list))

#474 of the 500 genes are in my data -- a good start

#log transform these counts
glioma_GBM_gene_list_log2 = log2(glioma_GBM_gene_list + 1)

#make PCA object
glioma_GBM_gene_list_log2_pca = fast.prcomp(glioma_GBM_gene_list_log2, center = T)

#plot
pdf(paste("results/glioma_gene_list_pca.pdf"), width = 4, height = 4)
print(
ggbiplot(glioma_GBM_gene_list_log2_pca, var.axes = F)
)
dev.off()
#what about a heatmap
pdf(paste("results/glioma_gene_list_heatmap.pdf"))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
print(
pheatmap(glioma_GBM_gene_list_log2, show_rownames = F, show_colnames = F)
)
setHook("grid.newpage", NULL, "replace")
grid.text("Genes", y=-0.03, x = 0.5, gp=gpar(fontsize=16)) #y axis
grid.text("Tumor Samples", x=-0.07, rot=90, gp=gpar(fontsize=16)) #x axis
dev.off()


#what about if i just subset down to GBM samples instead of using all glioma samples
#need to read in original tumor clinical data for this
tumor_clinical = fread("Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv")
GBM_ped = tumor_clinical[tumor_clinical$disease == "glioblastoma multiforme" & tumor_clinical$age_at_dx <=18,]
GBM_IDs = GBM_ped$th_sampleid

#subset from log2 gene list data
GBM_only_gene_list_log2 = glioma_GBM_gene_list_log2[GBM_IDs,]

#make PCA object
GBM_only_gene_list_log2_pca = fast.prcomp(GBM_only_gene_list_log2, center = T)

#plot
pdf(paste("results/GBMonly_gene_list_pca.pdf"), width = 4, height = 4)
print(
ggbiplot(GBM_only_gene_list_log2_pca, var.axes = F)
)
dev.off()

#what about a heatmap
pdf(paste("results/GBMonly_gene_list_heatmap.pdf"))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
print(
pheatmap(GBM_only_gene_list_log2, show_rownames = F, show_colnames = F)
)
setHook("grid.newpage", NULL, "replace")
grid.text("Genes", y=-0.03, x = 0.5, gp=gpar(fontsize=16)) #y axis
grid.text("Tumor Samples", x=-0.07, rot=90, gp=gpar(fontsize=16)) #x axis
dev.off()




###NB subtyping?
#Original paper
#https://www.nature.com/articles/s43018-020-00145-w#Sec37
#Made gene list by taking “target SYMBOL” column from Supplementary table 3 and pasting into a new Excel file without the header, and saving as CSV

NB_gene_list = fread("Raw_Data/NB/NB_gene_list.csv", header = F)
NB_gene_list = unique(NB_gene_list)

#run only if nb_expected_counts.csv not previously made -- faster just to load the NB counts once they exist in their own file
all_counts = fread("Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/TumorCompendium_v11_PolyA_ensembl_expected_count_58581genes_2020-04-09.tsv")
nb_samples = tumorMeta$th_sampleid[tumorMeta$comparison_disease == "neuroblastoma"]
nb = dplyr::select(all_counts, Gene, all_of(nb_samples))
write.csv(nb, file = "Intermediate_Data/nb_expected_counts.csv", row.names = F)

#read in neuroblastoma counts
nb = fread("Intermediate_Data/nb_expected_counts.csv")
nb = tibble::column_to_rownames(nb, "Gene")
nb = as.data.frame(t(nb))

#round counts
nb = nb %>% dplyr::mutate(across(everything(),round))

#remove genes with no variance across all samples
nb_vars = apply(nb, 2, var)
ordered_nb_vars = nb_vars[order(nb_vars, decreasing = T)]
non_zero_var_nb_genes = names(ordered_nb_vars[ordered_nb_vars>0])
nb = dplyr::select(nb, all_of(non_zero_var_nb_genes))



#select top 5000 most var genes for a PCA
nb_IQRs = apply(nb, 2, IQR)
ordered_NB_IQRs = nb_IQRs[order(nb_IQRs, decreasing = T)]
top_nb_genes = names(ordered_NB_IQRs)[1:5000]
nb_top = nb[,top_nb_genes]

#log transform for PCA
nb_expected_counts_log2 = log2(nb_top+1)

#make PCA object
nb_pca = fast.prcomp(nb_expected_counts_log2, center = TRUE)

#plot PCA
pdf(paste("results/nb_counts_top5k_pca.pdf"), width = 4, height = 4)
print(
  ggbiplot(nb_pca, var.axes = F)
)
dev.off()

#what about a heatmap
pdf(paste("results/nb_counts_top5k_heatmap.pdf"))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
print(
  pheatmap(nb_expected_counts_log2, show_rownames = F, show_colnames = F)
)
setHook("grid.newpage", NULL, "replace")
grid.text("Genes", y=-0.03, x = 0.5, gp=gpar(fontsize=16)) #y axis
grid.text("Tumor Samples", x=-0.07, rot=90, gp=gpar(fontsize=16)) #x axis
dev.off()

#now what if we use gene list from paper?
#first need to convert gene names to Hugo symbols

nb_names = colnames(nb)
trimmed_nb_names = stringr::str_extract(nb_names, "[^\\.]*")
mart = biomaRt::useMart('ensembl', "hsapiens_gene_ensembl")
nb_key = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values = trimmed_nb_names, mart = mart)
nb_key = nb_key[!nb_key$hgnc_symbol == "",]

nb_key$hgnc_symbol[duplicated(nb_key$hgnc_symbol)]
#removing duplicated symbols
#POLR2J4, LINC01238, PINX1, SIGLEC5, RMRP, DUXAP8, GOLGA8M, SNORA17B, ITFG2-AS1, SNORA50A
nb_key = nb_key[!nb_key$hgnc_symbol %in% nb_key$hgnc_symbol[duplicated(nb_key$hgnc_symbol)],]

nb_key$ensembl_gene_id[duplicated(nb_key$ensembl_gene_id)]
#same with ensembl ID
nb_key = nb_key[!nb_key$ensembl_gene_id %in% nb_key$ensembl_gene_id[duplicated(nb_key$ensembl_gene_id)],]

#double check for duplicates
nb_key$hgnc_symbol[duplicated(nb_key$hgnc_symbol)]
nb_key$ensembl_gene_id[duplicated(nb_key$ensembl_gene_id)]

#with no duplicates, trimming nb expected counts to just these genes
nb = as.data.frame(t(nb))
nb = tibble::rownames_to_column(nb, "ensembl_gene_id")
nb$ensembl_gene_id = stringr::str_extract(nb$ensembl_gene_id, "[^\\.]*")
nb = left_join(nb_key, nb, by = "ensembl_gene_id")
nb = dplyr::select(nb, -ensembl_gene_id)

#trimming further to just gene set genes
colnames(NB_gene_list) = "hgnc_symbol"
nb_only_gene_list = nb[nb$hgnc_symbol %in% NB_gene_list$hgnc_symbol,]
rownames(nb_only_gene_list) = c()
nb_only_gene_list = tibble::column_to_rownames(nb_only_gene_list, "hgnc_symbol")
nb_only_gene_list = as.data.frame(t(nb_only_gene_list))

#got 1248 genes out of the 1476 in the list, not bad

#log transform these
nb_only_gene_list_log2 = log2(nb_only_gene_list + 1)

#make PCA object
nb_only_gene_list_log2_pca = fast.prcomp(nb_only_gene_list_log2, center = T)

#plot
pdf(paste("results/nb_gene_list_pca.pdf"), width = 4, height = 4)
print(
  ggbiplot(nb_only_gene_list_log2_pca, var.axes = F)
)
dev.off()

#what about a heatmap
pdf(paste("results/nb_gene_list_heatmap.pdf"))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
print(
  pheatmap(nb_only_gene_list_log2, show_rownames = F, show_colnames = F)
)
setHook("grid.newpage", NULL, "replace")
grid.text("Genes", y=-0.03, x = 0.5, gp=gpar(fontsize=16)) #y axis
grid.text("Tumor Samples", x=-0.07, rot=90, gp=gpar(fontsize=16)) #x axis
dev.off()
