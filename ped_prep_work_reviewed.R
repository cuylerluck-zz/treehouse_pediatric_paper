library(data.table) #version 1.14.0
library(pheatmap) #version 1.0.12
library(gplots) #version 3.1.1
library(ggplot2) #version 3.3.5
library(pheatmap) #version 1.0.12
library(viridis) #version 0.6.1
library(edgeR) #version 3.30.3
library(dplyr) #version 1.0.7
library(tidyr) #version 1.1.3
library(ggbiplot) #version 0.55
library(gmodels) #version 2.18.1
library(umap) #version 0.2.7.0
library(corrplot) #version 0.90
library(grDevices) #version 4.0.3
library(sva) #version 3.36.0


#Set working directory as needed
setwd("/Volumes/Cuyler_T5/UCSF Sirota Lab Rotation")


###Loading in datasets

#Tumor Compendium v11 Public PolyA (April 2020) Clinical Data
tCompendium_v11_Clinical = fread(file = "Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv")

#Tumor Compendium v11 Public PolyA (April 2020) TPM Expression
tCompendium_v11_TPM = fread(file = "Raw_Data/Tumor_Compendium_v11_Public_PolyA_April2020/TumorCompendium_v11_PolyA_hugo_log2tpm_58581genes_2020-04-09.tsv")

#Cell Line Compendium v2 (December 2019) Clinical Data
cCompendium_v2_Clinical = fread(file = "Raw_Data/Cell_Line_Compendium_v2_December2019/clinical_TreehouseCellLineCompendium_v2_2019-12-02.tsv")

#Cell Line Compendium v2 (December 2019) TPM Expression
cCompendium_v2_TPM = fread(file= "Raw_Data/Cell_Line_Compendium_v2_December2019/TreehouseCellLineCompendium_v2_hugo_log2tpm_58581genes_2019-12-02.tsv")

#Tumor Compendium Metadata -- select diseases (>=30 samples) and age at dx <= 18 with comparison_disease that I defined
tumorMeta = fread(file = "Intermediate_Data/tCompendiumv11_pediatric_selectDiseases_forComparison_meta.csv", header = TRUE)
tumorMeta = select(tumorMeta, th_sampleid, comparison_disease, Database)

#Cell Compendium Metadata v2 -- cell lines from Treehouse compendium with TCGA acronyms taken from CCLE annotation file
#Excluding those with no acronym or unclassified acronym
cellMetaV2 = fread(file = "Intermediate_Data/cellMetaV2.csv", header = TRUE)
cellMetaV2 = select(cellMetaV2, th_sampleid, source_sample_ID, tcga_code)

#Select Pediatric Tumor TPM Data
pedData = fread("Intermediate_Data/tCompendium_v11_selectPedSamples_TPM.csv")

#Select Cell Line TPM Data
cellData = fread("Intermediate_Data/cellData_TPM.csv")

#CCLE Annotations
ccle = fread(file = "Raw_Data/ccle/Cell_lines_annotations_20181226.txt")




############

###Generating cellMetaV2

#Checking cell line disease annotations in Treehouse sheet vs. CCLE data

ccle_subset = select(ccle, c("CCLE_ID", "Name","Disease","tcga_code"))

th_originalCellClinical_subset = select(cCompendium_v2_Clinical, c("th_sampleid", "Source sample ID", "CCLE ID","TCGA Acronym","disease"))

#changing some column names so I know what variables are from what
colnames(ccle_subset) = sapply(colnames(ccle_subset), function(x){paste("CCLEanno ",x)})

colnames(th_originalCellClinical_subset) = sapply(colnames(th_originalCellClinical_subset), function(x){paste("TH Original ", x)})

#but also need to make shared column names for Source sample ID and CCLE ID across some datasets
#so that I can use these columns for left_joins

colnames(ccle_subset)[1] = "CCLE_ID"

colnames(th_originalCellClinical_subset)[2] = "Source sample ID"
colnames(th_originalCellClinical_subset)[3] = "CCLE_ID"

mergedMeta = left_join(th_originalCellClinical_subset, ccle_subset, by = "CCLE_ID")

View(table(mergedMeta$`CCLEanno  tcga_code`))

misfits = mergedMeta[is.na(mergedMeta$`CCLEanno  tcga_code`) | mergedMeta$`CCLEanno  tcga_code` == "UNABLE TO CLASSIFY",]

#dropping cell lines with either no TCGA code or unable to classify tcga code

mergedMeta = mergedMeta[!(is.na(mergedMeta$`CCLEanno  tcga_code`) | mergedMeta$`CCLEanno  tcga_code` == "UNABLE TO CLASSIFY"),]

cellMetaV2 = select(mergedMeta, c("TH Original  th_sampleid", "Source sample ID", "CCLEanno  tcga_code"))

colnames(cellMetaV2) = c("th_sampleid", "source_sample_ID", "tcga_code")

#writing updated cellMetaV2 sheet to the Intermediate_Data folder
write.csv(cellMetaV2, file = "Intermediate_Data/cellMetaV2.csv")

#visualizing the tumor types quickly
pdf(paste("results/cell_line_types.pdf"))
print(
ggplot(cellMetaV2, aes(x = tcga_code)) + geom_bar() + ggtitle("Cell Lines Grouped by Tumor Code") + theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", color = "black"))
)
dev.off()
tbl = table(cellMetaV2$tcga_code)

############







###Generating cellData

###Subsetting the cell line TPM data down to just cell lines with tumor codes (i.e. in cellMetaV2) and saving as .csv

cellData = select(cCompendium_v2_TPM,c("Gene",cellMetaV2$th_sampleid))

write.csv(cellData, "Intermediate_Data/cellData_TPM.csv", row.names = FALSE)

############









###Generating tumorMeta
###Read in tCompendium_v11_Clinical at top of script

###Looking at characteristics (e.g. child-young adult vs. adult? type of disease?) of the two datasets

#sorting out only pediatric samples. I am defining this as 18 years or younger at diagnosis, not using the pedaya status
#given in the dataset which was age < 30. I think 18 years is probably more representative of peds.
tCompendium_v11_Clinical_peds = tCompendium_v11_Clinical[tCompendium_v11_Clinical$age_at_dx<=18,]

#Looking at what tumor types are represented in the pediatric set
pedTumorTypes = dplyr::count(tCompendium_v11_Clinical_peds, disease)

#Visualizing counts of tumor types with a bar chart
ggplot(pedTumorTypes, aes(x=disease, y=n)) + geom_bar(stat="identity") + theme(axis.text.x = element_text (angle = 90)) + coord_flip()


#A rule of thumb I have heard is n = 30 is a lower bound for satisfactory sample size. There are many tumor types with less than 30 samples
#so I am going to filter only those that have 30 or more.

pedTumorTypes30Plus = filter(pedTumorTypes, n>=30)

#Visualizing filtered counts
ggplot(pedTumorTypes30Plus, aes(x=disease, y=n)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90)) + coord_flip()


#Now filtering the full pediatric clinical data to only contain samples for these 12 tumor types
tCompendium_v11_Clinical_peds_30Plus = filter(tCompendium_v11_Clinical_peds, disease %in% pedTumorTypes30Plus$disease)

#Verifying this worked properly
unique(tCompendium_v11_Clinical_peds_30Plus$disease)

#I am interested in which samples for each tumor type correspond to which database. I will add a variable to the df called Database which uses
#the sample ID to figure out which database each sample is from.

#see here for a full list of studies corresponding to site_ids: https://treehousegenomics.soe.ucsc.edu/public-data/dataset-accessions-legend.html

assignDB = function(sampleID){
  if(grepl("^TH\\d",sampleID) ==TRUE){
    return("Treehouse")
  }
  else if(grepl("^TCGA",sampleID)==TRUE){
    return("TCGA")
  }
  else if(grepl("^TARGET",sampleID)==TRUE){
    return("TARGET")
  }
  else if(grepl("^THR",sampleID)==TRUE){
    return("Other")
  }
  else{
    return("Unidentifiable")
  }
}

tCompendium_v11_Clinical_peds_30Plus_DB = mutate(tCompendium_v11_Clinical_peds_30Plus, Database = sapply(th_sampleid,assignDB))

#Now looking again at number of samples for these tumor types, but stacked bars colored for database
ggplot(tCompendium_v11_Clinical_peds_30Plus_DB, aes(fill=Database, x=disease)) + geom_histogram(stat="count") + theme(axis.text.x=element_text(angle=90)) + coord_flip()

#I am now going to create a new variable for the tumor samples called comparison_disease where I am re defining the disease those samples
#are to be compared within. This will be the same as the original disease for almost all samples, but some are subtypes that are being
#grouped into more general categories because they lack the number of samples to make their own group.
#This re-grouping must be done before the 30 sample cutoff, which I then have to re-do after this re-grouping.

tumorSuperGroup = function(disease){
  
  if(disease == "glioblastoma multiforme"){
    return("glioma")
  }
  else if(disease == "gliomatosis cerebri"){
    return("glioma")
  }
  else{
    return(tolower(disease))
  }
  
}

tCompendium_v11_Clinical_peds_comparingDiseases = mutate(tCompendium_v11_Clinical_peds, comparison_disease = sapply(disease,tumorSuperGroup))


regroupedDiseaseCounts = dplyr::count(tCompendium_v11_Clinical_peds_comparingDiseases, comparison_disease)

pedTumorTypes30PlusRegrouped = filter(regroupedDiseaseCounts, n>=30)

#the number of gliomas is now 197, while it was 168 before. This reflects the 28 GBM and 1 gliomatosis cerebri that were regrouped.

tCompendium_v11_Clinical_peds_comparingDiseases_30Plus = filter(tCompendium_v11_Clinical_peds_comparingDiseases, comparison_disease %in% pedTumorTypes30PlusRegrouped$comparison_disease)

tCompendium_v11_Clinical_peds_comparingDiseases_30PlusDB = mutate(tCompendium_v11_Clinical_peds_comparingDiseases_30Plus, Database = sapply(th_sampleid,assignDB))

#Now let's re plot the tumor diseases with # samples faceted by database to see again
pdf(paste("results/ped_tumor_sample_types.pdf"))
print(
ggplot(tCompendium_v11_Clinical_peds_comparingDiseases_30PlusDB, aes(fill=Database, x=comparison_disease)) + geom_histogram(stat="count") + theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", color = "black"))
)
dev.off()
tbl = table(tCompendium_v11_Clinical_peds_comparingDiseases_30PlusDB$comparison_disease)

###Creating simplified metadata sheets for the pediatric tumor samples that just contain sample id, comparison disease, and  database
###These tumor samples are pediatric & diseases which have >=30 samples

tumorMeta = select(tCompendium_v11_Clinical_peds_comparingDiseases_30PlusDB, th_sampleid, comparison_disease, Database)
write.csv(tumorMeta, file = "Intermediate_Data/tCompendiumv11_pediatric_selectDiseases_forComparison_meta.csv")


############


###Generating pedData

###Segmenting out the raw tumor data into just pediatric so I don't have to load all of it all of the time

#Read in tumorMeta & tCompendium_v11_TPM from the top of the script
#tumorMeta already only contains sampleIDs for pediatric (age <=18 at dx) samples and disease types that I have chosen (those with >=30 samples)

allTumorSamples = colnames(tCompendium_v11_TPM)

#identifying the column indices that match samples in the tumorMeta sheet
samplesToKeep = which(allTumorSamples %in% tumorMeta$th_sampleid)

#must add the first column too, this is the column containing Gene information
samplesToKeep = c(1,samplesToKeep)

#subsetting the tumor RNAseq TPM data to only keep Gene names & samples that are in tumorMeta
tCompendium_v11_selectPedSamples_TPM = tCompendium_v11_TPM[, ..samplesToKeep]

write.csv(tCompendium_v11_selectPedSamples_TPM, file = "Intermediate_Data/tCompendium_v11_selectPedSamples_TPM.csv", row.names = FALSE)

#testing reading it back in
# test = fread("Intermediate_Data/tCompendium_v11_selectPedSamples_TPM.csv")
# smallTest = test[1:10,1:10]
# smallBeforeWritten = tCompendium_v11_selectPedSamples_TPM[1:10,1:10]
#looks like it worked perfectly!

############


