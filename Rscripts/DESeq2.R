#### ~ Library loading ~ ####
#BiocManager::install("DESeq2")
library("DESeq2")
library(dplyr)
library("RColorBrewer")
library("gplots")
library( "genefilter" )
library("pheatmap")
library(tools)
library("ggplot2")
library(magrittr)
library("biomaRt")
library(goseq)
#BiocManager::install("biomaRt")
library(apeglm)
#BiocManager::install("apeglm")
library("genefilter")
library(tidyr)
#BiocManager::install("EnhancedVolcano")
library(extrafont)
library(GenomicFeatures)
loadfonts(device = "win")
windowsFonts()
#### ~ DECLARE YOUR VARIABLES HERE ~ ####
myspecies <- "C.mydas"
directory <- "C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/htseq_test/htseq/"
setwd("C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/htseq_test/htseq/")
testgroup <- "analysis02"
sample_table <- read.table(paste0("C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/Latest_table_as_st15102021.txt"), sep = "\t", header=T)
length(sample_table)
dim(sample_table)
tail(sample_table)
tail(sample_table <- sample_table[c(1:81),])

meta2 <- read.table("C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/RNA_metadata_151021.txt", header = T, sep = "\t")
dim(meta2)
(meta2 <- meta2[c(1:81),])
tail(meta2)
backup <- sample_table

sample_table = merge(meta2, sample_table, by.x="ID.Number", by.y="ID.Number")

source('../../Rscripts/myheatmap_plot.R')
source('../../Rscripts/mytopgene.R')
source('../../Rscripts/Myvolcano_plot.R')
source("../../Rscripts/mygoplots.R")

sampleFiles <- paste0(pull(sample_table, 'ID.Number'),'.tsv')
sampleFiles
all(file.exists(sampleFiles))

names(sampleFiles) <- pull(sample_table, 'ID.Number')
names(sampleFiles)

head(sample_table)

gtf <- "C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/NCBI/GCF_015237465.2_rCheMyd1.pri.v2_genomic.gtf"
txdb = makeTxDbFromGFF( gtf, format = "gtf")
txsBygene=transcriptsBy(txdb,"gene")
lengthdata=mean(width(txsBygene))

sample_table$Arsenic.Cat
Locality <- 'Locality.x'
Age <- 'Age.x'
Sex <- 'Sex.x'
Season <- 'Season.x'
PCV <- 'PCVG'
CCL <- 'CCL'
BCI <- 'BCI'
sample_table$TP..g.L..x
FP <- 'Fibropapillomatosis.x'
TP <- 'TPG'

#Cadmium
cad <- "Camium.Cat"

#Aluminium
ala <- "Aluminium.Cat"

#Molibdinum
Mold <- "Molybdenum.Cat"

#Selenium
seln <- "Selenium.Cat"

#Arsenic
arsenic <- 'Arsenic.Cat'

#Magnesium
Magnesium <- 'Magnesium.Cat'
sample_table$MgG

#Iron
Iron <- 'Iron.Cat'

#Cobalt
sample_table$CobaltG
Cobalt <- 'Cobalt.Cat'

#Manganese
sample_table$MnG
Manganese <- 'Manganese.Cat'

#Copper
sample_table$CopperG
Copper <- 'Copper.Cat'

#Zinc
sample_table$ZincG
Zinc <- 'Zin.Cat'

#Strontium
sample_table$StrontiumG
Strontium <- 'Strontium.Cat'

#CCL
sample_table$CCLG
ccl <- 'CCLG'
sample_table$BCS

#BCI
sample_table$BCSG

sampleLoc<-pull(sample_table, Locality)
sampleLoc
sampleAge <- pull(sample_table, Age)
sampleSex <- pull(sample_table, Sex)
sampleSeason <- pull(sample_table, Season)
samplefibro <- pull(sample_table, 'Fibropapillomatosis.y')
samplefibro
pcv <- 'PCVG'
samplePCV <- pull(sample_table, pcv)
samplePCV
sampleCCL <- pull(sample_table, ccl)
sampleCCL 
sampleBCI <- pull(sample_table, BCSG)
sampleBCI

#TP
sample_table$TPG
bci <- 'BCIG'
tpg <- 'TPG'
sampleTP <- pull(sample_table, tpg)
sampleTP

arsenic <- "Arsenic.Cat"
sampleArsenic <- pull(sample_table, arsenic)
sampleArsenic
length(sampleFiles)
sample_table$ID.Number
sampleFiles <- paste0(pull(sample_table, 'ID.Number'),'.tsv')

#Final table creation
sampleTable <- data.frame(sampleName = sample_table$ID.Number, 
                          fileName = sampleFiles,
                          Age = as.factor(sampleAge),
                          Sex = as.factor(sampleSex),
                          CCL = as.factor(sampleCCL),
                          TP = as.factor(sample_table$TPG),
                          Arsenic = as.factor(sample_table$Arsenic.Cat),
                          Magnesium = as.factor(sample_table$Magnesium.Cat),
                          Manganese = as.factor(sample_table$Magnesium.Cat),
                          Iron = as.factor(sample_table$Iron.Cat),
                          Cobalt = as.factor(sample_table$Cobalt.Cat),
                          Locality = as.factor(sampleLoc),
                          #Manganese = as.factor(sample_table$MnG),
                          Copper = as.factor(sample_table$Copper.Cat),
                          Zinc = as.factor(sample_table$Zinc.Cat),
                          Strontium = as.factor(sample_table$Strontium.Cat),
                          Selnium = as.factor(sample_table$Selenium.Cat),
                          #Cobalt = as.factor(sample_table$Cobalt.Cat),
                          Molybdinum = as.factor(sample_table$Molybdenum.Cat),
                          Cadmium = as.factor(sample_table$Cadmium.Cat),
                          Aluminium = as.factor(sample_table$Aluminium.Cat),
                          PCV = as.factor(samplePCV),
                          Season = as.factor(sampleSeason),
                          BCI = as.factor(sampleBCI),
                          FP = as.factor(samplefibro))
head(sampleTable)
resultsdir <- 'C:/Users/kar131/OneDrive - CSIRO/15.Chritable_green_turtle/analysis_07122021/'
dir.create(resultsdir, recursive = T)
write.csv(sample_table, file = paste0(resultsdir, 'fulltable.csv'), quote = F)
write.csv(sampleTable, file = paste0(resultsdir, 'table_matrix.csv'), quote = F)

#Design #1
design1 = '~ Season + Age + BCI + Sex + FP + PCV + TP'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex + FP + PCV + TP)
dds
dds$Season
dds$Age
dds$Sex
dds$CCL
dds$BCI
dds$FP
dds$Arsenic
dds$Manganese
dds$Magnesium
dds$TP
dds$Strontium
dds$Iron

keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP

##Which is the control #DE analysis basic data generation
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
nonnormalized_counts <- counts(dds)
write.table(nonnormalized_counts, file=paste(resultsdir, 
                                             testgroup, 
                                             "_cmydas_non_normalized_counts.txt", sep = ""), 
            sep="\t", quote=F, col.names=NA)
normalized_counts <- counts(ddsKEEP, normalized=TRUE)
write.table(normalized_counts, file=paste(resultsdir, testgroup, 
                                          "_cmydas_normalized_counts.txt", sep = ""), 
            sep="\t", quote=F, col.names=NA)

#### Design #1 - sex####
resultsNames(ddsKEEP)
contrast <- 'Sex_Male_vs_Female'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange<0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
source('../../Rscripts/myheatmap_plot.R')
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Sex, ddsKEEP = ddsKEEP)
mytopgene(res, 'Sex', contrast, ddsKEEP)
myvolcano(res)


#### Design #2 - PCV ####
(resultsNames(ddsKEEP))
contrast <- 'PCV_Low_vs_High'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]

#### GSE - PCV Low vs High ####
source("../../Rscripts/preparing_for_GSE.R")
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$PCV, ddsKEEP = ddsKEEP)
mytopgene(res, 'PCV', contrast, ddsKEEP)
myvolcano(res)

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)
#### Design #3 - TP ####
(resultsNames(ddsKEEP))
contrast <- 'TP_Low_vs_High'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$TP, ddsKEEP = ddsKEEP)
mytopgene(res, 'TP', contrast, ddsKEEP)
myvolcano(res)

#### Design #4 - FP ####
(resultsNames(ddsKEEP))
contrast <- 'FP_Y_vs_N'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$FP, ddsKEEP = ddsKEEP)
mytopgene(res, 'FP', contrast, ddsKEEP)
myvolcano(res)




#### Design #5 - BCI ####
(resultsNames(ddsKEEP))
contrast <- 'BCI_Normal_vs_Low'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$BCI, ddsKEEP = ddsKEEP)
mytopgene(res, 'BCI', contrast, ddsKEEP)
myvolcano(res)

#### GSE - BCI Normal vs Low ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$PCV, ddsKEEP = ddsKEEP)
mytopgene(res, 'PCV', contrast, ddsKEEP)
myvolcano(res)
#### Design #6 - Season Spring vs Autumn ####
(resultsNames(ddsKEEP))
contrast <- 'Season_Spring_vs_Autumn'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Season, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Season Spring vs Autumn ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #7 - Season Winter vs Autumn ####
(resultsNames(ddsKEEP))
contrast <- 'Season_Winter_vs_Autumn'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Season, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Season Winter vs Autumn ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #8 - Season Spring vs Winter ####
sub_sampleTable <- sampleTable[sampleTable$Season == c("Spring", "Winter"),]
sub_sampleTable
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Season)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP

##Which is the control #DE analysis basic data generation
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)


(resultsNames(ddsKEEP))
contrast <- 'Season_Winter_vs_Spring'
res <- results(ddsKEEP, name = contrast, alpha = 0.05)
res
table(is.na(res$padj))
sum(res$padj<0.05, na.rm = T)
resSig <- res[which(res$padj<0.05),]
sum(res$log2FoldChange>0.58, na.rm = T)
res[which(res$log2FoldChange>0.58),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Season, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Season Winter vs Spring ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 


mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #9 - Arsenic####
design2 = '~ Season + Age + BCI + Sex  + FP  + PCV   + TP + Arsenic'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Arsenic)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)


contrast <- 'Arsenic_High_vs_Control'
resArsenic <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resArsenic
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Arsenic, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Arsenic vs control ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #10 - Magnesium####
design3 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Magnesium'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Magnesium)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Magnesium_Low_vs_Control'
resMg <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resMg
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cadmium
          , ddsKEEP = ddsKEEP)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Magnesium, ddsKEEP = ddsKEEP)
myvolcano(res)

contrast <- 'Magnesium_High_vs_Control'
resMg <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resMg
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cadmium
          , ddsKEEP = ddsKEEP)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Magnesium, ddsKEEP = ddsKEEP)
myvolcano(res)


#### Design #11 - Iron####
design4 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Iron'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Iron)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- DESeq(ddsKEEP)
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Iron_High_vs_Control'
resIron <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resIron
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Iron, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Iron high vs Control ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #12 - Cobalt####
design5 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Cobalt'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Cobalt)
colData(dds)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP

ddsKEEP <- DESeq(ddsKEEP)

#ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Cobalt_High_vs_Control'
resCo <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resCo
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cobalt, ddsKEEP = ddsKEEP)
myvolcano(res)

#### Design #13 - Manganese####
design6 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Manganese'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Manganese)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Manganese_Low_vs_Control'
resMn <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resMn
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Manganese, ddsKEEP = ddsKEEP)

contrast <- 'Manganese_High_vs_Control'
resMn <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resMn
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Manganese, ddsKEEP = ddsKEEP)

#### Design #14 - Copper####
design7 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Copper'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Copper)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Copper_Low_vs_Control'
resCu <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resCu
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Copper, ddsKEEP = ddsKEEP)
myvolcano(res)


contrast <- 'Copper_High_vs_Control'
resCu <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resCu
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Copper, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Copper High vs Control ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #15 - Zinc####
design8 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Zinc'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Zinc)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Zinc_Low_vs_Control'
resZn <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resZn
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Zinc, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Zinc High vs Control ####
contrast <- 'Zinc_High_vs_Control'
resZn <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resZn
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Zinc, ddsKEEP = ddsKEEP)
myvolcano(res)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #16 - Strontium####
design9 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Strontium'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Strontium)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- 'Strontium_Low_vs_Control'
resStrontium <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resStrontium
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Strontium, ddsKEEP = ddsKEEP)
myvolcano(res)

contrast <- 'Strontium_High_vs_Control'
resStrontium <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resStrontium
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Strontium, ddsKEEP = ddsKEEP)
myvolcano(res)

#### Design #17 - Locality ####

#design10 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Locality'
sampleTable

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

ddsKEEP$Locality

contrast <- 'Locality_Facing.Island_vs_Central.South.Curtis'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)

contrast <- 'Locality_Wiggins_vs_Central.South.Curtis'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)

###
sub_sampleTable <- sampleTable[sampleTable$Locality == c("South Trees", "Wiggins"),]
dim(sub_sampleTable)
sub_sampleTable$Locality <- as.factor(sub_sampleTable$Locality)
sub_sampleTable$Locality
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ ~ Season + Age + BCI + Sex + FP + PCV + TP + Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

contrast <- 'Locality_Wiggins_vs_South.Trees'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)

###
sub_sampleTable <- sampleTable[sampleTable$Locality == c("South Trees", "Pelican Banks"),]
dim(sub_sampleTable)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ ~ Season + Age + BCI + Sex + FP + PCV + TP + Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

contrast <- 'Locality_South.Trees_vs_Pelican.Banks'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)

#### GSE - Locality Wiggins vs PB ####
sub_sampleTable <- sampleTable[sampleTable$Locality == c("Wiggins", "Pelican Banks"),]
dim(sub_sampleTable)
sub_sampleTable$Locality
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

contrast <- 'Locality_Wiggins_vs_Pelican.Banks'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)
###
sampleTable$Locality
sub_sampleTable <- sampleTable[sampleTable$Locality == c("South Trees", "MotB"),]
dim(sub_sampleTable)
sub_sampleTable$Locality
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)


#### GSE - Locality South tree vs MotB ####
contrast <- 'Locality_South.Trees_vs_MotB'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
(res <- resLocality)
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

###
sampleTable$Locality
sub_sampleTable <- sampleTable[sampleTable$Locality == c("Pelican Banks", "MotB"),]
dim(sub_sampleTable)
sub_sampleTable$Locality
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Locality Pelican Bay vs MotB ####
contrast <- 'Locality_Pelican.Banks_vs_MotB'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

###
sampleTable$Locality
sub_sampleTable <- sampleTable[sampleTable$Locality == c("Wiggins", "MotB"),]
dim(sub_sampleTable)
sub_sampleTable$Locality
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Locality)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Locality Wiggins vs MotB ####
contrast <- 'Locality_Wiggins_vs_MotB'
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myvolcano(res)
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Locality, ddsKEEP = ddsKEEP)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #18 - Aluminium ####
design10 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Aluminium'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Aluminium)
colData(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Aluminium Low vs High ####
contrast <- "Aluminium_Low_vs_High"
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Aluminium, ddsKEEP = ddsKEEP)
myvolcano(res)
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #19 - Cadmium ####
sampleTable$Molybdinum
design10 = '~ Season + Age + BCI + Sex + FP + PCV + TP + cad'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Cadmium)
colData(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- "Cadmium_Low_vs_High"
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cadmium, ddsKEEP = ddsKEEP)
myvolcano(res)

#### GSE - Cadmium Low vs High ####
assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #20 - Molybdenum ####
sampleTable$Molybdinum
design12 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Molybdinum'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Molybdinum)
colData(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)
contrast <- "Molybdinum_High_vs_Control"
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cadmium, ddsKEEP = ddsKEEP)
myvolcano(res)

#### Design #21 - Selenium ####
sampleTable$Selnium
design12 = '~ Season + Age + BCI + Sex + FP + PCV + TP + Selnium'
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Season + CCL + BCI + Sex  + FP  + PCV   + TP + Selnium)
colData(dds)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Selenium High vs Control ####
contrast <- "Selnium_High_vs_Control"
resLocality <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resLocality
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Cadmium, ddsKEEP = ddsKEEP)
myvolcano(res)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

#### Design #22 - Age ####
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ Age)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Age Juvenile vs Adult ####
contrast <- "Age_Juvenile_vs_Adult"
resAge <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resAge
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Age, ddsKEEP = ddsKEEP)
myvolcano(res)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

contrast <- "Age_Subadult_vs_Adult"
resAge <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resAge
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Age, ddsKEEP = ddsKEEP)
myvolcano(res)

###
sampleTable$Age
sub_sampleTable <- sampleTable[sampleTable$Age == c("Subadult", "Juvenile"),]
sub_sampleTable$Age
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ Age)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - Age Sub-adult vs Juvenile ####
contrast <- "Age_Subadult_vs_Juvenile"
resAge <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resAge
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$Age, ddsKEEP = ddsKEEP)
myvolcano(res)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)

sampleTable$fileName<- as.factor(sampleTable$fileName)
sampleTable$fileName

sub_sampleTable <- sampleTable[which(sampleTable$fileName %in% c("697.tsv", "755.tsv","801.tsv")),]
sub_sampleTable <- sampleTable[which(sampleTable$fileName %in% c("697.tsv", "755.tsv", "801.tsv", "803.tsv", 
                                                         "724.tsv", "772.tsv", "776.tsv", "842.tsv")),]
sub_sampleTable
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sub_sampleTable,
                                  directory = directory,
                                  design= ~ FP)
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP
ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP, maxit=500)
ddsKEEP <- ddsKEEP[mcols(ddsKEEP)$betaConv,]
colData(ddsKEEP)
resultsNames(ddsKEEP)

#### GSE - FP Y vs N ####
contrast <- "FP_Y_vs_N"
resAge <- results(ddsKEEP, alpha = 0.05, name = contrast)
res <- resAge
resSig <- res[which(res$padj<0.05),]
write.csv( as.data.frame(resSig), file=paste(resultsdir, contrast, "_DEGS_cmydas.csv", sep = "") )
myheatmap(res = res, contrast = contrast, intgroup = ddsKEEP$FP, ddsKEEP = ddsKEEP)
myvolcano(res)

assayed.genes <- rownames(res)
de.genes <- rownames(resSig)
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

G2GO <- gene2go[which(gene2go$symbol %in% names(gene.vector)),]
dim(G2GO)

geneID2GO <- by(G2GO$symbol,
                G2GO$V14,
                function(x) as.data.frame(x))

head(geneID2GO,15)
head(lengthdata)
where <- match(names(gene.vector), names(lengthdata))
lengthdata <- lengthdata[where]

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

goResults

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", contrast, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", contrast, "_enriched_go_", myspecies, ".csv"))

goResults = goResults[goResults$over_represented_pvalue<0.05,]
goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 
goBPRes

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term)) + 
  geom_bar(stat = "identity", colour = "green", fill = "forestgreen") +
  theme(text = element_text(size = 12, family = "serif")) +
  labs(x="Hits (%)", y=NULL, colour="p value", size="Count") 

mygoplot(contrast = contrast, pwf = pwf, goResults = goResults, myspecies = myspecies
         , resultsdir = resultsdir, method = "Wallenius", repcnt = 2000, geneID2GO = geneID2GO)
