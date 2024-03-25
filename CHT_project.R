# Load the necessary packages
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
BiocManager::install("GEOquery")
library(GEOquery)

BiocManager::install("DESeq2")
library(DESeq2)

files <- list.files('~/Downloads/GSE222992_RAW/New Folder With Items/')
expr_data <- list()
for(i in files){
  name <- strsplit(i,'_')[[1]][1]
  a <- read.csv(paste0('~/Downloads/GSE222992_RAW/New Folder With Items/',i),sep='\t',header = FALSE)
  expr_data[[name]] = a$V2
}

expr_data <- data.frame(expr_data)
rownames(expr_data) <- a$V1

geo <- getGEO('GSE222992')[[1]]
pheno_data <- pData(phenoData(geo))
cell_lines1 <- pheno_data$`cell line:ch1`
tissue <- pheno_data$`cell type:ch1`
cell_lines <- cbind(cell_lines1, tissue)

# exporting data to the local to further do machine learning prediction using python
write.csv(cell_lines, '~/Desktop/cell_line.csv')

write.csv(t(expr_data), '~/Desktop/EXP_DATA.csv')

#Reading in the files from GEO data base after differential expression using limma in GEO
cell_1_T1 <- read.csv("cell_1_T1.tsv", sep = "\t")
cell_1_T2 <- read.csv("cell_1_T2.tsv", sep = "\t")
cell_2_T1 <- read.csv("cell_2_T1.tsv", sep = "\t")
cell_2_T2 <- read.csv("cell_2_T2.tsv", sep = "\t")

install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
library(ggplot2)
library(ggrepel)
library(dplyr)

#seeing the gene expression of the diseased of teh data set 
install.packages("limma")
library(limma)

#Differential expression visualization
cell_1_T1$diffexpressed <- "NO"

cell_1_T1$diffexpressed[cell_1_T1$log2FoldChange>0.58 & cell_1_T1$padj<0.01] <- "UP"
cell_1_T1$diffexpressed[cell_1_T1$log2FoldChange<0.58 & cell_1_T1$padj<0.01] <- "DOWN"

cell_1_T1$delabel<- NA

ggplot(data = cell_1_T1, aes(x=LogFC, y=-log10(Pval), col=diffexpressed,label=delabel))+
    geom_point()+
    theme_minimal()+
    geom_text_repel()+
    scale_colour_manual(values=c("blue","black","red"))+
    theme(text=element_text(size = 20))


cell_1_T2$diffexpressed <- "NO"

cell_1_T2$diffexpressed[cell_1_T2$log2FoldChange>0.58 & cell_1_T2$padj<0.01] <- "UP"
cell_1_T2$diffexpressed[cell_1_T2$log2FoldChange<0.58 & cell_1_T2$padj<0.01] <- "DOWN"

cell_1_T1$delabel<- NA

ggplot(data = cell_1_T2, aes(x=LogFC, y=-log10(Pval), col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_colour_manual(values=c("blue","black","red"))+
  theme(text=element_text(size = 20))

cell_2_T1$diffexpressed <- "NO"

cell_2_T1$diffexpressed[cell_2_T1$log2FoldChange>0.58 & cell_2_T1$padj<0.01] <- "UP"
cell_2_T1$diffexpressed[cell_2_T1$log2FoldChange<0.58 & cell_2_T1$padj<0.01] <- "DOWN"

cell_2_T1$delabel<- NA

ggplot(data = cell_2_T1, aes(x=LogFC, y=-log10(Pval), col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_colour_manual(values=c("blue","black","red"))+
  theme(text=element_text(size = 20))


cell_2_T2$diffexpressed <- "NO"

cell_2_T2$diffexpressed[cell_2_T2$log2FoldChange>0.58 & cell_2_T2$padj<0.01] <- "UP"
cell_2_T2$diffexpressed[cell_2_T2$log2FoldChange<0.58 & cell_2_T2$padj<0.01] <- "DOWN"

cell_2_T1$delabel<- NA

ggplot(data = cell_2_T1, aes(x=LogFC, y=-log10(Pval), col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_colour_manual(values=c("blue","black","red"))+
  theme(text=element_text(size = 20))


#Filtering out teh genes according to the p value 
D1_cell_1_T1 <- subset(cell_1_T1, padj < 0.05 & abs(log2FoldChange) > 1)
D2_cell_1_T2 <- subset(cell_1_T2, padj < 0.05 & abs(log2FoldChange) > 1)
D3_cell_2_T1 <- subset(cell_2_T1, padj < 0.05 & abs(log2FoldChange) > 1)
D4_cell_2_T2 <- subset(cell_2_T2, padj < 0.05 & abs(log2FoldChange) > 1)


# Define the gene sets to analyze (e.g. differentially expressed genes)
de_genes1 <- D1_cell_1_T1$Symbol
de_genes2 <- D2_cell_1_T2$Symbol
de_genes3 <- D3_cell_2_T1$Symbol
de_genes4 <- D4_cell_2_T2$Symbol


#saved the files to do the pathway analysis in shiny 
write.csv(de_genes1, "de_genes1.csv")

write.csv(de_genes2, "de_genes2.csv")

write.csv(de_genes3, "de_genes3.csv")

write.csv(de_genes4, "de_genes4.csv")

enrich_C1_T1 <- read.csv("c1_t1_enrichment_all.csv")
enrich_C2_T2  <- read.csv("c1_t2_enrichment_all.csv")


