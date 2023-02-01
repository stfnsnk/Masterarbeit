#!/usr/bin/env R

## RNA-seq analysis using DESEQ2
# code template from 
# - Stephen Turner (https://gist.github.com/stephenturner/f60c1934405c127f09a6)
# - Hugo Varet (https://github.com/PF2-pasteur-fr/SARTools/blob/master/R/loadCountData.R)

sample_info_file <- "./ONT_results/Sample_Info.txt"
feature_count_files <- "./ONT_results/featureCounts_output/"

#starting log file
line <- paste("sample_info_file: ", sample_info_file,"count files: ", feature_count_files, sep = "\n")
write(paste("DESEQ2 - log \n",line,sep = "\n"), file="DESEQ2-log.txt")

# reading in count data
Sample_info <- read.table(sample_info_file, header = TRUE); Sample_info
Sample_info <- Sample_info[order(Sample_info$condition),]; Sample_info

# find all count files in featureCounts output directory
file_paths <- list.files(path = feature_count_files, pattern = "*.tsv$", full.names = TRUE); file_paths

# find the first file of the Sample info file
first_file <- file.path(feature_count_files, paste(Sample_info$samplename[1],".counts.tsv",sep = "")); first_file

#prepare raw count table with first file
first_file <- read.table(file = first_file , sep="\t", header=TRUE)
raw_counts <- first_file[,c(1,8)]
colnames(raw_counts) <- c("geneID", Sample_info$samplename[1]); colnames(raw_counts)
cat("File stats: \n",Sample_info$samplename[1],": ",length(raw_counts[,Sample_info$samplename[1]])," rows and ",sum(raw_counts[,Sample_info$samplename[1]]==0)," null count(s)\n",sep="",file = "DESEQ2-log.txt", append = TRUE)

#create gene set for later
annotation <- first_file[,c(1,7)]
colnames(annotation) <- c("geneID", "gene_name")

#append all files according to Sample_info.txt
for (i in 2:length(file_paths)){
  tmp <- read.table(file.path(feature_count_files, paste(Sample_info$samplename[i],".counts.tsv",sep = "")), sep="\t", header=TRUE)
  tmp <- tmp[,c(1,8)]
  colnames(tmp) <- c("geneID", Sample_info$samplename[i])
  if (any(duplicated(tmp$GeneID))){
    stop("Duplicated feature names in ", file_paths[i], ": ", 
         paste(unique(tmp$GeneID[duplicated(tmp$GeneId)]), collapse=", "))
  }
  raw_counts <- merge(raw_counts, tmp, by="geneID", all=TRUE)
  cat(Sample_info$samplename[i],": ",length(tmp[,Sample_info$samplename[i]])," rows and ",sum(tmp[,Sample_info$samplename[i]]==0)," null count(s)\n",sep="",file = "DESEQ2-log.txt", append = TRUE)
}

#genereate a count matrix: one column per sample, one row per feature
counts <- as.matrix(raw_counts[,-1])
rownames(counts) <- raw_counts[,1]
counts <- counts[order(rownames(counts)),]
#rownames(counts)

# check that input counts are integers to fit edgeR and DESeq2 requirements
if (any(counts %% 1 != 0)) print("Input counts are not integer values as required by DESeq2")


#==========================================================================================#
#=======================================DESEQ2=============================================#
#==========================================================================================#

library(DESeq2)

#create design matrix from Sample_info
col_data <- as.data.frame(Sample_info)
rownames(col_data) <-Sample_info[,1]
col_data$condition <- factor(col_data$condition)               

#generate DESEQ2 object (rowData = raw_counts[,1] wirklich notwendig?)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
#counts(dds)
#colData(dds)
#mcols(dds)

#adding gene names to DESEQ2 Dataset 

#checking if gene names in the annotation is in the same order as in the DESEQ2 Dataset
all(rownames(dds) == annotation$geneID)

#reorder the annotation according to dds
annotation <- annotation[match(rownames(dds), annotation$geneID),]

# !!!check before adding -> Should be ZERO!!!!
sum(is.na(annotation$Geneid)) == 0
sum(duplicated(annotation$Geneid)) == 0

#adding new gene annotation to dds metadata
mcols(dds)$gene_name <- annotation[,2]


# pre-filtering low count genes, not necessary but reduces file size and processing time
# maybe introduces a skewed p-value distribution?!
#dds <- dds[ rowMean(counts(dds)) > 10, ]; dds

# versuch um low counts zu entfernen
dds <- dds[ rowMeans(counts(dds)) > 4, ]; dds


#perform differential gene expression analysis (including normalisation, estimating dispersion and statistical testing)
dds <- DESeq(dds)

#see all comparisons
resultsNames(dds)

#gene expression table
#at this step independent filtering is applied by default to remove low count genes
res <- results(dds, name=resultsNames(dds)[2], alpha = 0.05,)
res$gene_name <- mcols(dds)$gene_name
res_nona <- na.omit(res)
MAplot.results <- plotMA(res_nona, ylim=c(-7,7))
cat("\nResult summary of DESEQ2 output: ",file = "DESEQ2-log.txt", append = TRUE)
capture.output(summary(res_nona),file = "DESEQ2-log.txt", append = TRUE)

#shrinking of LFC estimates for visualisation and ranking of genes
lfcs_res <- lfcShrink(dds, res = res, contrast = "condition", type = "ashr")
plotMA(lfcs_res, ylim=c(-7,7))
summary(lfcs_res)

#identifying dots in the MA plot
#idx <- identify(lfcs_res$baseMean, lfcs_res$log2FoldChange)
#rownames( lfcs_res)[idx]

#Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method)
res <- res[order(res$padj),]
lfcs_res[order(lfcs_res$padj),]

colnames(res)
#res$gene_name
write.csv(as.data.frame(res), file=paste(resultsNames(dds)[2],"_results.csv"))
#-------------------------------------
#------------VISUALISATION------------
#-------------------------------------

#BiocManager::install("regionReport")
#BiocManager::install("pheatmap")
library("regionReport")

#setwd("./reports")
#dir.create("DESeq2Report", showWarnings = FALSE, recursive = TRUE)
report <- DESeq2Report(dds,output = "DESEQ-Report",project = "DESEQ2-Exploration", intgroup = "condition", browse = TRUE)



