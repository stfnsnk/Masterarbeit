#!/usr/bin/env R

## RNA-seq analysis using DESEQ2
# code template from 
# - Stephen Turner (https://gist.github.com/stephenturner/f60c1934405c127f09a6)
# - Hugo Varet (https://github.com/PF2-pasteur-fr/SARTools/blob/master/R/loadCountData.R)

args <- commandArgs(TRUE)
sample_info_file <- args[1]
feature_count_files <- args[2]

#==========================================================================================#
#==============================reading in count data=======================================#
#==========================================================================================#

Sample_info <- read.table(sample_info_file, header = TRUE); Sample_info
Sample_info <- Sample_info[order(Sample_info$condition),]; Sample_info

# find all count files in featureCounts output directory
file_paths <- list.files(path = feature_count_files, pattern = "*.tsv$", full.names = TRUE); file_paths

# find the first file of the Sample info file
first_file <- file.path(feature_count_files, paste(Sample_info$samplename[1],".counts.tsv",sep = "")); first_file

#prepare raw count table with first file
raw_counts <- read.table(file = first_file , sep="\t", header=TRUE)
raw_counts <- raw_counts[,c(1,7)]
colnames(raw_counts) <- c("geneID", Sample_info$samplename[1]); colnames(raw_counts)

#append all files according to Sample_Info.txt
for (i in 2:length(file_paths)){
  tmp <- read.table(file.path(feature_count_files, paste(Sample_info$samplename[i],".counts.tsv",sep = "")), sep="\t", header=TRUE)
  tmp <- tmp[,c(1, 7)]
  colnames(tmp) <- c("geneID", Sample_info$samplename[i])
  if (any(duplicated(tmp$GeneID))){
    stop("Duplicated feature names in ", file_paths[i], ": ", 
         paste(unique(tmp$GeneID[duplicated(tmp$GeneId)]), collapse=", "))
  }
  raw_counts <- merge(raw_counts, tmp, by="geneID", all=TRUE)
  cat(Sample_info$samplename[i],": ",length(tmp[,Sample_info$samplename[i]])," rows and ",sum(tmp[,Sample_info$samplename[i]]==0)," null count(s)\n",sep="")
}

#genereate a count matrix: one column per sample, one row per feature
counts <- as.matrix(raw_counts[,-1])
rownames(counts) <- raw_counts[,1]
counts <- counts[order(rownames(counts)),]


#==========================================================================================#
#=======================================DESEQ2=============================================#
#==========================================================================================#

library(DESeq2)

#create design matrix from Sample_info
col_data <- as.data.frame(Sample_info)
rownames(col_data) <- Sample_info$samplename
col_data$condition <- factor(col_data$condition)

#generate DESEQ2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)

#adding gene colum to DESEQ2 Dataset for visualisation
featureData <- data.frame(gene=rownames(counts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# pre-filtering low count genes, not necessary but reduces file size and processing time
dds <- dds[ rowMeans(counts(dds)) > 4, ]; dds

#perform differential gene expression analysis (including normalisation, estimating dispersion and statistical testing)
dds <- DESeq(dds)


#generate resulttable with given p-value cutoff (alpha)
  #Quote: "Note that the results function automatically performs independent filtering 
  #       based on the mean of normalized counts for each gene, optimizing the number 
  #       of genes which will have an adjusted p value below a given FDR cutoff, alpha (Default = 0.1"

res05 <- results(dds, alpha = 0.05)

#Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method)
res05[order(res05$padj),]
summary(res05)

#export gene expression table
write.csv(as.data.frame(res05), file=paste(resultsNames(dds)[2],"_results.csv"))

#==========================================================================================#
#================================Exploration Report========================================#
#==========================================================================================#


library("regionReport")

report <- DESeq2Report(dds,output ="DESEQ-Report", project = "DESEQ2-Exploration", intgroup = "condition", browse = FALSE)


#==========================================================================================#
#==============================interactive html plots======================================#
#==========================================================================================#

#library(Glimma)
