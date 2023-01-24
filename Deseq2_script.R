#!/usr/bin/env R

## RNA-seq analysis using DESEQ2
# code inspiration from 
# Stephen Turner (https://gist.github.com/stephenturner/f60c1934405c127f09a6)
#       
# merging count data into one data.frame - @author Marie-Agnes Dillies and Hugo Varet 
#                                           (https://github.com/PF2-pasteur-fr/SARTools/blob/master/R/loadCountData.R)


args <- commandArgs(TRUE)
sample_info_file <- args[1]
feature_count_files <- args[2:length(args)]

#==========================================================================================#
#==============================reading in count data=======================================#
#==========================================================================================#

Sample_info <- read.table(sample_info_file, header = TRUE); 
Sample_info <- Sample_info[order(Sample_info$condition),]; 
paste("Sample Info File: ", sample_info_file)
Sample_info

# find the first file of the Sample info file
first_file <- paste(Sample_info$samplename[1],".counts.tsv",sep = "")

#prepare raw count table with first file
first_file <- read.table(file = first_file , sep="\t", header=TRUE)
raw_counts <- first_file[,c(1,8)]
colnames(raw_counts) <- c("geneID", Sample_info$samplename[1]); colnames(raw_counts)

#creating gene annotation for later
annotation <- first_file[,c(1,7)]
colnames(annotation) <- c("geneID", "gene_name")

#append all files according to Sample_Info.txt
for (i in 2:length(Sample_info$samplename)){
  next_file = match(paste(Sample_info$samplename[i],".counts.tsv",sep = ""), feature_count_files);next_file #find index in feature_count_files
  tmp <- read.table(file = feature_count_files[next_file], sep="\t", header=TRUE); 
  tmp <- tmp[,c(1, 8)]
  colnames(tmp) <- c("geneID", Sample_info$samplename[i])
  if (any(duplicated(tmp$GeneID))){
    stop("Duplicated feature names in ", Sample_info$samplename[i], ": ", 
         paste(unique(tmp$GeneID[duplicated(tmp$GeneId)]), collapse=", "))
  }
  raw_counts <- merge(raw_counts, tmp, by="geneID", all=TRUE)
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

#adding gene annotation to DESEQ2 results 
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
dds <- dds[ rowMeans(counts(dds)) > 4, ]; dds

#perform differential gene expression analysis (including normalisation, estimating dispersion and statistical testing)
dds <- DESeq(dds)


#generate resulttable with given p-value cutoff (alpha)
  #Quote: "Note that the results function automatically performs independent filtering 
  #       based on the mean of normalized counts for each gene, optimizing the number 
  #       of genes which will have an adjusted p value below a given FDR cutoff, alpha (Default = 0.1"

res05 <- results(dds, alpha = 0.05)
res05$gene_name <- mcols(dds)$gene_name
#Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method)
res05 <- res05[order(res05$padj),]

#export gene expression table
write.csv(as.data.frame(res05), file=paste(resultsNames(dds)[2],"_results.csv"))

#==========================================================================================#
#================================Exploration Report========================================#
#==========================================================================================#

library("regionReport")

report <- DESeq2Report(dds,output ="DESEQ-Report", project = "DESEQ2-Exploration", intgroup = "condition", browse = FALSE)


