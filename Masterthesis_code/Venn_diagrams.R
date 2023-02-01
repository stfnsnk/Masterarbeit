
## Code for creating the Venn Diagrams 
## Author: Stefan Senk - BIOINF22/23
## Masterthesis

ONT_file_path <- "./ONT_results/condition_wildtype_vs_knockout _results.csv"
Illumina_file_path <- "./Illumina_results/condition_wildtype_vs_knockout _results.csv"

ONT_results <- read.delim(ONT_file_path, header = TRUE, sep = "," )
ONT_results_padj <- ONT_results[order(ONT_results$padj),]

#ONT_results_pval <- ONT_results[order(ONT_results$pvalue),]

Illumina_results <- read.delim(Illumina_file_path, header = TRUE, sep = "," )
Illumina_results_padj <- Illumina_results[order(Illumina_results$padj),]

Illumina_nona <- na.omit(Illumina_results)
ONT_nona <- na.omit(ONT_results)

#BiocManager::install("ggvenn")
# load ggvenn package
library("ggvenn")

##Venn diagram for all genes
LIST_allgenes <- list('ONT'=ONT_nona$X,'Illumina'=Illumina_nona$X)
ggvenn(LIST_allgenes)


`####Venn diagram for significant genes < 0.05
Illumina_sig_nona <- Illumina_nona[Illumina_nona$padj<0.05,]
ONT_sig_nona <- ONT_nona[ONT_nona$padj<0.05,]
LIST_sig_nona <- list('ONT'=ONT_sig_nona$gene_name,'Illumina'=Illumina_sig_nona$gene_name)
ggvenn(LIST_sig_nona)
