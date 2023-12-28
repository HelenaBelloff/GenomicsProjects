library(readxl)
library(stringr)

psoriasis <- read_excel("Downloads/Supplemental_Table_2.xlsx")

DEG_Group <- function(i) {
  if(i >= 0.00000000) group <- "up"
  if(i < 0.00000000) group <- "down"
  return(group)
}

psoriasis$group <- sapply(psoriasis$`Asterand log(ratio)`, DEG_Group)

sig_group <- function(i) {
  if(i >= 0.05) sig <- "unch"
  if(i < 0.05) sig <- "sig"
  return(sig)
}

psoriasis$sig <- sapply(psoriasis$`Asterand p-value`, sig_group)

psoriasis$group <- ifelse(psoriasis$sig == "unch", psoriasis$group == "unch", psoriasis$group)
psoriasis$group <- str_replace(psoriasis$group, "FALSE", "unch")
colnames(psoriasis)[2] <- "Gene"


# Start GSEA
library(GOtest)
library(msigdb)

options(stringsAsFactors = FALSE)

# 3 input
GOsets = c('C5.BP', 'C5.CC', 'C5.MF')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species = 'human', return.data.frame = T)

# Our sample of genes
universe = psoriasis$Gene

result = GOtest(x=psoriasis[,c('Gene', 'group')], go=gosets_genes, query.population = universe, 
                background='query', name.x = 'significant_DEG', name.go='GOsets', method = 'hypergeometric')

View(result)


# Heatmap
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyr)


result_subset <- result[1:50,]
# p.adj on -log scale
# sep df for groups
# DEG_Diabetes_Genes[,c('Gene', 'group')]

colfunc <- colorRampPalette(c("white", "orangered"))
aka3 <- list(DEG_Group = c(up = "red", down = "darkblue", unch = "grey"))

result_subset <- subset(result_subset, select  = c(GOsets, significant_DEG, Pvalue))
result_subset <- data.frame(result_subset)
result_subset <- spread(result_subset, significant_DEG, Pvalue)
row.names(result_subset) = result_subset$GOsets
result_subset[1] <- NULL


keyDriversHeat = pheatmap(mat = result_subset, color=colfunc(10),
                          cluster_cols = FALSE, cluster_rows = FALSE, 
                          show_rownames = TRUE, main = "Top 50 Significantly Enriched Pathways Using DEG For Psoriasis") 
                          annotation_row = Top20_DEG_Psoriasis_Group, annotation_colors = aka3)


test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)
