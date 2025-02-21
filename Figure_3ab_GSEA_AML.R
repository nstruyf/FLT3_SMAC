
# Figure 3a-b: GSEA using AML cell type-specific genesets comparing FLT3 mutant midostaurin responders and non-responders

# load libraries

library(dplyr)
library(clusterProfiler)
library(ggridges)
library(ggplot2)
library(readxl)


# load protein data and annotation files

setwd()
protein_data = read.table("protein.txt", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "",sep = "\t")
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE) 


# filter out NA values from protein table

data <- na.omit(protein_data)
rownames(data) <- data$Name
data$Name <- NULL  


# assign patients to respective groups and perform multiple t tests

responder <- annotation$ID[annotation$Label == "responder_mut"]
non_responder <- annotation$ID[annotation$Label == "non_responder_mut"]

t_test <- function(Name) {
  data_responder <- as.numeric(na.omit(data[Name, colnames(data) %in% responder]))
  data_non_responder <- as.numeric(na.omit(data[Name, colnames(data) %in% non_responder]))
  p_value <- t.test(data_responder, data_non_responder)$p.value
  mean1 <- mean(data_responder, na.rm = TRUE)
  mean2 <- mean(data_non_responder, na.rm = TRUE)
  
  data.frame(p = p_value, mean1, mean2)
}


# merge t test results into a single table

result <- do.call(rbind, lapply(rownames(data), t_test)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    diff = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q),
    Protein = rownames(data)
  ) %>%
  arrange(p) %>%
  as.data.frame()


# create sorted gene list for GSEA

gene_list <- result$diff
names(gene_list) <- result$Protein 
gene_list <- sort(gene_list, decreasing = T) 


# load gene set and run GSEA. Genesets are available at https://github.com/andygxzeng/AMLHierarchies and https://doi.org/10.1016/j.cell.2019.01.031 (Table S3 - Tumor-derived, combined)

gene_list_LPSC<- read.gmt("AMLCellType_Genesets.gmt")
LSPC_GSEA_res <- GSEA(geneList = gene_list,
                      pvalueCutoff = 10,
                      pAdjustMethod = 'fdr',
                      TERM2GENE = gene_list_LPSC[c('term', 'gene')])

gene_list_VG<- read.table("vangalen_set3.txt", header = TRUE, row.names = 1, sep = "\t")
VG_GSEA_res <- GSEA(geneList = gene_list,
                    pvalueCutoff = 10,
                    pAdjustMethod = 'fdr',
                    TERM2GENE = gene_list_VG[c('gs_name', 'gene_symbol')],
                    minGSSize = 1)


# filter and organize GSEA results

LSPC_GSEA_df <- LSPC_GSEA_res@result
filt_LSPC_GSEA_df <- LSPC_GSEA_df %>% dplyr::filter(rank > 40)
filt_LPSC_GSEA_res <- LSPC_GSEA_res
filt_LPSC_GSEA_res@result <- filt_LSPC_GSEA_df


# plot data

options(enrichplot.colours = c("#B2182B","#2166AC"))

ridgeplot(filt_LPSC_GSEA_res) +
  labs(x = "Enrichment distribution", fill = expression("Adjusted" ~ italic(p))) + 
  geom_vline(xintercept=0, col="black") + 
  theme_bw() +
  scale_fill_gradient(low = "#B2182B", high = "#2166AC", limits = c(0.0000000001, 0.05)) +
  theme(text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12))

ggsave(file.path(folder_path,"fig_3a.png"), ridgeplot_LPSC, width = 5.5, height = 5.5)



ridgeplot(VG_GSEA_res) +
  labs(x = "Enrichment distribution") + 
  theme_bw() +
  geom_vline(xintercept=0, col="black") + 
  scale_fill_gradient(low = "#B2182B", high = "#2166AC", limits = c(0.0000000001, 0.05)) +
  xlim(-1,1) + 
  theme(text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none")

ggsave(file.path(folder_path,"fig_3b.png"), ridgeplot_VG, width = 3.3, height = 2.2)
