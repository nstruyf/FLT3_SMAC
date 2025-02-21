
# Figure 7e-g: GSEA using AML cell type-specific genesets comparing FLT3 wild type and mutant midostaurin responders and non-responders

# load libraries

library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ggridges)
library(ggplot2)
library(readxl)


# load protein data and annotation files

setwd()
protein_data = read.table("protein.txt", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "",sep = "\t")
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE) 


# filter out NA values from protein table

if (any(duplicated(protein_data$Name))) {
  protein_data <- protein_data[!duplicated(protein_data$Name), ]
}
protein_data_filtered <- na.omit(protein_data)
rownames(protein_data_filtered) <- protein_data_filtered$Name
protein_data_filtered$Name <- NULL  


# assign patients to respective groups and perform multiple t tests

wt <- annotation$ID[annotation$Group == "wt"]
mut <- annotation$ID[annotation$Group == "mut"]
wt_responder <- annotation$ID[annotation$Label == "responder_wt"]
wt_non_responder <- annotation$ID[annotation$Label == "non_responder_wt"]
mut_responder <- annotation$ID[annotation$Label == "responder_mut"]
mut_non_responder <- annotation$ID[annotation$Label == "non_responder_mut"]

t_test1 <- function(Name) {
  data_wt <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% wt]))
  data_mut <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% mut]))
  p_value <- t.test(data_wt, data_mut)$p.value
  mean1 <- mean(data_wt, na.rm = TRUE)
  mean2 <- mean(data_mut, na.rm = TRUE)
  data.frame(p = p_value, mean1, mean2)
}

t_test2 <- function(Name) {
  data_responder <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% wt_responder]))
  data_non_responder <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% wt_non_responder]))
  p_value <- t.test(data_responder, data_non_responder)$p.value
  mean1 <- mean(data_responder, na.rm = TRUE)
  mean2 <- mean(data_non_responder, na.rm = TRUE)
  data.frame(p = p_value, mean1, mean2)
}

t_test3 <- function(Name) {
  data_responder <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% mut_responder]))
  data_non_responder <- as.numeric(na.omit(protein_data_filtered[Name, colnames(protein_data_filtered) %in% mut_non_responder]))
  p_value <- t.test(data_responder, data_non_responder)$p.value
  mean1 <- mean(data_responder, na.rm = TRUE)
  mean2 <- mean(data_non_responder, na.rm = TRUE)
  data.frame(p = p_value, mean1, mean2)
}


# merge t test results into a single table

result1 <- do.call(rbind, lapply(rownames(protein_data_filtered), t_test1)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    diff = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q),
    Protein = rownames(protein_data_filtered)
  ) %>%
  arrange(p) %>%
  as.data.frame()

result2 <- do.call(rbind, lapply(rownames(protein_data_filtered), t_test2)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    diff = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q),
    Protein = rownames(protein_data_filtered)
  ) %>%
  arrange(p) %>%
  as.data.frame()

result3 <- do.call(rbind, lapply(rownames(protein_data_filtered), t_test3)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    diff = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q),
    Protein = rownames(protein_data_filtered)
  ) %>%
  arrange(p) %>%
  as.data.frame()


# create sorted gene list for GSEA

gene_list_1 <- result1$diff
names(gene_list_1) <- result1$Protein 
gene_list_1 <- sort(gene_list_1, decreasing = T) 

gene_list_2 <- result2$diff
names(gene_list_2) <- result2$Protein 
gene_list_2 <- sort(gene_list_1, decreasing = T) 

gene_list_3 <- result3$diff
names(gene_list_3) <- result3$Protein 
gene_list_3 <- sort(gene_list_3, decreasing = T) 


# load gene set and run GSEA. Genesets are available at https://github.com/andygxzeng/AMLHierarchies and https://doi.org/10.1016/j.cell.2019.01.031 (Table S3 - Tumor-derived, combined)

gene_list_LPSC<- read.gmt("AMLCellType_Genesets.gmt")
LSPC_GSEA_res1 <- GSEA(geneList = gene_list_1,
                      pvalueCutoff = 10,
                      pAdjustMethod = 'fdr',
                      TERM2GENE = gene_list_LPSC[c('term', 'gene')])
LSPC_GSEA_res2 <- GSEA(geneList = gene_list_2,
                      pvalueCutoff = 10,
                      pAdjustMethod = 'fdr',
                      TERM2GENE = gene_list_LPSC[c('term', 'gene')])
LSPC_GSEA_res3 <- GSEA(geneList = gene_list_3,
                      pvalueCutoff = 10,
                      pAdjustMethod = 'fdr',
                      TERM2GENE = gene_list_LPSC[c('term', 'gene')])

gene_list_VG<- read.table("vangalen_set3.txt", header = TRUE, row.names = 1, sep = "\t")
VG_GSEA_res1 <- GSEA(geneList = gene_list_1,
                    pvalueCutoff = 10,
                    pAdjustMethod = 'fdr',
                    TERM2GENE = gene_list_VG[c('gs_name', 'gene_symbol')],
                    minGSSize = 1)


# filter and organize GSEA results

LSPC_GSEA_df <- LSPC_GSEA_res1@result
filt_LSPC_GSEA_df <- LSPC_GSEA_df %>% dplyr::filter(rank > 40)
filt_LPSC_GSEA_res <- LSPC_GSEA_res1
filt_LPSC_GSEA_res@result <- filt_LSPC_GSEA_df

LSPC_GSEA_wt <- LSPC_GSEA_res2@result
LSPC_GSEA_wt <- LSPC_GSEA_wt %>% dplyr::filter(rank > 40)

LSPC_GSEA_mut <- LSPC_GSEA_res3@result
LSPC_GSEA_mut <- LSPC_GSEA_mut %>% dplyr::filter(rank > 40)


# plot wt vs mut data

options(enrichplot.colours = c("#B2182B","#2166AC"))

ridgeplot(VG_GSEA_res1) +
  labs(x = "Enrichment distribution") + 
  theme_bw() +
  geom_vline(xintercept=0, col="black") + 
  scale_fill_gradient(low = "#B2182B", high = "#2166AC", limits = c(0.0001, 0.05)) +
  xlim(-1,1) + 
  theme(text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none")

ggsave("fig_7e.png", width = 3.3, height = 2.2)

ridgeplot(filt_LPSC_GSEA_res) +
  labs(x = "Enrichment distribution", fill = expression("Adjusted" ~ italic(p))) + 
  geom_vline(xintercept=0, col="black") + 
  theme_bw() +
  scale_fill_gradient(low = "#B2182B", high = "#2166AC", limits = c(0.0001, 0.05)) +
  theme(text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12))

ggsave("fig_7f.png", width = 5.5, height = 5.5)


# merge GSEA results

LSPC_comparison1 <- LSPC_GSEA_wt %>%
  dplyr::select(ID, NES, p.adjust) %>%
  dplyr::rename(NES_comparison1 = NES, p.adjust_comparison1 = p.adjust)
LSPC_comparison2 <- LSPC_GSEA_mut %>%
  dplyr::select(ID, NES, p.adjust) %>%
  dplyr::rename(NES_comparison2 = NES, p.adjust_comparison2 = p.adjust)

merged_results <- full_join(LSPC_comparison1, LSPC_comparison2, by = "ID") 

long_data <- merged_results %>%
  pivot_longer(cols = starts_with("NES"), names_to = "Comparison", values_to = "NES") %>%
  pivot_longer(cols = starts_with("p.adjust"), names_to = "p.adjust_name", values_to = "p.adjust") %>%
  dplyr::filter(substr(Comparison, 5, 20) == substr(p.adjust_name, 10, 25)) %>%
  mutate(Comparison = case_when(
    Comparison == "NES_comparison1" ~ "FLT3wt",
    Comparison == "NES_comparison2" ~ "FLT3mut"))

long_data <- long_data %>%
  mutate(log10_p_adjust = -log10(p.adjust)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) 


# plot merged data

ggplot(long_data, aes(x = abs(NES), y = ID, size = -log10(p.adjust), color = NES > 0)) +   
  geom_point(alpha = 0.8) +
  facet_wrap(~Comparison) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(
    values = c("TRUE" = "#870052", "FALSE" = "#4DB5BC"),
    name = "",
    labels = c("TRUE" = "Upregulated in non-responders", "FALSE" = "Upregulated in responders"), guide = "none") +
  theme_bw() +
  theme(legend.title = element_text(family = "Arial", size = 12, color = "black"),
        text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.x = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        legend.text = element_text(family = "Arial", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "Normalize enrichment score (absolute value)", y = "", size = expression(-log[10]~italic(p)~value)) 

ggsave("fig_7g.png", width = 7.1, height = 5.4)
