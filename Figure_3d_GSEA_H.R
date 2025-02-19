
# Figure 3d: GSEA using MsigDB Hallmarks genesets comparing FLT3 mutant midostaurin responders and non-responders

# load libraries
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
library(DESeq2)
library(stringr)
library(plyr)

# load protein/RNA data and annotation files

setwd()
protein_data = read.table("protein.txt", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "",sep = "\t")
count <- read.table("RNA.txt", header = TRUE)
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

result_p <- do.call(rbind, lapply(rownames(data), t_test)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    diff = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q),
    Protein = rownames(data)
  ) %>%
  arrange(p) %>%
  as.data.frame()


# Process RNAseq data

colData <- annotation %>%
  filter(Label %in% c("responder_mut", "non_responder_mut")) %>%
  mutate(Group = ifelse(Label == "responder_mut", "1", "2")) %>%
  select(SampleID = ID, Group)
count$Name <- make.unique(count$Name, sep = "_")
rownames(count) <- count$Name
count <- round(count[, -1])
colData <- colData[colData$SampleID %in% colnames(count), ]
colData <- colData[match(colnames(count), colData$SampleID), ]

dds <- DESeqDataSetFromMatrix(countData = count, colData = colData, design = ~ Group)
dds <- DESeq(dds)
result_r <- results(dds, alpha = 0.1)


# create sorted gene list for GSEA

gene_list_p <- result_p$diff
names(gene_list_p) <- result_p$Protein 
gene_list_p <- sort(gene_list_p, decreasing = T) 


gene_list_r <- result_r$log2FoldChange 
names(gene_list_r) <- rownames(result_r) 
gene_list_r <- sort(gene_list_r, decreasing = T)


# load gene set and run GSEA. 

h_gene_sets <- msigdbr(species = 'Homo sapiens', category = 'H') 

msigdb_GSEA_res_p <- GSEA(geneList = gene_list_p,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        TERM2GENE = h_gene_sets[c('gs_name', 'gene_symbol')],
                        minGSSize = 10)



msigdb_GSEA_res_r <- GSEA(geneList = gene_list_r,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'fdr',
                        TERM2GENE = h_gene_sets[c('gs_name', 'gene_symbol')])


# Filter and annotate GSEA results table

GSEA_df_p <- msigdb_GSEA_res_p@result
GSEA_df_p <- GSEA_df_p %>%
  mutate(
    count = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = (count / setSize) * 100,
    ID = str_to_title(str_replace_all(str_remove(ID, "HALLMARK_"), "_", " ")),
    status = ifelse(NES > 0, "Upregulated", "Downregulated")
  ) %>%
  arrange (-p.adjust)
GSEA_df_p$ID <- factor(GSEA_df_p$ID, levels = GSEA_df_p$ID)
new_id <- c("Uv Response Dn" = "UV Response down",
            "Il2 Stat5 Signaling" = "IL2 STAT5 Signaling",
            "Kras Signaling Dn" = "KRAS Signaling down",
            "Myc Targets V2" = "MYC Targets V2")
GSEA_df_p$ID <- revalue(GSEA_df_p$ID, new_id)


GSEA_df_r <- msigdb_GSEA_res_r@result
GSEA_df_r <- GSEA_df_r %>%
  mutate(
    count = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = (count / setSize) * 100,
    ID = str_to_title(str_replace_all(str_remove(ID, "HALLMARK_"), "_", " ")),
    status = ifelse(NES > 0, "Upregulated", "Downregulated")
  ) %>%
  arrange (-p.adjust)
GSEA_df_r$ID <- factor(GSEA_df_r$ID, levels = GSEA_df_r$ID)
new_id <- c("E2f Targets" = "E2F Targets",
            "Myc Targets V1" = "MYC Targets V1",
            "G2m Checkpoint" = "G2M Checkpoint")
GSEA_df_r$ID <- revalue(GSEA_df_r$ID, new_id)


# plot data

ggplot(GSEA_df_p, aes(x = -log10(p.adjust), y = ID, size = GeneRatio, color = status)) + 
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0, 10)) +
  scale_color_manual(values = c("Upregulated" = "#870052", "Downregulated" = "#4DB5BC"),
                     labels = c("Upregulated" = "Upregulated in \nnon-responders", "Downregulated" = "Downregulated in \nnon-responders"),
                     guide = "none") +
  labs(size = "GeneRatio", color = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(),
        text = element_text(family = "sans", size = 12, color = "black"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = c(1.3, 0.675),
        plot.margin = margin(5,100,5,5)) +
  xlab(expression(-log[10]~(adjusted~italic(p)~value)))

ggsave("fig_3d_protein.png", width = 6.1, height = 3.8)



dotplot_r <- ggplot(GSEA_df_r, aes(x = -log10(p.adjust), y = ID, size = GeneRatio, color = status)) + 
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0, 10)) +
  scale_color_manual(values = c("Upregulated" = "#870052", "Downregulated" = "#4DB5BC"),
                     labels = c("Upregulated" = "Upregulated in \nnon-responders", "Downregulated" = "Downregulated in \nnon-responders"),
                     guide = "none") +
  labs(size = "GeneRatio", color = element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(),
        text = element_text(family = "sans", size = 12, color = "black"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  xlab(expression(-log[10]~(adjusted~italic(p)~value)))
dotplot_r

ggsave("fig_3d_rna.png", width = 3.5, height = 2)

