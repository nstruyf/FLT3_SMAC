# Figure 3g: Linear regression and visualisation of midostaurin response in FLT3 mutant AML patients in relation to gene expression

# load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# load sDSS data and annotation files
setwd()
RNA_data <- read.table("/Users/nonastruyf/Projects/FLT3_project/RNA/FLT3_RNA_new.txt", header = TRUE)
annotation <- read.csv("/Users/nonastruyf/Projects/FLT3_project/data_sharing/files/annotation.csv", stringsAsFactors = FALSE) 

# filter out NA values from RNA table
RNA_data <- na.omit(RNA_data)
rownames(RNA_data) <- RNA_data$Name
RNA_data$Name <- NULL  
RNA_data <- as.matrix(RNA_data)
RNA_data <- t(RNA_data)
RNA_data <- as.data.frame(RNA_data)

# load midostuarin sDSS and add to table
drug_scores <- annotation %>%
  filter(Label %in% c("responder_mut", "non_responder_mut")) %>%
  select(SampleID = ID, Midostaurin)
drug_scores <- drug_scores[drug_scores$SampleID %in% rownames(RNA_data), ]
drug_scores <- drug_scores[match(rownames(RNA_data), drug_scores$SampleID), ]
RNA_data$sdss <- drug_scores$Midostaurin

# initlialize empty lists and vectors to store results
regression_coeficients <- list()
r_squared_values <- numeric()
p_values <- numeric()
spearman_correlations <- numeric()
spearman_p_values <- numeric()
gene_names <- numeric()

# loop through each gene
for (gene_name in colnames(RNA_data)) {
  if (all(RNA_data[, gene_name] == 0)) {
    next
  }
  formula <- as.formula(paste("sdss ~", paste0("`", gene_name, "`")))
  model <- lm(formula, data = RNA_data)
  r_squared <- summary(model)$r.squared
  f <- summary(model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  spearman_result <- cor.test(RNA_data[[gene_name]], RNA_data$sdss, method = "spearman")
  spearman_corr <- spearman_result$estimate
  spearman_p <- spearman_result$p.value
  gene_names <- c(gene_names, gene_name)
  p_values = c(p_values, p)
  r_squared_values <- c(r_squared_values, r_squared)
  spearman_correlations <- c(spearman_correlations, spearman_corr)
  spearman_p_values <- c(spearman_p_values, spearman_p)
}

# combine results in a single table
results_table_RNA <- data.frame(
  Gene = gene_names,
  R_squared = r_squared_values,
  P_value = p_values,
  Spearman = spearman_correlations,
  Spearman_p_value = spearman_p_values)

# select the top genes from the results table
top_genes <- results_table_RNA %>%
  arrange(Spearman_p_value) %>%
  head(1500)

# create a ranked results table for all genes
results_table_RNA$log10_p_value <- -log10(results_table_RNA$P_value)
results_table_RNA$rank_last <- rank(results_table_RNA$Spearman, ties.method = "first")
results_table_RNA_ranked <- results_table_RNA[order(results_table_RNA$rank),]

# choose genes to display on plot
genes <- c("CD28","CD14",'ITGAX','IL10RB','XIAP',"TET2","BIRC2", "BIRC3", "TLR6", "CD164", 
           "BCL2L1", "PIR", "DIABLO", "PIR", "DIABLO", "MIR155HG", "CSF2", "IL7")

# create rankplot
cutoff <- abs(results_table_RNA_ranked$Spearman) > 0.6
ggplot(results_table_RNA_ranked, aes(x = rank_last, y = Spearman, color = Spearman)) +
  geom_point(size = 2, alpha = 1) +
  geom_hline(yintercept=c(0), col="black", linetype = 'longdash') + 
  geom_text_repel(data = subset(results_table_RNA_ranked, Gene %in% genes), aes(label = as.character(Gene)), 
                  point.padding = 1.5,
                  box.padding = 0.3,
                  color = "black",
                  max.overlaps = 50,
                  segment.color = "black",
                  segment.alpha = 0.6) +
  scale_color_gradient2(low = "#870052", mid = "#F1F1F1", high = "#4DB5BC") +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.x = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        legend.text = element_text(family = "Arial", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        legend.position = "none") + 
  labs(y = "Spearman 's rank \ncorrelation coefficient", x = "Gene rank")

ggsave("fig_3g.png", height = 3, width = 5)

# choose gene to visualize in XY plot
chosen_gene <- "TET2" # replace with gene of choice
chosen_gene_data <- RNA_data[, chosen_gene]

# calculate rho and p value
spearman_test <- cor.test(RNA_data[[chosen_gene]], RNA_data$sdss, method = "spearman")
rho <- spearman_test$estimate  
spearman_p_value <- spearman_test$p.value  

# create XY plot
ggplot(RNA_data, aes(x = chosen_gene_data, y = sdss)) +
  geom_point(size = 3, color = "#FFC66D") +
  theme_minimal(base_family = "Arial") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = chosen_gene, y = "Midostaurin sDSS") +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.x = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
        legend.text = element_text(family = "Arial", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")  +
  ylim(-5,15) +
  annotate("text", x = Inf, y = Inf, 
           label = bquote(atop(rho == .(round(rho, 4)), italic(p) == .(format(spearman_p_value, scientific = TRUE)))),
           hjust = 1.05, vjust = 4.7)

ggsave(paste0("fig_3g_", chosen_gene,".png" ), height = 2.5, width = 2.5)

# create sorted gene list for GSEA
gene_list <- results_table_RNA$Spearman
names(gene_list) <- results_table_RNA$Gene 
gene_list <- sort(gene_list, decreasing = T) 

# load gene set and run GSEA. 
h_gene_sets <- msigdbr(species = 'Homo sapiens', category = 'H')
msigdb_GSEA_res <- GSEA(geneList = gene_list,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        TERM2GENE = h_gene_sets[c('gs_name', 'gene_symbol')],
                        minGSSize = 1)

# Filter and annotate GSEA results table
GSEA_df <- msigdb_GSEA_res@result
GSEA_df <- GSEA_df %>%
  mutate(
    count = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = (count / setSize) * 100,
    ID = str_to_title(str_replace_all(str_remove(ID, "HALLMARK_"), "_", " ")),
    status = ifelse(NES > 0, "Upregulated", "Downregulated")
  ) %>%
  arrange (-p.adjust)
GSEA_df$ID <- factor(GSEA_df$ID, levels = GSEA_df$ID)

# plot data
ggplot(GSEA_df, aes(x = -log10(p.adjust), y = ID, size = GeneRatio, color = status)) + 
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0, 10)) +
  scale_color_manual(values = c("Upregulated" = "#4DB5BC", "Downregulated" = "#870052"),
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
        legend.position = "none",
        plot.margin = margin(5,100,5,5)) +
  xlab(expression(-log[10]~(adjusted~italic(p)~value)))

ggsave("supp_fig_3c.png", width = 6.5, height = 3.8)

# heatmap of top 40 most significantly correlated genes (Supplementary fig 3a)
# sort and filter top positive and negative associations
results_table_RNA$log10_p_value <- -log10(results_table_RNA$P_value)
alpha <- 0.05
sig_results <- results_table_RNA %>% 
  dplyr::filter(Spearman_p_value < alpha)
top_n <- 40 
top_pos <- sig_results %>%
  plyr::arrange(plyr::desc(Spearman)) %>%
  slice(1:top_n)
top_neg <- sig_results %>%
  plyr::arrange(Spearman) %>%
  slice(1:top_n)
top_genes <- bind_rows(top_pos, top_neg)
expr <- RNA_data
expr <- na.omit(expr)
rownames(expr) <- expr$Name
expr$Name <- NULL  
expr_sub <- expr[rownames(expr) %in% top_genes$Gene, ]

# add sDSS annotation and set heatmap parameters
annotation_sorted <- annotation %>% plyr::arrange(Midostaurin)
annotation_filtered <- annotation_sorted[annotation_sorted$ID %in% colnames(expr_sub), ]
expr_sub <- expr_sub[, annotation_filtered$ID]
expr_scaled <- t(scale(t(expr_sub)))
annotation_hm2 <- data.frame(Midostaurin = annotation_filtered$Midostaurin)
dot_colors <- ifelse(annotation_hm2$Midostaurin > 3.46,
                     "#4DB5BC",  # high
                     "#870052")  # low
ylim_fixed <- c(0, 15)
axis_ticks <- seq(0, 15, by = 5)
colAnn2 <- HeatmapAnnotation(
  Midostaurin_Dots = anno_points(
    annotation_hm2$Midostaurin,
    size = unit(3, "mm"),
    gp = gpar(col = dot_colors),
    ylim = ylim_fixed,
    axis_param = list(
      at = axis_ticks,
      labels = axis_ticks,
      gp = gpar(fontsize = 12)
    )
  ),
  which = "col",
  annotation_height = unit(2, "cm"),
  annotation_width = unit(1, "cm"),
  show_annotation_name = FALSE
)

# draw and save heatmap
HM <- Heatmap(
  expr_scaled,
  name = "Z-score",
  top_annotation = colAnn2,
  cluster_rows = TRUE,
  cluster_columns = FALSE,   
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = colorRamp2(c(-2,-1,0,1,2), c("#053061","#4393C3","#F7F7F7", "#D6604D", "#67001F")),
  row_title = "Top correlated genes",
  column_title = "",
  heatmap_legend_param = list(legend_direction = "vertical",
                              title_gp = gpar(fontsize = 12,fontface = "bold"),
                              row_names_gp = gpar(fontsize = 12),
                              labels_gp = gpar(fontsize = 12)))
draw(HM)
HM_drawn <- draw(HM, merge_legend = TRUE)
png("supp_fig_3a.png", width = 16, height = 30, units = "cm", res = 300)
HM_drawn
dev.off()
