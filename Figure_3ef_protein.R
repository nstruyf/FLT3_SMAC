# Figure 3e-f: Linear regression and visualisation of midostaurin response in FLT3 mutant AML patients in relation to protein expression 

# load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# load sDSS data and annotation files
setwd()
protein_data = read.table("protein.txt", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "",sep = "\t")
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE) 

# filter out NA values from protein table
protein_data <- protein_data[complete.cases(protein_data), ]
rownames(protein_data) <- protein_data$Name
protein_data$Name <- NULL  
protein_data <- t(protein_data)
protein_data <- as.data.frame(protein_data)

# load midostuarin sDSS and add to table
drug_scores <- annotation %>%
  filter(Label %in% c("responder_mut", "non_responder_mut")) %>%
  select(SampleID = ID, Midostaurin)
drug_scores <- drug_scores[drug_scores$SampleID %in% rownames(protein_data), ]
drug_scores <- drug_scores[match(rownames(protein_data), drug_scores$SampleID), ]
protein_data$sdss <- drug_scores$Midostaurin

# initlialize empty lists and vectors to store results
regression_coeficients <- list()
r_squared_values <- numeric()
p_values <- numeric()
spearman_correlations <- numeric()
spearman_p_values <- numeric()
protein_names <- numeric()

# loop through each protein
for (protein_name in colnames(protein_data)) {
  formula <- as.formula(paste("sdss ~", paste0("`", protein_name, "`")))
  model <- lm(formula, data = protein_data)
  r_squared <- summary(model)$r.squared
  f <- summary(model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  spearman_result <- cor.test(protein_data[[protein_name]], protein_data$sdss, method = "spearman")
  spearman_corr <- spearman_result$estimate
  spearman_p <- spearman_result$p.value
  protein_names <- c(protein_names, protein_name)
  p_values = c(p_values, p)
  r_squared_values <- c(r_squared_values, r_squared)
  spearman_correlations <- c(spearman_correlations, spearman_corr)
  spearman_p_values <- c(spearman_p_values, spearman_p)
}

# combine results in a single table
results_table_protein <- data.frame(
  Protein = protein_names,
  R_squared = r_squared_values,
  P_value = p_values,
  Spearman_correlation = spearman_correlations,
  Spearman_p_value = spearman_p_values)

# select top 25 positive and negative proteins for heatmap
results_table_protein$log10_p_value <- -log10(results_table_protein$Spearman_p_value)
alpha <- 0.05
sig_results <- results_table_protein %>% 
  dplyr::filter(Spearman_p_value < alpha)
top_n <- 25 
top_pos <- sig_results %>%
  plyr::arrange(plyr::desc(Spearman_correlation)) %>%
  dplyr::slice(1:top_n)
top_neg <- sig_results %>%
  plyr::arrange(Spearman_correlation) %>%
  dplyr::slice(1:top_n)
top_proteins <- bind_rows(top_pos, top_neg)
expr <- protein_data
expr <- na.omit(expr)
rownames(expr) <- expr$Name
expr$Name <- NULL  
expr_sub <- expr[rownames(expr) %in% top_proteins$Protein, ]

# prepare heatmap annotation and parameters
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

# generate and draw heatmap
HM <- Heatmap(
  expr_scaled,
  name = "Z-score",
  top_annotation = colAnn2,
  cluster_rows = TRUE,
  cluster_columns = FALSE,   
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = colorRamp2(c(-2,-1,0,1,2), c("#053061","#4393C3","#F7F7F7", "#D6604D", "#67001F")),
  row_title = "Top correlated proteins",
  column_title = "",
  heatmap_legend_param = list(legend_direction = "vertical",
                              title_gp = gpar(fontsize = 12,fontface = "bold"),
                              row_names_gp = gpar(fontsize = 12),
                              labels_gp = gpar(fontsize = 12)))
draw(HM)
HM_drawn <- draw(HM, merge_legend = TRUE)
png("fig_3e.png", width = 14, height = 22, units = "cm", res = 300)
HM_drawn
dev.off()

# select the top significant proteins from the results table
top_proteins <- results_table_protein %>%
  arrange(P_value) %>%
  head(450)

# create a ranked results table for all proteins
results_table_protein$log10_p_value <- -log10(results_table_protein$P_value)
results_table_protein$rank_last <- rank(results_table_protein$Spearman_correlation, ties.method = "first")
results_table_protein_ranked <- results_table_protein[order(results_table_protein$rank),]

# choose proteins to display on plot
proteins <- c("CCDC88A","PML","TNFAIP8","KMT2A","BRAF",'AKT1','AKT2','DIABLO','TET2',
              'MAP2K1','MAP2K2','MAP2K3','MAP2K4','PIK3R4', "ATF2", "PYCR2", "BCL6", 
              "NFKB1","FAM3A","TACO1","CCDC90B","NUP58","MEIS1","EI24")

# create rankplot
cutoff <- abs(results_table_protein_ranked$Spearman_correlation) > 0.6
ggplot(results_table_protein_ranked, aes(x = rank_last, y = Spearman_correlation, color = Spearman_correlation)) +
  geom_point(size = 2, alpha = 1) +
  geom_hline(yintercept=c(0), col="black", linetype = 'longdash') + 
  geom_text_repel(data = subset(results_table_protein_ranked, Protein %in% proteins), aes(label = as.character(Protein)), 
                  point.padding = 1, 
                  color = "black",
                  max.overlaps = 60,
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
  labs(y = "Spearman 's rank \ncorrelation coefficient", x = "Protein rank")

ggsave("fig_3f.png", height = 3, width = 5)

# choose protein to visualize in XY plot
chosen_protein <- "TET2"  # replace with protein of choice
chosen_protein_data <- protein_data[, chosen_protein]

# calculate rho and p value
spearman_test <- cor.test(protein_data[[chosen_protein]], protein_data$sdss, method = "spearman")
rho <- spearman_test$estimate  
spearman_p_value <- spearman_test$p.value  

# create XY plot
ggplot(protein_data, aes(x = chosen_protein_data, y = sdss)) +
  geom_point(size = 3, color = "#094334") +
  theme_minimal() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = chosen_protein, y = "Midostaurin sDSS") +
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
        legend.position = "none")  +
  ylim(-5,15) +
  annotate("text", x = Inf, y = Inf, 
           label = bquote(atop(rho == .(round(rho, 4)), italic(p) == .(format(spearman_p_value, scientific = TRUE)))),
           hjust = 1.05, vjust = 4.7)

ggsave(paste0("fig_3f_", chosen_protein,".png" ), height = 2.5, width = 2.5)

# create sorted gene list for GSEA
gene_list <- results_table_protein$Spearman_correlation
names(gene_list) <- results_table_protein$Protein 
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
        legend.position = "right") +
  xlab(expression(-log[10]~(adjusted~italic(p)~value)))

ggsave("supp_fig_3b.png", width = 6.25, height = 4)
