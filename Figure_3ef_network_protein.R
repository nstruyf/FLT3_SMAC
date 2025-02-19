
# Figure 3e-f: Linear regression and visualisation of midostaurin response in FLT3 mutant AML patients in relation to protein expression 

# load libraries

library(dplyr)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(ggrepel)


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


# select the top significant proteins from the results table

top_proteins <- results_table_protein %>%
  arrange(P_value) %>%
  head(450)


# load the STRINGdb package for protein network analysis

string_db <- STRINGdb$new(version = "11", species = 9606, score_threshold = 400, input_directory = "")


# map proteins to STRING IDs and retrieve protein-protein interactions

mapped_proteins <- string_db$map(top_proteins, "Protein", removeUnmappedRows = TRUE, takeFirst = TRUE)
interaction_network <- string_db$get_interactions(mapped_proteins$STRING_id)


# convert interaction data into a graph object

interaction_graph <- graph_from_data_frame(d = interaction_network, directed = FALSE)


# compute centrality measures for each protein in the network

centrality_measures <- data.frame(
  protein = V(interaction_graph)$name,
  degree = degree(interaction_graph),
  betweenness = betweenness(interaction_graph),
  closeness = closeness(interaction_graph))


# merge mapped protein information and create a named vector to associate STRING IDs with original protein names

merged_results <- merge(mapped_proteins, centrality_measures, by.x = "STRING_id", by.y = "protein")
protein_names <- setNames(mapped_proteins$Protein, mapped_proteins$STRING_id)
V(interaction_graph)$Protein <- protein_names[V(interaction_graph)$name]


# assign colors and set color mapping parameters

set.seed(42)
color_palette <- colorRampPalette(c("#BA6E9B", "#EDDBE4","#F1F1F1","#CCEBED","#4DB5BC"))(100)
binned_rsq <- cut(top_proteins$Spearman_correlation, breaks = 100, labels = FALSE)
vertex_colors <- color_palette[binned_rsq]
protein_colors <- setNames(vertex_colors, top_proteins$Protein)
V(interaction_graph)$color <- protein_colors[V(interaction_graph)$Protein]


# function to generate graph

generate_layout <- function(interaction_graph){
  set.seed(2)
  layout <- layout_with_fr(interaction_graph)
  return(layout)
}


# compute layout 

layout <- layout_with_fr(interaction_graph)


# create network visualization without labels

output_file <- "fig_3e.png"
png(filename = output_file, width = 30, height = 30, units = "cm", res = 300, bg = "transparent")
plot(interaction_graph, 
     layout = layout,
     vertex.size = 3.5, 
     vertex.size2 = 2.5,
     vertex.shape = "rectangle",
     vertex.frame.width=0.5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.family = "arial",
     vertex.label.cex = 0.6, 
     edge.arrow.size = 0.5,
     edge.color = "grey50",
     edge.width = 0.1,
     edge.curved = 0.4,
     asp = 0.6)
dev.off()


# create network visualization with labels

output_file <- "fig_3e2.png"
png(filename = output_file, width = 50, height = 30, units = "cm", res = 600)
plot(0, type="n", ann = FALSE, axes = FALSE, xlim=extendrange(layout[,1]), ylim=extendrange(layout[,2]))
plot(interaction_graph, 
     layout = layout,
     rescale = FALSE,
     add = TRUE,
     vertex.size = (strwidth(V(interaction_graph)$Protein) + strwidth("oo")) * 25, 
     vertex.size2 = strheight("I") * 2 * 25,
     vertex.shape = "rectangle",
     vertex.frame.width=0.2,
     vertex.label = V(interaction_graph)$Protein, 
     vertex.label.color = "black",
     vertex.label.cex = 0.25, 
     edge.arrow.size = 0.3,
     edge.color = "grey80",
     edge.curved = 0.4,
     asp = 0.2)
dev.off()


# create a ranked results table for all proteins

results_table_protein$log10_p_value <- -log10(results_table_protein$P_value)
results_table_protein$rank_last <- rank(results_table_protein$Spearman_correlation, ties.method = "first")
results_table_protein_ranked <- results_table_protein[order(results_table_protein$rank),]


# choose proteins to display on plot

proteins <- c("CCDC88A","PML","TNFAIP8","KMT2A","BRAF",'AKT1','AKT2','DIABLO','TET2','MAP2K1','MAP2K2','MAP2K3','MAP2K4','PIK3R4', "ATF2", "PYCR2", "BCL6", "NFKB1",
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
