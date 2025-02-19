
# Figure 4a-e: spatial single cell proteomics comparing FLT3 mutant midostaurin responders and non-responders. 

# load libraries

library(SeuratObject)
library(Seurat)
library(pixelatorR)
library(tibble)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)


# read .pxl files and merge them into a single Seurat object

setwd("/Pixelgen/")
data_files <-
  c(S1 = "FU-198-AML-094.dataset.pxl",
    S3 = "FU-198-AML-112-2.dataset.pxl",
    S4 = "FU-198-AML-116-1.dataset.pxl",
    S6 = "FU-198-AML-123.dataset.pxl",
    S7 = "FU-198-AML-154.dataset.pxl",
    S8 = "FU-198-AML-162.dataset.pxl")

pg_data <- lapply(data_files, ReadMPX_Seurat)
pg_data_combined <- merge (pg_data[[1]], y = pg_data[-1], add.cell.ids = names(pg_data))


# assign sample IDs as metadata and filter based on cell quality

pg_data_combined <- pg_data_combined %>%
  AddMetaData(metadata = str_remove(rownames(pg_data_combined[[]]), "_.*"), col.name = "sample") %>%
  subset(edges >= 2000 & tau_type == "normal")  


# merge all mpxCells layers into one

pg_data_combined[["mpxCells"]] <- JoinLayers(pg_data_combined[["mpxCells"]])


# identify proteins with a median expression of at least 5 while excluding QC markers

filtered_proteins <- LayerData(pg_data_combined, layer = "counts") %>%
  apply(MARGIN = 1, median) %>%
  enframe("marker", "median") %>%
  arrange(-median) %>% 
  filter(median >= 5, !marker %in% c("ACTB", "mIgG1", "mIgG2a", "mIgG2b")) 

pg_data_combined_processed <- 
  pg_data_combined %>%
  subset(features = filtered_proteins$marker)
pg_data_combined_processed <-
  pg_data_combined_processed %>%
  NormalizeData(normalization.method = "CLR",
                margin = 2)

# normalize data using CLR normalization, identify variable features, scale the data, and run PCA & UMAP for visualization

pg_data_combined_processed <- pg_data_combined %>%
  subset(features = filtered_proteins$marker) %>%
  NormalizeData(normalization.method = "CLR", margin = 2) %>%
  FindVariableFeatures(nfeatures = 67) %>%
  ScaleData(do.scale = TRUE, do.center = TRUE, verbose = FALSE) %>%
  RunPCA(npcs = 30) %>%
  RunUMAP(dims = 1:10)


# assign condition to samples based on midostaurin response

sensitive_samples <- c("S1", "S3", "S6")
pg_data_combined_processed$condition <- sapply(pg_data_combined_processed$sample, function(ita) ifelse (ita %in% sensitive_samples, "Responders", "Non-responders"))


# identify nearest neighbors based on PCA dimensions and perform clustering to group similar samples together

pg_data_combined_processed <-
  pg_data_combined_processed %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(random.seed = 1)


# run UMAP using the first 10 PCs to reduce dimensions for visualization

pg_data_combined_processed <- RunUMAP(pg_data_combined_processed, dims = 1:10, n.neighbors = 30, min.dist = 0.3)


# annotate clusters based on marker expression

DoHeatmap(subset(pg_data_combined_processed, downsample = 100), size = 2, 
          features = c("CD44","CD38","CD33","HLA-DR","CD11c","CD14","CD64","CD7","CD71","CD41","CD20","CD4", "CD8"))

cell_annotation <-
  c("0" = "Stem cell-like",
    "1" = "Myelomonocytic",
    "2" = "Myelomonocytic",
    "3" = "Myeloblastic",
    "4" = "Promonocytic",
    "5" = "T cells",
    "6" = "Monocytic, CD86+",
    "7" = "Stem cell-like",
    "8" = "Myeloblastic",
    "9" = "Monocytic, CD86-",
    "10" = "Monocytes",
    "11" = "Stem cell-like",
    "12" = "Stem cell-like",
    "13" = "B cells")

pg_data_combined_processed <- pg_data_combined_processed %>%
  AddMetaData(setNames(cell_annotation[pg_data_combined_processed$seurat_clusters],
                       nm = colnames(pg_data_combined_processed)), "cell_type")


# create annotated UMAP

plot1 <- pg_data_combined_processed %>%
  DimPlot(group.by = "cell_type", pt.size = 1,
          cols = c("#CCEBED","#CAB2D6","#B84145", "#FF876F","#F59A00", "#FFC66D","#6A3D9A","#54B986","#4DB5BC")) + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),                     
    axis.text = element_text(size = 12),                      
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),                   
    plot.title = element_text(size = 12, hjust = 0.5),        
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  NoLegend() +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)
plot1 <- LabelClusters(plot = plot1, id = "cell_type", repel = TRUE, position = "median", box.padding = 0.4, segment.size = 0)
plot1

ggsave("fig_4a.png", height = 3.2, width = 3.2)


# remove T and B cell clusters 

Idents(pg_data_combined_processed) <- "seurat_clusters"
removed_clusters <- WhichCells(pg_data_combined_processed, idents = c(5, 13))
selected_clusters <- subset(pg_data_combined_processed, cells = setdiff(Cells(pg_data_combined_processed), removed_clusters))


# re-cluster myeloid cells

selected_clusters <- selected_clusters %>%
  RunPCA(npcs = 30) %>%
  RunUMAP(dims = 1:10, n.neighbors = 30, min.dist = 0.3)


# create annotated UMAP without lymphocytes

plot2 <- selected_clusters %>%
  DimPlot(group.by = "cell_type", pt.size = 1,
          cols = c("#CAB2D6","#B84145", "#FF876F","#F59A00", "#FFC66D","#6A3D9A","#54B986")) +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),                     
    axis.text = element_text(size = 12),                      
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),                   
    plot.title = element_text(size = 12, hjust = 0.5),        
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  NoLegend() +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)
plot2 <- LabelClusters(plot = plot2, id = "cell_type", repel = TRUE, position = "median", box.padding = 0.5, segment.size = 0)
plot2

ggsave("fig_4a2.png", height = 3.2, width = 3.2)


# group samples by condition and perform differential expression analysis

Idents(selected_clusters) <- "condition"
sens_vs_res <- FindMarkers(selected_clusters, ident.1 = "Non-responders", ident.2 = "Responders", verbose = FALSE)
sens_vs_res$gene <- rownames(sens_vs_res)


# plot differential expression

cutoff <- sens_vs_res$p_val_adj < 0.05 &  abs(sens_vs_res$avg_log2FC) > 0.4

ggplot(sens_vs_res, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(size = 3, color = "grey70") +
  geom_point(data = sens_vs_res[cutoff,], size = 3, color = ifelse(sens_vs_res[cutoff,"avg_log2FC"] > 0, "#870052","#4DB5BC")) +
  geom_text_repel(aes(label = ifelse(cutoff, as.character(gene), ""))) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed') + 
  geom_vline(xintercept=c(-0.4, 0.4), col="black", linetype = 'dashed') +
  ylab(expression(-log[10]~(italic(q)~value))) + 
  scale_x_continuous(breaks = seq(-1,1, by = 1), labels = c(-1, 0, 1), limits = c(-1,1))  +
  xlab("Differential marker expression") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "Arial", size = 12,color = "black"),         
        axis.title = element_text(size = 12,color = "black"),                     
        axis.text = element_text(size = 12,color = "black"),                      
        legend.title = element_text(size = 12,color = "black"),                  
        legend.text = element_text(size = 12, color = "black"))

ggsave(file.path(folder_path,"fig_4b.png"), plot = p, height = 2.5, width = 2.8)


# create featureplots for CD20 and CD45RA

Idents(selected_clusters) <- "seurat_clusters"


fp1 <- FeaturePlot(selected_clusters, pt.size = 0.5, max.cutoff = 0.5, features = c("CD200"), reduction = "umap", coord.fixed = FALSE) +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),                     
    axis.text = element_text(size = 12),                      
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),                   
    plot.title = element_blank(),     
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(color = "Z-score", x = "UMAP 1", y = "UMAP 2") + NoLegend() +
  xlim(-8.5,13) +
  scale_color_gradientn(colors = c("#C7ECDC", "#54B986","#094334"), breaks = c(0,0.5), labels = c("min", "max"))

fp2 <- FeaturePlot(selected_clusters, pt.size = 0.5, max.cutoff = 1, features = c("CD45RA"), reduction = "umap", coord.fixed = FALSE) +
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),                     
    axis.text = element_text(size = 12), 
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  xlim(-8.5,13) +
  labs(color = "Z-score",x = "UMAP 1", y = "UMAP 2")  + NoLegend() +
  scale_color_gradientn(colors = c("#C7ECDC", "#54B986","#094334"), breaks = c(0,0.5), labels = c("min", "max"))

comb <- (fp1 | fp2) + plot_layout(ncol=2)

ggsave("fig_4c.png", height = 3, width = 6)


# run differential polarity analysis

Idents(selected_clusters) <- "condition"
pol_test <- RunDPA(selected_clusters, 
                   target = "Non-responders",
                   reference = "Responders", 
                   contrast_column = "condition")


# plot differential polarity

PolarizationScores(selected_clusters, meta_data_columns = "condition") %>%
  filter(condition == "Non-responders") %>%
  group_by(marker) %>%
  summarize(pol_mean = mean(morans_z)) %>%
  left_join(pol_test, by = "marker") %>%
  ggplot(aes(estimate, -log10(p_adj), label = marker, color = pol_mean)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = 70, size = 4, color = "black",
                  point.padding = 0.5,
                  box.padding = 0.2) +
  scale_color_gradient2(low = "#C7ECDC", mid = "#54B986", high = "#094334", midpoint = 0.5, breaks = c(0,0.5,1)) +
  xlim(-1,1) +
  labs(x = "Median difference", y = expression(-log[10]~(adjusted~italic(p)~value)), color = "Mean polarity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "Arial", size = 12, color = "black"),         
        axis.title = element_text(size = 12, color = "black"),                     
        axis.text = element_text(size = 12, color = "black"),                      
        legend.title = element_text(size = 12, color = "black"),                  
        legend.text = element_text(size = 12, color = "black"),
        legend.position = c(0.3,0.75), legend.direction = "horizontal") +
  ylim(0,29) +
  xlim(-0.4, 0.6) + NoLegend() +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed')

ggsave("fig_4d.png", height = 2.5, width = 2.6)


# run differential colocalization analysis

progressr::with_progress({
  colocalization_test <- RunDCA(selected_clusters,
                                target = "Non-responders",
                                reference = "Responders",
                                contrast_column = "condition")
})


# plot differential colocalization

ColocalizationScores(selected_clusters, meta_data_columns = "condition") %>%
  filter(condition == "Non-responders") %>%
  group_by(marker_1, marker_2) %>%
  summarize(coloc_mean = mean(pearson_z), .groups = "drop") %>%
  left_join(colocalization_test, by = c("marker_1", "marker_2")) %>%
  mutate(label = ifelse(abs(estimate) > 1, paste(marker_1, marker_2, sep = "/"), "")) %>%
  ggplot(aes(estimate, -log10(p_adj), label = label, color = coloc_mean)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = label), max.overlaps = 20, size = 4, box.padding = 0.6, point.padding = 0.3, segment.size = 0.2, force = 2) +
  scale_color_gradient2(low = "#023858", mid = "#3690C0", high = "#D0D1E6", midpoint = -5) +
  scale_x_continuous(expand = expansion(0.1)) +
  scale_y_continuous(expand = expansion(0.1)) +
  theme_bw() +
  guides(color = guide_colorbar(title.position = "top")) +
  labs(x = "Median difference",
       y = expression(-log[10]~(adj.~p-value)), 
       color = "Mean colocalization score") +
  theme(legend.position = c(0.18,0.8), legend.direction = "horizontal", panel.grid.minor = element_blank())

ggsave("supp_figure_4b.png", width = 6, height = 6)


# calculate colocalization scores based on the condition and remove same marker pairs

DefaultAssay(selected_clusters) <- "mpxCells"
colocalization_scores <- ColocalizationScores(selected_clusters, meta_data_columns = "condition") %>%
  filter(marker_1 != marker_2)


# create new contrast variable and select marker pairs for visualization

plot_components <- colocalization_scores %>%
  unite("contrast", marker_1:marker_2, sep = "/") %>%
  filter(contrast == "CD45/CD82") 


# arrange and filter colocalization scores 

plot_components <- plot_components %>%
  arrange(-pearson_z) %>%
  filter(abs(pearson_z) < 10) %>%
  group_by(condition) %>%
  filter(row_number() %in% 1:10 | rev(row_number()) %in% 1:50) %>%
  ungroup()

# select relevant columns and convert component names

plot_components <- 
  dplyr::select(plot_components, contrast, pearson_z, component, condition) %>%
  mutate(component_number = str_extract(component, "[0-9]+") %>% as.numeric())


# order data by colocalization values

plot_components %>%
  ungroup() %>% 
  arrange(pearson_z) 


# select components for visualization

plot_components <- plot_components %>%
  filter(row_number() %in% c(1,118))


# load graph layout

selected_clusters <- selected_clusters %>%
  LoadCellGraphs(cells = plot_components$component) %>%
  ComputeLayout(layout_method = "fr", seed = 1)


# plot 2D graph based on computed layout

Plot2DGraph(selected_clusters, cells = plot_components$component, layout_method = "fr") +
  plot_layout(ncol = 4) &
  theme(plot.title = element_text(size = 8))


# plot and customize 2D graph

plot_markers <- c("CD45", "CD82")
plot_titles <- paste(plot_components$condition, plot_components$component_number, sep = "\n") %>%
  setNames(nm = plot_components$component)

Plot2DGraphM(selected_clusters,
             layout_method = "fr",
             node_size = 2,
             cells = plot_components$component,
             markers = plot_markers,
             colors = viridis::inferno(n = 11),
             titles = plot_titles) +
  plot_layout(heights = c(1, rep(2, 2)))

ggsave("fig_4e.png", height = 8, width = 6)
