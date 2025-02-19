
# Figure 2g: Volcano plot of all tested drugs in FLT3 mutant midostaurin responders and non-responders

# load libraries

library(dplyr)
library(ggplot2)
library(ggrepel)

# load sDSS data and annotation files

setwd()
sDSS_data <- as.matrix(read.csv("sDSS.csv", header = TRUE, row.names = 1))
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE) 
drug_annotation <- read.csv("drug_class.csv", stringsAsFactors = FALSE) 


# assign patients to respective groups and perform multiple t tests

responder <- annotation$ID[annotation$Label == "responder_mut"]
non_responder <- annotation$ID[annotation$Label == "non_responder_mut"]

t_test <- function(Name) {
  data_responder <- as.numeric(sDSS_data[Name, responder])
  data_non_responder <- as.numeric(sDSS_data[Name, non_responder])
  t_test_result <- t.test(data_responder, data_non_responder, na.action = na.omit)
  
  data.frame(
    Drug = Name,
    p = t_test_result$p.value,
    mean1 = mean(data_responder, na.rm = TRUE),
    mean2 = mean(data_non_responder, na.rm = TRUE)
  )
}


# merge t test results into a single table

result <- do.call(rbind, lapply(rownames(sDSS_data), t_test)) %>%
  mutate(
    q = p.adjust(p, method = "fdr"),
    dsDSS = mean2 - mean1,
    log10p = -log10(p),
    log10q = -log10(q)
  ) %>%
  arrange(p) %>%
  as.data.frame()

result_sorted <- result[order(result$Drug), ]


# add drug annotations to result table

rownames(drug_annotation) <- drug_annotation$drug 
drug_annotation <- drug_annotation[order(drug_annotation$drug), ]
result_sorted$class <- drug_annotation$class
result_sorted$subclass <-drug_annotation$subclass
result_sorted$class <- gsub("_", " ", sub("^[^_]*_", "", result_sorted$class))
result_sorted$subclass <- gsub("_", " ", gsub("inhibitor", "i", result_sorted$subclass))


# plot data

cutoff <- result_sorted$p < 0.05 & abs(result_sorted$dsDSS) > 6
result_sorted$alpha <- ifelse(result_sorted$log10p > 1.3, 1, 0.5)
class_color <- c("Conventional chemotherapy" = "#66C2A5",
                 "Kinase inhibitor" = "#FC8D62",
                 "Differentiating/epigenetic modifier" = "#8DA0CB",
                 "Apoptotic modulator" = "#E78AC3",
                 "Immunomodulatory" = "#A6D854",
                 "Hormone therapy" = "#FFD92F",
                 "Metabolic modifier" = "#E5C494",
                 "Other"= "#B3B3B3")

ggplot(result_sorted, aes(x = dsDSS, y = log10p, color = class)) +
  geom_point(aes(alpha = alpha), size = 2) + 
  geom_vline(xintercept = c(-6, 6), col = "black", alpha = 0.5, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, col = "black", alpha = 0.5, linetype = 'dashed') +
  geom_text_repel(data = result_sorted[cutoff, ], aes(label = as.character(Drug)), 
                  point.padding = 0.15, max.overlaps = 50, color = "black", size = 4,
                  vjust = 0.1,
                  min.segment.length = 0.3,
                  segment.size = 0.5,
                  box.padding = 0.3) + 
  scale_color_manual(values = class_color, name = "Drug class") +
  scale_alpha(range = c(0.5, 1), guide = "none") + 
  labs( x = "Differential sDSS", y = expression(-log[10](italic(p)~value))) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
    axis.title = element_text(family = "Arial", size = 12, color = "black"),
    axis.text = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
    axis.text.x = element_text(family = "Arial", size = 12, color = "black", face = "plain"),
    legend.text = element_text(family = "Arial", size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid = element_blank(),
  ) +
  ylim(0, 5.6) +
  xlim(-10, 10)

ggsave("Figure_2g.png", width = 7, height = 3.5)