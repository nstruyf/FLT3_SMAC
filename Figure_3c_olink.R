
# Figure 3c: Volcano plot of soluble proteins measured with olink in FLT3 mutant midostaurin responders and non-responders

# load libraries

library(dplyr)
library(ggplot2)
library(ggrepel)


# load sDSS data and annotation files

setwd()
data <- as.matrix(read.csv("olink_NPX.csv", header = TRUE, row.names = 1))
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE) 


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


# plot data

cutoff <- result$p < 0.05

ggplot(result, aes(x = diff, y = log10p)) +
  geom_point(size = ifelse(cutoff, 1.5, 2),
             color = ifelse(cutoff, "grey50", "grey70")) +
  geom_point(data = result[cutoff,], size = 3,
             color = ifelse(result[cutoff,"diff"] > 0, "#870052", "#4DB5BC")) +
  geom_text_repel(aes(label = ifelse(cutoff, as.character(Protein), ""))) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed') +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = 'dashed') +
  ylab(expression(-log[10]~(italic(p)~value))) +
  xlab("Differential soluble \nprotein expression") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans", size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

ggsave("fig_3c.png", width = 3.3, height = 3)


