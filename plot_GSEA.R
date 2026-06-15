#### GSEA Plot


library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

# read file
gsea <- read_csv("/path/to/PDX_analysis_R/GSEA_lowESR1_nonprolif_vs_highESR1.csv")

# clean and prepare
gsea2 <- gsea %>%
  mutate(
    pathway_short = str_wrap(pathway, width = 45),
    direction = ifelse(NES >= 0, "Enriched in low_ESR1_nonprolif", "Enriched in high_ESR1")
  ) %>%
  arrange(desc(abs(NES)))

# choose top pathways to avoid overcrowding
top_n_plot <- 55

gsea_top <- gsea2 %>%
  slice_max(order_by = abs(NES), n = top_n_plot) %>%
  mutate(pathway_short = fct_reorder(pathway_short, NES))

p1 <- ggplot(gsea_top, aes(x = NES, y = pathway_short)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = `-log10(padj)`, color = direction), alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adj.~p))) +
  scale_color_manual(values = c(
    "Enriched in low_ESR1_nonprolif" = "#D95F02",
    "Enriched in high_ESR1" = "#1B9E77"
  )) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Top enriched GSEA signatures",
    subtitle = paste("Top", top_n_plot, "signatures ranked by |NES|")
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

p1

