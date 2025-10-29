library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm) 

meta <- read.delim(file = "GSE231489_cellMetadata.tsv")

# Calculate DMSO means
dmso_means <- meta %>%
  filter(Treatment == "DMSO") %>%
  summarise(across(c("MES.like", "AC.like", "NPC.like", "OPC.like"), \(x) mean(x, na.rm = TRUE)))

# Create new columns for normalized scores (the "shift"), keeping original scores
df_with_norm_wide <- meta %>%
  mutate(across(c("MES.like", "AC.like", "NPC.like", "OPC.like"), 
                ~ if_else(Treatment == "Alisertib", .x - as.numeric(dmso_means[cur_column()]), NA_real_),
                .names = "{.col}_norm")) # Creates MES.like_norm, AC.like_norm, etc.

# Pivot the *normalized* data to long format for plotting
df_long <- df_with_norm_wide %>%
  pivot_longer(cols = c("MES.like_norm", "AC.like_norm", "NPC.like_norm", "OPC.like_norm"), 
               names_to = "Subtype", 
               values_to = "Score") %>%
  filter(Treatment == "Alisertib") %>%  # Only plot normalized Alisertib scores
  mutate(Subtype = gsub("_norm$", "", Subtype)) # Remove "_norm" suffix to get clean subtype names

# Calculate p-values manually for each subtype
pvals <- df_long %>%
  group_by(Subtype) %>%
  summarise(n = n(), 
            p.value = wilcox.test(Score, mu = 0)$p.value,
            mean_score = mean(Score, na.rm = TRUE),
            max_score = max(Score, na.rm = TRUE), # Get max score for positioning
            .groups = "drop")

# Format p-values
pvals <- pvals %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         label = paste0("FDR =", signif(p.adj, 2)),
         y.position_bar = mean_score + 0.1 * max(abs(df_long$Score), na.rm = TRUE),
         y.position_violin = max_score + 0.05 * max(df_long$Score, na.rm = TRUE)
  )

# Save Statistics CSV ---
write.csv(pvals, file = "figure_3i_v2_statistics.csv", row.names = FALSE)


# Save Source Data CSV ---
write.csv(df_long, file = "figure_3i_v2_source_data.csv", row.names = FALSE)

pdf(file = "figure_3i_v2.pdf")
ggplot(df_long, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.3) + # Show distribution shape
  # geom_quasirandom(aes(fill = Subtype), color = "black", shape = 21, groupOnX = T, alpha = 0.2, size = 0.3, stroke = 0.2) + # Original quasirandom
  geom_boxplot(aes(fill = Subtype), width = 0.5, 
               # fill = "white", 
               outlier.shape = NA, alpha = 0.7) + # Wider boxplots
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Reference line
  # Add p-value labels
  geom_text(data = pvals, aes(x = Subtype, y = y.position_violin, label = label), 
            inherit.aes = FALSE, vjust = 0, size = 3.5, color = "black") +
  labs(# title = "Violin Plot of Normalized Subtype Scores (Alisertib vs DMSO)",
    # subtitle = "Shows distribution, median/quartiles, and all individual points",
    x = "State",
    y = "Normalized Enrichment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") # Hide legend
dev.off()

