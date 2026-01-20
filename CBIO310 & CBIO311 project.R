# Installing the project's core libraries
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("VennDiagram")


# Recalling core libraries:
library(tidyverse)
library(ggrepel)
library(VennDiagram)
library(tidyverse)
library(tidyverse)
library(pheatmap)


# Download data:
covid_df <- read.csv(file.choose())
alz_df <- read.csv(file.choose())


# Cleaning up Alzheimer's data:
clean_alz <- alz_df %>%
  mutate(Gene = str_extract(SPOT_ID, "(?<=// )[^ ]+")) %>%
  filter(!is.na(Gene)) %>%
  select(Gene, logFC, adj.P.Val)

# COVID file cleanup: Change the ORF column name to Gene for standardization
clean_covid <- covid_df %>%
  rename(Gene = ORF) %>%
  select(Gene, logFC, adj.P.Val)

# Filtering application: P-value < 0.05 and LogFC > 0.5
degs_covid <- clean_covid %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)
degs_alz <- clean_alz %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)

# Finding the genes common to both diseases
common_genes <- intersect(degs_covid$Gene, degs_alz$Gene)

# Print results:
cat("Number of influential COVID genes: ", nrow(degs_covid), "\n")
cat("Number of genes involved in Alzheimer's: ", nrow(degs_alz), "\n")
cat("Number of genes they share: ", length(common_genes), "\n")


# Drawing a Venn Diagram (to calculate the intersection)
grid.newpage()
venn_plot <- draw.pairwise.venn(
  area1 = nrow(degs_covid),
  area2 = nrow(degs_alz),
  cross.area = length(common_genes),
  category = c("COVID-19", "Alzheimer's"),
  fill = c("#00AFBB", "#E7B800"),
  lty = "blank"
)


# Volcano Plot drawing of COVID data (as an example of analysis)
ggplot(clean_covid, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 0.5), alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: SARS-CoV-2",
       subtitle = "Red points represent significant DEGs",
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value")

# Linkage Matrix for Shared Genes:
comparison_df <- inner_join(
  degs_covid %>% select(Gene, logFC_covid = logFC),
  degs_alz %>% select(Gene, logFC_alz = logFC),
  by = "Gene"
)

ggplot(comparison_df, aes(x = logFC_covid, y = logFC_alz)) +
  geom_point(color = "purple") +
  geom_smooth(method = "lm", color = "darkblue", fill = "lightblue") +
  geom_text_repel(aes(label = Gene)) +
  theme_bw() +
  labs(title = "Common Genes Expression Correlation",
       x = "LogFC in COVID-19", y = "LogFC in Alzheimer's")


# Distribution of Gene Expression (Alzheimer’s):
ggplot(clean_alz, aes(x = logFC)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Gene Expression in Alzheimer’s",
       x = "Log Fold Change",
       y = "Frequency")

# Distribution of Gene Expression (SARS‑CoV‑2):
ggplot(clean_covid, aes(x = logFC)) +
  geom_histogram(bins = 50, fill = "darkred") +
  theme_minimal() +
  labs(title = "Distribution of Gene Expression in SARS‑CoV‑2",
       x = "Log Fold Change",
       y = "Frequency")


head(covid_df)
head(alz_df)

comparison_df <- inner_join(
  degs_covid %>% select(Gene, logFC_covid = logFC),
  degs_alz %>% select(Gene, logFC_alz = logFC),
  by = "Gene"
)

top_shared <- comparison_df %>%
  mutate(mean_logFC = (abs(logFC_covid) + abs(logFC_alz)) / 2) %>%
  arrange(desc(mean_logFC)) %>%
  slice(1:10)

# Bar Plot for the strongest shared genes (Top Shared DEGs)
# We select the top 10 most shared genes based on the average |logFC|:
ggplot(top_shared, aes(x = reorder(Gene, mean_logFC), y = mean_logFC)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Shared Differentially Expressed Genes",
       x = "Gene",
       y = "Mean Absolute LogFC")


# Density Plot for comparing distributions:
combined_df <- bind_rows(
  clean_covid %>% mutate(Disease = "COVID-19"),
  clean_alz %>% mutate(Disease = "Alzheimer")
)

ggplot(combined_df, aes(x = logFC, fill = Disease)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Density Plot of LogFC Distributions",
       x = "Log Fold Change",
       y = "Density")


# Converting data to a suitable heatmap format
heatmap_df <- comparison_df %>%
  column_to_rownames("Gene")

# Heatmap:
pheatmap(
  heatmap_df,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of Shared Differentially Expressed Genes"
)