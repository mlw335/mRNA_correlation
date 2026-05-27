library(ggplot2)
library(uwot)

# read in files
cor_mat <- readRDS("../data/correlation_matrix_all_genes.rds")
uniprot_to_function <- read_tsv("../data/uniprotkb_proteome_UP000000625_2026_03_19.tsv")

#compute distances
dist_mat <- 1 - cor_mat

set.seed(1)

# filter NAs
keep <- complete.cases(dist_mat) & apply(dist_mat, 1, function(x) all(!is.na(x)))
dist_mat_clean <- dist_mat[keep, keep]
cor_mat_clean <- cor_mat[keep, keep]
umap_coords <- uwot::umap(cor_mat_clean)

# generate umap df
umap_df <- data.frame(
  Gene = rownames(cor_mat_clean),
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2]
)

write.csv(umap_df, file = "../data/UMAP_df.csv")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_classic()
