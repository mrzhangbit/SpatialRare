rm(list = ls())  # Clear environment variables

# Load required libraries
library(cluster)
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
library(scran)
library(DT)
library(Matrix)

# Set working directory (relative path)
setwd("../data/")
getwd()

# Load custom functions
source("../R_tools/slice.R")
source("../R_tools/cci.R")

# Load dataset
data <- read.csv("../data/expression_matrix.csv", row.names = 1, header = TRUE)

# Data Preprocessing: Check and clean missing values
if (sum(is.na(data)) > 0) {
  cat("Warning: Missing values detected. Imputing with mean values.\n")
  for (col in colnames(data)) {
    data[is.na(data[, col]), col] <- mean(data[, col], na.rm = TRUE)
  }
}

# Remove zero-variance features
data <- data[, apply(data, 2, var) > 1e-10]

# Transpose the dataset
transposed_data <- t(data)

# Apply log transformation to stabilize variance
transformed_data <- log2(transposed_data + 1)

# K-means clustering
set.seed(123)  # Ensure reproducibility
kmeans_result <- kmeans(data, centers = 13)
my_meta <- data.frame(row.names = rownames(data), Classification = as.factor(kmeans_result$cluster))

# PCA visualization
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
pca_components <- pca_result$x[, 1:2]
pca_df <- data.frame(x = pca_components[, 1], y = pca_components[, 2], meta = factor(my_meta$Classification))

# Plot PCA
p <- ggplot(pca_df, aes(x = x, y = y, color = meta)) +
  geom_point(size = 5) +
  ggtitle("K-means Classification") +
  scale_color_discrete(name = "K-means Classification") +
  theme_minimal()
p

# Save PCA plot as SVG
ggsave("../results/kmeans_Classification.svg", plot = p, width = 8, height = 4)

# Convert data to sparse matrix format
data_matrix <- as.matrix(transformed_data)
data.input <- as(transformed_data, "dgCMatrix")

# Ensure metadata consistency
if (!all(rownames(my_meta) %in% colnames(data_matrix))) {
  stop("Error: Mismatch between metadata and expression matrix.")
}

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = my_meta, group.by = "Classification")
cellchat <- addMeta(cellchat, meta = my_meta)
cellchat <- setIdent(cellchat, ident.use = "Classification")
groupSize <- as.numeric(table(cellchat@idents))  # Count cells per group

# Load CellChat database (for mouse)
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

# Data preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

# Compute cell-cell communication probability
gc()
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualize interaction network
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weights")

# Visualize individual cell interactions
interaction_matrix <- cellchat@net$count
par(mfrow = c(3,3), xpd = TRUE)

for (i in 1:nrow(interaction_matrix)) {
  temp_matrix <- matrix(0, nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix),
                        dimnames = dimnames(interaction_matrix))
  temp_matrix[i,] <- interaction_matrix[i,]
  netVisual_circle(temp_matrix, vertex.weight = groupSize, weight.scale = TRUE, 
                   arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(interaction_matrix), 
                   title.name = rownames(interaction_matrix)[i])
}

# Extract inferred interactions
df.net <- subsetCommunication(cellchat)

# Perform CCI analysis
gc()
my_cci <- get_all_cci(data, my_meta, df.net)

# Save analysis results

write.csv(df.net, file = "../results/df_net.csv")
write.csv(my_cci, file = "../results/cci.csv")

