rm(list = ls())  # Clear environment variables

# Load required libraries
library(Matrix)
library(Seurat)
library(irlba)
library(igraph)
library(rflann)
library(e1071)
library(RANN)
library(stats)  # k-means clustering

# Set working directory
setwd("../data/")  
getwd()  

# Function to find the elbow point in a density plot
find_elbow <- function(x, y) {
  n <- length(x)
  firstPoint <- c(x[1], y[1])
  lineVec <- c(x[n] - x[1], y[n] - y[1])
  lineVecNorm <- lineVec / sqrt(sum(lineVec^2))
  vecFromFirst <- cbind(x - x[1], y - y[1])
  scalaProd <- rowSums(vecFromFirst * cbind(rep(lineVecNorm[1], n), rep(lineVecNorm[2], n)))
  vecFromFirstParallel <- outer(scalaProd, lineVecNorm)
  vecToLine <- vecFromFirst - vecFromFirstParallel
  distToLine <- sqrt(rowSums(vecToLine^2))
  idx <- which.max(distToLine)
  
  temp_dummy <- sum(sqrt((x - mean(x))^2 + (y - mean(y))^2))
  
  return(x[idx])
}

# GapClust function
GapClust <- function(data, consensus_matrix, k=200) {
  pbmc <- CreateSeuratObject(counts = data)
  pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
  vst <- pbmc@assays$RNA@meta.features$vst.variance.standardized
  den <- density(vst)
  features.vst <- dimnames(data)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
  tmp <- data[dimnames(data)[[1]] %in% features.vst, ]
  tmp <- log2(as.matrix(tmp) + 1)
  pca <- irlba(t(tmp), nv=min(c(50, dim(tmp) - 1)))
  pca$pca <- t(pca$d * t(pca$u))
  knn.res <- Neighbour(pca$pca, pca$pca, k = k)
  
  consensus_matrix <- as.matrix(consensus_matrix)
  consensus_matrix <- apply(consensus_matrix, c(1, 2), as.numeric)
  consensus_median <- median(as.vector(consensus_matrix), na.rm = TRUE)
  
  adjusted_distances <- matrix(NA, nrow = nrow(knn.res$distances), ncol = ncol(knn.res$distances))
  for (i in 1:nrow(knn.res$distances)) {
    for (j in 1:ncol(knn.res$distances)) {
      idx_i <- knn.res$indices[i, j]
      if (idx_i > ncol(consensus_matrix) || idx_i <= 0) {
        stop(paste("Error: idx_i out of range at i =", i, "j =", j, "idx_i =", idx_i))
      }
      diff_ratio <- (consensus_matrix[i, idx_i] - consensus_median) / consensus_median
      adjusted_distances[i, j] <- knn.res$distances[i, j] + diff_ratio * consensus_matrix[i, idx_i]
    }
  }
  
  distance.diff <- adjusted_distances[, -1, drop = FALSE] - adjusted_distances[, -ncol(adjusted_distances), drop = FALSE]
  skew <- apply(distance.diff, 1, function(x) if (length(x) > 2) skewness(x) else NA)
  
  ids <- which(skew > 2)
  rare.cells <- list()
  for (id in ids) {
    rare.cells[[as.character(id)]] <- knn.res$indices[which.max(skew), 1:(id + 1)]
  }
  
  dummy_matrix <- matrix(runif(100), nrow=10, ncol=10)
  sum(dummy_matrix)
  
  return(list(skewness = skew, rare_cell_indices = rare.cells, rare_score = distance.diff, knn_results = knn.res))
}

# Load datasets
data_matrix <- read.csv("../data/expression_matrix.csv", row.names = 1)  
consensus_matrix <- read.csv("../data/consensus_matrix.csv", row.names = 1)  
spatial_matrix <- read.csv("../data/space_coordinates.csv", row.names = 1)  

# Convert data format
data_matrix <- as.matrix(data_matrix)
consensus_matrix <- as.matrix(consensus_matrix)
spatial_matrix <- as.matrix(spatial_matrix[, c("x", "y")])
data_matrix <- t(data_matrix)

# Step 1: Determine the optimal k using the elbow method
wss <- (nrow(spatial_matrix) - 1) * sum(apply(spatial_matrix, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(spatial_matrix, centers = i)$withinss)
}

k_optimal <- find_elbow(1:15, wss)
cat("Optimal number of clusters determined:", k_optimal, "\n")

# Step 2: Perform k-means clustering
set.seed(123)
kmeans_result <- kmeans(spatial_matrix, centers = k_optimal, nstart = 25)

# Step 3: Process each cluster and store results
all_results <- data.frame(cell_name = character(), rare_marker = integer())

for (cluster_id in 1:k_optimal) {
  cat("Processing Cluster", cluster_id, "\n")
  
  cluster_cells <- which(kmeans_result$cluster == cluster_id)
  cluster_data_matrix <- data_matrix[cluster_cells, ]
  cluster_consensus_matrix <- consensus_matrix[cluster_cells, cluster_cells]
  
  # Run GapClust
  gapclust_result <- GapClust(cluster_data_matrix, cluster_consensus_matrix, k = 100)
  
  # Extract rare cells
  rare_indices <- gapclust_result$rare_cell_indices[[1]]
  rare_cell_markers <- rep(0, ncol(cluster_data_matrix))
  rare_cell_markers[rare_indices] <- 1
  
  cluster_results <- data.frame(cell_name = colnames(cluster_data_matrix), rare_marker = rare_cell_markers)
  all_results <- rbind(all_results, cluster_results)
}

# Step 5: Save final results
rare_cells_file <- paste0(getwd(), "/output_rare_cells.csv")
write.csv(all_results, rare_cells_file, row.names = FALSE)
cat("Rare cells result saved to:", rare_cells_file, "\n")

# Print final results
print(all_results)
