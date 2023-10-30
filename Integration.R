#' Construct a KNN (K nearest neighbors) network
#' @param data The gene expression matrix in log scale. Each row corresponds to the gene and column corresponds to the cell
#' @param k_num The number of nearest neighbors to construct KNN network
#'
reticulate::py_install(packages ='umap-learn')
library(umap)
library(parallelDist)
library(irlba)
library(loe)
MakeKNN <- function(data = in_data, k_num = in_k){
  edist <- parDist(data, method = 'euclidean')
  edist <- as.matrix(edist)
  return(edist)
}

set.seed(1)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
a=inData1[sampled_genes,]
umap_res1 = umap(t(inData1[sampled_genes,]), method ='umap-learn')
umap_knn1 <- MakeKNN(data = umap_res1$layout, k_num = 10)
umap_knn1 <- umap_knn1 + t(umap_knn1)

set.seed(2)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
a=inData1[sampled_genes,]
umap_res2 = umap(t(inData1[sampled_genes,]), method ='umap-learn')
umap_knn2 <- MakeKNN(data = umap_res2$layout, k_num = 10)
umap_knn2 <- umap_knn2 + t(umap_knn2)

set.seed(3)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
a=inData1[sampled_genes,]
umap_res3 = umap(t(inData1[sampled_genes,]), method ='umap-learn')
umap_knn3 <- MakeKNN(data = umap_res3$layout, k_num = 10)
umap_knn3 <- umap_knn3 + t(umap_knn3)

set.seed(4)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
a=inData1[sampled_genes,]
umap_res4 = umap(t(inData1[sampled_genes,]), method ='umap-learn')
umap_knn4 <- MakeKNN(data = umap_res4$layout, k_num = 10)
umap_knn4 <- umap_knn4 + t(umap_knn4)

set.seed(5)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
a=inData1[sampled_genes,]
umap_res5 = umap(t(inData1[sampled_genes,]), method ='umap-learn')
umap_knn5 <- MakeKNN(data = umap_res5$layout, k_num = 10)
umap_knn5 <- umap_knn5 + t(umap_knn5)

umap=1/5*(umap_knn1+umap_knn2+umap_knn3+umap_knn4+umap_knn5)




set.seed(1)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
my_pca_log <- prcomp_irlba(t(inData1[sampled_genes,]), n = 10)
PCs <- (my_pca_log$x)
cor_dist <- cor(t(PCs), method ='pearson')
cor_knn1 <- make.kNNG(1-cor_dist, k = 10, symm = F)
cor_knn1 <- cor_knn1 + t(cor_knn1)

set.seed(2)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
my_pca_log <- prcomp_irlba(t(inData1[sampled_genes,]), n = 10)
PCs <- (my_pca_log$x)
cor_dist <- cor(t(PCs), method ='pearson')
cor_knn2 <- make.kNNG(1-cor_dist, k = 10, symm = F)
cor_knn2 <- cor_knn2 + t(cor_knn2)

set.seed(3)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
my_pca_log <- prcomp_irlba(t(inData1[sampled_genes,]), n = 10)
PCs <- (my_pca_log$x)
cor_dist <- cor(t(PCs), method ='pearson')
cor_knn3 <- make.kNNG(1-cor_dist, k = 10, symm = F)
cor_knn3 <- cor_knn3 + t(cor_knn3)

set.seed(4)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
my_pca_log <- prcomp_irlba(t(inData1[sampled_genes,]), n = 10)
PCs <- (my_pca_log$x)
cor_dist <- cor(t(PCs), method ='pearson')
cor_knn4 <- make.kNNG(1-cor_dist, k = 10, symm = F)
cor_knn4 <- cor_knn4 + t(cor_knn4)

set.seed(5)
sample_rng <- round(runif(n=1, min = 0.6, max = 1), digit = 2)
sampled_genes <- sample(features, size = ceiling(sample_rng*length(features)), replace = F)
my_pca_log <- prcomp_irlba(t(inData1[sampled_genes,]), n = 10)
PCs <- (my_pca_log$x)
cor_dist <- cor(t(PCs), method ='pearson')
cor_knn5 <- make.kNNG(1-cor_dist, k = 10, symm = F)
cor_knn5 <- cor_knn5 + t(cor_knn5)

PCA=1/5*(cor_knn1+cor_knn2+cor_knn3+cor_knn4+cor_knn5)


library(ANF)
K = 10

distL<- list(1-PCAKNN, 1-UMAPKNN)
#distL=as.numeric(distL)
affinityL = lapply(distL, function(x) affinity_matrix(x, K)) 
knn_PCAKNN=affinityL[[1]]
knn_UMAPKNN=affinityL[[2]]
Integrated_graph=1/2*(PCAKNN+UMAPKNN)


