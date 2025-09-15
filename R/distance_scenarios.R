library(terra)

jaccard_matrix <- function(raster, threshold = 0, na.rm = TRUE) {
  # Efficient binarization
  bin_raster <- classify(raster, 
                         rcl = matrix(c(-Inf, threshold, 0,
                                        threshold, Inf, 1), 
                                      ncol = 3, byrow = TRUE))
  
  # Conversion to compact binary matrix
  mat <- values(bin_raster, mat = TRUE)
  if (na.rm) mat[is.na(mat)] <- 0
  
  # Intersection calculation (A ∩ B)
  a <- crossprod(mat)  # Much more efficient than t(mat) %*% mat
  
  # Marginal sums calculation (|A|, |B|)
  b <- colSums(mat)
  
  # Union calculation (|A ∪ B| = |A| + |B| - |A ∩ B|)
  union_mat <- outer(b, b, "+") - a
  
  # Jaccard similarity calculation with numerical stabilization
  jaccard_sim <- a / (union_mat + .Machine$double.eps)
  jaccard_sim[union_mat == 0] <- 1  # Special case: both empty sets
  
  # Conversion to distance
  jaccard_dist <- 1 - jaccard_sim
  
  # Adding dimension names
  dimnames(jaccard_dist) <- list(names(raster), names(raster))
  
  # Symmetry enforcement (for floating-point precision)
  return(as.dist((jaccard_dist + t(jaccard_dist)) / 2))
}

## Usage example:
# r <- rast(ncols=100, nrows=100, nl=5)
# values(r) <- runif(ncell(r)*5)
# j_mat <- jaccard_matrix(r, threshold=0.5)
# print(j_mat)