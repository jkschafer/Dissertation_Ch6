library(MCMCglmm)

# Function of getting phylogenetic correlation
phylo_corr <- function(covar_XY, var_X, var_Y) {
  post.mode <- posterior.mode(covar_XY/sqrt(var_X * var_Y))
  return(post.mode)
}