library(MCMCglmm)

# Function of getting phylogenetic correlation, HPDinterval, plot
phylo_corr <- function(covar_XY, var_X, var_Y) {
  post.mode <- posterior.mode(as.mcmc(covar_XY/sqrt(var_X * var_Y)))
  post.hpd <- HPDinterval(as.mcmc(covar_XY/sqrt(var_X * var_Y)))
  post.plot <- plot(as.mcmc(covar_XY/sqrt(var_X * var_Y)))
  return(list(post.mode, post.hpd, post.plot))
}