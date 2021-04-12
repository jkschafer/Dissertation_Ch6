
# Functions for RPANDA analysis from 
# http://www.phytools.org/***SanJuan2016/ex/14/Complex-diversification-models.html

lambda.cst <- function(x,y){y}
lambda.var <- function(x,y){y[1]*exp(y[2]*x)}
mu.cst <- function(x,y){y}
mu.var <- function(x,y){y[1]*exp(y[2]*x)}

bd <- function(x){
  if(class(x)!="birthdeath") stop("x should be an object of class 'birthdeath'")
  b <- x$para[2]/(1-x$para[1])
  d <- b-x$para[2]
  setNames(c(b, d),c("b", "d"))
}

fit.multi.rpanda <- function(tree, par)
{
  bcstdcst <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.cst, 
                     f.mu = mu.cst, 
                     lamb_par = par[[1]][1],
                     mu_par = par[[1]][2],
                     cst.lamb = TRUE,
                     cst.mu = TRUE,
                     cond = "crown",
                     f = 99/301, # n of species in phylogeny/total in 10Ktree
                     dt = 1e-3)
  bvardcst <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.cst, 
                     lamb_par = par[[2]][c(1,2)],
                     mu_par = par[[2]][3],
                     expo.lamb = TRUE,
                     cst.mu = TRUE,
                     cond = "crown",
                     f = 99/301,
                     dt = 1e-3)
  bcstdvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.cst, 
                     f.mu = mu.var, 
                     lamb_par = par[[3]][1],
                     mu_par = par[[3]][c(2,3)],
                     cst.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 99/301,
                     dt = 1e-3)
  bvardvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.var, 
                     lamb_par = par[[4]][c(1, 2)],
                     mu_par = par[[4]][c(3, 4)],
                     expo.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 99/301,
                     dt = 1e-3)
  return(list("bcstdcst" = bcstdcst,
              "bvardcst" = bvardcst,
              "bcstdvar" = bcstdvar,
              "bvardvar" = bvardvar))
}


plotdtt <- function (fit.bd, tot_time, 
                     N0, col = 1, 
                     add = FALSE, 
                     div.time, 
                     xlim, ylim)
{
  if (!inherits(fit.bd, "fit.bd"))
    stop("object \"fit.bd\" is not of class \"fit.bd\"")
  t <- seq(tot_time-div.time, tot_time, 0.01)
  if ("f.mu" %in% attributes(fit.bd)$names) {
    r <- function(t) {
      -fit.bd$f.lamb(t) + fit.bd$f.mu(t)
    }
    R <- function(s) {
      RPANDA:::.Integrate(Vectorize(r), 0, s)
    }
    N <- N0 * exp(Vectorize(R)(t))
    #dev.new()
    if(add==FALSE)
    {
      plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
           main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
    }
    else
    {
      lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
            main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
    }
  }
  else {
    r <- function(t) {
      -fit.bd$f.lamb(t)
    }
    R <- function(s) {
      RPANDA:::.Integrate(Vectorize(r), 0, s)
    }
    N <- N0 * exp(Vectorize(R)(t))
    #dev.new()
    if(add==FALSE)
    {
      plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
           main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
    }
    else
    {
      lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
            main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
    }
  }
}

