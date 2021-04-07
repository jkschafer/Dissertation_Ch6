#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Model (PGLMM)   #
#   Results, Diagnostics, Analysis, and Plots             #
#---------------------------------------------------------#

# Loading packages
list_of_packages <- c("tidyverse", "ape", 
                      "MCMCglmm", "phytools",
                      "ggpmisc", "ggpubr",
                      "broom.mixed", "reshape2",
                      "diversitree", "RPANDA",
                      "plyr")
lapply(list_of_packages, library, character.only = TRUE)

# Load data
load("./Data/reduced_data.Rdata")
load("./Data/10ktree.Rdata")

# Make tree bifurcating
tree2 <- multi2di(tree)

ovsign <- reduced_data$Ovulation_Signs

ovsign <- mapvalues(ovsign,
                    from = c("0", "1", "2", "3"), 
                    to = c("1", "2", "3", "4")) 
ovsign <- as.numeric(as.character(ovsign))
names(ovsign) <- reduced_data$Species

  
ov_musse <- make.musse(tree2, ovsign, 4, sampling.f=NULL, strict=TRUE,
                       control=list())
starting.point.musse(tree2, 4, q.div = 5, yule = FALSE)

argnames(ov_musse)
diversitree:::default.argnames.musse(4)

col <- c("blue", "orange", "red", "green")

# RPANDA method
plot(tree2, cex = 0.35)
nodelabels(cex = 0.5)

lambda.cst <- function(x,y){y}
lambda.var <- function(x,y){y[1]*exp(y[2]*x)}
mu.cst <- function(x,y){y}
mu.var <- function(x,y){y[1]*exp(y[2]*x)}


#Use these node numbers to extract clades.
apes.tree  <- extract.clade(tree2, 153)
papionins.tree <- extract.clade(tree2, 115)
platyrrhines.tree <- extract.clade(tree2, 170)
colobines.tree <- extract.clade(tree2, 139)
cercopiths.tree <- extract.clade(tree2, 104)

par(mfcol = c(3, 2))
plot(apes.tree, main = "Apes")
plot(papionins.tree, main = "Papionins") 
plot(platyrrhines.tree, main = "Platyrrhines")
plot(colobines.tree, main = "Colobines")
plot(cercopiths.tree, main = "Other Cercopiths")

# We also need to capture the rest of the cetaceans 
# that do not fall into these clades

fit_bd(apes.tree, 
       max(branching.times(apes.tree)), 
       f.lamb=lambda.cst, 
       f.mu=mu.cst, 
       lamb_par=0.4, 
       mu_par=0,
       cst.lamb=TRUE,
       cst.mu=TRUE,
       cond="crown",
       f=87/89,
       dt=1e-3)

bd <- function(x){
  if(class(x)!="birthdeath") stop("x should be an object of class 'birthdeath'")
  b <- x$para[2]/(1-x$para[1])
  d <- b-x$para[2]
  setNames(c(b, d),c("b", "d"))
}

fit.bd <- birthdeath(apes.tree)
bd(fit.bd)

fit_bd(apes.tree, 
       max(branching.times(apes.tree)), 
       f.lamb=lambda.var, 
       f.mu=mu.cst, 
       lamb_par= c(0.4,-0.05),
       mu_par=0,
       expo.lamb=TRUE,
       cst.mu=TRUE,
       cond="crown",
       f=87/89,dt=1e-3)


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
                     f = 87/89,
                     dt = 1e-3)
  bvardcst <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.cst, 
                     lamb_par = par[[2]][c(1,2)],
                     mu_par = par[[2]][3],
                     expo.lamb = TRUE,
                     cst.mu = TRUE,
                     cond = "crown",
                     f = 87/89,
                     dt = 1e-3)
  bcstdvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.cst, 
                     f.mu = mu.var, 
                     lamb_par = par[[3]][1],
                     mu_par = par[[3]][c(2,3)],
                     cst.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 87/89,
                     dt = 1e-3)
  bvardvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.var, 
                     lamb_par = par[[4]][c(1, 2)],
                     mu_par = par[[4]][c(3, 4)],
                     expo.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 87/89,
                     dt = 1e-3)
  return(list("bcstdcst" = bcstdcst,
              "bvardcst" = bvardcst,
              "bcstdvar" = bcstdvar,
              "bvardvar" = bvardvar))
}

#we need to supply starting parameter values for optimization to RPANDA
primates.par <- list(c(0.4, 0),
                     c(0.4, -0.05, 0),
                     c(0.4, 0.1, 0.05),
                     c(0.4, -0.05, 0.1, 0.05)) 

results <- list("apes.res" = fit.multi.rpanda(apes.tree,
                                              primates.par),
                "papionins.res" = fit.multi.rpanda(papionins.tree,
                                                   primates.par),
                "platyrrhines.res" = fit.multi.rpanda(platyrrhines.tree,
                                                      primates.par),
                "colobines.res" = fit.multi.rpanda(colobines.tree,
                                                   primates.par),
                "cercopiths.res"= fit.multi.rpanda(cercopiths.tree,
                                                   primates.par))


aic.table <- matrix(nrow=4,ncol=5,NA)
for(i in 1:5)
{
  for(j in 1:4)
  {
    aic.table[j,i] <- results[[i]][[j]]$aicc
  }
}
colnames(aic.table) <- c("Apes", "Papionins",
                         "Platyrrhines", "Colobines",
                         "Other Cercopithecines")
rownames(aic.table) <- c("bcstdcst", "bvardcst",
                         "bcstdvar", "bvardvar")
aic.table

par.table <- data.frame("Apes"=c(results[[1]]$bcstdcst$lamb_par[1:2],
                                            results[[1]]$bcstdcst$mu_par[1:2]),
                        "Papionins"=c(results[[2]]$bvardcst$lamb_par[1:2],
                                        results[[2]]$bvardcst$mu_par[1:2]),
                        "Platyrrhines"=c(results[[3]]$bcstdcst$lamb_par[1:2],
                                        results[[3]]$bcstdcst$mu_par[1:2]),
                        "Colobines"=c(results[[4]]$bcstdcst$lamb_par[1:2],
                                     results[[4]]$bcstdcst$mu_par[1:2]),
                        "Other Cercopithecines"=c(results[[5]]$bcstdvar$lamb_par[1:2],
                                            results[[5]]$bcstdvar$mu_par[1:2]))
par.table


# Function to calculate species richness in a given point in time
div.times <- c(max(branching.times(apes.tree)), 
               max(branching.times(papionins.tree)),
               max(branching.times(platyrrhines.tree)),
               max(branching.times(colobines.tree)),
               max(branching.times(cercopiths.tree)))

# Function modified from plot_dtt from RPANDA package
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

plotdtt(results$apes.res$bcstdcst, 
        div.times[1],
        N0 = Ntip(apes.tree),
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time=div.times[1])
plotdtt(results$papionins.res$bvardvar, 
        div.times[2],
        N0 = Ntip(papionins.tree),
        col = 6,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim=c(0, 2300),
        div.time = div.times[2])
plotdtt(results$platyrrhines.res$bcstdcst,
        div.times[3],
        N0 = Ntip(platyrrhines.tree),
        col = "goldenrod", 
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[3])
plotdtt(results$colobines.res$bcstdcst, 
        div.times[4],
        N0 = Ntip(colobines.tree),
        col = 4,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[4])
plotdtt(results$cercopiths.res$bvardcst,
        div.times[5],
        N0 = Ntip(cercopiths.tree),
        col = "darkred",
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[5])
legend("topleft",legend=c("Apes",
                          "Papionins",
                          "Platyrrhines",
                          "Colobines",
                          "Other Cercopithecines"),
       text.col = c(1, 6, 
                    "goldenrod", 4,
                    "darkred"))



plotdtt(results$apes.res$bvardvar, 
        div.times[1],
        N0 = Ntip(apes.tree),
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time=div.times[1])
plotdtt(results$papionins.res$bvardvar, 
        div.times[2],
        N0 = Ntip(papionins.tree),
        col = 6,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim=c(0, 2300),
        div.time = div.times[2])
plotdtt(results$platyrrhines.res$bvardvar,
        div.times[3],
        N0 = Ntip(platyrrhines.tree),
        col = "goldenrod", 
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[3])
plotdtt(results$colobines.res$bvardvar, 
        div.times[4],
        N0 = Ntip(colobines.tree),
        col = 4,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[4])
plotdtt(results$cercopiths.res$bvardvar,
        div.times[5],
        N0 = Ntip(cercopiths.tree),
        col = "darkred",
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 2300),
        div.time = div.times[5])
legend("topleft",legend=c("Apes",
                          "Papionins",
                          "Platyrrhines",
                          "Colobines",
                          "Other Cercopithecines"),
       text.col = c(1, 6, 
                    "goldenrod", 4,
                    "darkred"))



