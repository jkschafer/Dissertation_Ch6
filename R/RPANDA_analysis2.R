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

source("./R/Functions/RPANDAS_Revell.R")

# Load data
load("./Data/reduced_data.Rdata")
load("./Data/10ktree.Rdata")

# Make tree bifurcating
tree2 <- multi2di(tree)
plot(tree2, cex = 0.35)
nodelabels(cex = 0.5)

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
                        "Papionins"=c(results[[2]]$bcstdcst$lamb_par[1:2],
                                        results[[2]]$bcstdcst$mu_par[1:2]),
                        "Platyrrhines"=c(results[[3]]$bcstdcst$lamb_par[1:2],
                                        results[[3]]$bcstdcst$mu_par[1:2]),
                        "Colobines"=c(results[[4]]$bcstdcst$lamb_par[1:2],
                                     results[[4]]$bcstdcst$mu_par[1:2]),
                        "Other Cercopithecines"=c(results[[5]]$bcstdcst$lamb_par[1:2],
                                            results[[5]]$bcstdcst$mu_par[1:2]))
par.table

# Function to calculate species richness in a given point in time
div.times <- c(max(branching.times(apes.tree)), 
               max(branching.times(papionins.tree)),
               max(branching.times(platyrrhines.tree)),
               max(branching.times(colobines.tree)),
               max(branching.times(cercopiths.tree)))

# Plotting diversity through time for different clades
plotdtt(results$apes.res$bcstdcst, 
        div.times[1],
        N0 = Ntip(apes.tree),
        xlim = c(-max(div.times), 0),
        ylim = c(0, 30),
        div.time=div.times[1])
plotdtt(results$papionins.res$bcstdcst, 
        div.times[2],
        N0 = Ntip(papionins.tree),
        col = 6,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim=c(0, 30),
        div.time = div.times[2])
plotdtt(results$platyrrhines.res$bcstdcst,
        div.times[3],
        N0 = Ntip(platyrrhines.tree),
        col = "goldenrod", 
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 30),
        div.time = div.times[3])
plotdtt(results$colobines.res$bcstdcst, 
        div.times[4],
        N0 = Ntip(colobines.tree),
        col = 4,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 30),
        div.time = div.times[4])
plotdtt(results$cercopiths.res$bcstdcst,
        div.times[5],
        N0 = Ntip(cercopiths.tree),
        col = "darkred",
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 30),
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



