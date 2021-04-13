


full_tree <- read.nexus("./Data/consensusTree_10kTrees_Primates_Version3.nex")
plot(full_tree, cex = 0.35)
nodelabels(cex = 0.5)


anth_tree  <- extract.clade(full_tree, 303)
anth_tree <- multi2di(anth_tree)
plot(anth_tree, cex = 0.5)
nodelabels(cex = 0.5)

# Extract clades
platyrrhines_tree <- extract.clade(anth_tree, 339)
apes_tree <- extract.clade(anth_tree, 315)
colobus_tree <- extract.clade(anth_tree, 277)
papio_tree <- extract.clade(anth_tree, 239)
cerc_tree <- extract.clade(anth_tree, 209)

par(mfcol = c(3, 2))
plot(apes_tree, main = "Apes")
plot(papio_tree, main = "Papionins") 
plot(platyrrhines_tree, main = "Platyrrhines")
plot(colobus_tree, main = "Colobines")
plot(cerc_tree, main = "Other Cercopithocines")
plot(anth_tree, main = "All Anthropoids")


primates.par <- list(c(0.4, 0),
                     c(0.4, -0.05, 0),
                     c(0.4, 0.1, 0.05),
                     c(0.4, -0.05, 0.1, 0.05)) 

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
                     f = 204/301, # n of species in phylogeny/total in 10Ktree
                     dt = 1e-3)
  bvardcst <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.cst, 
                     lamb_par = par[[2]][c(1,2)],
                     mu_par = par[[2]][3],
                     expo.lamb = TRUE,
                     cst.mu = TRUE,
                     cond = "crown",
                     f = 204/301,
                     dt = 1e-3)
  bcstdvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.cst, 
                     f.mu = mu.var, 
                     lamb_par = par[[3]][1],
                     mu_par = par[[3]][c(2,3)],
                     cst.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 204/301,
                     dt = 1e-3)
  bvardvar <- fit_bd(tree, max(branching.times(tree)), 
                     f.lamb = lambda.var, 
                     f.mu = mu.var, 
                     lamb_par = par[[4]][c(1, 2)],
                     mu_par = par[[4]][c(3, 4)],
                     expo.lamb = TRUE,
                     expo.mu = TRUE,
                     cond = "crown",
                     f = 204/301,
                     dt = 1e-3)
  return(list("bcstdcst" = bcstdcst,
              "bvardcst" = bvardcst,
              "bcstdvar" = bcstdvar,
              "bvardvar" = bvardvar))
}

results <- list("apes.res" = fit.multi.rpanda(apes_tree,
                                              primates.par),
                "papionins.res" = fit.multi.rpanda(papio_tree,
                                                   primates.par),
                "platyrrhines.res" = fit.multi.rpanda(platyrrhines_tree,
                                                      primates.par),
                "colobines.res" = fit.multi.rpanda(colobus_tree,
                                                   primates.par),
                "cercopiths.res"= fit.multi.rpanda(cerc_tree,
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
div.times <- c(max(branching.times(apes_tree)), 
               max(branching.times(papio_tree)),
               max(branching.times(platyrrhines_tree)),
               max(branching.times(colobus_tree)),
               max(branching.times(cerc_tree)))

# Plotting diversity through time for different clades
plotdtt(results$apes.res$bcstdcst, 
        div.times[1],
        N0 = Ntip(apes_tree),
        xlim = c(-max(div.times), 0),
        ylim = c(0, 80),
        div.time=div.times[1])
plotdtt(results$papionins.res$bcstdcst, 
        div.times[2],
        N0 = Ntip(papio_tree),
        col = 6,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim=c(0, 80),
        div.time = div.times[2])
plotdtt(results$platyrrhines.res$bcstdcst,
        div.times[3],
        N0 = Ntip(platyrrhines_tree),
        col = "goldenrod", 
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 80),
        div.time = div.times[3])
plotdtt(results$colobines.res$bcstdcst, 
        div.times[4],
        N0 = Ntip(colobus_tree),
        col = 4,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 80),
        div.time = div.times[4])
plotdtt(results$cercopiths.res$bcstdcst,
        div.times[5],
        N0 = Ntip(cerc_tree),
        col = "darkred",
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim = c(0, 80),
        div.time = div.times[5])
legend("topleft",legend=c("Apes",
                          "Papionins",
                          "Platyrrhines",
                          "Colobines",
                          "Other Cercopithecines"),
       text.col = c(1, 6, 
                    "goldenrod", 4,
                    "darkred"))

