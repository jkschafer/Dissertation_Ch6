

other_prim_bd <- fit_bd(other.primates.tree, 
                        max(branching.times(other.primates.tree)), 
                        f.lamb=lambda.cst, 
                        f.mu=mu.cst, 
                        lamb_par=0.4, 
                        mu_par=0,
                        cst.lamb=TRUE,
                        cst.mu=TRUE,
                        cond="crown",
                        f=97/99,
                        dt=1e-3)


pap_bd <- fit_bd(papionins.tree, 
                        max(branching.times(papionins.tree)), 
                        f.lamb=lambda.cst, 
                        f.mu=mu.cst, 
                        lamb_par=0.4, 
                        mu_par=0,
                        cst.lamb=TRUE,
                        cst.mu=TRUE,
                        cond="crown",
                        f=97/99,
                        dt=1e-3)

plot_dtt(pap_bd, N0 = Ntip(papionins.tree))
plot_dtt(other_prim_bd, N0 = Ntip(other.primates.tree))

papionins.tree <- extract.clade(tree2, 115)
other.primates.tree <- drop.tip(tree2, papionins.tree$tip.label)

#we need to supply starting parameter values for optimization to RPANDA
primates.par <- list(c(0.4, 0),
                     c(0.4, -0.05, 0),
                     c(0.4, 0.1, 0.05),
                     c(0.4, -0.05, 0.1, 0.05)) 

results2 <- list("primates.res" = fit.multi.rpanda(other.primates.tree,
                                                   primates.par),
                 "papionins.res" = fit.multi.rpanda(papionins.tree,
                                                    primates.par))


aic.table <- matrix(nrow = 4, ncol = 2, NA)
for(i in 1:2)
{
  for(j in 1:4)
  {
    aic.table[j,i] <- results2[[i]][[j]]$aicc
  }
}
colnames(aic.table) <- c("Other Primates", "Papionins")
rownames(aic.table) <- c("bcstdcst", "bvardcst",
                         "bcstdvar", "bvardvar")
aic.table

par.table <- data.frame("Other Primates"=c(results2[[1]]$bvardvar$lamb_par[1:2],
                                           results2[[1]]$bvardvar$mu_par[1:2]),
                        "Papionins"=c(results2[[2]]$bvardvar$lamb_par[1:2],
                                      results2[[2]]$bvardvar$mu_par[1:2]))
par.table


div.times <- c(max(branching.times(other.primates.tree)), 
               max(branching.times(papionins.tree)))

plotdtt(results2$primates.res$bvardvar, 
        div.times[1],
        N0 = Ntip(other.primates.tree),
        xlim = c(-max(div.times), 0),
        ylim = c(0, 230000),
        div.time = div.times[1])
plotdtt(results$papionins.res$bvardvar, 
        div.times[2],
        N0 = Ntip(papionins.tree),
        col = 6,
        add = TRUE,
        xlim = c(-max(div.times), 0),
        ylim=c(0, 230000),
        div.time = div.times[2])
legend("topleft",legend=c("Other Primates",
                          "Papionins"),
       text.col = c(1, 6, 
                    "goldenrod", 4,
                    "darkred"))
