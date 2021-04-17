library(phytools)

load("./Data/reduced_data.Rdata")
load("./Data/10Ktree.Rdata")

# Plotting tree and characters
row.names(reduced_data) <- reduced_data$Species
plotTree(tree, ftype = "i", fsize = 0.5, lwd = 1) # Plot of tree
ov_signs <- setNames(reduced_data[, 13], rownames(reduced_data))
cols <- setNames(c("red", "blue", "green", "yellow"), 
                 levels(ov_signs))
tiplabels(pie = to.matrix(ov_signs[tree$tip.label], levels(ov_signs)),
          piecol = cols, cex = 0.3)
add.simmap.legend(colors = cols, 
                  prompt = FALSE, 
                  x = 0.9 * par()$usr[1],
                  y = 0.8 * par()$usr[3], fsize = 0.8)

mtrees <- make.simmap(tree, 
                      ov_signs, 
                      model = "ER", 
                      nsim = 100)
pd <- summary(mtrees)
plot(pd, fsize = 0.6, 
     ftype = "i",
     colors = cols, 
     ylim = c(-2, Ntip(tree)))
add.simmap.legend(colors = cols[4:1], 
                  prompt = FALSE, 
                  x = 0, y = -4, 
                  vertical = FALSE)

#---------- Plot of ovulation signals on tree ---------------#
rownames(reduced_data) <- reduced_data$Species

ovsig <- setNames(reduced_data$Ovulation_Signs, 
                  rownames(reduced_data))

cols <- setNames(c("#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7"),
                 unique(as.vector(sapply(ovsig, levels))))
cols

dotTree(tree, 
        ovsig, 
        length = 10, 
        fsize = 0.5, 
        lwd = 2,
        ftype = "i",
        legend = F,
        colors = cols)
legend("bottomleft",
       legend = c("absent",
                  "absent/slight",
                  "slight",
                  "present"),
       col = cols,
       cex = 1,
       pch = 19)

# ---------------------- Contour map of SSD ----------------#
ssd <- log(reduced_data$SSD)
names(ssd) <- reduced_data$Species
ssd_cont <- contMap(tree, 
                    ssd,
                    plot = FALSE)
ssd_cont_map <- setMap(ssd_cont,
                       colors = cols)
plot(ssd_cont_map,
     lwd = 5,
     fsize = 0.5,
     xlim = c(0, 53))

# ---------------------- Contour map of VTD ----------------#
vtd <- log(reduced_data$VTDwoSD + 1)
names(vtd) <- reduced_data$Species
vtd_cont <- contMap(tree, 
                    vtd,
                    plot = FALSE)
vtd_cont_map <- setMap(vtd_cont,
                       colors = cols)
plot(vtd_cont_map,
     lwd = 5,
     fsize = 0.5,
     xlim = c(0, 53))

#------------- SSD and VTD ---------------#
is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]

layout(matrix(1:3, 1, 3),
       widths=c(0.39, 0.2, 0.39))
plot(vtd_cont_map,
     lwd = 6,
     outline = TRUE,
     fsize=c(0, 1.2),
     ftype = "off",
     xlim = c(0, 53),
     leg.txt = "log(VTD + 1)")
plot.new()
plot.window(xlim = c(-0.1, 0.1),
            ylim = get("last_plot.phylo",
                       envir = .PlotPhyloEnv)$y.lim)
par(cex = 0.6)
text(rep(0, length(tree$tip.label[ordered_tips])),
     1:Ntip(tree),
     gsub("_", " ",
          rev(tree$tip.label[ordered_tips])), font = 3)
plot(ssd_cont_map,
     lwd = 6,
     outline = TRUE,
     fsize = c(0, 1.2),
     xlim = c(0, 53),
     direction = "leftwards",
     ftype = "off",
     leg.txt = "log(SSD)")


#---------- Plot of ovulation signals on tree ---------------#
rownames(reduced_data) <- reduced_data$Species

ovsig <- setNames(reduced_data$Ovulation_Signs, 
                  rownames(reduced_data))

cols <- setNames(c("#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7"),
                 unique(as.vector(sapply(ovsig, levels))))
cols

mtrees <- make.simmap(tree, 
                      ovsig, 
                      model = "ER", 
                      nsim = 100)


densityTree(mtrees, method = "plotSimmap",
            lwd = c(8), nodes = "intermediate",
            colors = cols, ylim = c(0, 99),
            compute.consensus = FALSE,
            use.gradient = FALSE,
            use.edge.length = T,
            fsize = 0.5,
            ftype = "i",
            show.axis = F)
legend("bottomleft",
       legend = c("absent",
                  "absent/slight",
                  "slight",
                  "present"),
       col = cols,
       cex = 1,
       pch = 19)


#-------------------- Starting Over --------------------------#
library(phytools)

load("./Data/reduced_data.Rdata")
load("./Data/10Ktree.Rdata")

# Plotting tree and characters
tree <- multi2di(tree)
row.names(reduced_data) <- reduced_data$Species
plotTree(tree, ftype = "i", 
         fsize = 0.5, lwd = 1) # Plot of tree
rtree <- rotateNodes(tree, "all")
plotTree(rtree, ftype = "i", 
         fsize = 0.5, lwd = 1)
rtree <- force.ultrametric(rtree)

ssd <- log(reduced_data$SSD)
names(ssd) <- reduced_data$Species

vtd <- log(reduced_data$VTDwoSD + 1)
names(vtd) <- reduced_data$Species

ovsig <- setNames(reduced_data$Ovulation_Signs, 
                  rownames(reduced_data))
cols <- setNames(c("#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7"),
                 unique(as.vector(sapply(ovsig, levels))))
cols


mtrees <- make.simmap(rtree, 
                      ovsig, 
                      model = "ER", 
                      nsim = 100)

tips <- getStates(mtrees, "tips")
tip.cols <- cols[tips]


barcols <- setNames(sapply(tips,
                           function(x, y) y[which(names(y) == x)],
                           y = cols),
                    names(tips))

par(mar=c(5.1, 4.1, 4.1, 20), xpd=TRUE)
plotTree.barplot(rtree,
                 ssd,
                 args.plotTree = list(ftype = "off"),
                 args.barplot = list(col = barcols,
                                     border = barcols, xlab = "",
                                   xlim = c(-0.12, 1)),
                 args.axis = list(at = seq(-0.12, 1.35,
                                           by = 0.2)))
legend(x = "topright",
       legend = c("absent",
                  "absent/slight",
                  "slight",
                  "present"),
       col = cols,
       pch = 22,
       pt.cex=2,
       pt.bg = cols,
       box.col="transparent",
       horiz = F,
       inset=c(-0.2,0))
mtext("log(SSD)", 1, at = 0.5, line = 2.5)
## this is how I go back to my tree frame
par(mfg=c(1,1))
plot(tree, cols, ftype = 'off', mar=c(5.1, 1.1, 2.1, 0))
axis(1, at = seq(0, 40, by = 40))
mtext("time", 1 , at = 25, line = 2)


