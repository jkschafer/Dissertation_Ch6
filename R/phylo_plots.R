list_of_packages <- c("tidyverse", "ggtree", 
                      "ape", "phytools",
                      "ggpmisc", "ggpubr",
                      "ggstance", "rcartocolor")
lapply(list_of_packages, library, character.only = TRUE)


library(ggtree)
library(ggstance)
library(ape)
library(phytools)

# load data
load("./Data/reduced_data.Rdata")
load("./Data/10Ktree.Rdata")

# Plot tree, make bifurcating and find node labels
plotTree(tree, ftype = "i", 
         fsize = 0.5, lwd = 1) # Plot of tree
tree <- multi2di(tree)
nodelabels(cex = 0.5)

# Plotting tree with clade labels
p <- ggtree(tree) + 
        theme_tree2() +
        geom_cladelabel(node = 104, 
                        label = "Cercopithecini",
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 115, 
                        label = "Papioini", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 139, 
                        label = "Colobinae", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 154, 
                        label = "Hylobatidae", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 164, 
                        label = "Hominidae", 
                        #align = TRUE, 
                        angle = 270,
                        offset = 2,
                        hjust = 0.5,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 172, 
                        label = "Atelidae", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 180, 
                        label = "Callitrichidae", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") +
        geom_cladelabel(node = 194, 
                        label = "Cebidae", 
                        #align = TRUE, 
                        angle = 270,
                        hjust = 0.5,
                        offset = 2,
                        offset.text = 1,
                        fontsize = 3.5,
                        color = "#000000") 
p1 <- revts(p)

df <- data.frame(id = reduced_data$Species,
                 SSD = log(reduced_data$SSD),
                 VTD = log(reduced_data$VTDwSD + 1),
                 category = reduced_data$Ovulation_Signs)

df$category2 <- ifelse(df$category == "3", "present",
                       ifelse(df$category == "2", "slight",
                              ifelse(df$category == "1", "absent/slight",
                                     "absent")))

p2 <- facet_plot(p1, panel = 'log(SSD)', 
                 data = df, 
                 geom = geom_barh,
                 mapping = aes(x = SSD, 
                               fill = as.factor(category2)),
                 stat = 'identity') +
        theme_tree2()
p2 + theme(strip.text.x = element_text(colour = "black", 
                                       face = "bold", 
                                       size = 15))

df$category2 <- factor(df$category2, 
                       levels = c("absent", 
                                  "absent/slight", 
                                  "slight",
                                  "present"))

# Color palette for plot
mycolors <- rcartocolor::carto_pal(4, "Earth")
# Reverse first two colors of palette manually
mycolors <- c("#d6bd8d", "#A16928",
              "#b5c8b8", "#2887a1")
# Final plot
p3 <- facet_plot(p2, panel = 'log(VTD + 1)', 
                 data = df, 
                 geom = geom_barh,
                 mapping = aes(x = VTD, 
                               fill = as.factor(category2)),
                 stat = 'identity') +
        scale_fill_manual(values = mycolors,
                          limits = c("absent", 
                                     "absent/slight", 
                                     "slight",
                                     "present"),
                          name = "Ovulation Signal") +
        theme_tree2(legend.position = c(0.05, 0.90)) 
p3 + theme(strip.text.x = element_text(colour = "black", 
                                     face = "bold", 
                                     size = 12),
           legend.text = element_text(size = 12),
           legend.title = element_text(size = 12))
