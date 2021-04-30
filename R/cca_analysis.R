#---------- Species specific rates -----------------#
load("./Results/Data/tip_rates_ah.Rdata")

not_in_spdat <- tip_rates[!tip_rates$Species %in% 
                            df_bf_coefs$Species,]
unique(not_in_spdat$Species) # many
in_spdata <- tip_rates[tip_rates$Species %in%
                         df_bf_coefs$Species,]
unique(in_spdata)

df_bf_coefs_divrates <- merge(df_bf_coefs, tip_rates,
                              by.x = "Species",
                              by.y = "Species")

pca_df <- df_bf_coefs_divrates %>%
  dplyr::select(traitOvulation_Signs,
                traitSSD, traitVTDwSD, 
                SpRate, ExRate)
rownames(pca_df) <- df_bf_coefs_divrates$Species

# PCA analysis of blups and diversification rates
res.pca <- prcomp(pca_df, scale = TRUE)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 55))

# Variable contributions
fviz_contrib(res.pca, choice = "var", axes = 1, top = 5)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 5)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 5)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 5)
fviz_contrib(res.pca, choice = "var", axes = 5, top = 5)

# Individual contributions
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

# Correlation plot
corrplot(res.pca$x, is.corr=FALSE)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(df_bf_coefs_divrates$clade), # color by groups
             palette = c("#56B4E9", "#009E73", "#F0E442", 
                         "#0072B2", "#D55E00", "#CC79A7"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups") + 
  theme_classic(base_size = 15)


#-------------------- Canonical Correlation Analysis ----------------------#
pca_df
traits <- pca_df[, 1:3]
rates <- pca_df[, 4:5]
res_cc <- cancor(traits, rates)
res_cc_2 <- cca(traits, rates)
res_cc_3 <- cc(traits, rates)
F.test.cca(res_cc_2)

CC1_X <- as.matrix(traits) %*% res_cc$xcoef[, 1]
CC1_Y <- as.matrix(rates) %*% res_cc$ycoef[, 1]
CC2_X <- as.matrix(traits) %*% res_cc$xcoef[, 2]
CC2_Y <- as.matrix(rates) %*% res_cc$ycoef[, 2]

cca_df <- pca_df %>% 
  mutate(CC1_X = CC1_X,
         CC1_Y = CC1_Y,
         CC2_X = CC2_X,
         CC2_Y = CC2_Y)

cca_df$Species <- rownames(cca_df)
cca_df$clade <- NULL
cca_df$clade <- ifelse(cca_df$Species %in% gr_ape,
                       paste0("Hominidae"),
                       ifelse(cca_df$Species %in% l_ape,
                              paste0("Hylobatidae"),
                              ifelse(cca_df$Species %in% papios,
                                     paste0("Papioini"), 
                                     ifelse(cca_df$Species %in% platys,
                                            paste0("Platyrrhini"),
                                            ifelse(cca_df$Species %in% colobs,
                                                   paste0("Colobinae"),
                                                   paste0("Cercopithecini"))))))


# Plotting raw correlation matrices
ggpairs(cca_df, columns = 1:3, # trait blups
        ggplot2::aes(colour = clade)) + 
  theme_classic(base_size = 15)

ggpairs(cca_df, columns = 4:5, # diversification rates
        ggplot2::aes(colour = clade)) + 
  theme_classic(base_size = 15)

cca_df %>% 
  ggplot(aes(x = CC1_X,
             y = CC1_Y,
             color = clade))+
  geom_point() +
  theme_classic(base_size = 15)

cca_df %>% 
  ggplot(aes(x = CC2_X,
             y = CC2_Y,
             color = clade))+
  geom_point() +
  theme_classic(base_size = 15)

cca_df %>% 
  ggplot(aes(x = clade,
             y = CC1_X)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.15) +
  theme(legend.position = "none") +
  theme_classic(base_size = 15)

cca_df %>% 
  ggplot(aes(x = clade,
             y = CC1_Y)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.15) +
  theme(legend.position = "none") +
  theme_classic(base_size = 15)

comput(traits, rates, res_cc_3)
rho <- res_cc_3$cor
n <- dim(traits)[1]
p <- length(traits)
q <- length(rates)
p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Hotelling")
p.asym(rho, n, p, q, tstat = "Pillai")
p.asym(rho, n, p, q, tstat = "Roy")