#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Model (PGLMM)   #
#   Results, Diagnostics, Analysis, and Plots             #
#---------------------------------------------------------#

# Loading packages
list_of_packages <- c("tidyverse", "ape", 
                      "MCMCglmm", "phytools",
                      "ggpmisc", "ggpubr",
                      "broom.mixed", "reshape2")
lapply(list_of_packages, library, character.only = TRUE)

# load function for phylogenetic correlation
source("./R/Functions/PhyloCorr.R")

# Loading data and tree from data cleaning
load("./Data/reduced_data.Rdata")
load("./Data/10ktree.Rdata")

# Load outputs from macroevolutionary MCMC models
load("./Results/Data/TrivMacro_Model_4levResp.Rdata")
load("./Results/Data/TrivMacro_Model_2levResp.Rdata")
load("./Results/Data/TrivMacro_Model_Reduced_4levResp.Rdata")


#------ Diagnostics to check posteriors-------#
# Model with 4 level Ovulation signs
summary(Mod1)
plot(Mod1$VCV)
autocorr.plot(Mod1$VCV)
autocorr.plot(Mod1$Sol)
heidel.diag(Mod1$VCV)
heidel.diag(Mod1$Sol)

# Model with binary level Ovulation signs
summary(Mod2)
plot(Mod2$VCV)
autocorr.plot(Mod2$VCV)
autocorr.plot(Mod2$Sol)
heidel.diag(Mod2$VCV)
heidel.diag(Mod2$Sol)

# RAM Model with 4 level Ovulation signs
summary(Mod3)
plot(Mod3$VCV)
autocorr.plot(Mod3$VCV)
autocorr.plot(Mod3$Sol)
heidel.diag(Mod3$VCV)
heidel.diag(Mod3$Sol)

# Table of model results
Mod1_Res <- tidy(Mod1, 
                 effects = c("fixed", "ran_pars"),
                 conf.int = TRUE, 
                 conf.method = "HPDinterval",
                 conf.level = 0.95)

Mod2_Res <- tidy(Mod2, 
                 effects = c("fixed", "ran_pars"),
                 conf.int = TRUE, 
                 conf.method = "HPDinterval",
                 conf.level = 0.95)

Mod3_Res <- tidy(Mod3, 
                 effects = c("fixed", "ran_pars"),
                 conf.int = TRUE, 
                 conf.method = "HPDinterval",
                 conf.level = 0.95)


#-----------------------------------------#
#      Extracting model parameters        #
#      Correlations, Blups                # 
#      Model 1 w/ 4 level OvSig           #
#-----------------------------------------#
# Data frame of phylogenetic covariance and variance components
m1_vcv <- data.frame(
  "vtd_os" = c(as.mcmc(
    Mod1$VCV[, "traitOvulation_Signs:traitVTDwSD.Species"])),
  "ssd_os" = c(as.mcmc(
    Mod1$VCV[, "traitOvulation_Signs:traitSSD.Species"])),
  "vtd_ssd" = c(as.mcmc(
    Mod1$VCV[, "traitSSD:traitVTDwSD.Species"])),
  "vtd" = c(as.mcmc(
    Mod1$VCV[, "traitVTDwSD:traitVTDwSD.Species"])),
  "ssd" = c(as.mcmc(
    Mod1$VCV[, "traitSSD:traitSSD.Species"])),
  "os" = c(as.mcmc(
    Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"]))
  )

m2_vcv <- data.frame(
  "vtd_os" = c(as.mcmc(
    Mod2$VCV[, "traitOvulation_Signs_bin:traitVTDwSD.Species"])),
  "ssd_os" = c(as.mcmc(
    Mod2$VCV[, "traitOvulation_Signs_bin:traitSSD.Species"])),
  "vtd_ssd" = c(as.mcmc(
    Mod2$VCV[, "traitSSD:traitVTDwSD.Species"])),
  "vtd" = c(as.mcmc(
    Mod2$VCV[, "traitVTDwSD:traitVTDwSD.Species"])),
  "ssd" = c(as.mcmc(
    Mod2$VCV[, "traitSSD:traitSSD.Species"])),
  "os" = c(as.mcmc(
    Mod2$VCV[, "traitOvulation_Signs_bin:traitOvulation_Signs_bin.Species"]))
)

m3_vcv <- data.frame(
  "vtd_os" = c(as.mcmc(
    Mod3$VCV[, "traitOvulation_Signs:traitVTDwSD.Species"])),
  "ssd_os" = c(as.mcmc(
    Mod3$VCV[, "traitOvulation_Signs:traitSSD.Species"])),
  "vtd_ssd" = c(as.mcmc(
    Mod3$VCV[, "traitSSD:traitVTDwSD.Species"])),
  "vtd" = c(as.mcmc(
    Mod3$VCV[, "traitVTDwSD:traitVTDwSD.Species"])),
  "ssd" = c(as.mcmc(
    Mod3$VCV[, "traitSSD:traitSSD.Species"])),
  "os" = c(as.mcmc(
    Mod3$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"]))
)


# Phylogenetic correlations - 4 level OS Model with Lambda estimated
m1_vtd_os_corr <- phylo_corr(covar_XY = m1_vcv$vtd_os,
                             var_X = m1_vcv$vtd,
                             var_Y = m1_vcv$os)

m1_ssd_os_corr <- phylo_corr(covar_XY = m1_vcv$ssd_os,
                             var_X = m1_vcv$ssd,
                             var_Y = m1_vcv$os)

m1_ssd_vtd_corr <- phylo_corr(covar_XY = m1_vcv$vtd_ssd,
                              var_X = m1_vcv$ssd,
                              var_Y = m1_vcv$vtd)

# Phylogenetic correlations - 2 level OS Model with Lambda estimated
m2_vtd_os_corr <- phylo_corr(covar_XY = m2_vcv$vtd_os,
                             var_X = m2_vcv$vtd,
                             var_Y = m2_vcv$os)

m2_ssd_os_corr <- phylo_corr(covar_XY = m2_vcv$ssd_os,
                             var_X = m2_vcv$ssd,
                             var_Y = m2_vcv$os)

m2_ssd_vtd_corr <- phylo_corr(covar_XY = m2_vcv$vtd_ssd,
                              var_X = m2_vcv$ssd,
                              var_Y = m2_vcv$vtd)

# Phylogenetic correlations - 4 level OS RAM Model with Lambda fixed
m3_vtd_os_corr <- phylo_corr(covar_XY = m3_vcv$vtd_os,
                             var_X = m3_vcv$vtd,
                             var_Y = m3_vcv$os)

m3_ssd_os_corr <- phylo_corr(covar_XY = m3_vcv$ssd_os,
                             var_X = m3_vcv$ssd,
                             var_Y = m3_vcv$os)

m3_ssd_vtd_corr <- phylo_corr(covar_XY = m3_vcv$vtd_ssd,
                              var_X = m3_vcv$ssd,
                              var_Y = m3_vcv$vtd)

#----------- Ploting posterior correlations ---------------#
# 4 level OS sign with lambda estimated
m1_post_cor_vtd_os <- m1_vcv$vtd_os/
  sqrt(m1_vcv$os * m1_vcv$vtd)

m1_post_cor_ssd_os <- m1_vcv$ssd_os/
  sqrt(m1_vcv$os * m1_vcv$ssd)

m1_post_cor_ssd_vtd <- m1_vcv$vtd_ssd/
  sqrt(m1_vcv$vtd * m1_vcv$ssd)

mod1dens <- data.frame(OS_VTD = c(m1_post_cor_vtd_os),
                       OS_SSD = c(m1_post_cor_ssd_os),
                       SSD_VTD = c(m1_post_cor_ssd_vtd))

# Melt data for density plots
mod1dens_melted <- melt(mod1dens)

# Phylogenetic model plot
p1 <- ggplot(mod1dens_melted, 
       aes(x = value, 
           fill = variable)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0,
             color = "#000000",
             linetype = "dashed") +
  xlim(-1, 1) +
  labs(x = "Posterior",
       y = "Density") +
  scale_fill_manual(name = "Phylogenetic Correlation",
                      labels = c("Ovulation Signs & VTD", 
                                 "Ovulation Signs & SSD", 
                                 "VTD & SSD"),
                      values = c("#D55E00", 
                                 "#0072B2", 
                                 "#009E73")) +
  theme_classic(base_size = 15)
p1 + theme(legend.position = c(0.2, 0.9),
          legend.direction = "vertical",
          legend.background = element_rect(fill = "darkgray"))

# 4 level OS sign with lambda fixed
m3_post_cor_vtd_os <- m3_vcv$vtd_os/
  sqrt(m3_vcv$os * m3_vcv$vtd)

m3_post_cor_ssd_os <- m3_vcv$ssd_os/
  sqrt(m3_vcv$os * m3_vcv$ssd)

m3_post_cor_ssd_vtd <- m3_vcv$vtd_ssd/
  sqrt(m3_vcv$vtd * m3_vcv$ssd)

mod3dens <- data.frame(OS_VTD = c(m3_post_cor_vtd_os),
                       OS_SSD = c(m3_post_cor_ssd_os),
                       SSD_VTD = c(m3_post_cor_ssd_vtd))

# Melt data for density plots
mod3dens_melted <- melt(mod3dens)

# Phylogenetic model plot
p2 <- ggplot(mod3dens_melted, 
       aes(x = value, 
           fill = variable)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0,
             color = "#000000",
             linetype = "dashed") +
  xlim(-0.1, 1) +
  labs(x = "Posterior",
       y = "Density") +
  scale_fill_manual(name = "Phylogenetic Correlation",
                    labels = c("Ovulation Signs & VTD", 
                               "Ovulation Signs & SSD", 
                               "VTD & SSD"),
                    values = c("#D55E00", 
                               "#0072B2", 
                               "#009E73")) +
  theme_classic(base_size = 15)
p2 + theme(legend.position = c(0.3, 0.9),
           legend.direction = "vertical",
           legend.background = element_rect(fill = "darkgray"))

#--------------------------------------------------------------------#
#    Different plots for visualizing magnitude of correlations       #
#--------------------------------------------------------------------#
# Extracting posterior correlations and HPD intervals
m3_corr_tbl <- data.frame(Traits = c("Ovulation Signs & VTD",
                                     "Ovulation Signs & SSD",
                                     "VTD & SSD"),
                          Estimate = c(m3_vtd_os_corr[[1]],
                                       m3_ssd_os_corr[[1]],
                                       m3_ssd_vtd_corr[[1]]),
                          Upper = c(m3_vtd_os_corr[[2]][, "upper"],
                                    m3_ssd_os_corr[[2]][, "upper"],
                                    m3_ssd_vtd_corr[[2]][, "upper"]),
                          Lower = c(m3_vtd_os_corr[[2]][, "lower"],
                                    m3_ssd_os_corr[[2]][, "lower"],
                                    m3_ssd_vtd_corr[[2]][, "lower"]))

p3 <- ggplot(m3_corr_tbl, 
             aes(x = Traits, 
                 y = Estimate)) + 
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) + 
  #geom_hline(yintercept = 0, 
  #           linetype = "dashed",
  #           alpha = 1) +
  scale_x_discrete(limits = c("Ovulation Signs & VTD",
                              "Ovulation Signs & SSD",
                              "VTD & SSD")) +
  labs(x = "Trait combinations",
       y = "Correlation (Estimate +/- 95% HPD)") +
  #ylim(0, 1) +
  coord_flip() + 
  theme_classic(base_size = 15)
p3

#-------------------------- BLUPs analysis --------------------------#
# Table of BLUPs (aka "ancestral states" in PGLMM)
df_bf_coefs <- tibble(Trait = attr(colMeans(Mod3$Sol), "names"), 
                      Value = colMeans(Mod3$Sol)) %>%
  separate(Trait, c("Trait","Type","Species"), 
           sep = "\\.", fill = "right") %>% 
  filter(Type == "Species") %>%
  filter(Trait %in% c("traitVTDwSD", 
                      "traitSSD", 
                      "traitOvulation_Signs")) %>% 
  select(-Type) %>%
  spread(Trait, Value)

# Table of BLUPs with HPD intervals for estimates
df_bf_coefs_error <- tibble(Trait = attr(colMeans(Mod3$Sol), "names"), 
                            Value = colMeans(Mod3$Sol),
                            L_HPD = HPDinterval(Mod3$Sol)[,"lower"],
                            U_HPD = HPDinterval(Mod3$Sol)[,"upper"]) %>%
  separate(Trait, c("Trait","Type","Species"), 
           sep = "\\.", fill = "right") %>% 
  filter(Type == "Species") %>%
  filter(Trait %in% c("traitVTDwSD", 
                      "traitSSD", 
                      "traitOvulation_Signs")) %>% 
  select(-Type)


#--------- Plot of posterior BLUPs ------------------------#
# Slope of relationship between OS and VTD
os_vtd_slope <- Mod3$VCV[,"traitOvulation_Signs:traitVTDwSD.Species"]/
  Mod3$VCV[,"traitVTDwSD:traitVTDwSD.Species"]
mean(os_vtd_slope); HPDinterval(os_vtd_slope)

# Plot for ovulation signals and VTD
p4 <- ggplot(df_bf_coefs, 
       aes(x = traitVTDwSD, 
           y = traitOvulation_Signs, 
           group = Species)) + 
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, 
              slope = mean(os_vtd_slope)) +
  labs(x = "VTD (BLUP)",
       y = "Ovulation Signals (BLUP)") + 
  theme_classic(base_size = 15)
p4

# Slope of relationship between OS and SSD
os_ssd_slope <- Mod3$VCV[,"traitOvulation_Signs:traitSSD.Species"]/
  Mod3$VCV[,"traitSSD:traitSSD.Species"]
mean(os_ssd_slope); HPDinterval(os_ssd_slope)

# Plot for ovulation signals and SSD
p5 <- ggplot(df_bf_coefs, 
             aes(x = traitSSD, 
                 y = traitOvulation_Signs, 
                 group = Species)) + 
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, 
              slope = mean(os_ssd_slope)) +
  labs(x = "SSD (BLUP)",
       y = "Ovulation Signals (BLUP)") + 
  theme_classic(base_size = 15)
p5

# Slope of relationship between OS and SSD
vtd_ssd_slope <- Mod3$VCV[,"traitVTDwSD:traitSSD.Species"]/
  Mod3$VCV[,"traitSSD:traitSSD.Species"]
mean(vtd_ssd_slope); HPDinterval(vtd_ssd_slope)

# Plot for VTD and SSD
p6 <- ggplot(df_bf_coefs, 
             aes(x = traitSSD, 
                 y = traitVTDwSD, 
                 group = Species)) + 
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, 
              slope = mean(vtd_ssd_slope)) +
  labs(x = "SSD (BLUP)",
       y = "VTD (BLUP)") + 
  theme_classic(base_size = 15)
p6


# Looking at clade specific patterns 
tree2 <- multi2di(tree)
plot(tree2, cex = 0.5)
nodelabels(cex = 0.5)

ape_tree  <- extract.clade(tree2, 153)
papio_tree <- extract.clade(tree2, 115)
platy_tree <- extract.clade(tree2, 170)
colob_tree <- extract.clade(tree2, 139)
cercop_tree <- extract.clade(tree2, 104)

apes <- ape_tree$tip.label
papios <- papio_tree$tip.label
platys <- platy_tree$tip.label
colobs <- colob_tree$tip.label
cercs <- cercop_tree$tip.label

df_bf_coefs$clade <- NULL
df_bf_coefs$clade <- ifelse(df_bf_coefs$Species %in% apes,
                            paste0("apes"), 
                            ifelse(df_bf_coefs$Species %in% papios,
                            paste0("papionins"), 
                            ifelse(df_bf_coefs$Species %in% platys,
                                   paste0("platyrrhines"),
                            ifelse(df_bf_coefs$Species %in% colobs,
                                   paste0("colobines"),
                            paste0("guenons")))))

ggplot(data = df_bf_coefs,
       aes(x = traitVTDwSD,
           y = traitOvulation_Signs,
           color = clade)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + 
  geom_abline(intercept = 0, 
              slope = mean(os_vtd_slope), 
              size = 1) +
  labs(y = "Ovulation Signals (BLUP)",
       x = "VTD (BLUP)") +
  theme_classic(base_size = 15)

ggplot(data = df_bf_coefs,
       aes(x = traitSSD,
           y = traitOvulation_Signs,
           color = clade)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + 
  geom_abline(intercept = 0, 
              slope = mean(os_ssd_slope),
              size = 1) +
  labs(y = "Ovulation Signals (BLUP)",
       x = "SSD (BLUP)") +
  theme_classic(base_size = 15)

ggplot(data = df_bf_coefs,
       aes(x = traitSSD,
           y = traitVTDwSD,
           color = clade)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + 
  geom_abline(intercept = 0, 
              slope = mean(vtd_ssd_slope),
              size = 1) +
  labs(y = "VTD (BLUP)",
       x = "SSD (BLUP)") +
  theme_classic(base_size = 15)


# Correlation between trait correlations and speciation rate
cor_ape_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "apes")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "apes")])

cor_pap_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "papionins")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "papionins")])

cor_plat_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "platyrrhines")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "platyrrhines")])

cor_colob_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "colobines")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "colobines")])

cor_guen_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "guenons")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "guenons")])
# OS and SSD
cor_ape_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "apes")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "apes")])

cor_pap_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "papionins")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "papionins")])

cor_plat_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "platyrrhines")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "platyrrhines")])

cor_colob_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "colobines")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "colobines")])

cor_guen_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "guenons")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "guenons")])

df_corrs <- data.frame("Clade" = c("Apes",
                                   "Papionins",
                                   "Platyrrhines",
                                   "Colobines",
                                   "Guenons"),
                       "Speciation" = c(0.4301551,
                                        2.828000e-01,
                                        0.3465664,
                                        0.3008649,
                                        2.741188e-01),
                       "Corr_OSVTD" = c(cor_ape_os_vtd$estimate,
                                        cor_pap_os_vtd$estimate,
                                        cor_plat_os_vtd$estimate,
                                        cor_colob_os_vtd$estimate,
                                        cor_guen_os_vtd$estimate),
                       "Corr_OSSSD" = c(cor_ape_os_ssd$estimate,
                                        cor_pap_os_ssd$estimate,
                                        cor_plat_os_ssd$estimate,
                                        cor_colob_os_ssd$estimate,
                                        cor_guen_os_ssd$estimate))
