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
# Data frame of covariance and variance components
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
post_cor_vtd_os <- m1_vcv$vtd_os/
  sqrt(m1_vcv$os * m1_vcv$vtd)

post_cor_ssd_os <- m1_vcv$ssd_os/
  sqrt(m1_vcv$os * m1_vcv$ssd)

post_cor_ssd_vtd <- m1_vcv$vtd_ssd/
  sqrt(m1_vcv$vtd * m1_vcv$ssd)

mod1dens <- data.frame(OS_VTD = c(post_cor_vtd_os),
                       OS_SSD = c(post_cor_ssd_os),
                       SSD_VTD = c(post_cor_ssd_vtd))

# Melt data for density plots
mod1dens_melted <- melt(mod1dens)

# Phylogenetic model plot
ggplot(mod1dens_melted, 
       aes(x = value, 
           fill = variable)) + 
  geom_density(alpha = 0.25) +
  geom_vline(xintercept = 0,
             color = "#000000",
             linetype = "dashed") +
  xlim(-1, 1) +
  labs(x = "Posterior",
       y = "Density") +
  scale_fill_discrete(name = "Posterior Correlation") +
  theme_classic(base_size = 15)



# Table of BLUPs (aka "ancestral states" in PGLMM)
df_bf_coefs <- tibble(Trait = attr(colMeans(Mod1$Sol), "names"), 
                      Value = colMeans(Mod1$Sol)) %>%
  separate(Trait, c("Trait","Type","Species"), sep = "\\.", fill = "right") %>% 
  filter(Type == "Species") %>%
  filter(Trait %in% c("traitVTDwSD", "traitSSD", "traitOvulation_Signs")) %>% 
  select(-Type) %>%
  spread(Trait, Value)

# Table of BLUPs with HPD intervals for estimates
df_bf_coefs_error <- tibble(Trait = attr(colMeans(Mod1$Sol), "names"), 
                            Value = colMeans(Mod1$Sol),
                            L_HPD = HPDinterval(Mod1$Sol)[,"lower"],
                            U_HPD = HPDinterval(Mod1$Sol)[,"upper"]) %>%
  separate(Trait, c("Trait","Type","Species"), sep = "\\.", fill = "right") %>% 
  filter(Type == "Species") %>%
  filter(Trait %in% c("traitVTDwSD", "traitSSD", "traitOvulation_Signs")) %>% 
  select(-Type)
