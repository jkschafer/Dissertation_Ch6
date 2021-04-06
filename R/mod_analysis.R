#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Model (PGLMM)   #
#   Results, Diagnostics, Analysis, and Plots             #
#---------------------------------------------------------#

# Loading packages
list_of_packages <- c("tidyverse", "ape", 
                      "MCMCglmm", "phytools",
                      "ggpmisc", "ggpubr",
                      "broom.mixed")
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
#      Lambda, Correlations, Blups        #
#      Model 1 w/ 4 level OvSig           #
#-----------------------------------------#

# Phylogenetic signal (Pagel's lambda)
phyloSig_OvSigs <- Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"]/
  (Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"] + 
     Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.units"])
posterior.mode(phyloSig_OvSigs); HPDinterval(phyloSig_OvSigs); plot(phyloSig_OvSigs)

phyloSig_VTD <- Mod1$VCV[, "traitVTDwSD:traitVTDwSD.Species"]/
  (Mod1$VCV[, "traitVTDwSD:traitVTDwSD.Species"] + 
     Mod1$VCV[, "traitVTDwSD:traitVTDwSD.units"])
posterior.mode(phyloSig_VTD); HPDinterval(phyloSig_VTD); plot(phyloSig_VTD)

phyloSig_SSD <- Mod1$VCV[, "traitSSD:traitSSD.Species"]/
  (Mod1$VCV[, "traitSSD:traitSSD.Species"] + 
     Mod1$VCV[, "traitSSD:traitSSD.units"])
posterior.mode(phyloSig_SSD); HPDinterval(phyloSig_SSD); plot(phyloSig_SSD)

# Phylogenetic correlations
phyloCorr_VTD_OS <- Mod1$VCV[, "traitOvulation_Signs:traitVTDwSD.Species"]/
  (sqrt(Mod1$VCV[, "traitVTDwSD:traitVTDwSD.Species"] * 
          Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"]))
mean(phyloCorr_VTD_OS); HPDinterval(phyloCorr_VTD_OS); plot(phyloCorr_VTD_OS)

phyloCorr_VTD_SSD <- Mod1$VCV[, "traitSSD:traitVTDwSD.Species"]/
  (sqrt(Mod1$VCV[, "traitVTDwSD:traitVTDwSD.Species"] * 
          Mod1$VCV[, "traitSSD:traitSSD.Species"]))
mean(phyloCorr_VTD_SSD); HPDinterval(phyloCorr_VTD_SSD); plot(phyloCorr_VTD_SSD)

phyloCorr_SSD_OS <- Mod1$VCV[, "traitOvulation_Signs:traitSSD.Species"]/
  (sqrt(Mod1$VCV[, "traitSSD:traitSSD.Species"] * 
          Mod1$VCV[, "traitOvulation_Signs:traitOvulation_Signs.Species"]))
mean(phyloCorr_SSD_OS); HPDinterval(phyloCorr_SSD_OS); plot(phyloCorr_SSD_OS)

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
