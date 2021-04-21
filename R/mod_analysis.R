#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Model (PGLMM)   #
#   Results, Diagnostics, Analysis, and Plots             #
#---------------------------------------------------------#

# Loading packages
list_of_packages <- c("tidyverse", "ape", 
                      "MCMCglmm", "phytools",
                      "ggpmisc", "ggpubr",
                      "broom.mixed", "reshape2",
                      "factoextra", "corrplot",
                      "CCA", "yacca", "CCP")
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
load("./Results/Data/UniMacro_Model_OvSign.Rdata")
load("./Results/Data/UniMacro_Model_VTD.Rdata")
load("./Results/Data/UniMacro_Model_SSD.Rdata")

#--------------------------------------------------------------# 
#    Calculating Pagel's lambda from univariate analyses
#--------------------------------------------------------------# 
# Ovulation signals
OvSign_lambda <- Mod6$VCV[, "Species"]/
  (Mod6$VCV[, "Species"] + 
     Mod6$VCV[, "units"])
mean(OvSign_lambda) 
HPDinterval(OvSign_lambda) 
plot(OvSign_lambda)

# VTD
VTD_lambda <- Mod5$VCV[, "Species"]/
  (Mod5$VCV[, "Species"] + 
     Mod5$VCV[, "units"])
mean(VTD_lambda) 
HPDinterval(VTD_lambda) 
plot(VTD_lambda)

# SSD
SSD_lambda <- Mod4$VCV[, "Species"]/
  (Mod4$VCV[, "Species"] + 
     Mod4$VCV[, "units"])
mean(SSD_lambda) 
HPDinterval(SSD_lambda) 
plot(SSD_lambda)

lambda_tbl <- data.frame(Trait = c("Ovulation_Signs",
                                   "VTD",
                                   "SSD"),
                         Lambda = c(mean(OvSign_lambda),
                                    mean(VTD_lambda),
                                    mean(SSD_lambda)),
                         Lower = c(HPDinterval(OvSign_lambda)[, "lower"],
                                   HPDinterval(VTD_lambda)[, "lower"],
                                   HPDinterval(SSD_lambda)[, "lower"]),
                         Upper = c(HPDinterval(OvSign_lambda)[, "upper"],
                                   HPDinterval(VTD_lambda)[, "upper"],
                                   HPDinterval(SSD_lambda)[, "upper"]))
stargazer(lambda_tbl, type = "latex", summary = F)
# Plot lambdas
p_lamb <- ggplot(lambda_tbl, 
             aes(x = Trait, 
                 y = Lambda)) + 
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) + 
  geom_hline(yintercept = 0, 
             linetype = "solid",
             alpha = 1) +
  scale_x_discrete(limits = c("Ovulation_Signs",
                              "VTD",
                              "SSD"),
                   labels = c("Ovulation Signs",
                              "VTD",
                              "SSD")) +
  labs(x = "",
       y = expression("Pagel's" ~ lambda ~ "+/- 95% HPD")) +
  ylim(0, 1) +
  coord_flip() + 
  theme_classic(base_size = 15)
p_lamb + theme(axis.line.y = element_line(linetype = "blank"),
               axis.text.y = element_text(angle = 45))
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
stargazer(m3_corr_tbl, type = "latex", summary = F)

p3 <- ggplot(m3_corr_tbl, 
             aes(x = Traits, 
                 y = Estimate)) + 
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) + 
  geom_hline(yintercept = 0, 
             linetype = "solid",
             alpha = 1) +
  scale_x_discrete(limits = c("Ovulation Signs & VTD",
                              "Ovulation Signs & SSD",
                              "VTD & SSD")) +
  labs(x = "",
       y = "Phylogenetic Correlation (Estimate +/- 95% HPD)") +
  ylim(0, 1) +
  coord_flip() + 
  theme_classic(base_size = 15)
p3 + theme(axis.line.y = element_line(linetype = "blank"),
           axis.text.y = element_text(angle = 45))

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

gr_ape_tree  <- extract.clade(tree2, 164)
l_ape_tree <- extract.clade(tree2, 154)
papio_tree <- extract.clade(tree2, 115)
colob_tree <- extract.clade(tree2, 139)
cercop_tree <- extract.clade(tree2, 104)

platy_tree <- extract.clade(tree2, 170)
callit_tree <- extract.clade(tree2, 181)
cebid_tree <- extract.clade(tree2, 194)


gr_ape <- gr_ape_tree$tip.label
l_ape <- l_ape_tree$tip.label
papios <- papio_tree$tip.label
platys <- platy_tree$tip.label
colobs <- colob_tree$tip.label
cercs <- cercop_tree$tip.label
#pithec <- tree2$tip.label[59]

df_bf_coefs$clade <- NULL
df_bf_coefs$clade <- ifelse(df_bf_coefs$Species %in% gr_ape,
                            paste0("Hominidae"),
                            ifelse(df_bf_coefs$Species %in% l_ape,
                                   paste0("Hylobatidae"),
                            ifelse(df_bf_coefs$Species %in% papios,
                            paste0("Papioini"), 
                            ifelse(df_bf_coefs$Species %in% platys,
                                   paste0("Platyrrhini"),
                            ifelse(df_bf_coefs$Species %in% colobs,
                                   paste0("Colobinae"),
                            paste0("Cercopithecini"))))))

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
cor_g_ape_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Hominidae")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Hominidae")])

cor_l_ape_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Hylobatidae")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Hylobatidae")])

cor_pap_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Papioini")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Papioini")])

cor_plat_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Platyrrhini")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Platyrrhini")])

cor_colob_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Colobinae")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Colobinae")])

cor_guen_os_vtd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Cercopithecini")], df_bf_coefs$traitVTDwSD[which(
    df_bf_coefs$clade == "Cercopithecini")])
# OS and SSD
cor_g_ape_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Hominidae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Hominidae")])

cor_l_ape_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Hylobatidae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Hylobatidae")])

cor_pap_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Papioini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Papioini")])

cor_plat_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Platyrrhini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Platyrrhini")])

cor_colob_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Colobinae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Colobinae")])

cor_guen_os_ssd <- cor.test(df_bf_coefs$traitOvulation_Signs[which(
  df_bf_coefs$clade == "Cercopithecini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Cercopithecini")])

# VTD and SSD
cor_g_ape_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Hominidae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Hominidae")])

cor_l_ape_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Hylobatidae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Hylobatidae")])

cor_pap_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Papioini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Papioini")])

cor_plat_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Platyrrhini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Platyrrhini")])

cor_colob_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Colobinae")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Colobinae")])

cor_guen_vtd_ssd <- cor.test(df_bf_coefs$traitVTDwSD[which(
  df_bf_coefs$clade == "Cercopithecini")], df_bf_coefs$traitSSD[which(
    df_bf_coefs$clade == "Cercopithecini")])


df_corrs <- data.frame(Clade = c("Hominidae",
                                 "Hylobatidae",
                                 "Papioini",
                                 "Platyrrhini",
                                 "Colobinae",
                                 "Cercopithecini"),
                       Corr_OSVTD = c(cor_g_ape_os_vtd$estimate,
                                      cor_l_ape_os_vtd$estimate,
                                      cor_pap_os_vtd$estimate,
                                      cor_plat_os_vtd$estimate,
                                      cor_colob_os_vtd$estimate,
                                      cor_guen_os_vtd$estimate),
                       LCI_Corr_OSVTD = c(cor_g_ape_os_vtd$conf.int[1],
                                          cor_l_ape_os_vtd$conf.int[1],
                                          cor_pap_os_vtd$conf.int[1],
                                          cor_plat_os_vtd$conf.int[1],
                                          cor_colob_os_vtd$conf.int[1],
                                          cor_guen_os_vtd$conf.int[1]),
                       UCI_Corr_OSVTD = c(cor_g_ape_os_vtd$conf.int[2],
                                          cor_l_ape_os_vtd$conf.int[2],
                                          cor_pap_os_vtd$conf.int[2],
                                          cor_plat_os_vtd$conf.int[2],
                                          cor_colob_os_vtd$conf.int[2],
                                          cor_guen_os_vtd$conf.int[2]),
                       Corr_OSSSD = c(cor_g_ape_os_ssd$estimate,
                                      cor_l_ape_os_ssd$estimate,
                                      cor_pap_os_ssd$estimate,
                                      cor_plat_os_ssd$estimate,
                                      cor_colob_os_ssd$estimate,
                                      cor_guen_os_ssd$estimate),
                       LCI_Corr_OSSSD = c(cor_g_ape_os_ssd$conf.int[1],
                                          cor_l_ape_os_ssd$conf.int[1],
                                          cor_pap_os_ssd$conf.int[1],
                                          cor_plat_os_ssd$conf.int[1],
                                          cor_colob_os_ssd$conf.int[1],
                                          cor_guen_os_ssd$conf.int[1]),
                       UCI_Corr_OSSSD = c(cor_g_ape_os_ssd$conf.int[2],
                                          cor_l_ape_os_ssd$conf.int[2],
                                          cor_pap_os_ssd$conf.int[2],
                                          cor_plat_os_ssd$conf.int[2],
                                          cor_colob_os_ssd$conf.int[2],
                                          cor_guen_os_ssd$conf.int[2]),
                       Corr_VTDSSD = c(cor_g_ape_vtd_ssd$estimate,
                                       cor_l_ape_vtd_ssd$estimate,
                                       cor_pap_vtd_ssd$estimate,
                                       cor_plat_vtd_ssd$estimate,
                                       cor_colob_vtd_ssd$estimate,
                                       cor_guen_vtd_ssd$estimate),
                       LCI_Corr_VTDSSD = c(cor_g_ape_vtd_ssd$conf.int[1],
                                           cor_l_ape_vtd_ssd$conf.int[1],
                                           cor_pap_vtd_ssd$conf.int[1],
                                           cor_plat_vtd_ssd$conf.int[1],
                                           cor_colob_vtd_ssd$conf.int[1],
                                           cor_guen_vtd_ssd$conf.int[1]),
                       UCI_Corr_VTDSSD = c(cor_g_ape_vtd_ssd$conf.int[2],
                                           cor_l_ape_vtd_ssd$conf.int[2],
                                           cor_pap_vtd_ssd$conf.int[2],
                                           cor_plat_vtd_ssd$conf.int[2],
                                           cor_colob_vtd_ssd$conf.int[2],
                                           cor_guen_vtd_ssd$conf.int[2]))


# Merge clade rates and correlations
load("./Results/Data/clade_rates_as.Rdata")

corr_rates_df <- merge(x = df_corrs, 
                       y = clade_rates, 
                       by = "Clade")

ggplot(data = corr_rates_df,
       aes(x = Corr_OSVTD,
           y = Sp_Rate,
           color = Clade)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  geom_errorbar(aes(ymin = L_Sp_Rate,
                    ymax = U_Sp_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_OSVTD,
                     xmax = UCI_Corr_OSVTD)) +
  labs(y = "Speciation Rate +/- 95% CI",
       x = "Ovulation Signal-VTD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)

ggplot(data = corr_rates_df,
       aes(x = Corr_OSSSD,
           y = Sp_Rate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  geom_errorbar(aes(ymin = L_Sp_Rate,
                    ymax = U_Sp_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_OSSSD,
                     xmax = UCI_Corr_OSSSD)) +
  labs(y = "Speciation Rate +/- 95% CI",
       x = "Ovulation Signal-SSD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)

ggplot(data = corr_rates_df,
       aes(x = Corr_VTDSSD,
           y = Sp_Rate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  geom_errorbar(aes(ymin = L_Sp_Rate,
                    ymax = U_Sp_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_VTDSSD,
                     xmax = UCI_Corr_VTDSSD)) +
  labs(y = "Speciation Rate +/- 95% CI",
       x = "VTD-SSD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)

# Extinction rates
ggplot(data = corr_rates_df,
       aes(x = Corr_OSVTD,
           y = Ex_Rate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T) +
  geom_errorbar(aes(ymin = L_Ex_Rate,
                    ymax = U_Ex_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_OSVTD,
                     xmax = UCI_Corr_OSVTD)) +
  labs(y = "Extinction Rate +/- 95% CI",
       x = "Ovulation Signal-VTD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)

ggplot(data = corr_rates_df,
       aes(x = Corr_OSSSD,
           y = Ex_Rate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T) +
  geom_errorbar(aes(ymin = L_Ex_Rate,
                    ymax = U_Ex_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_OSSSD,
                     xmax = UCI_Corr_OSSSD)) +
  labs(y = "Extinction Rate +/- 95% CI",
       x = "Ovulation Signal-SSD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)

ggplot(data = corr_rates_df,
       aes(x = Corr_VTDSSD,
           y = Ex_Rate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = T) +
  geom_errorbar(aes(ymin = L_Ex_Rate,
                    ymax = U_Ex_Rate)) + 
  geom_errorbarh(aes(xmin = LCI_Corr_VTDSSD,
                     xmax = UCI_Corr_VTDSSD)) +
  labs(y = "Extinction Rate +/- 95% CI",
       x = "VTD-SSD Correlation +/- 95% CI") +
  theme_classic(base_size = 15)



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
  select(traitOvulation_Signs,
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
# Attempt to extract species MCMC samples for correlations
# instead of correlation of blup means
dftest <- Mod3$Sol[, grep(pattern = "Species*", 
                          x = colnames(Mod3$Sol))]



