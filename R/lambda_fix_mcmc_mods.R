# MCMCglmmRAM not on CRAN - must download from site
install.packages("https://jarrod.bio.ed.ac.uk/MCMCglmmRAM_2.24.tar.gz", 
                 repos = NULL, type = "source")

# Loading packages
list_of_packages <- c("ape", 
                      "MCMCglmmRAM")
lapply(list_of_packages, 
       library, 
       character.only = TRUE)

# Loading data and tree from data cleaning
load("./Data/reduced_data.Rdata")
load("./Data/10ktree.Rdata")

# -------------------------------------------------------------------#
#  Models where lambda is fixed to 1 (i.e. reduced animal model)     #
# -------------------------------------------------------------------#

# Create new column with "animal" variable for RAM
reduced_data$animal <- factor(reduced_data$Species)

PEprior_Reduced <- list(R = list(V = diag(3) * 1e-15, 
                                 nu = 0.002, 
                                 fix = 3), 
                        G = list(G1 = list(V = diag(3), 
                                           nu = 0.002, 
                                           fix = 3)))

set.seed(8675309) # Calls Jenny for a good time

# Run parameters for posterior sample size of N = 2000
nsamp <- 2000
BURN <- 500000; THIN <- 1000; (NITT <- BURN + THIN*nsamp)


Mod3 <- MCMCglmmRAM::MCMCglmm(cbind(log(VTDwSD+1), 
                                    log(SSD), 
                                    Ovulation_Signs) ~ 
                                trait - 1,
                              random = ~ us(trait):Species, 
                              rcov = ~ us(trait):units, 
                              pedigree = tree, 
                              family = c("gaussian", 
                                         "gaussian", 
                                         "threshold"), 
                              data = reduced_data,
                              prior = PEprior_Reduced,
                              pr = TRUE, 
                              pl = TRUE,
                              reduced = TRUE,
                              nitt = NITT,
                              burnin = BURN,
                              thin = THIN)
save(Mod3, file = "TrivMacro_Model_Reduced_4levResp.Rdata")

set.seed(8675309) # Calls Jenny for a good time
Mod7 <- MCMCglmmRAM::MCMCglmm(cbind(log(VTDwSD+1), 
                                    log(SSD), 
                                    Ovulation_Signs_bin) ~ 
                                trait - 1,
                              random = ~ us(trait):Species, 
                              rcov = ~ us(trait):units, 
                              pedigree = tree, 
                              family = c("gaussian", 
                                         "gaussian", 
                                         "threshold"), 
                              data = reduced_data,
                              prior = PEprior_Reduced,
                              pr = TRUE, 
                              pl = TRUE,
                              reduced = TRUE,
                              nitt = NITT,
                              burnin = BURN,
                              thin = THIN)
save(Mod7, file = "BivMacro_Model_Reduced_2levResp.Rdata")

# Reduced models where VTD and SSD lambda are also fixed at 1
PEprior_Reduced1 <- list(R = list(V = diag(3) * 1e-15, 
                                  nu = 0.002, 
                                  fix = 1), 
                         G = list(G1 = list(V = diag(3), 
                                            nu = 0.002)))
set.seed(8675309) # Calls Jenny for a good time

Mod3.1 <- MCMCglmmRAM::MCMCglmm(cbind(log(VTDwSD+1), 
                                      log(SSD), 
                                      Ovulation_Signs) ~ 
                                  trait - 1,
                                random = ~ corg(trait):Species, 
                                rcov = ~ corg(trait):units, 
                                pedigree = tree, 
                                family = c("gaussian", 
                                           "gaussian", 
                                           "threshold"), 
                                data = reduced_data,
                                prior = PEprior_Reduced1,
                                pr = TRUE, 
                                pl = TRUE,
                                reduced = TRUE,
                                nitt = NITT,
                                burnin = BURN,
                                thin = THIN)
save(Mod3.1, file = "TrivMacro_Model_Reduced_4levRespFixAll.Rdata")


set.seed(8675309) # Calls Jenny for a good time

Mod7.1 <- MCMCglmmRAM::MCMCglmm(cbind(log(VTDwSD+1), 
                                      log(SSD), 
                                      Ovulation_Signs_bin) ~ 
                                  trait - 1,
                                random = ~ corg(trait):Species, 
                                rcov = ~ corg(trait):units, 
                                pedigree = tree, 
                                family = c("gaussian", 
                                           "gaussian", 
                                           "threshold"), 
                                data = reduced_data,
                                prior = PEprior_Reduced1,
                                pr = TRUE, 
                                pl = TRUE,
                                reduced = TRUE,
                                nitt = NITT,
                                burnin = BURN,
                                thin = THIN)
save(Mod7.1, file = "TrivMacro_Model_Reduced_BinRespFixAll.Rdata")
