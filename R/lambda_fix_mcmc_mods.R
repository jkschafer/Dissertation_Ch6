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
                              nitt = 13000,
                              burnin = 3000,
                              thin = 10)
save(mod3, file = "TrivMacro_Model_Reduced_4levResp.Rdata")
