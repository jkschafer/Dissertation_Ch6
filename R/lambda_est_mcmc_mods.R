#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Models (PGLMM)  #
#   Constructing and running models                       #
#---------------------------------------------------------#

# Loading packages
list_of_packages <- c("ape", 
                      "MCMCglmm")
lapply(list_of_packages, 
       library, 
       character.only = TRUE)

# Loading data and tree from data cleaning
load("./Data/reduced_data.Rdata")
load("./Data/10ktree.Rdata")

# --------------- Phylogenetic mixed models -------------------------#
Ainv_phylo <- inverseA(tree, nodes="ALL", scale = FALSE)$Ainv

#Parameter expanded prior
PEprior <- list(R = list(V = diag(3), nu = 3, fix = 3),
                G = list(G1 = list(V = diag(3) * 0.02, 
                                   nu = 4, 
                                   alpha.mu = rep(0, 3), 
                                   alpha.V = diag(3) * 1000)))

set.seed(8675309) # Calls Jenny for a good time

# Run parameters for posterior sample size of N = 2000
nsamp <- 2000
BURN <- 500000; THIN <- 1000; (NITT <- BURN + THIN*nsamp)

# -------------------------------------------------------------------#
#       Models where lambda is estimated instead of fixed            #
# -------------------------------------------------------------------#

# Trivariate model with 4 level ovulation signs response
Mod1 <- MCMCglmm:MCMCglmm(cbind(log(VTDwSD + 1), 
                                log(SSD), 
                                Ovulation_Signs) ~ 
                            trait - 1,
                          random = ~ us(trait):Species, 
                          rcov = ~ us(trait):units, 
                          ginverse = list(Species = Ainv_phylo), 
                          family = c("gaussian", 
                                     "gaussian", 
                                     "threshold"), 
                          data = reduced_data,
                          prior = PEprior,
                          pr = TRUE, 
                          pl = TRUE,
                          slice = TRUE, 
                          trunc = TRUE, 
                          DIC = TRUE,
                          nitt = NITT,
                          burnin = BURN,
                          thin = THIN)
save(Mod1, file = "TrivMacro_Model_4levResp.Rdata")

# Model with binary ovulatory sign variable -> 0 = absent; 1 = present
Mod2 <- MCMCglmm:MCMCglmm(cbind(log(VTDwSD + 1), 
                                log(SSD), 
                                Ovulation_Signs_bin) ~ 
                            trait - 1,
                          random = ~ us(trait):Species, 
                          rcov = ~ us(trait):units, 
                          ginverse = list(Species = Ainv_phylo), 
                          family = c("gaussian", 
                                     "gaussian", 
                                     "threshold"), 
                          data = reduced_data,
                          prior = PEprior,
                          pr = TRUE, 
                          pl = TRUE,
                          slice = TRUE,
                          trunc = TRUE,
                          DIC = TRUE,
                          nitt = NITT,
                          burnin = BURN,
                          thin = THIN)
save(Mod2, file = "TrivMacro_Model_2levResp.Rdata")
