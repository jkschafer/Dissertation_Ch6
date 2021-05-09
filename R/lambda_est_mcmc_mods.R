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
Mod1 <- MCMCglmm(cbind(log(VTDwSD + 1), 
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
Mod2 <- MCMCglmm(cbind(log(VTDwSD + 1), 
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

# Univariate models for each variable to calculate lambda
PEprior2 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(1) * 0.02, 
                                    nu = 1, 
                                    alpha.mu = rep(0, 1), 
                                    alpha.V = diag(1) * 1000)))


# Run parameters for posterior sample size of N = 2000
reduced_data$animal <- reduced_data$Species
nsamp <- 2000
BURN <- 500000; THIN <- 1000; (NITT <- BURN + THIN*nsamp)

set.seed(8675309) # Calls Jenny for a good time
Mod4 <- MCMCglmm(log(SSD) ~ 1,
                 random = ~ animal,
                 rcov = ~ units,
                 pedigree = tree, 
                 family = "gaussian", 
                 data = reduced_data,
                 prior = PEprior2,
                 pr = TRUE, 
                 DIC = TRUE,
                 nitt = NITT,
                 burnin = BURN,
                 thin = THIN)
save(Mod4, file = "UniMacro_Model_SSD.Rdata")

set.seed(8675309) # Calls Jenny for a good time
Mod5 <- MCMCglmm(log(VTDwSD + 1) ~ 1,
                 random = ~ animal,
                 rcov = ~ units,
                 pedigree = tree, 
                 family = "gaussian", 
                 data = reduced_data,
                 prior = PEprior2,
                 pr = TRUE, 
                 DIC = TRUE,
                 nitt = NITT,
                 burnin = BURN,
                 thin = THIN)
save(Mod5, file = "UniMacro_Model_VTD.Rdata")


PEprior3 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(1) * 0.02, 
                                    nu = 1, 
                                    alpha.mu = rep(0, 1), 
                                    alpha.V = diag(1) * 1000)))
set.seed(8675309) # Calls Jenny for a good time
Mod6 <- MCMCglmm(Ovulation_Signs ~ 1,
                 random = ~ animal,
                 rcov = ~ units,
                 pedigree = tree, 
                 family = "threshold", 
                 data = reduced_data,
                 prior = PEprior3,
                 pr = TRUE,
                 pl = TRUE,
                 trunc = TRUE,
                 DIC = TRUE,
                 nitt = NITT,
                 burnin = BURN,
                 thin = THIN)
save(Mod6, file = "UniMacro_Model_OvSign.Rdata")


vtd <- log(reduced_data$VTDwSD + 1)
names(vtd) <- reduced_data$Species
sig_lamb_vtd <- phylosig(tree, vtd, 
                         method = "lambda", 
                         test = TRUE)
sig_lamb_vtd

sig_k_vtd <- phylosig(tree, vtd, 
                      method = "K", 
                      test = TRUE, 
                      nsim = 10000)
sig_k_vtd

ssd <- log(reduced_data$SSD)
names(ssd) <- reduced_data$Species
sig_lamb_ssd <- phylosig(tree, ssd, 
                         method = "lambda", 
                         test = TRUE)
sig_lamb_ssd

sig_k_ssd <- phylosig(tree, ssd, 
                      method = "K", 
                      test = TRUE, 
                      nsim = 10000)
sig_k_ssd

