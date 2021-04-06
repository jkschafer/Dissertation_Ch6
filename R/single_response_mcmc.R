#---------------------------------------------------------#
#   Phylogenetic Generalized Linear Mixed Models (PGLMM)  #
#   Constructing and running single-response models       #
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
PEprior <- list(R = list(V = diag(1), fix = 1),
                G = list(G1 = list(V = diag(1), 
                                   nu = 0.002)))

set.seed(8675309) # Calls Jenny for a good time

# Run parameters for posterior sample size of N = 2000
nsamp <- 2000
BURN <- 500000; THIN <- 1000; (NITT <- BURN + THIN*nsamp)

Mod4 <- MCMCglmm(Ovulation_Signs ~ log(SSD) + log(VTDwSD + 1),
                 random = ~ Species,
                 rcov = ~ units,
                 prior = PEprior,
                 data = reduced_data,
                 ginverse = list(Species = Ainv_phylo),
                 family = "threshold",
                 trunc = TRUE,
                 pr = TRUE,
                 pl = TRUE,
                 nitt = NITT,
                 thin = THIN,
                 burnin = BURN)
save(Mod4, file = "UnivMacro_Model_4levResp.Rdata")
