library(BAMMtools)
library(ape)
library(phytools)
library(coda)

# Load data
tree <- read.tree("./Results/Data/BAMM_out/primatetree.tre")
edata <- getEventData(tree, 
                      eventdata = "./Results/Data/Arbour_Santana_Data/BAMM filesevent_data.txt", 
                      burnin=0.1)
mcmcout <- read.csv("./Results/Data/BAMM_out/mcmc_out.txt", 
                    header=T)

# Assess convergence
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs['X'] / post_probs['Y']

plotPrior(mcmcout, expectedNumberOfShifts=1)


plot.bammdata(edata, lwd=2)
plot.bammdata(edata, lwd=2, legend=T)

index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2)
addBAMMshifts(e2, cex=2)

css <- credibleShiftSet(edata, 
                        expectedNumberOfShifts=1, 
                        threshold=5, 
                        set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

best <- getBestShiftConfiguration(edata, 
                                  expectedNumberOfShifts=1, 
                                  threshold=5)
plot.bammdata(best, lwd=1.25)
addBAMMshifts(best, cex=2)

allrates <- getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

plot.phylo(tree, cex = 0.5)
nodelabels(cex = 0.5)

# OWM
plot(extract.clade(tree, 426), cex = 0.5) #Hominidae
plot(extract.clade(tree, 414), cex = 0.5) #Hylobatidae
plot(extract.clade(tree, 306), cex = 0.5) #Cercopithecidae

plot(extract.clade(tree, 308), cex = 0.5) #Cercopithecini
plot(extract.clade(tree, 337), cex = 0.5) #Papionini
plot(extract.clade(tree, 375), cex = 0.5) #Colobinae
# NWMs
plot(extract.clade(tree, 437), cex = 0.5) # Platyrrhines
plot(extract.clade(tree, 462), cex = 0.5) #Callitrichids
plot(extract.clade(tree, 494), cex = 0.5) #Pitheciidae
plot(extract.clade(tree, 486), cex = 0.5) #Cebidae
plot(extract.clade(tree, 486), cex = 0.5) #Cebidae
plot(extract.clade(tree, 453), cex = 0.5) #Aotidae
plot(extract.clade(tree, 439), cex = 0.5) #Atelidae

#------------ Speciation rates -----------------------#
# Papionini tribe
papio_rates <- getCladeRates(edata, node = 337)
mean(papio_rates$lambda)
quantile(papio_rates$lambda, c(0.05, 0.95))

# Hominidae
homini_rates <- getCladeRates(edata, node = 426)
mean(homini_rates$lambda)
quantile(homini_rates$lambda, c(0.05, 0.95))

# Hylobatidae
hylob_rates <- getCladeRates(edata, node = 414)
mean(hylob_rates$lambda)
quantile(hylob_rates$lambda, c(0.05, 0.95))

# Cercopithecini
cerco_rates <- getCladeRates(edata, node = 308)
mean(cerco_rates$lambda)
quantile(cerco_rates$lambda, c(0.05, 0.95))

# All primates except papioini
nonpapiorate <- getCladeRates(edata, node = 337, nodetype = "exclude")
mean(nonpapiorate$lambda)
quantile(nonpapiorate$lambda, c(0.05, 0.95))

plotRateThroughTime(edata, ratetype="speciation")

plotRateThroughTime(edata, node = 337, 
                    nodetype="include")
plotRateThroughTime(edata, node = 426, 
                    nodetype="include")


cmat <- getCohortMatrix(edata)
cohorts(cmat, edata)

# Tip speciation rates
rates <- data.frame(Species = c(tree$tip.label),
                    SpRate = c(edata$meanTipLambda),
                    ExRate = c(edata$meanTipMu))





