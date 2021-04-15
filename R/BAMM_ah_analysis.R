library(BAMMtools)
library(ape)
library(phytools)
library(coda)

# Load data
tree <- read.tree("./Results/Data/BAMM_out/AH_phy/springer_primates_as_tree.phy")
edata <- getEventData(tree, 
                      eventdata = "./Results/Data/BAMM_out/AH_phy/event_data.txt", 
                      burnin = 0.1)
mcmcout <- read.csv("./Results/Data/BAMM_out/AH_phy/mcmc_out.txt", 
                    header = TRUE)

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

bamm.primates <- plot.bammdata(edata, 
                               lwd=2, 
                               labels = T, 
                               cex = 0.35,
                               legend = T)
addBAMMshifts(edata, cex = 2)


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
par(mfrow = c(4, 3))
plot.credibleshiftset(css)

best <- getBestShiftConfiguration(edata, 
                                  expectedNumberOfShifts=1, 
                                  threshold=5)
plot.bammdata(best, lwd=1.25)
addBAMMshifts(best, cex=2)

allrates <- getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

# Finding clade lables for taxa
plot.phylo(tree, cex = 0.5)
nodelabels(cex = 0.5)

# OWM
plot(extract.clade(tree, 548), cex = 0.5) #Hominoidea
plot(extract.clade(tree, 549), cex = 0.5) #Hominidae
plot(extract.clade(tree, 555), cex = 0.5) #Hylobatidae
plot(extract.clade(tree, 570), cex = 0.5) #Cercopithecidae

plot(extract.clade(tree, 609), cex = 0.5) #Cercopithecini
plot(extract.clade(tree, 572), cex = 0.5) #Papionini
plot(extract.clade(tree, 633), cex = 0.5) #Colobinae
# NWMs
plot(extract.clade(tree, 461), cex = 0.5) # Platyrrhines
plot(extract.clade(tree, 512), cex = 0.5) #Callitrichids
plot(extract.clade(tree, 462), cex = 0.5) #Pitheciidae
plot(extract.clade(tree, 501), cex = 0.5) #Cebidae
plot(extract.clade(tree, 541), cex = 0.5) #Aotidae
plot(extract.clade(tree, 480), cex = 0.5) #Atelidae

#------------ Speciation rates -----------------------#
#Cercopithecidae
cercd_rates <- getCladeRates(edata, node = 570)
mean(cercd_rates$lambda)
quantile(cercd_rates$lambda, c(0.05, 0.95))

#Platyrrhines
platy_rates <- getCladeRates(edata, node = 461)
mean(platy_rates$lambda)
quantile(platy_rates$lambda, c(0.05, 0.95))

# Papionini tribe
papio_rates <- getCladeRates(edata, node = 572)
mean(papio_rates$lambda)
quantile(papio_rates$lambda, c(0.05, 0.95))

# Cercopithecini tribe
cerci_rates <- getCladeRates(edata, node = 609)
mean(cerci_rates$lambda)
quantile(cerci_rates$lambda, c(0.05, 0.95))

# Colobinae 
colob_rates <- getCladeRates(edata, node = 633)
mean(colob_rates$lambda)
quantile(colob_rates$lambda, c(0.05, 0.95))

# Hominidae
homini_rates <- getCladeRates(edata, node = 549)
mean(homini_rates$lambda)
quantile(homini_rates$lambda, c(0.05, 0.95))

# Hylobatidae
hylob_rates <- getCladeRates(edata, node = 555)
mean(hylob_rates$lambda)
quantile(hylob_rates$lambda, c(0.05, 0.95))

# All primates except papioini
nonpapiorate <- getCladeRates(edata, node = 572, nodetype = "exclude")
mean(nonpapiorate$lambda)
quantile(nonpapiorate$lambda, c(0.05, 0.95))

# Plot of speciation rate through time for all primates
plotRateThroughTime(edata, ratetype="speciation")

plotRateThroughTime(edata, node = 572, 
                    nodetype="include")
plotRateThroughTime(edata, node = 426, 
                    nodetype="include")

# Table of speciation rates and 95%CI
clade_rates <- data.frame(Clade = c("Great_Apes",
                                    "Lesser_Apes",
                                    "Papionins",
                                    "Platyrrhines",
                                    "Colobines",
                                    "Guenons"),
                          Sp_Rate = c(mean(homini_rates$lambda),
                                      mean(hylob_rates$lambda),
                                      mean(papio_rates$lambda),
                                      mean(platy_rates$lambda),
                                      mean(colob_rates$lambda),
                                      mean(cerci_rates$lambda)),
                          L_Sp_Rate = c(quantile(homini_rates$lambda, c(0.05)),
                                        quantile(hylob_rates$lambda, c(0.05)),
                                        quantile(papio_rates$lambda, c(0.05)),
                                        quantile(platy_rates$lambda, c(0.05)),
                                        quantile(colob_rates$lambda, c(0.05)),
                                        quantile(cerci_rates$lambda, c(0.05))),
                          U_Sp_Rate = c(quantile(homini_rates$lambda, c(0.95)),
                                        quantile(hylob_rates$lambda, c(0.95)),
                                        quantile(papio_rates$lambda, c(0.95)),
                                        quantile(platy_rates$lambda, c(0.95)),
                                        quantile(colob_rates$lambda, c(0.95)),
                                        quantile(cerci_rates$lambda, c(0.95))),
                          Ex_Rate = c(mean(homini_rates$mu),
                                      mean(hylob_rates$mu),
                                      mean(papio_rates$mu),
                                      mean(platy_rates$mu),
                                      mean(colob_rates$mu),
                                      mean(cerci_rates$mu)),
                          L_Ex_Rate = c(quantile(homini_rates$mu, c(0.05)),
                                        quantile(hylob_rates$mu, c(0.05)),
                                        quantile(papio_rates$mu, c(0.05)),
                                        quantile(platy_rates$mu, c(0.05)),
                                        quantile(colob_rates$mu, c(0.05)),
                                        quantile(cerci_rates$mu, c(0.05))),
                          U_Ex_Rate = c(quantile(homini_rates$mu, c(0.95)),
                                        quantile(hylob_rates$mu, c(0.95)),
                                        quantile(papio_rates$mu, c(0.95)),
                                        quantile(platy_rates$mu, c(0.95)),
                                        quantile(colob_rates$mu, c(0.95)),
                                        quantile(cerci_rates$mu, c(0.95))))
save(clade_rates, file = "clade_rates.Rdata")

# Cohort matrix
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata)

# Tip speciation rates
tip_rates <- data.frame(Species = c(tree$tip.label),
                        SpRate = c(edata$meanTipLambda),
                        ExRate = c(edata$meanTipMu))
save(tip_rates, file = "tip_rates.Rdata")

full_tip_data <- getTipRates(edata)

ggplot(corr_rates_df,
       aes(x = Clade,
           y = Sp_Rate)) + 
  geom_boxplot(aes(fill = Clade)) +
  geom_errorbar(aes(ymin = L_Sp_Rate,
                    ymax = U_Sp_Rate))

boxplot(full_tip_data$lambda.avg[120:137],
        full_tip_data$lambda.avg[138:158],
        full_tip_data$lambda.avg[159:169],
        full_tip_data$lambda.avg[170:199],
        full_tip_data$lambda.avg[200:206],
        full_tip_data$lambda.avg[207:213],
        full_tip_data$lambda.avg[214:229],
        full_tip_data$lambda.avg[230:267], 
        full_tip_data$lambda.avg[268:292],
        full_tip_data$lambda.avg[293:340],
        #col = c("red", "blue"), 
        names = c("Pitheciidae",
                  "Atelidae",
                  "Cebidae",
                  "Callitrichidae",
                  "Aotidae",
                  "Hominidae",
                  "Hylobatidae",
                  "Papionini",
                  "Cercopithecini",
                  "Colobinae"),
        ylim = c(0.26, 0.32))

# Papionins vs other anthropoids
boxplot(c(full_tip_data$lambda.avg[-c(1:119)],
          full_tip_data$lambda.avg[-c(230:267)]),
        full_tip_data$lambda.avg[230:267],
        col = c("red", "blue"),
        names = c("Other Anthropoids",
                  "Papioini"))





