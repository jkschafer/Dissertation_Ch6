library(tidyverse)
library(ape)

mainData <- read.csv("./Data/Grueter_et_al_2015_Rooker&Gavrilets2018_combined_data.csv",
                     stringsAsFactors = FALSE)

names(mainData)[names(mainData) == "Visual.trait.dimorphism.score.without..sexual.dichromatism"] <- "VTDwoSD"
names(mainData)[names(mainData) == "Visual.trait.dimorphism.score.value.with.sexual.dichromatism"] <- "VTDwSD"
mainData$Species <- as.factor(mainData$Species)
mainData$Habitat <- as.factor(mainData$Habitat)
mainData$Social.organization[mainData$Social.organization == "unk"] <- NA
mainData$Social.organization[mainData$Social.organization == ""] <- NA
mainData$Social.organization <- as.factor(mainData$Social.organization)
mainData$Fission.fusion <- as.factor(mainData$Fission.fusion)
mainData$Ovulation_Signs[mainData$Ovulation_Signs == ""] <- NA

# round VTD scores so that they are integers -> can use poisson distribution
mainData$VTDwSD_roundedUP <- ceiling(mainData$VTDwSD)
mainData$VTDwoSD_roundedUP <- ceiling(mainData$VTDwoSD)
mainData$VTDwSD_roundedDN <- floor(mainData$VTDwSD)
mainData$VTDwoSD_roundedDN <- floor(mainData$VTDwoSD)

# Creating binary variable for ovulation signals -> absent, slight, absent_slight = 0; present = 1
mainData$Ovulation_Signs_bin <- ifelse(mainData$Ovulation_Signs == "Present",
                                       "1", "0")
mainData$Ovulation_Signs_bin <- as.factor(mainData$Ovulation_Signs_bin)
# Creating dummy variables for ovulation signals categories
mainData$Ovulation_Signs <- ifelse(mainData$Ovulation_Signs == "Absent", 
                                   "0", 
                                   mainData$Ovulation_Signs)
mainData$Ovulation_Signs <- ifelse(mainData$Ovulation_Signs == "Absent_Slight", 
                                   "1", 
                                   mainData$Ovulation_Signs)
mainData$Ovulation_Signs <- ifelse(mainData$Ovulation_Signs == "Slight", 
                                   "2", 
                                   mainData$Ovulation_Signs)
mainData$Ovulation_Signs <- ifelse(mainData$Ovulation_Signs == "Present", 
                                   "3", 
                                   mainData$Ovulation_Signs)
# Make ovulation signals factor
mainData$Ovulation_Signs <- as.factor(mainData$Ovulation_Signs)

mainData$Species <- gsub(" ", "_", mainData$Species)

names(mainData)[names(mainData) == "Body.mass.dimorphism..g."] <- "SSD"
reduced_data <- mainData[complete.cases(mainData[,c("Habitat", "Social.organization",
                                                    "Group.size", "Ovulation_Signs",
                                                    "VTDwSD", "SSD")]),]

#------------------ Phylogenetic tree -----------------#
tenKtree <- read.nexus("./Data/consensusTree_10kTrees_Primates_Version3.nex")

not_in_tree_reduced <- reduced_data[!reduced_data$Species %in% tenKtree$tip.label,]
unique(not_in_tree_reduced$Species) # many

reduced_data$Species <- gsub("Gorilla_gorilla_", 
                             "Gorilla_gorilla_gorilla", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Callithrix_kuhlii", 
                             "Callithrix_kuhli", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Hylobates_hoolock", 
                             "Bunopithecus_hoolock", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Pan_troglodytes", 
                             "Pan_troglodytes_troglodytes", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Alouatta_palliala", 
                             "Alouatta_palliata", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Cercopithecus/Allochrocebus_lhoseti", 
                             "Cercopithecus_lhoesti", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Cercopithecus_diana_", 
                             "Cercopithecus_diana", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Cebuella_pygmaea", 
                             "Callithrix_pygmaea", 
                             reduced_data$Species)
reduced_data$Species <- gsub("Mico_argentatus", 
                             "Callithrix_argentata", 
                             reduced_data$Species)

# Run it again
not_in_tree_reduced <- reduced_data[!reduced_data$Species %in% tenKtree$tip.label,]
unique(not_in_tree_reduced$Species) # NONE!

tree <- drop.tip(tenKtree, 
                 tenKtree$tip.label[-na.omit(match(reduced_data$Species, 
                                                   tenKtree$tip.label))])
tree <- force.ultrametric(tree)

# Saving as Rdata files for analyses
save(reduced_data, file = "reduced_data.Rdata")
save(tree, file = "10Ktree.Rdata")
