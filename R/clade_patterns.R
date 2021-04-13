

# Looking at clade specific patterns 
dfsol <- data.frame(Mod3$Sol)
dfsol_melt <- melt(dfsol)

papio_patterns <- c("Macaca", "Mandrillus", 
                    "Papio", "Cercocebus",
                    "Lophocebus")

papdf <- filter(dfsol_melt,
                grepl(paste(papio_patterns,
                            collapse = "|"),
                      variable)) %>%
  separate(variable, c("trait","Type","Species"), 
           sep = "\\.", fill = "right") %>%
  select(-Type) %>%
  filter(trait %in% c("traitVTDwSD", 
                      "traitSSD", 
                      "traitOvulation_Signs")) 

plot(papdf$value[which(papdf$trait == "traitVTDwSD")], 
     papdf$value[which(papdf$trait == "traitOvulation_Signs")])

ggplot() +
  geom_hex(aes(x = papdf$value[which(papdf$trait == "traitVTDwSD")], 
               y = papdf$value[which(papdf$trait == "traitOvulation_Signs")]),
           bins = 100) +
  theme_classic()


papdf2 <- papdf %>% 
  dplyr::group_by(Species) %>% 
  dplyr::mutate(grouped_id = row_number())

papdf2 <- papdf2 %>%
  spread(trait, value) %>%
  select(-grouped_id)

ggplot(data = papdf2,
       aes(x = traitSSD,
           y = traitOvulation_Signs)) +
  geom_hex(bins = 50)



