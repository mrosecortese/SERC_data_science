# Load packages:
library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2)
library(rfishbase)
library(gridExtra)
library(lemon)
library(multcompView)
library(rcompanion)
library(cowplot)
library(FSA)
library(scales)
library(cowplot)
library(effects)
library(emmeans)
library(lme4)
library(lmtest)
library(lubridate)
library(nlme)
library(piecewiseSEM)
library(ggpubr)
library(vegan)
library(leaflet)
library(viridis)
library(scales)

# Import the survey data
RLS <- read_csv("data/CBC_RLS_2015-2019_taxonomic.csv")
CBCcoords <- read.csv("data/RLS metadata.csv")
biomass_coefs <- read_csv("data/CBC_biomass_coefs.csv")


se<-function(x)sqrt(var(x)/length(x))
pd <- position_dodge(width = 0.4)


RLS <- RLS %>% dplyr::rename("sort_order" = "...1") %>%
  select(-c(total_count, valid_name))

depths <- RLS %>% group_by(year, location_name) %>%
  summarize(depth = mean(depth_m))

RLSLong<-melt(RLS, id.vars =c("sort_order","year","habitat","location_name",
                              "transect_decimal_latitude","transect_decimal_longitude","orig_lat",
                              "orig_lon","event","sample_collection_date","method","block",
                              "direction","depth_m","visibility_m","diver","buddy",
                              "kingdom","phylum","class","order","family","genus","rank",
                              "taxon_id","orig_scientific_name","scientific_name"))

RLSLong <- RLSLong %>% rename("size_class" = "variable", "size_count" = "value")

RLSLong <- RLSLong %>% subset(size_count > 0 & method == 1) %>% subset(phylum == "Chordata")

#remove Blueground and Tobacco seagrass; each only surveyed once
RLSLong <- RLSLong %>% subset(location_name != "Blueground Seagrass" &
                                location_name != "Tobacco Seagrass")

###########################
#Refine taxonomy so that lower-resolution IDs aren't counted as separate species

RLSLong$habitat <- recode(RLSLong$habitat, "forereef" = "Forereef",
                             "patch reef" = "Patch Reef",
                             "mangrove" = "Mangrove",
                             "sand" = "Sand",
                             "seagrass" = "Seagrass")

RLSLong$habitat <- factor(RLSLong$habitat, levels = c("Forereef", "Patch Reef", "Mangrove", "Sand", "Seagrass"))



RLSLong <- RLSLong %>%
  mutate(resolution = case_when(
    order == 0 ~ "class",
    family == 0 ~ "order",
    genus == 0 ~ "family",
    rank == "Genus" ~ "genus",
    TRUE ~ as.character(rank)))

RLSLong$resolution <- recode(RLSLong$resolution, "Species" = "species")

res_df <- RLSLong %>% group_by(resolution) %>%
  summarize(n = n())

#####################################
###Number of fish species identified
###Borrowed code from taxonomy_aggregated to sort by taxon rank
RLSLong$habitat <- as.factor(RLSLong$habitat)
habitats <- levels(RLSLong$habitat)

# Only keep observations at the finest level of identification
# This is a conservative sorting..
# if an organism was id at the genus level and other species level ID exists elsewhere with that genus, it is dropped

# Only keep observations at the finest level of identification
# This is a conservative sorting..
# if an organism was id at the genus level and other species level ID exists elsewhere with that genus, it is dropped

df <- data.frame()

id_df <- RLSLong %>%
  filter(year > 2010)

species.df <- id_df %>%
  filter(resolution == "species") %>%
  mutate(identification = scientific_name) %>%
  left_join(id_df)

genus.df <- id_df %>%
  filter(resolution == "genus") %>%
  mutate(identification = genus) %>%
  left_join(id_df)

family.df <- id_df %>%
  filter(resolution == "family") %>%
  mutate(identification = family) %>%
  left_join(id_df)

order.df <- id_df %>%
  filter(resolution == "order") %>%
  mutate(identification = order) %>%
  left_join(id_df)

class.df <- id_df %>%
  filter(resolution == "class") %>%
  mutate(identification = class) %>%
  left_join(id_df)

phylum.df <- id_df %>%
  filter(resolution == "phylum") %>%
  mutate(identification = phylum) %>%
  left_join(id_df)


unique_species <- species.df

unique_species <- genus.df %>%
  filter(!(genus %in% unique_species$genus)) %>%
  bind_rows(unique_species)
unique_species <- family.df %>%
  filter(!(family %in% unique_species$family)) %>%
  bind_rows(unique_species)
unique_species <- order.df %>%
  filter(!(order %in% unique_species$order)) %>%
  bind_rows(unique_species)
unique_species <- class.df %>%
  filter(!(class %in% unique_species$class)) %>%
  bind_rows(unique_species)
unique_species <- phylum.df %>%
  filter(!(phylum %in% unique_species$phylum)) %>%
  bind_rows(unique_species)

df <- df %>%
  bind_rows(unique_species)

fish_clean <- df

####presence/absence#######################
#PLOTS

library(betapart)

##create presence/absence matrix of all species to get total dissimilarity, nestedness, and turnover

RLScast<-dcast(fish_clean, year + habitat + location_name + event~scientific_name,
               value.var = 'size_count',sum)

RLSmat = RLScast[,4:ncol(RLScast)]
RLSmat <- RLSmat %>% remove_rownames %>% column_to_rownames(var="event")
RLSmat[RLSmat > 0] <- 1
RLSmat <- as.matrix(RLSmat)
set.seed(123456)


beta.multi(RLSmat, index.family="sorensen")

##create separate matrices for each year

fish_clean$year <- as.factor(fish_clean$year)
years <- levels(fish_clean$year)

RLSmatList = list()

for(current_year in years) {

  yearmat <- fish_clean %>% subset(year == current_year) %>%
    group_by(year, habitat, location_name, scientific_name) %>%
    summarize(count = sum(size_count)) %>%
    pivot_wider(id_cols = c(year, habitat, location_name), names_from = scientific_name,
                values_from = count) %>%
    replace(is.na(.), 0)

  yearmat = yearmat[,3:ncol(yearmat)]
  yearmat <- yearmat %>% remove_rownames %>% column_to_rownames(var="location_name")
  yearmat[yearmat > 0] <- 1
  yearmat <- as.matrix(yearmat)

  RLSmatList[[length(RLSmatList) + 1]] <- yearmat

}

names(RLSmatList) <- c("RLSmat15", "RLSmat16", "RLSmat17", "RLSmat18", "RLSmat19")


#calculate overall diversity indices within each year

df <- data.frame()

bmult_df <- do.call(rbind, lapply(names(RLSmatList), function(i) {

  mat <- RLSmatList[[i]]

  bmult <- beta.multi(mat, index.family = "sorensen")

  df <- bind_cols(year = as.numeric(paste0("20", gsub("RLSmat([0-9]+)", "\\1", i))),
                          total_diss=bmult$beta.SOR,nest=bmult$beta.SNE,turn=bmult$beta.SIM)

} ) )


#calculate pairwise diversity indices for each site pairing within each year


df <- data.frame()

dissdata <- do.call(rbind, lapply(names(RLSmatList), function(i) {

  mat <- RLSmatList[[i]]

  bp <- beta.pair(mat, index.family = "sorensen")

  df1 <- melt(as.matrix(bp$beta.sor), varnames = c("row", "col")) %>%
    rename("site1" = row, "site2" = col, "total_diss" = value)


  df1$site1 <- as.character(df1$site1)
  df1$site2 <- as.character(df1$site2)
  df1 <- df1 %>%
    rowwise() %>%      # for each row
    mutate(comparison = paste(sort(c(site2, site1)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
    ungroup()


  df2 <- melt(as.matrix(bp$beta.sim), varnames = c("row", "col")) %>%
    rename("site1" = row, "site2" = col, "turnover" = value)
  df2$site1 <- as.character(df2$site1)
  df2$site2 <- as.character(df2$site2)
  df2 <- df2 %>%
    rowwise() %>%      # for each row
    mutate(comparison = paste(sort(c(site2, site1)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
    ungroup()


  df3 <- melt(as.matrix(bp$beta.sne), varnames = c("row", "col")) %>%
    rename("site1" = row, "site2" = col, "nestedness" = value)
  df3$site1 <- as.character(df3$site1)
  df3$site2 <- as.character(df3$site2)
  df3 <- df3 %>%
    rowwise() %>%      # for each row
    mutate(comparison = paste(sort(c(site2, site1)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
    ungroup()


  dissdf <- left_join(df1,df2, by = "comparison")
  dissdf <- select(dissdf, -c(site1.y, site2.y))

  dissdf <- left_join(dissdf, df3, by = "comparison")
  dissdf <- select(dissdf, c(site1, site2, comparison, nestedness, turnover, total_diss)) %>%
    subset(site1 != site2) %>%
    mutate(year = as.numeric(paste0("20", gsub("RLSmat([0-9]+)", "\\1", i))))


  df <- df %>%
    bind_rows(dissdf)

} ) )


dissdata <- select(dissdata, c(year, site1, site2, comparison, nestedness, turnover, total_diss))

dissdata$habitat1 <- dissdata$site1
dissdata$habitat1 <- recode(dissdata$habitat1, "Blueground Mangrove" = "Mangrove",
                            "Twin Cays Mangrove" = "Mangrove",
                            "Tobacco Mangrove" = "Mangrove",
                            "Blueground Seagrass" = "Seagrass",
                            "CBC Seagrass" = "Seagrass",
                            "Twin Cays Seagrass" = "Seagrass",
                            "Curlew Seagrass" = "Seagrass",
                            "Tobacco Seagrass" = "Seagrass",
                            "CBC Sand" = "Sand",
                            "Curlew Sand" = "Sand",
                            "Tobacco Sand" = "Sand",
                            "CBC House Reef" = "Patch Reef",
                            "CBC Lagoon Reef" = "Patch Reef",
                            "Curlew Patch Reef" = "Patch Reef",
                            "CBC Reef Central" = "Forereef",
                            "South Reef Central" = "Forereef",
                            "Tobacco Reef" = "Forereef")

dissdata$habitat2 <- dissdata$site2
dissdata$habitat2 <- recode(dissdata$habitat2, "Blueground Mangrove" = "Mangrove",
                            "Twin Cays Mangrove" = "Mangrove",
                            "Tobacco Mangrove" = "Mangrove",
                            "Blueground Seagrass" = "Seagrass",
                            "CBC Seagrass" = "Seagrass",
                            "Twin Cays Seagrass" = "Seagrass",
                            "Curlew Seagrass" = "Seagrass",
                            "Tobacco Seagrass" = "Seagrass",
                            "CBC Sand" = "Sand",
                            "Curlew Sand" = "Sand",
                            "Tobacco Sand" = "Sand",
                            "CBC House Reef" = "Patch Reef",
                            "CBC Lagoon Reef" = "Patch Reef",
                            "Curlew Patch Reef" = "Patch Reef",
                            "CBC Reef Central" = "Forereef",
                            "South Reef Central" = "Forereef",
                            "Tobacco Reef" = "Forereef")


dissdata <- dissdata %>%
  rowwise() %>%
  mutate(hab_comparison = paste(sort(c(habitat2, habitat1)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
  ungroup()


diss <- dissdata %>% group_by(year, comparison, nestedness, turnover, total_diss) %>%
  summarize(site1 = first(site1), site2 = first(site2),
            habitat1 = first(habitat1), habitat2 = first(habitat2),
            hab_comparison = first(hab_comparison))


disssum2 <- diss %>% group_by(habitat1, habitat2) %>%
  summarize(meand = mean(total_diss), sed = se(total_diss),
            meant = mean(turnover), set = se(turnover),
            meann = mean(nestedness), sen = se(nestedness))

#Want habitats to be in a certain order.  Cannot figure out why this only works if I relevel habitat1 twice,
#but this is my solution for now...
disssum2$habitat1 <- as.factor(disssum2$habitat1)
disssum2$habitat2 <- as.factor(disssum2$habitat2)

disssum2 <- disssum2 %>% mutate(habitat2 = fct_relevel(habitat2,
                                 "Forereef", "Patch Reef", "Mangrove", "Sand", "Seagrass"))
disssum2 <- disssum2 %>% mutate(habitat1 = fct_relevel(habitat1,
                                                       "Forereef", "Patch Reef", "Mangrove", "Sand", "Seagrass"))
disssum2 <- disssum2 %>% mutate(habitat1 = fct_relevel(habitat1,
                                                       "Forereef", "Patch Reef", "Mangrove", "Sand", "Seagrass"))

#get significance letters for turnover
diss$hab_comparison <- gsub("\\ - ", "", diss$hab_comparison)

model <- glm(turnover ~ hab_comparison*year, data = diss)

hist(resid(model))

car::Anova(model)

rsquared(model)

out <- emmeans(model, "hab_comparison")

library(sjPlot)
plot_model(model, type = "pred", terms = "year")

plot(out)


out2 <- as.data.frame(pairs(out))

turnletters <- cldList(p.value ~ contrast,
                       data = out2,
                       threshold = 0.05)

turnletterdf <- turnletters %>%
  separate(Group, c("hab0", "hab1", "hab2", "hab3", "hab4"), sep="(?=[A-Z])") %>%
  mutate(hab1 = case_when(
    hab1 == "Patch" ~ "Patch Reef",
    hab1 != "Patch" ~ hab1)) %>%
  mutate(hab2 = case_when(
    hab2 == "Patch" ~ "Patch Reef",
    hab2 != "Patch" ~ hab2)) %>%
  mutate(hab3 = case_when(
    hab3 == "Patch" ~ "Patch Reef",
    hab3 != "Patch" ~ hab3)) %>%
  mutate(hab2 = case_when(
    hab2 == "Reef" ~ NA_character_,
    hab2 != "Reef" ~ hab2)) %>%
  mutate(hab3 = case_when(
    hab3 == "Reef" ~ NA_character_,
    hab3 != "Reef" ~ hab3)) %>%
  mutate(hab2 = case_when(
    is.na(hab2) ~ hab3,
    !(is.na(hab2)) ~ hab2)) %>%
  rename("habitat2" = "hab1", "habitat1" = "hab2") %>%
  select(habitat1, habitat2, Letter) %>%
  mutate(habitat1 = case_when(
    Letter == "b" ~ "Mangrove",
    Letter != "b" ~ habitat1)) %>%
  mutate(habitat2 = case_when(
    Letter == "b" ~ "Patch Reef",
    Letter != "b" ~ habitat2)) %>%
  arrange(habitat1)



disssum2$habitat1_a <- disssum2$habitat1
disssum2$habitat2_a <- disssum2$habitat2
disssum2 <- disssum2 %>%
  unite_("Hab_Comp", c("habitat1_a", "habitat2_a"), sep = "_")

turnletterdf$habitat1_a <- turnletterdf$habitat1
turnletterdf$habitat2_a <- turnletterdf$habitat2
turnletterdf <- turnletterdf %>%
  unite_("Hab_Comp", c("habitat1_a", "habitat2_a"), sep = "_")

turndf <- left_join(disssum2, turnletterdf, by = "Hab_Comp") %>%
  select(habitat1.x:Hab_Comp, Letter)

colnames(turndf) <- gsub("\\.x","",colnames(turndf))

turndf <- turndf %>%
  mutate(color = "black") %>%
  mutate(color = case_when(
    meant > 0.63 ~ "white",
    meant < 0.63 ~ color))

turnmap <- ggplot(turndf, aes(habitat1, habitat2, fill= meant)) +
  geom_tile(color = "black") +
  scale_fill_viridis("Mean Repl", discrete = FALSE, na.value = "white", direction = -1) +
  scale_y_discrete("") +
  ggtitle("Species Replacement") +
  geom_text(aes(habitat1, habitat2, label = Letter, color = color), size = 2.5) +
  scale_color_identity() +
  labs(tag = "A") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8, angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8, margin = margin(b = 0.2, unit = 'cm')),
        strip.placement = "bottom",
        strip.text.x = element_text(size = 8),
        plot.margin = unit(c(.5,.5,.5,0), "cm"))

plot(turnmap)

#get significance letters for nestedness
model <- glm(nestedness ~ hab_comparison*year, data = diss)

car::Anova(model)

summary(model)$r.squared

rsquared(model)

out <- emmeans(model, "hab_comparison")

plot(out)

plot_model(model, type = "pred", terms = "year")

out2 <- as.data.frame(pairs(out))

nestletters <- cldList(p.value ~ contrast,
                       data = out2,
                       threshold = 0.05)

nestletterdf <- nestletters %>%
  separate(Group, c("hab0", "hab1", "hab2", "hab3", "hab4"), sep="(?=[A-Z])") %>%
  mutate(hab1 = case_when(
    hab1 == "Patch" ~ "Patch Reef",
    hab1 != "Patch" ~ hab1)) %>%
  mutate(hab2 = case_when(
    hab2 == "Patch" ~ "Patch Reef",
    hab2 != "Patch" ~ hab2)) %>%
  mutate(hab3 = case_when(
    hab3 == "Patch" ~ "Patch Reef",
    hab3 != "Patch" ~ hab3)) %>%
  mutate(hab2 = case_when(
    hab2 == "Reef" ~ NA_character_,
    hab2 != "Reef" ~ hab2)) %>%
  mutate(hab3 = case_when(
    hab3 == "Reef" ~ NA_character_,
    hab3 != "Reef" ~ hab3)) %>%
  mutate(hab2 = case_when(
    is.na(hab2) ~ hab3,
    !(is.na(hab2)) ~ hab2)) %>%
  rename("habitat2" = "hab1", "habitat1" = "hab2") %>%
  select(habitat1, habitat2, Letter) %>%
  mutate(marker = "leave") %>%
  mutate(marker = case_when(
    habitat1 == "Patch Reef" & habitat2 == "Mangrove" ~ "change",
    !(habitat1 == "Patch Reef" & habitat2 == "Mangrove") ~ marker)) %>%
  mutate(habitat1 = case_when(
    habitat1 == "Patch Reef" & habitat2 == "Mangrove" ~ "Mangrove",
    !(habitat1 == "Patch Reef" & habitat2 == "Mangrove") ~ habitat1)) %>%
  mutate(habitat2 = case_when(
    habitat1 == "Mangrove" & habitat2 == "Mangrove" & marker == "change" ~ "Patch Reef",
    !(habitat1 == "Patch Reef" & habitat2 == "Mangrove" & marker == "change") ~ habitat2)) %>%
  arrange(habitat1) %>%
  select(habitat1, habitat2, Letter)

nestletterdf$habitat1_a <- nestletterdf$habitat1
nestletterdf$habitat2_a <- nestletterdf$habitat2
nestletterdf <- nestletterdf %>%
  unite_("Hab_Comp", c("habitat1_a", "habitat2_a"), sep = "_")

nestdf <- left_join(disssum2, nestletterdf, by = "Hab_Comp") %>%
  select(habitat1.x:Hab_Comp, Letter)

colnames(nestdf) <- gsub("\\.x","",colnames(nestdf))

nestdf <- nestdf %>%
  mutate(color = "black") %>%
  mutate(color = case_when(
    meann > 0.17 ~ "white",
    meann < 0.17 ~ color))


nestmap <- ggplot(nestdf, aes(habitat1, habitat2, fill= meann)) +
  geom_tile(color = "black") +
  scale_fill_viridis("Mean RichDiff", discrete = FALSE, na.value = "white", direction = -1) +
  scale_y_discrete("") +
  ggtitle("Richness Difference") +
  labs(tag = "B") +
  geom_text(aes(habitat1, habitat2, label = Letter, color = color), size = 2.5) +
  scale_color_identity() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8, angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8, margin = margin(b = 0.2, unit = 'cm')),
        strip.placement = "bottom",
        strip.text.x = element_text(size = 8),
        plot.margin = unit(c(.5,.5,.5,0), "cm"))

plot(nestmap)

tiff("TurnNestMap.tif",width = 8, height = 4, units = "in", res = 800)
grid.arrange(turnmap, nestmap, nrow = 1)
dev.off()

#model total dissimilarity
modeldiss <- glm(total_diss ~ hab_comparison*year, data = diss)

hist(resid(modeldiss))

car::Anova(modeldiss)

rsquared(modeldiss)

out <- emmeans(modeldiss, "hab_comparison")

plot(out)

plot_model(modeldiss, type = "pred", terms = "year")

###species uniqueness tables####
#species uniqueness tables

fish_clean$year <- as.factor(fish_clean$year)

years <- levels(fish_clean$year)

#create artificial observation to prevent empty dataframe for seagrass 2019
fish_clean_narrow <- fish_clean %>% select(year, habitat, scientific_name, resolution)
artificial_obs <- data.frame("2019", "Seagrass", "Artificial Observation", "species")
names(artificial_obs) <- c("year", "habitat", "scientific_name", "resolution")

fish_clean_narrow <- bind_rows(fish_clean_narrow, artificial_obs)

df <- data.frame()

for(current_year in years) {

  splistall <- fish_clean_narrow %>% subset(year == current_year) %>%
    group_by(habitat) %>%
    distinct(habitat, scientific_name, resolution)


  forereefsp <- splistall %>% subset(habitat == "Forereef")
  patchreefsp <- splistall %>% subset(habitat == "Patch Reef")
  mangrovesp <- splistall %>% subset(habitat == "Mangrove")
  seagrasssp <- splistall %>% subset(habitat == "Seagrass")
  sandsp <- splistall %>% subset(habitat == "Sand")


  forereefunique <- forereefsp %>% subset((!(scientific_name %in% patchreefsp$scientific_name)) &
                                            (!(scientific_name %in% mangrovesp$scientific_name)) &
                                            (!(scientific_name %in% seagrasssp$scientific_name)) &
                                            (!(scientific_name %in% sandsp$scientific_name)))

  patchreefunique <- patchreefsp %>% subset((!(scientific_name %in% forereefsp$scientific_name)) &
                                              (!(scientific_name %in% mangrovesp$scientific_name)) &
                                              (!(scientific_name %in% seagrasssp$scientific_name)) &
                                              (!(scientific_name %in% sandsp$scientific_name)))

  mangroveunique <- mangrovesp %>% subset((!(scientific_name %in% patchreefsp$scientific_name)) &
                                            (!(scientific_name %in% forereefsp$scientific_name)) &
                                            (!(scientific_name %in% seagrasssp$scientific_name)) &
                                            (!(scientific_name %in% sandsp$scientific_name)))

  seagrassunique <- seagrasssp %>% subset((!(scientific_name %in% patchreefsp$scientific_name)) &
                                            (!(scientific_name %in% mangrovesp$scientific_name)) &
                                            (!(scientific_name %in% forereefsp$scientific_name)) &
                                            (!(scientific_name %in% sandsp$scientific_name)))

  sandunique <- sandsp %>% subset((!(scientific_name %in% patchreefsp$scientific_name)) &
                                    (!(scientific_name %in% mangrovesp$scientific_name)) &
                                    (!(scientific_name %in% seagrasssp$scientific_name)) &
                                    (!(scientific_name %in% forereefsp$scientific_name)))


  patchun <- patchreefunique %>% summarize(nunique = n())
  patchtot <- patchreefsp %>% ungroup %>% summarize(ntot = n())
  patchreefshare <- patchun %>% mutate(percent = (patchun$nunique/patchtot$ntot) * 100) %>% cbind(patchtot)


  foreun <- forereefunique %>% summarize(nunique = n())
  foretot <- forereefsp %>% ungroup %>% summarize(ntot = n())
  forereefshare <- foreun %>% mutate(percent = (foreun$nunique/foretot$ntot) * 100) %>% cbind(foretot)


  mangroveun <- mangroveunique %>% summarize(nunique = n())
  mangrovetot <- mangrovesp %>% ungroup %>% summarize(ntot = n())
  mangroveshare <- mangroveun %>% mutate(percent = (mangroveun$nunique/mangrovetot$ntot) * 100) %>% cbind(mangrovetot)


  seagrassun <- seagrassunique %>% summarize(nunique = n())
  seagrasstot <- seagrasssp %>% ungroup %>% summarize(ntot = n())
  seagrassshare <- seagrassun %>% mutate(percent = (seagrassun$nunique/seagrasstot$ntot) * 100) %>% cbind(seagrasstot)


  sandun <- sandunique %>% summarize(nunique = n())
  sandtot <- sandsp %>% ungroup %>% summarize(ntot = n())
  sandshare <- sandun %>% mutate(percent = (sandun$nunique/sandtot$ntot) * 100) %>% cbind(sandtot)

  df1 <- bind_rows(forereefshare, patchreefshare, mangroveshare, seagrassshare, sandshare) %>%
    mutate("year" = current_year)

  df <- df %>% bind_rows(df1)

}

#correct for artificial observation
uniquesp_df <- df %>% mutate(nunique = ifelse(year == "2019" & habitat == "Seagrass",
                                               0, nunique)) %>%
  mutate(ntot = ifelse(year == "2019" & habitat == "Seagrass",
                       ntot-1, ntot)) %>%
  mutate(percent = ifelse(year == "2019" & habitat == "Seagrass",
                       0, percent))

uniquesp_df <- uniquesp_df %>% mutate(habitat = fct_relevel(habitat,
                                                       "Forereef", "Patch Reef", "Mangrove", "Sand", "Seagrass"))



unmod <- lm(nunique ~ habitat, data = uniquesp_df)
Anova(unmod)
rsquared(unmod)

pairwise <- as.data.frame(pairs(emmeans(unmod, ~ habitat)))
letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))
letters$Group <- recode(letters$Group, "PatchReef" = "Patch Reef")

unsum <- uniquesp_df %>% group_by(habitat) %>%
  summarize(mean = mean(nunique), se = se(nunique))

unsum <- left_join(unsum, letters, by = c("habitat" = "Group"))
unsum <- unsum %>% mutate(habitat = fct_relevel(habitat,
                                                          "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))


#compute log odds (probability of finding unique species relative to forereef)
exp(coefficients(unmod)[2])

unique <- ggplot(unsum, aes(x = habitat, y = mean, fill = habitat)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = pd, width = 0.2) +
  geom_text(aes(x = habitat, y = mean + se, label = Letter), nudge_y = 0.5, size = 5) +
  geom_text(label = "ANOVA p<0.001", x = 5, y = 15) +
  labs(tag = "C") +
  scale_y_continuous("# Unique Species", expand = c(0,0), limits = c(0, 15.5)) +
  scale_x_discrete("") +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  ggtitle("Number of Unique Species") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 0.5, size = 14),
        axis.title = element_text(size = 15),
        legend.position = "none")

unique

permod <- lm(percent ~ habitat, data = uniquesp_df)
Anova(permod)
rsquared(permod)

pairwise <- as.data.frame(pairs(emmeans(permod, ~ habitat)))
letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))
letters$Group <- recode(letters$Group, "PatchReef" = "Patch Reef")

persum <- uniquesp_df %>% group_by(habitat) %>%
  summarize(mean = mean(percent), se = se(percent))

persum <- left_join(persum, letters, by = c("habitat" = "Group"))

persum <- persum %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

percent <- ggplot(persum, aes(x = habitat, y = mean, fill = habitat)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = pd, width = 0.2) +
  geom_text(aes(x = habitat, y = mean + se, label = Letter), nudge_y = 2, size = 5) +
  geom_text(label = "ANOVA p<0.001", x = 5, y = 51) +
  scale_y_continuous("% Unique", expand = c(0,0), limits = c(0,53)) +
  scale_x_discrete("") +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  ggtitle("Percent Unique Species") +
  labs(tag = "D") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 0.5, size = 14),
        axis.title = element_text(size = 15),
        legend.position = "none")

percent

richmod <- lm(ntot ~ habitat, data = uniquesp_df)
Anova(richmod)
rsquared(richmod)

pairwise <- as.data.frame(pairs(emmeans(richmod, ~ habitat)))
letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))
letters$Group <- recode(letters$Group, "PatchReef" = "Patch Reef")

richsum <- uniquesp_df %>% group_by(habitat) %>%
  summarize(mean = mean(ntot), se = se(ntot))

richsum <- left_join(richsum, letters, by = c("habitat" = "Group"))

richsum <- richsum %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

richness <- ggplot(richsum, aes(x = habitat, y = mean, fill = habitat)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = pd, width = 0.2) +
  geom_text(aes(x = habitat, y = mean + se, label = Letter), nudge_y = 2, size = 5) +
  geom_text(label = "ANOVA p<0.001", x = 5, y = 58) +
  scale_y_continuous("# Species", expand = c(0,0), limits = c(0,60)) +
  scale_x_discrete("") +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  ggtitle("Species Richness") +
  labs(tag = "B") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 0.5, size = 14),
        axis.title = element_text(size = 15),
        legend.position = "none")
richness

abundance <- fish_clean %>%
  group_by(habitat, year, location_name) %>%
  summarize(abundance = sum(size_count)) %>%
  mutate(logab = log10(abundance))

histogram(abundance$logab)

abmod <- lm(logab ~ habitat, data = abundance)
Anova(abmod)
rsquared(abmod)

pairwise <- as.data.frame(pairs(emmeans(abmod, ~ habitat)))
letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

letters$Group <- recode(letters$Group, "PatchReef" = "Patch Reef")

absum <- abundance %>%
  group_by(habitat) %>%
  summarize(mean = mean(logab), se = se(logab))

absum <- left_join(absum, letters, by = c("habitat" = "Group"))

absum <- absum %>% mutate(habitat = fct_relevel(habitat,
                                                    "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

abund <- ggplot(absum, aes(x = habitat, y = mean, fill = habitat)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = pd, width = 0.2) +
  scale_y_continuous("Fish Abundance (log10)", expand = c(0,0), limits = c(0,3.5)) +
  scale_x_discrete("") +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  geom_text(aes(x = habitat, y = mean + se, label = Letter), nudge_y = 0.1, size = 5) +
  geom_text(label = "ANOVA p<0.001", x = 5, y = 3.4) +
  ggtitle("Fish Abundance") +
  labs(tag = "A") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 0.5, size = 14),
        axis.title = element_text(size = 15),
        legend.position = "none")

abund

tiff("FishBarplots.tif",width = 12, height = 10, units = "in", res = 800)
grid.arrange(abund, richness, unique, percent, nrow = 2)
dev.off()

#########
#PRESENCE-ABSENCE
#########
#calculations per Legendre 2014
library(adespatial)

####################
#beta.div.comp <- function(mat, coef="J", quant=FALSE, save.abc=FALSE)


Fishcast<-dcast(fish_clean, habitat + location_name~scientific_name,
               value.var = 'size_count',sum)

Fishmat = Fishcast[,2:ncol(Fishcast)]
Fishmat <- Fishmat %>% remove_rownames %>% column_to_rownames(var="location_name")
Fishmat[Fishmat > 0] <- 1
Fishmat <- as.matrix(Fishmat)
Fishmat <- sqrt(Fishmat)
set.seed(123456)


div1 <- adespatial::beta.div.comp(Fishmat, coef = "S", quant = FALSE, save.abc = FALSE)
rich_mat <- as.matrix(div1$rich)
repl_mat <- as.matrix(div1$repl)
D <- as.matrix(div1$D)


LCBD <- LCBD.comp(D, sqrt.D = TRUE)
rich <- LCBD.comp(rich_mat, sqrt.D = TRUE)
repl <- LCBD.comp(repl_mat, sqrt.D = TRUE)


LCBD_col <- as.data.frame(LCBD$LCBD)
rich_col <- as.data.frame(rich$LCBD)
repl_col <- as.data.frame(repl$LCBD)

sites <- as.data.frame(rownames(D))

df <- cbind(sites,LCBD_col,rich_col,repl_col) %>%
  rename("location_name" = "rownames(D)",
         "LCBD" = "LCBD$LCBD",
         "richdiff" = "rich$LCBD",
         "repl" = "repl$LCBD")

LCBD_indices <- left_join(df, CBCcoords, by = "location_name") %>%
  select(-c("X")) %>%
  rename("lon" = "transect_decimal_longitude",
         "lat" = "transect_decimal_latitude")

LCBD_indices$habitat <- recode(LCBD_indices$habitat, "forereef" = "Forereef",
                               "patch reef" = "Patch Reef",
                               "mangrove" = "Mangrove",
                               "seagrass" = "Seagrass",
                               "sand" = "Sand")


#SCBD plot

#presence absence for SCBD

div3 <- beta.div(Fishmat, sqrt.D = TRUE)
LCBD_df <- as.data.frame(div3$LCBD)
SCBD_df <- as.data.frame(div3$SCBD) %>%
  rename("SCBD" = "div3$SCBD") %>%
  arrange(-(SCBD))

sites <- as.data.frame(rownames(D))


SCBD_df <- SCBD_df %>%
  slice(1:30) %>% rownames_to_column(var="scientific_name")

top20 <- fish_clean %>% subset(scientific_name %in% SCBD_df$scientific_name)

top20$scientific_name <- as.factor(top20$scientific_name)
levels(top20$scientific_name)

top20sum <- top20 %>% group_by(scientific_name, habitat) %>%
  summarize(n = sum(size_count))


spsum <- top20 %>% group_by(scientific_name) %>%
  summarize(nsp = sum(size_count))

top20per <- left_join(top20sum, spsum, by = "scientific_name")

top20prop <- top20per %>% mutate(Proportion = n/nsp) %>%
    arrange(scientific_name, -Proportion) %>%
    slice(1)

top20_fin <- left_join(top20prop, SCBD_df, by = "scientific_name") %>%
  arrange(-(SCBD)) %>%
  mutate(hab = 0.001)


top20_fin <- top20_fin %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

scbd_pa <- ggplot() +
  geom_bar(data = top20_fin,aes(x = reorder(scientific_name, - SCBD), y = SCBD, alpha = Proportion), fill = "black", stat = "identity", position = "dodge", color = "black") +
  geom_bar(data = top20_fin,aes(x=reorder(scientific_name, - SCBD),y=hab,fill = habitat),stat = "identity", position = "dodge", color = "black") +
  scale_y_continuous("SCBD", expand = c(0,0), limits = c(0,0.016)) +
  scale_x_discrete("") +
  scale_fill_viridis("Dominant Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5))

#tiff("SCBDFish.tif",width = 9, height = 6, units = "in", res = 600)
tiff("SCBDFish_PA.tif",width = 9, height = 6, units = "in", res = 800)
scbd_pa
dev.off()

c <- fish_clean %>% subset(habitat == "Mangrove" & scientific_name == "Sphyraena barracuda")

top20pa <- fish_clean %>% subset(scientific_name %in% SCBD_df$scientific_name) %>%
  group_by(year, habitat, location_name, scientific_name) %>%
  summarize(sum = sum(size_count)) %>%
  ungroup %>% complete(year,location_name, scientific_name) %>%
  left_join(CBCcoords, by = "location_name") %>%
  select(year:scientific_name, sum, habitat.y) %>%
  rename(habitat = habitat.y) %>%
  mutate(sum = case_when(
    is.na(sum) ~ 0,
    sum > 0 ~ 1))

top20pa$year <- as.character(top20pa$year)
top20pa$year <- as.numeric(top20pa$year)

topspmod <- glm(sum ~ habitat*year,
               family = "binomial",
               # control = glmerControl(opt = "optim"),
               data = top20pa)


hist(resid(topspmod))
qqnorm(residuals(topspmod, type = "deviance"))
abline(a=0,b=1)

Anova(topspmod, Type = "II")
rsquared(topspmod)

library(emmeans)

out <- emmeans(topspmod, "habitat")

out

plot(out)

out2 <- as.data.frame(pairs(out))

cldList(p.value ~ contrast,
        data = out2,
        threshold = 0.05)



#################################
##random matrix testing
#################################

# Create fake species-by-site matrix
mat <- matrix(rpois(50*30, 1), ncol = 30, dimnames = list(paste0("species", 1:50),paste0("site", 1:30)))

# Transpose
Fishmat_t <- t(Fishmat)

##################################################
#beta.div.comp version

# Create simulation
niter <- 10000 # number of iterations, 10,000 is usually good but takes a while
out <- do.call(rbind, lapply(1:niter, function(i) {

  # Work over each column and randomly assign species
  newmat <- apply(Fishmat_t, 2, function(x) {

    # Determine richness
    nspecies <- sum(x)

    # Randomly sample observed richness from regional species pool
    samp <- sample(rownames(Fishmat_t), nspecies, replace = F)

    # Set all abundances to zero
    x[x > 0] <- 0

    # Repopulate abundances for new "randomized" species set
    x[names(x) %in% samp] <- 1

    return(x)

  } )


  div1 <- adespatial::beta.div.comp(t(newmat), coef = "S", quant = FALSE, save.abc = FALSE)
  repl <- as.matrix(div1$repl)
  repl <- LCBD.comp(repl, sqrt.D = TRUE)


  data.frame(
    habitat = Fishcast$habitat,
    location_name = Fishcast$location_name,
    #iteration = i,
    repl = repl$LCBD
  ) %>%
    group_by(habitat) %>%
    summarize(mean_repl = mean(repl))

} ) )
out # ta-da!

sim_summary <- out %>% group_by(habitat) %>%
  summarize(mean = mean(mean_repl), se = se(mean_repl)) %>%
  mutate(method = "Simulation")

real_summary <- LCBD_indices %>%
  group_by(habitat) %>%
  summarize(mean = mean(repl), se = se(repl)) %>%
  mutate(method = "Observed")

sim_v_real <- rbind(sim_summary, real_summary) %>%
  mutate(habitat = fct_reorder(habitat, desc(mean)))


repl_simvrealp <- ggplot(sim_v_real, aes(x = habitat, y = mean, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = position_dodge(0.9), width = 0.2) +
  scale_y_continuous("LCBD: Species Replacement", expand = c(0,0), limits = c(0,0.12)) +
  scale_x_discrete("") +
  scale_fill_manual("", values = c("gray36", "gray80")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 11))

tiff("replsimFish.tif",width = 8, height = 6, units = "in", res = 600)
repl_simvrealp
dev.off()


##################################3
#LCBD by year


names(RLSmatList) <- c("RLSmat15", "RLSmat16", "RLSmat17", "RLSmat18", "RLSmat19")

lcbd_year <- do.call(rbind, lapply(names(RLSmatList), function(i) {

  mat <- RLSmatList[[i]]

  div1 <- adespatial::beta.div.comp(mat, coef = "S", quant = FALSE, save.abc = FALSE)
  rich_mat <- as.matrix(div1$rich)
  repl_mat <- as.matrix(div1$repl)
  D <- as.matrix(div1$D)


  LCBD <- LCBD.comp(D, sqrt.D = TRUE)
  rich <- LCBD.comp(rich_mat, sqrt.D = TRUE)
  repl <- LCBD.comp(repl_mat, sqrt.D = TRUE)

#
# LCBD_col <- as.data.frame(LCBD$LCBD)
# rich_col <- as.data.frame(rich$LCBD)
# repl_col <- as.data.frame(repl$LCBD)


  out <- cbind.data.frame(location_name=rownames(RLSmatList$RLSmat15),
                          year = as.numeric(paste0("20", gsub("RLSmat([0-9]+)", "\\1", i))),
                          richness = rowSums(mat),
                         LCBD=LCBD$LCBD,rich=rich$LCBD,repl=repl$LCBD)

} ) )


lcbd_year <- left_join(lcbd_year, LCBD_indices[, c("location_name", "habitat")], by = "location_name")


lcbd_year <- lcbd_year %>% mutate(habitat = fct_relevel(habitat,
                                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))
LCBD_year_sum <- lcbd_year %>% group_by(habitat, year) %>%
  summarize(meanrepl = mean(repl), serepl = se(repl))


LCBDyearplot <- ggplot(LCBD_year_sum, aes(x = year, y = meanrepl, colour = habitat)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymax = meanrepl + serepl, ymin = meanrepl - serepl),
                position = pd, width = 0.2) +
  scale_y_continuous("LCBD-replacement", expand = c(0,0)) +
  scale_x_continuous("", breaks = c(2015,2016,2017,2018,2019)) +
  scale_color_viridis("Habitat", begin = 0.1, end = 0.85, option = "magma", discrete = TRUE) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))

tiff("replyearFish.tif",width = 8, height = 6, units = "in", res = 600)
LCBDyearplot
dev.off()

model <- lm(repl ~ habitat*year, data = lcbd_year)

car::Anova(model)

summary(model)$r.squared


out <- emmeans(model, "habitat")

out

plot(out)

out2 <- as.data.frame(pairs(out))

cldList(p.value ~ contrast,
        data = out2,
        threshold = 0.05)


model2 <- lm(richness ~ habitat*year, data = lcbd_year)

car::Anova(model2)

summary(model2)$r.squared


out <- emmeans(model2, "habitat")

out

plot(out)

out2 <- as.data.frame(pairs(out))

cldList(p.value ~ contrast,
        data = out2,
        threshold = 0.001)


########################
#NMDS
########################
#Moved to end to capture top SCBD species
fish_clean <- fish_clean %>% subset(phylum == "Chordata") #gets rid of reef squid...necessary?

##Reduce to top species
top <- fish_clean %>%
  group_by(habitat, scientific_name) %>% summarize(abundance = sum(size_count))

top20_fin$scientific_name <- as.factor(top20_fin$scientific_name)
fishsplist <- levels(top20_fin$scientific_name)


RLScast<-dcast(fish_clean, year + habitat + location_name + event~scientific_name,
               value.var = 'size_count',sum)

RLSmat = RLScast[,4:ncol(RLScast)]
RLSmat <- RLSmat %>% remove_rownames %>% column_to_rownames(var="event")
RLSmat[RLSmat > 0] <- 1
RLSmat <- as.matrix(RLSmat)
RLSmat <- sqrt(RLSmat)
set.seed(123456)


RLSNMDS1 <-
  metaMDS(RLSmat,
          distance = "bray",
          k = 2,
          maxit = 999,
          trymax = 500,
          wascores = TRUE)

goodness(RLSNMDS1)
stressplot(RLSNMDS1)
plot(RLSNMDS1, type = "t")

#set up grouping variables
data.scores1 = as.data.frame(scores(RLSNMDS1)$sites)

data.scores1$year = RLScast$year
data.scores1$habitat = RLScast$habitat
data.scores1$location_name = RLScast$location_name

data.scores1$year <- as.factor(data.scores1$year)

data.scores1 <- data.scores1 %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

species.scores1 <- as.data.frame(scores(RLSNMDS1, "species"))
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores
species.scores1 <- species.scores1 %>% subset(species %in% fishsplist)



forereef1 <- data.scores1[data.scores1$habitat == "Forereef", ][chull(data.scores1[data.scores1$habitat ==
                                                                                     "Forereef", c("NMDS1", "NMDS2")]), ]
patchreef1 <- data.scores1[data.scores1$habitat == "Patch Reef", ][chull(data.scores1[data.scores1$habitat ==
                                                                                        "Patch Reef", c("NMDS1", "NMDS2")]), ]
mangrove1 <- data.scores1[data.scores1$habitat == "Mangrove", ][chull(data.scores1[data.scores1$habitat ==
                                                                                     "Mangrove", c("NMDS1", "NMDS2")]), ]
sand1 <- data.scores1[data.scores1$habitat == "Sand", ][chull(data.scores1[data.scores1$habitat ==
                                                                             "Sand", c("NMDS1", "NMDS2")]), ]
seagrass1 <- data.scores1[data.scores1$habitat == "Seagrass", ][chull(data.scores1[data.scores1$habitat ==
                                                                                     "Seagrass", c("NMDS1", "NMDS2")]), ]
hull.data1 <- rbind(forereef1, patchreef1, mangrove1, seagrass1, sand1) %>%
  mutate(habitat = fct_relevel(habitat,
 "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))


hull.data1

library(ggrepel)

tiff("FishNMDS_SCBDsp2.tif",width = 10, height = 10, units = "in", res = 800)
ggplot() +
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=habitat,group=habitat),alpha=0.30, color = "black") +
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,shape=as.factor(year),colour = habitat, fill = habitat),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 4, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 4, min.segment.length = 0.5, max.overlaps = 18) +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Year",values=c("2015" = 21, "2016" = 22, "2017" = 23, "2018" = 24,
                                     "2019" = 25)) +
  xlim(-1,2) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


RLScast$year <- as.numeric(RLScast$year)
#PERMANOVA using adonis function
adonis2(formula = RLSmat ~ habitat + year + year*habitat, data = RLScast, permutations = 10000)


############
#Abundance-weighted SCBD & LCBD
############

Fishcast<-dcast(fish_clean, habitat + location_name~scientific_name,
                value.var = 'size_count',sum)

FishmatAb = Fishcast[,2:ncol(Fishcast)]
FishmatAb <- FishmatAb %>% remove_rownames %>% column_to_rownames(var="location_name")
FishmatAb <- as.matrix(FishmatAb)
FishmatAb <- sqrt(FishmatAb)
set.seed(123456)


div1 <- adespatial::beta.div.comp(FishmatAb, coef = "S", quant = TRUE, save.abc = FALSE)
rich_mat <- as.matrix(div1$rich)
repl_mat <- as.matrix(div1$repl)
D <- as.matrix(div1$D)


LCBD <- LCBD.comp(D, sqrt.D = TRUE)
rich <- LCBD.comp(rich_mat, sqrt.D = TRUE)
repl <- LCBD.comp(repl_mat, sqrt.D = TRUE)


LCBD_col <- as.data.frame(LCBD$LCBD)
rich_col <- as.data.frame(rich$LCBD)
repl_col <- as.data.frame(repl$LCBD)

sites <- as.data.frame(rownames(D))

df <- cbind(sites,LCBD_col,rich_col,repl_col) %>%
  rename("location_name" = "rownames(D)",
         "LCBD" = "LCBD$LCBD",
         "richdiff" = "rich$LCBD",
         "repl" = "repl$LCBD")

LCBD_indices_Ab <- left_join(df, CBCcoords, by = "location_name") %>%
  select(-c("X")) %>%
  rename("lon" = "transect_decimal_longitude",
         "lat" = "transect_decimal_latitude")

LCBD_indices_Ab$habitat <- recode(LCBD_indices_Ab$habitat, "forereef" = "Forereef",
                                  "patch reef" = "Patch Reef",
                                  "mangrove" = "Mangrove",
                                  "seagrass" = "Seagrass",
                                  "sand" = "Sand")
div3 <- beta.div(FishmatAb, sqrt.D = TRUE)
LCBD_df <- as.data.frame(div3$LCBD)
SCBD_df <- as.data.frame(div3$SCBD) %>%
  rename("SCBD" = "div3$SCBD") %>%
  arrange(-(SCBD))

sites <- as.data.frame(rownames(D))


SCBD_df <- SCBD_df %>%
  slice(1:30) %>% rownames_to_column(var="scientific_name")

top30Ab <- fish_clean %>% subset(scientific_name %in% SCBD_df$scientific_name)

top30Ab$scientific_name <- as.factor(top30Ab$scientific_name)
fishsplistAb <- levels(top30Ab$scientific_name)

top30Absum <- top30Ab %>% group_by(scientific_name, habitat) %>%
  summarize(n = sum(size_count))


spsum <- top30Ab %>% group_by(scientific_name) %>%
  summarize(nsp = sum(size_count))

top30Abper <- left_join(top30Absum, spsum, by = "scientific_name")

top30Abprop <- top30Abper %>% mutate(Proportion = n/nsp) %>%
  arrange(scientific_name, -Proportion) %>%
  slice(1)

top30Ab_fin <- left_join(top30Abprop, SCBD_df, by = "scientific_name") %>%
  arrange(-(SCBD)) %>%
  mutate(hab = 0.003)


top30Ab_fin <- top30Ab_fin %>%
  mutate(habitat = fct_relevel(habitat,
                               "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

scbd_Ab <- ggplot() +
  geom_bar(data = top30Ab_fin,aes(x = reorder(scientific_name, - SCBD), y = SCBD, alpha = Proportion), fill = "black", stat = "identity", position = "dodge", color = "black") +
  geom_bar(data = top30Ab_fin,aes(x=reorder(scientific_name, - SCBD),y=hab,fill = habitat),stat = "identity", position = "dodge", color = "black") +
  scale_y_continuous("SCBD", expand = c(0,0)) +
  scale_x_discrete("") +
  scale_fill_viridis("Dominant Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5))

#tiff("SCBDFish.tif",width = 9, height = 6, units = "in", res = 600)
tiff("SCBDFishAbundance.tif",width = 9, height = 6, units = "in", res = 600)
scbd_Ab
dev.off()



real_summaryAb <- LCBD_indices_Ab %>%
  group_by(habitat) %>%
  summarize(mean = mean(repl), se = se(repl)) %>%
  mutate(method = "real")



LCBD_Ab <- ggplot(real_summaryAb, aes(x = reorder(habitat, - mean), y = mean)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = position_dodge(0.9), width = 0.2) +
  scale_y_continuous("LCBD: Species Replacement", expand = c(0,0), limits = c(0,0.12)) +
  scale_x_discrete("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 11))

tiff("LCBD_Ab.tif",width = 8, height = 6, units = "in", res = 600)
LCBD_Ab
dev.off()

##Abundance NMDS
Abcast <-dcast(fish_clean, habitat + year + location_name + event~scientific_name,
                value.var = 'size_count',sum)

FishmatAb = Abcast[,4:ncol(Abcast)]
FishmatAb <- FishmatAb %>% remove_rownames %>% column_to_rownames(var="event")
FishmatAb <- as.matrix(FishmatAb)
FishmatAb <- sqrt(FishmatAb)
set.seed(123456)


RLSNMDSAb <-
  metaMDS(FishmatAb,
          distance = "bray",
          k = 2,
          maxit = 999,
          trymax = 500,
          wascores = TRUE)

goodness(RLSNMDSAb)
stressplot(RLSNMDSAb)
plot(RLSNMDSAb, type = "t")

#set up grouping variables
data.scoresAb = as.data.frame(scores(RLSNMDSAb)$sites)

data.scoresAb$year = Abcast$year
data.scoresAb$habitat = Abcast$habitat
data.scoresAb$location_name = Abcast$location_name

data.scoresAb$year <- as.factor(data.scoresAb$year)

data.scoresAb <- data.scoresAb %>% mutate(habitat = fct_relevel(habitat,
                                                              "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

species.scoresAb <- as.data.frame(scores(RLSNMDSAb, "species"))
species.scoresAb$species <- rownames(species.scoresAb)  # create a column of species, from the rownames of species.scores
species.scoresAb <- species.scoresAb %>% subset(species %in% fishsplistAb)



forereefAb <- data.scoresAb[data.scoresAb$habitat == "Forereef", ][chull(data.scoresAb[data.scoresAb$habitat ==
                                                                                     "Forereef", c("NMDS1", "NMDS2")]), ]
patchreefAb <- data.scoresAb[data.scoresAb$habitat == "Patch Reef", ][chull(data.scoresAb[data.scoresAb$habitat ==
                                                                                        "Patch Reef", c("NMDS1", "NMDS2")]), ]
mangroveAb <- data.scoresAb[data.scoresAb$habitat == "Mangrove", ][chull(data.scoresAb[data.scoresAb$habitat ==
                                                                                     "Mangrove", c("NMDS1", "NMDS2")]), ]
sandAb <- data.scoresAb[data.scoresAb$habitat == "Sand", ][chull(data.scoresAb[data.scoresAb$habitat ==
                                                                             "Sand", c("NMDS1", "NMDS2")]), ]
seagrassAb <- data.scoresAb[data.scoresAb$habitat == "Seagrass", ][chull(data.scoresAb[data.scoresAb$habitat ==
                                                                                     "Seagrass", c("NMDS1", "NMDS2")]), ]
hull.dataAb <- rbind(forereefAb, patchreefAb, mangroveAb, seagrassAb, sandAb) %>%
  mutate(habitat = fct_relevel(habitat,
                               "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))


hull.dataAb

library(ggrepel)

tiff("FishNMDS_Ab.tif",width = 10, height = 10, units = "in", res = 800)
ggplot() +
  geom_polygon(data=hull.dataAb,aes(x=NMDS1,y=NMDS2,fill=habitat,group=habitat),alpha=0.30, color = "black") +
  geom_point(data=data.scoresAb,aes(x=NMDS1,y=NMDS2,shape=as.factor(year),colour = habitat, fill = habitat),size=4, alpha = 0.8) +
  geom_point(data=species.scoresAb,aes(x=NMDS1,y=NMDS2), pch = 4, color = "black", size=2) +
  geom_text_repel(data=species.scoresAb,aes(x=NMDS1,y=NMDS2,label=species),size = 4, min.segment.length = 0.5, max.overlaps = 18) +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Year",values=c("2015" = 21, "2016" = 22, "2017" = 23, "2018" = 24,
                                     "2019" = 25)) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


Abcast$year <- as.numeric(Abcast$year)
#PERMANOVA using adonis function
adonis2(formula = FishmatAb ~ habitat + year + year*habitat, data = Abcast, permutations = 10000)


################################
#Supplemental Table 2: Unique species across all years
################################
traits <- read.csv("Traits_all-species_LHedit.csv")

traits$CURRENT_TAXONOMIC_NAME <- recode(traits$CURRENT_TAXONOMIC_NAME,
                                               "Strongylura notata" = "Strongylura notata notata",
                                               "Dasyatis centroura" = "Bathytoshia centroura")


splistall2 <- fish_clean %>% group_by(habitat) %>%
  distinct(habitat, scientific_name, resolution)


forereefsp2 <- splistall2 %>% subset(habitat == "Forereef")
patchreefsp2 <- splistall2 %>% subset(habitat == "Patch Reef")
mangrovesp2 <- splistall2 %>% subset(habitat == "Mangrove")
seagrasssp2 <- splistall2 %>% subset(habitat == "Seagrass")
sandsp2 <- splistall2 %>% subset(habitat == "Sand")




forereefunique2 <- forereefsp2 %>% subset((!(scientific_name %in% patchreefsp2$scientific_name)) &
                                          (!(scientific_name %in% mangrovesp2$scientific_name)) &
                                          (!(scientific_name %in% seagrasssp2$scientific_name)) &
                                          (!(scientific_name %in% sandsp2$scientific_name)))

patchreefunique2 <- patchreefsp2 %>% subset((!(scientific_name %in% forereefsp2$scientific_name)) &
                                            (!(scientific_name %in% mangrovesp2$scientific_name)) &
                                            (!(scientific_name %in% seagrasssp2$scientific_name)) &
                                            (!(scientific_name %in% sandsp2$scientific_name)))

mangroveunique2 <- mangrovesp2 %>% subset((!(scientific_name %in% patchreefsp2$scientific_name)) &
                                          (!(scientific_name %in% forereefsp2$scientific_name)) &
                                          (!(scientific_name %in% seagrasssp2$scientific_name)) &
                                          (!(scientific_name %in% sandsp2$scientific_name)))

seagrassunique2 <- seagrasssp2 %>% subset((!(scientific_name %in% patchreefsp2$scientific_name)) &
                                          (!(scientific_name %in% mangrovesp2$scientific_name)) &
                                          (!(scientific_name %in% forereefsp2$scientific_name)) &
                                          (!(scientific_name %in% sandsp2$scientific_name)))

sandunique2 <- sandsp2 %>% subset((!(scientific_name %in% patchreefsp2$scientific_name)) &
                                  (!(scientific_name %in% mangrovesp2$scientific_name)) &
                                  (!(scientific_name %in% seagrasssp2$scientific_name)) &
                                  (!(scientific_name %in% forereefsp2$scientific_name)))


#mean abundances of unique species
forereefmeanab <- fish_clean %>% subset(scientific_name %in% forereefunique2$scientific_name) %>%
  group_by(habitat, year, location_name) %>% summarise(abundance = sum(size_count)) %>%
  group_by(habitat) %>% summarize(mean = mean(abundance))

patchreefmeanab <- fish_clean %>% subset(scientific_name %in% patchreefunique2$scientific_name) %>%
  group_by(habitat, year, location_name) %>% summarise(abundance = sum(size_count)) %>%
  group_by(habitat) %>% summarize(mean = mean(abundance))

mangrovemeanab <- fish_clean %>% subset(scientific_name %in% mangroveunique2$scientific_name) %>%
  group_by(habitat, year, location_name) %>% summarise(abundance = sum(size_count)) %>%
  group_by(habitat) %>% summarize(mean = mean(abundance))

seagrassmeanab <- fish_clean %>% subset(scientific_name %in% seagrassunique2$scientific_name) %>%
  group_by(habitat, year, location_name) %>% summarise(abundance = sum(size_count)) %>%
  group_by(habitat) %>% summarize(mean = mean(abundance))

sandmeanab <- fish_clean %>% subset(scientific_name %in% sandunique2$scientific_name) %>%
  group_by(habitat, year, location_name) %>% summarise(abundance = sum(size_count)) %>%
  group_by(habitat) %>% summarize(mean = mean(abundance))

meanabunique <- rbind(forereefmeanab, patchreefmeanab, mangrovemeanab, sandmeanab, seagrassmeanab) %>%
  select(mean)

abundancesum <- fish_clean %>%
  group_by(habitat, year, location_name) %>%
  summarize(abundance = sum(size_count)) %>%
  group_by(habitat) %>%
  summarize(meantotal = mean(abundance)) %>%
  mutate(logtotal = log10(meantotal))

table <- cbind(meanabunique, abundancesum) %>%
  select(habitat, mean, meantotal) %>%
  mutate(proportion = mean/meantotal)

#total abundances of each unique species
forereeflist <- fish_clean %>% subset(scientific_name %in% forereefunique2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(abundance = sum(size_count))

patchreeflist <- fish_clean %>% subset(scientific_name %in% patchreefunique2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(abundance = sum(size_count))

mangrovelist <- fish_clean %>% subset(scientific_name %in% mangroveunique2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(abundance = sum(size_count))

seagrasslist <- fish_clean %>% subset(scientific_name %in% seagrassunique2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(abundance = sum(size_count))

sandlist <- fish_clean %>% subset(scientific_name %in% sandunique2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(abundance = sum(size_count))


uniquespp <- rbind(forereeflist,patchreeflist,mangrovelist,sandlist,seagrasslist)


write.csv(uniquespp, "unique_by_hab.csv")

totalab <- fish_clean %>% group_by(habitat) %>%
  summarize(abundance = sum(size_count))

fish_traits <- left_join(uniquespp, traits, by = c("scientific_name" = "CURRENT_TAXONOMIC_NAME"))

ft <- fish_traits %>% select(habitat, scientific_name, MaxLength:Habitat) %>% distinct()

summary <- ft %>% group_by(habitat, Trophic.group) %>%
  summarize(n = n())

SCBDfish_traits <- left_join(top20_fin, traits, by = c("scientific_name" = "CURRENT_TAXONOMIC_NAME"))

ft2 <- SCBDfish_traits %>% select(habitat, scientific_name, MaxLength:Habitat) %>% distinct() %>%
  arrange(habitat)

summary2 <- ft2 %>% group_by(habitat, Trophic.group) %>%
  summarize(n = n())


sandlisttotal <- fish_clean %>%
  subset(scientific_name %in% sandsp2$scientific_name) %>%
  group_by(habitat, scientific_name) %>% summarise(n = sum(size_count)) %>%
  subset(habitat == "Sand")


##############
#Map Sites
##############
CBCcoords <- CBCcoords %>% subset(location_name != "Tobacco Seagrass" &
                                    location_name != "Blueground Seagrass")

CBCcoords$habitat <- recode(CBCcoords$habitat, "forereef" = "Forereef",
                            "patch reef" = "Patch Reef",
                            "mangrove" = "Mangrove",
                            "seagrass" = "Seagrass",
                            "sand" = "Sand")

CBCcoords <- CBCcoords %>% mutate(habitat = fct_relevel(habitat,
                                                         "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

library(rgeos)
library(ggmap)

habColor <- colorFactor(c("purple4", "deeppink", "magenta3", "khaki", "salmon1"),
                        domain = c("Forereef", "Patch Reef", "Mangrove",
                                   "Seagrass", "Sand"))
# A map of all CBC sites
CBCsites <- leaflet(CBCcoords) %>%
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  addCircleMarkers(lng = ~transect_decimal_longitude, lat = ~transect_decimal_latitude, weight = 1,
                   radius = 6, popup = ~location_name, fillColor = ~habColor(habitat),
                   stroke = TRUE, color = "black", fillOpacity = 1.0) %>%
  addLegend("bottomright", habColor, values = ~habitat,
            title = "Habitats",
            opacity = 1) %>%
  addScaleBar("topright", options = scaleBarOptions(imperial = FALSE))

CBCsites

inset <- leaflet(CBCcoords) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  addCircleMarkers(lng = ~transect_decimal_longitude, lat = ~transect_decimal_latitude, weight = 1,
                   radius = 6, popup = ~location_name, fillColor = ~habColor(habitat),
                   stroke = TRUE, color = "black", fillOpacity = 1.0) %>%
  addLegend("bottomright", habColor, values = ~habitat,
            title = "Habitats",
            opacity = 1)

inset




###################################################################################
# FISH BIOMASS                                                                    #
###################################################################################
biomass_coefs <- read_csv("CBC_biomass_coefs.csv")

#Figure out what constants are in the RLS teams database--should use theirs
##whenever possible efs$CURRENT_TAXONOMIC_NAME)]
nasp <- biomass_coefs %>% subset(is.na(a))

biomass_coefs$CURRENT_TAXONOMIC_NAME <- recode(biomass_coefs$CURRENT_TAXONOMIC_NAME,
                                               "Strongylura notata" = "Strongylura notata notata",
                                               "Dasyatis centroura" = "Bathytoshia centroura")


biomass_coefs <- select(biomass_coefs, "CURRENT_TAXONOMIC_NAME":"SGFGu")

spp_biomass <- left_join(fish_clean, biomass_coefs, by = c("scientific_name" = "CURRENT_TAXONOMIC_NAME"))

check <- spp_biomass %>% subset(is.na(a) & phylum == "Chordata")#should be empty

### CALCULATE BIOMASS:
##Convert size classes to diver-bias corrected values
##provided by Rick Stuart-Smith
fishspp_biomass <- spp_biomass %>% subset(phylum == "Chordata") %>% droplevels()

fishspp_biomass$orig_size_class <- fishspp_biomass$size_class
fishspp_biomass$size_class <- recode(fishspp_biomass$size_class, "invert_count" = "0",
                                     "2.5" = "3",
                                     "5" = "6",
                                     "7.5" = "9",
                                     "10" = "12",
                                     "12.5" = "15",
                                     "15" = "18",
                                     "20" = "23.1",
                                     "25" = "27",
                                     "30" = "30.5",
                                     "35" = "33.8",
                                     "40" = "36.8",
                                     "50" = "57.5",
                                     "62.5" = "69",
                                     "75" = "80.5",
                                     "87.5" = "92",
                                     "100" = "104",
                                     "112.5" = "115",
                                     "125" = "127",
                                     "137.5" = "138",
                                     "150" = "150",
                                     "162.5" = "161",
                                     "175" = "173",
                                     "187.5" = "184",
                                     "200" = "230",
                                     "250" = "276",
                                     "300" = "322",
                                     "350" = "368")

###calculate biomass using equation provided by Rick Stuart-Smith:
###EXP(LN(a)+b*LN(size class 1*CF))*abundance in size class 1
fishspp_biomass$size_class <- as.numeric(levels(fishspp_biomass$size_class))[fishspp_biomass$size_class]

fishspp_biomass2 <- fishspp_biomass %>% mutate(mass =
                                                 exp(log(a)+b*log(size_class*cf))*size_count)


fish_clean <- fishspp_biomass2

###########
#Biomass-weighted SCBD and LCBD
###########

FishcastMass <-dcast(fish_clean, habitat + location_name~scientific_name,
                     value.var = 'mass',sum)

FishmatMass = FishcastMass[,2:ncol(FishcastMass)]
FishmatMass <- FishmatMass %>% remove_rownames %>% column_to_rownames(var="location_name")
FishmatMass <- as.matrix(FishmatMass)
FishmatMass <- sqrt(FishmatMass)
set.seed(123456)


div1 <- adespatial::beta.div.comp(FishmatMass, coef = "S", quant = TRUE, save.abc = FALSE)
rich_mat <- as.matrix(div1$rich)
repl_mat <- as.matrix(div1$repl)
D <- as.matrix(div1$D)


LCBD <- LCBD.comp(D, sqrt.D = TRUE)
rich <- LCBD.comp(rich_mat, sqrt.D = TRUE)
repl <- LCBD.comp(repl_mat, sqrt.D = TRUE)


LCBD_col <- as.data.frame(LCBD$LCBD)
rich_col <- as.data.frame(rich$LCBD)
repl_col <- as.data.frame(repl$LCBD)

sites <- as.data.frame(rownames(D))

df <- cbind(sites,LCBD_col,rich_col,repl_col) %>%
  rename("location_name" = "rownames(D)",
         "LCBD" = "LCBD$LCBD",
         "richdiff" = "rich$LCBD",
         "repl" = "repl$LCBD")

LCBD_indices_Mass <- left_join(df, CBCcoords, by = "location_name") %>%
  select(-c("X")) %>%
  rename("lon" = "transect_decimal_longitude",
         "lat" = "transect_decimal_latitude")

LCBD_indices_Mass$habitat <- recode(LCBD_indices_Mass$habitat, "forereef" = "Forereef",
                                    "patch reef" = "Patch Reef",
                                    "mangrove" = "Mangrove",
                                    "seagrass" = "Seagrass",
                                    "sand" = "Sand")
div3 <- beta.div(FishmatMass, sqrt.D = TRUE)
LCBD_df <- as.data.frame(div3$LCBD)
SCBD_df <- as.data.frame(div3$SCBD) %>%
  rename("SCBD" = "div3$SCBD") %>%
  arrange(-(SCBD))

sites <- as.data.frame(rownames(D))


SCBD_df <- SCBD_df %>%
  slice(1:30) %>% rownames_to_column(var="scientific_name")

top30Mass <- fish_clean %>% subset(scientific_name %in% SCBD_df$scientific_name)

top30Mass$scientific_name <- as.factor(top30Mass$scientific_name)
splistMass <-levels(top30Mass$scientific_name)

top30Masssum <- top30Mass %>% group_by(scientific_name, habitat) %>%
  summarize(n = sum(size_count))


spsum <- top30Mass %>% group_by(scientific_name) %>%
  summarize(nsp = sum(size_count))

top30Massper <- left_join(top30Masssum, spsum, by = "scientific_name")

top30Massprop <- top30Massper %>% mutate(Proportion = n/nsp) %>%
  arrange(scientific_name, -Proportion) %>%
  slice(1)

top30Mass_fin <- left_join(top30Massprop, SCBD_df, by = "scientific_name") %>%
  arrange(-(SCBD)) %>%
  mutate(hMass = 0.003)

top30Mass_fin <- top30Mass_fin %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))
scbd_Mass <- ggplot() +
  geom_bar(data = top30Mass_fin,aes(x = reorder(scientific_name, - SCBD), y = SCBD, alpha = Proportion), fill = "black", stat = "identity", position = "dodge", color = "black") +
  geom_bar(data = top30Mass_fin,aes(x=reorder(scientific_name, - SCBD),y=hMass,fill = habitat),stat = "identity", position = "dodge", color = "black") +
  scale_y_continuous("SCBD", expand = c(0,0)) +
  scale_x_discrete("") +
  scale_fill_viridis("Dominant Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1,1,1,1, "cm"))

tiff("SCBDFishMass.tif",width = 9, height = 6, units = "in", res = 800)
scbd_Mass
dev.off()


real_summaryMass <- LCBD_indices_Mass %>%
  group_by(habitat) %>%
  summarize(mean = mean(repl), se = se(repl)) %>%
  mutate(method = "real")

real_summaryMass <- real_summaryMass %>% mutate(habitat = fct_relevel(habitat,
                                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

LCBD_Mass <- ggplot(real_summaryMass, aes(x = reorder(habitat, - mean), y = mean)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = position_dodge(0.9), width = 0.2) +
  scale_y_continuous("LCBD: Species Replacement", expand = c(0,0), limits = c(0,0.13)) +
  scale_x_discrete("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5))

tiff("replMassFish.tif",width = 6, height = 4, units = "in", res = 600)
LCBD_Mass
dev.off()

##Biomass NMDS
MassCast<-dcast(fish_clean, habitat + year + location_name + event~scientific_name,
                value.var = 'mass',sum)

FishmatMass = MassCast[,4:ncol(MassCast)]
FishmatMass <- FishmatMass %>% remove_rownames %>% column_to_rownames(var="event")
FishmatMass <- as.matrix(FishmatMass)
FishmatMass <- sqrt(FishmatMass)
set.seed(123456)


RLSNMDSMass <-
  metaMDS(FishmatMass,
          distance = "bray",
          k = 2,
          maxit = 999,
          trymax = 500,
          wascores = TRUE)

goodness(RLSNMDSMass)
stressplot(RLSNMDSMass)
plot(RLSNMDSMass, type = "t")

#set up grouping variables
data.scoresMass = as.data.frame(scores(RLSNMDSMass)$sites)

data.scoresMass$year = MassCast$year
data.scoresMass$habitat = MassCast$habitat
data.scoresMass$location_name = MassCast$location_name

data.scoresMass$year <- as.factor(data.scoresMass$year)

data.scoresMass <- data.scoresMass %>% mutate(habitat = fct_relevel(habitat,
                                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

species.scoresMass <- as.data.frame(scores(RLSNMDSMass, "species"))
species.scoresMass$species <- rownames(species.scoresMass)  # create a column of species, from the rownames of species.scores
species.scoresMass <- species.scoresMass %>% subset(species %in% splistMass)



forereefMass <- data.scoresMass[data.scoresMass$habitat == "Forereef", ][chull(data.scoresMass[data.scoresMass$habitat ==
                                                                                         "Forereef", c("NMDS1", "NMDS2")]), ]
patchreefMass <- data.scoresMass[data.scoresMass$habitat == "Patch Reef", ][chull(data.scoresMass[data.scoresMass$habitat ==
                                                                                            "Patch Reef", c("NMDS1", "NMDS2")]), ]
mangroveMass <- data.scoresMass[data.scoresMass$habitat == "Mangrove", ][chull(data.scoresMass[data.scoresMass$habitat ==
                                                                                         "Mangrove", c("NMDS1", "NMDS2")]), ]
sandMass <- data.scoresMass[data.scoresMass$habitat == "Sand", ][chull(data.scoresMass[data.scoresMass$habitat ==
                                                                                 "Sand", c("NMDS1", "NMDS2")]), ]
seagrassMass <- data.scoresMass[data.scoresMass$habitat == "Seagrass", ][chull(data.scoresMass[data.scoresMass$habitat ==
                                                                                         "Seagrass", c("NMDS1", "NMDS2")]), ]
hull.dataMass <- rbind(forereefMass, patchreefMass, mangroveMass, seagrassMass, sandMass) %>%
  mutate(habitat = fct_relevel(habitat,
                               "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))


hull.dataMass

tiff("FishNMDS_Mass.tif",width = 10, height = 10, units = "in", res = 800)
ggplot() +
  geom_polygon(data=hull.dataMass,aes(x=NMDS1,y=NMDS2,fill=habitat,group=habitat),alpha=0.30, color = "black") +
  geom_point(data=data.scoresMass,aes(x=NMDS1,y=NMDS2,shape=as.factor(year),colour = habitat, fill = habitat),size=4, alpha = 0.8) +
  geom_point(data=species.scoresMass,aes(x=NMDS1,y=NMDS2), pch = 4, color = "black", size=2) +
  geom_text_repel(data=species.scoresMass,aes(x=NMDS1,y=NMDS2,label=species),size = 4, min.segment.length = 0.5, max.overlaps = 18) +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Year",values=c("2015" = 21, "2016" = 22, "2017" = 23, "2018" = 24,
                                     "2019" = 25)) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


MassCast$year <- as.numeric(MassCast$year)
#PERMANOVA using adonis function
adonis2(formula = FishmatMass ~ habitat + year + year*habitat, data = MassCast, permutations = 10000)

############################
#Species Accumulation Curves
############################

##create separate matrices for each habitat

fish_clean$habitat <- as.factor(fish_clean$habitat)
habitats <- levels(fish_clean$habitat)
habitats

RLSmatListHab = list()

for(current_habitat in habitats) {

  habmat <- fish_clean %>% subset(habitat == current_habitat) %>%
    group_by(year, habitat, location_name, event, scientific_name) %>%
    summarize(count = sum(size_count)) %>%
    pivot_wider(id_cols = c(year, habitat, location_name, event), names_from = scientific_name,
                values_from = count) %>%
    replace(is.na(.), 0)

  habmat = habmat[,4:ncol(habmat)]
  habmat <- habmat %>% remove_rownames %>% column_to_rownames(var="event")
  habmat[habmat > 0] <- 1
  habmat <- as.matrix(habmat)

  RLSmatListHab[[length(RLSmatListHab) + 1]] <- habmat

}

names(RLSmatListHab) <- c("Forereef", "PatchReef", "Mangrove", "Sand", "Seagrass")


accurveFore<-specaccum(RLSmatListHab$Forereef, method="random", permutations=100)
accurvePatch<-specaccum(RLSmatListHab$PatchReef, method="random", permutations=100)
accurveMang<-specaccum(RLSmatListHab$Mangrove, method="random", permutations=100)
accurveSea<-specaccum(RLSmatListHab$Seagrass, method="random", permutations=100)
accurveSand<-specaccum(RLSmatListHab$Sand, method="random", permutations=100)

plot(accurvePatch, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Sampling Events",
     ylab = "Species Richness")

habcolors <- c("Forereef" = "purple4", "Patch Reef" = "magenta3", "Mangrove" = "deeppink",
            "Seagrass" = "salmon1", "Sand" = "khaki")

AccumPlot <- ggplot() +
  geom_line(aes(x = accurvePatch$sites, y = accurvePatch$richness), color = "magenta3", size = 1) +
  geom_point(aes(x = accurvePatch$sites, y = accurvePatch$richness, fill = "Patch Reef"), pch = 21, color = "black", size = 3) +
  geom_line(aes(x = accurveFore$sites, y = accurveFore$richness), color = "purple4", size = 1) +
  geom_point(aes(x = accurveFore$sites, y = accurveFore$richness, fill = "Forereef"), pch = 21, color = "black", size = 3) +
  geom_line(aes(x = accurveMang$sites, y = accurveMang$richness), color = "deeppink", size = 1) +
  geom_point(aes(x = accurveMang$sites, y = accurveMang$richness, fill = "Mangrove"), pch = 21, color = "black", size = 3) +
  geom_line(aes(x = accurveSea$sites, y = accurveSea$richness), color = "salmon1", size = 1) +
  geom_point(aes(x = accurveSea$sites, y = accurveSea$richness, fill = "Seagrass"), pch = 21, color = "black", size = 3) +
  geom_line(aes(x = accurveSand$sites, y = accurveSand$richness), color = "khaki", size = 1) +
  geom_point(aes(x = accurveSand$sites, y = accurveSand$richness, fill = "Sand"), pch = 21, color = "black", size = 3) +
  scale_x_continuous(breaks = seq(0,15, by = 5)) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x = "Number of Sampling Events",
       y = "Species Richness",
       fill = "Habitat") +
  scale_fill_manual(values = c("Forereef" = "purple4", "Patch Reef" = "magenta3", "Mangrove" = "deeppink",
                                "Seagrass" = "salmon1", "Sand" = "khaki")) +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 12),
        legend.key = element_rect(fill = "white"))

AccumPlot


tiff("AccumCurves.tif",width = 7, height = 5, units = "in", res = 600)
AccumPlot
dev.off()

accurveALL<-specaccum(RLSmat, method="random", permutations=100)
AllAccum <- ggplot() +
  geom_line(aes(x = accurveALL$sites, y = accurveALL$richness), color = "gray50", size = 1) +
  geom_point(aes(x = accurveALL$sites, y = accurveALL$richness), pch = 21, fill = "gray50", color = "black", size = 3) +
  scale_x_continuous(breaks = seq(0,75, by = 15)) +
  labs(x = "Number of Sampling Events",
       y = "Species Richness",
       fill = "Habitat") +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 12),
        legend.key = element_rect(fill = "white"))
AllAccum

tiff("AccumCurvesALL.tif",width = 7, height = 5, units = "in", res = 600)
AllAccum
dev.off()

#biomass bar plot - not used in manuscript
biomass <- fish_clean %>%
  group_by(habitat, year, location_name) %>%
  summarize(mass = sum(mass)) %>%
  mutate(logmass = log10(mass))

histogram(biomass$logmass)

massmod <- lm(logmass ~ habitat, data = biomass)
Anova(massmod)
rsquared(massmod)

pairwise <- as.data.frame(pairs(emmeans(massmod, ~ habitat)))
letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

letters$Group <- recode(letters$Group, "PatchReef" = "Patch Reef")

masssum <- biomass %>%
  group_by(habitat) %>%
  summarize(mean = mean(logmass), se = se(logmass))

masssum <- left_join(masssum, letters, by = c("habitat" = "Group"))

masssum <- masssum %>% mutate(habitat = fct_relevel(habitat,
                                                "Forereef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

massp <- ggplot(masssum, aes(x = habitat, y = mean, fill = habitat)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se),
                position = pd, width = 0.2) +
  scale_y_continuous("Fish Biomass (log10)", expand = c(0,0), limits = c(0,4.5)) +
  scale_x_discrete("") +
  scale_fill_viridis("Habitat", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  geom_text(aes(x = habitat, y = mean + se, label = Letter), nudge_y = 0.2, size = 5) +
  geom_text(label = "ANOVA p<0.001", x = 5, y = 4) +
  ggtitle("Fish Biomass") +
  labs(tag = "A") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 0.5, size = 14),
        axis.title = element_text(size = 15),
        legend.position = "none")

massp
