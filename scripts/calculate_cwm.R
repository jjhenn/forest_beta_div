library(tidyverse)
library(readxl)
library(FD)

#Calculate CWM for each plot for each trait (both functional and metabolomic)
#Also could calculate FD for each quadrat. I have not done this yet.

#### load in community data ####
SERC_comp <- read.table("data/SERC_comp.txt") %>% 
  mutate(quadrat = row.names(.),
         site = "SERC") %>% 
  pivot_longer(cols = c(-quadrat, -site), names_to = "sp", values_to = "abu") %>% 
  filter(abu != 0)

TRCP_comp <- read.table("data/TRCP_comp.txt") %>% 
  mutate(quadrat = row.names(.),
         site = "TRCP") %>% 
  pivot_longer(cols = c(-quadrat, -site), names_to = "sp", values_to = "abu") %>% 
  filter(abu != 0)

WFDP_comp <- read.table("data/WFDP_comp.txt") %>% 
  mutate(quadrat = row.names(.),
         site = "WFDP") %>% 
  pivot_longer(cols = c(-quadrat, -site), names_to = "sp", values_to = "abu") %>% 
  filter(abu != 0)

all_comp <- bind_rows(WFDP_comp, TRCP_comp, SERC_comp) %>% 
  group_by(quadrat, site) %>% 
  mutate(rel_abund = abu/sum(abu)) %>% 
  separate(sp, into = c("genus", "species"), sep = "_", remove = F) %>% 
  mutate(code = paste(toupper(substr(genus, 1, 4)), toupper(substr( species, 1, 2)), sep = "")) %>% 
  select(-species) %>% 
  mutate(code = ifelse(code == "ACERSA", "ACERS2", code),
         code = ifelse(code == "CARYOV", "CARYO2", code),
         code = ifelse(code == "DIOSVI", "DIO2VI", code),
         code = ifelse(code == "QUERPA", "QUERP1", code))



#### load in trait data ####
SERC_trait <- read_excel("data/SERC_Traits_SpeciesMeans_20150120.xlsx", sheet = 2) %>% 
  select(SpCode, "sp" = Species, Family, GrowthForm, LeafSize, SLA, LWC, ChlorophyllContent, BarkThickness, WoodDensity, SeedMassKEW) %>% 
  pivot_longer(cols = c(-SpCode, -sp, -Family, -GrowthForm), names_to = "trait", values_to = "values") %>% 
  mutate(values = as.numeric(values),
         site = "SERC")

TRCP_trait <- read_excel("data/TRCP_Traits_SpeciesMeans_20160501.xlsx", sheet = 2, col_types = "text") %>% 
  select("SpCode" = spcode, "sp" = Latin, "Family" = family, GrowthForm, "LeafSize" = LeafArea, SLA, "ChlorophyllContent" = LeafChlorophyll, LWC, SeedMassKEW, BarkThickness, WoodDensity) %>% 
  pivot_longer(cols = c(-SpCode, -sp, -Family, -GrowthForm), names_to = "trait", values_to = "values") %>% 
  mutate(values = as.numeric(values),
         site = "TRCP")

WFDP_trait <- read_excel("data/WFDP_Traits_SpeciesMeans_20150120.xlsx", sheet = 2, col_types = "text") %>% 
  select(SpCode, "sp" = Species, Family, GrowthForm, LeafSize, SLA, LWC, ChlorophyllContent, BarkThickness, WoodDensity, SeedMassKEW) %>% 
  pivot_longer(cols = c(-SpCode, -sp, -Family, -GrowthForm), names_to = "trait", values_to = "values") %>% 
  mutate(values = as.numeric(values),
         site = "WFDP")

all_trait = bind_rows(SERC_trait, TRCP_trait, WFDP_trait) %>% 
  mutate(sp = gsub(" ", "_", sp),
         values = ifelse(trait %in% c("LeafSize", "SLA", "SeedMassKEW"), log(values), values))


### combine data ####
all_dat <- all_comp %>%
  left_join(all_trait)

# write.csv(all_dat, file = "output/all_dat.csv", row.names = F)

cwm_func <- all_dat %>% 
  group_by(site, quadrat, trait) %>% 
  filter(!is.na(values)) %>% 
  summarize(cwm = sum(values*rel_abund),
            cover = sum(rel_abund)) %>% 
  filter(trait != "ChlorophyllContent") %>% 
  mutate(type = "func")


cwm_func %>% 
  ggplot(aes(x = cwm)) +
  geom_histogram() +
  facet_grid(site ~ trait, scales = "free")

write.csv(cwm_func, file = "output/cwm_func.csv", row.names = F)


#### load in metabolomic trait data ####
meta_traits <- read.csv("data/alltraits.csv") %>% 
  select("code" = Row.names, alkaloids, heterocyclics, lipids, hydrocarbons, organic_acids, benzenoids, lignans) %>% 
  pivot_longer(cols = -code, names_to = "trait", values_to = "values") %>% 
  left_join(all_comp %>% ungroup() %>% select(code, sp)) %>% 
  distinct()

all_dat_meta <- all_comp %>% 
  left_join(meta_traits) %>% 
  filter(values > 0)

cwm_meta <- all_dat_meta %>% 
  group_by(site, quadrat, trait) %>% 
  filter(!is.na(values)) %>% 
  summarize(cwm = sum(values*rel_abund),
            cover = sum(rel_abund)) %>% 
  filter(!trait %in% c("lignans", "hydrocarbons")) %>% 
  mutate(type = "meta")

cwm_all <- bind_rows(cwm_func, cwm_meta)

write.csv(cwm_all, file = "output/cwm_all.csv")

#### calculate FD for each quadrat with all traits ####
#create wide species and trait data frames
fd_comp <- all_comp %>% 
  select(-genus, -rel_abund, -code) %>% 
  pivot_wider(names_from = sp, values_from = abu, values_fill = 0)

out1 <- data.frame(NULL)

for(i in unique(fd_comp$site)){

  comp <- fd_comp %>% 
    filter(site == i) %>% 
    ungroup() %>% 
    select(-site) %>% 
    select(where(is.numeric)) %>% 
    select(where(~sum(.) != 0))
  
  trait <- all_trait %>% 
    filter(site == i) %>% 
    select(sp, trait, values) %>% 
    filter(sp %in% colnames(comp)) %>% 
    pivot_wider(names_from = trait, values_from = values)
  
  rownames(trait) <- trait$sp
  
  comp <- comp[, unique(trait$sp)]
  
  test <- dbFD(trait, comp, corr= "none")
  
  out <- fd_comp %>% filter(site == i) %>% ungroup() %>%  select(quadrat) %>% bind_cols(., as.data.frame(test$FDis)) %>% mutate(site = i) %>% mutate(trait = "fdis_func")
  colnames(out) <- c("quadrat", "values", "site", "trait")
  
  out1 <- bind_rows(out1, out)
}


out2 <- data.frame(NULL)

for(i in unique(fd_comp$site)){
  
  comp <- fd_comp %>% 
    filter(site == i) %>% 
    ungroup() %>% 
    select(-site) %>% 
    select(where(is.numeric)) %>% 
    select(where(~sum(.) != 0))
  
  trait <- meta_traits %>% 
    filter(trait != "lignans", trait != "hydrocarbons") %>% 
    select(sp, trait, values) %>% 
    filter(sp %in% colnames(comp)) %>% 
    pivot_wider(names_from = trait, values_from = values)
  
  rownames(trait) <- trait$sp
  
  comp <- comp[, unique(trait$sp)]
  
  test <- dbFD(trait, comp, corr = "none")
  
  out <- fd_comp %>% filter(site == i) %>% ungroup() %>%  select(quadrat) %>% bind_cols(., as.data.frame(test$FDis)) %>% mutate(site = i) %>% mutate(trait = "fdis_meta")
  colnames(out) <- c("quadrat", "values", "site", "trait")
  
  out2 <- bind_rows(out2, out)
}

fdis_dat <- bind_rows(out1, out2)

write.csv(fdis_dat, file = "output/fdis.csv")
