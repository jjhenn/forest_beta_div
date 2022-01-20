library(tidyverse)
library(readxl)

#Calculate CWM for each plot for each trait

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
  mutate(rel_abund = abu/sum(abu))

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


#### combine data ####
all_dat <- all_comp %>% 
  left_join(all_trait) 

write.csv(all_dat, file = "output/all_dat.csv", row.names = F)

cwm <- all_dat %>% 
  group_by(site, quadrat, trait) %>% 
  filter(!is.na(values)) %>% 
  summarize(cwm = sum(values*rel_abund))


cwm %>% 
  ggplot(aes(x = cwm)) +
  geom_histogram() +
  facet_grid(site ~ trait, scales = "free")

write.csv(cwm, file = "output/cwm.csv", row.names = F)
