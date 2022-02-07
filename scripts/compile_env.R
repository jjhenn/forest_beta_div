library(tidyverse)
library(readxl)

#### load in env. data and organize it ####
SERC_env <- read_excel("data/SERC_Env_20x20_20141224.xlsx", sheet = 2) %>% 
  select(-pHWater, -NH4min, -TotNmin, -NO3res, -NO3min, -H, -X, -Y) %>% 
  mutate(site = "SERC") %>% 
  mutate(quadrat = as.character(quadrat))

TRCP_env <- read_excel("data/TRCP_Env_20x20_20141224.xlsx", sheet = 2) %>% 
  mutate_all(.funs = as.numeric) %>% 
  select_if(~all(!is.na(.))) %>% 
  rename("gx" = gx25ha, "gy" = gy25ha)

colnames(TRCP_env) <- gsub("_13", "", colnames(TRCP_env))

TRCP_env <- TRCP_env %>% 
  select(-pHWater, -NH4min, -TotNmin, -NO3, -NO3min) %>% 
  mutate(site = "TRCP") %>% 
  mutate(quadrat = as.character(quadrat))

WFDP_env <- read_excel("data/WFDP_Env_20x20_20141224.xlsx", sheet = 2) %>% 
  select(-pHBaCl2) %>% 
  mutate(site = "WFDP") %>% 
  mutate(quadrat = as.character(quadrat))

all_env <- bind_rows(SERC_env, TRCP_env, WFDP_env) %>% 
  pivot_longer(cols = -c(quadrat, site), names_to = "variable", values_to = "vals")


write.csv(all_env, file = "output/all_env.csv", row.names = F)
