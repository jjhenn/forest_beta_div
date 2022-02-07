library(vegan)
#library(tidyr)
library(FD)
library(adespatial)
library(tidyverse)

#pick site of interest
# sites = "TRCP"
# trait_s = "SLA"
# traits <- read.csv("output/cwm_all.csv") %>% 
#   distinct(trait)

#### loop for each trait ####
res <- data.frame(trait = NULL, site = NULL, env = NULL, both = NULL, space = NULL, residuals = NULL)
vars <- data.frame(NULL)

for(i in c("TRCP", "WFDP", "SERC")){
  for(j in c("BarkThickness", "LeafSize", "LWC", "SeedMassKEW", "SLA", "WoodDensity", "alkaloids", "benzenoids", "heterocyclics", "lipids",  "organic_acids")){
    
    sites = i
    trait_s = j
## first we will load all of our datasets
## 1. species composition data
comp = read.csv("output/all_dat.csv") %>% 
  select(quadrat, site, sp, rel_abund) %>% 
  distinct() %>% 
  pivot_wider(names_from = sp, values_from = rel_abund, values_fill = 0) %>% 
  mutate(code = paste(site, quadrat, sep = "")) %>% 
  filter(site == sites)

cwm <- read.csv("output/cwm_all.csv") %>% 
  select(-X, -type, -cover) %>% 
  pivot_wider(names_from = trait, values_from = cwm) %>% 
  filter(site == sites)

good_plot <- comp$code


## 2. environmental data
env = read.csv("output/all_env.csv") %>% 
  mutate(code = paste(site, quadrat, sep = "")) %>% 
  pivot_wider(names_from = variable, values_from = vals) %>% 
  distinct() %>% 
  filter(code %in% good_plot) %>% 
  filter(site == sites)


SPdata1 = data.frame(env$gx, env$gy)
row.names(SPdata1) = env$code

## and calculate our spatial eigenvalues and eigenvectors  

##calculate euclidean distances among plots
xyD1 = dist(SPdata1)
SNmem =  as.data.frame(dbmem(xyD1))
summary(SNmem)

################isolate soil variables
SN_soils = env[,6:25]


SNpred = cbind(SNmem, SN_soils)

######################  calculate dissimilarity ##############
## on taxonomic diversity
# SN_dis <- vegdist(comp %>% select(-quadrat, -site, -code), method="bray")
SN_dis <- dist(cwm %>% select(trait_s))


formula1 <- as.formula(paste("SN_dis ~ ", paste(names(SNpred), collapse= "+")))

SNmodel1 <- dbrda(formula1, SNpred, add = T)

## check the model
SNmodel1

anova.cca(SNmodel1, alpha=0.05, beta=0.01, step=100, perm.max=99)

#create varpart model
SNmodel3 <- capscale(SN_dis~ 1)
attributes(SNmodel3$CA)
SNmodel3_sc= scores(SNmodel3, choices = c(1), display = c("wa"), scaling =1)
formula1.rda2T <- as.formula(paste("SNmodel3_sc ~ ", paste(names(SNpred), collapse= "+")))
formula1.rda2T
RDA_mode2T <- rda(formula1.rda2T, SNpred)
attributes(RDA_mode2T$CCA)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode2T)
plot(RDA_mode2T)

# test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(RDA_mode2T, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
r <- RsquareAdj(RDA_mode2T)

## run forward selection model to get significant variables
For_Modev = forward.sel(SNmodel3_sc, SN_soils, R2thresh=r$r.squared, adjR2thresh=r$adj.r.squared, alpha = 0.7)
f_ev <- For_Modev

For_Modsp = forward.sel(SNmodel3_sc,SNmem, R2thresh=r$r.squared, adjR2thresh=r$adj.r.squared, alpha = 0.7)
f_sp <- For_Modsp



## calculate various fraction of variance explained using distance-based redundancy analysis (dbRDA)
# dbRDA takes distance matrix, then makes an ordination (i.e., principal coordinates analysis), 
# then a constrained ordination in which ordination axes are a linear function of environmental and spatial variables.
# function "capscale" (using a formula with only a constant on the left hand side) 
# takes the distance and creates a multivariate space that preserves those distances (i.e., performsvprincipal coordinates analysis);
# the main message: the result of principal coordinate analysis (the object that below is named "PCO_model") 
# is the response variable used for variance partitioning using function 'varpart'.
PCO_model <- capscale(SN_dis ~ 1)

# create two formulas, one for environment ("d") and  one for space ("e")
# in c()below enter the columns corresponding to the env factor that are sig in the forward sel model
d<-c(paste("unlist(SN_soils[,", as.character(f_ev$order), "])", sep=""))

# for this model include significant spatial columns
e<-c(paste("unlist(SNmem[,", f_sp$order, "])", sep=""))

## put it all together
formula.environment <- as.formula(paste("SNmodel3_sc ~", paste(d, collapse= "+")))
formula.environment
formula.space <- as.formula(paste("SNmodel3_sc ~", paste(e, collapse= "+")))
formula.space
variance_explained_CAP_model_final <- varpart(SNmodel3_sc, formula.environment, formula.space) 
variance_explained_CAP_model_final
# [a] = pure enviro
# [b] = E+S
# [c] = pur space
# {d} = Unexplained

##graph it
showvarparts(2, bg = c("red", "blue"))
plot(variance_explained_CAP_model_final, bg = c("red", "blue"), Xnames=c("env", "space"))


out <- data.frame(trait = trait_s, site = sites, r2 =  variance_explained_CAP_model_final$part$indfract$Adj.R.squared, part = c("env", "both", "sp", "res"))

vars_out <- bind_rows(f_ev, f_sp) %>% 
  mutate(site = sites,
         trait = trait_s)

res <- bind_rows(res, out)
vars <- bind_rows(vars, vars_out)
}}

res %>% mutate(part = factor(part, levels = c("res",  "sp", "both",  "env"))) %>% 
  mutate(trait = factor(trait, levels = c("BarkThickness", "LeafSize", "LWC", "SeedMassKEW", "SLA", "WoodDensity", "alkaloids", "benzenoids", "heterocyclics", "lipids", "organic_acids"))) %>% 
  ggplot(aes(x = trait, y = r2, fill = part)) +
  geom_bar(position = "fill", stat= "identity")+
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0)) +
  facet_wrap(~site)

write.csv(res, file = "output/varpart_res.csv", row.names = F)
write.csv(vars, file = "output/vars_varpart.csv", row.names = F)


#### varpart calculations for all morph and all meta traits together
res <- data.frame(trait = NULL, site = NULL, env = NULL, both = NULL, space = NULL, residuals = NULL)
vars <- data.frame(NULL)

for(i in c("TRCP", "WFDP", "SERC")){
  for(j in c("func", "meta")){
    
    sites = i
    trait_s = j
    ## first we will load all of our datasets
    ## 1. species composition data
    comp = read.csv("output/all_dat.csv") %>% 
      select(quadrat, site, sp, rel_abund) %>% 
      distinct() %>% 
      pivot_wider(names_from = sp, values_from = rel_abund, values_fill = 0) %>% 
      mutate(code = paste(site, quadrat, sep = "")) %>% 
      filter(site == sites)
    
    cwm <- read.csv("output/cwm_all.csv") %>% 
      select(-X, -cover) %>% 
      filter(site == sites,
             type == trait_s) %>% 
      pivot_wider(names_from = trait, values_from = cwm)
    
    good_plot <- comp$code
    
    
    ## 2. environmental data
    env = read.csv("output/all_env.csv") %>% 
      mutate(code = paste(site, quadrat, sep = "")) %>% 
      pivot_wider(names_from = variable, values_from = vals) %>% 
      distinct() %>% 
      filter(code %in% good_plot) %>% 
      filter(site == sites)
    
    
    SPdata1 = data.frame(env$gx, env$gy)
    row.names(SPdata1) = env$code
    
    ## and calculate our spatial eigenvalues and eigenvectors  
    
    ##calculate euclidean distances among plots
    xyD1 = dist(SPdata1)
    SNmem =  as.data.frame(dbmem(xyD1))
    summary(SNmem)
    
    ################isolate soil variables
    SN_soils = env[,6:25]
    
    
    SNpred = cbind(SNmem, SN_soils)
    
    ######################  calculate dissimilarity ##############
    ## on taxonomic diversity
    # SN_dis <- vegdist(comp %>% select(-quadrat, -site, -code), method="bray")
    SN_dis <- dist(cwm[,-c(1:3)])
    
    
    formula1 <- as.formula(paste("SN_dis ~ ", paste(names(SNpred), collapse= "+")))
    
    SNmodel1 <- dbrda(formula1, SNpred, add = T)
    
    ## check the model
    SNmodel1
    
    anova.cca(SNmodel1, alpha=0.05, beta=0.01, step=100, perm.max=99)
    
    #create varpart model
    SNmodel3 <- capscale(SN_dis~ 1)
    attributes(SNmodel3$CA)
    SNmodel3_sc= scores(SNmodel3, choices = c(1), display = c("wa"), scaling =1)
    formula1.rda2T <- as.formula(paste("SNmodel3_sc ~ ", paste(names(SNpred), collapse= "+")))
    formula1.rda2T
    RDA_mode2T <- rda(formula1.rda2T, SNpred)
    attributes(RDA_mode2T$CCA)
    # obtain adjusted R^2 for global model
    RsquareAdj(RDA_mode2T)
    plot(RDA_mode2T)
    
    # test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
    #Forward selection of explanatory variables. Ecology 89, 2623?2632.
    anova.cca(RDA_mode2T, alpha=0.05, beta=0.01, step=100, perm.max=999)
    # obtain adjusted R^2 for global model
    r <- RsquareAdj(RDA_mode2T)
    
    ## run forward selection model to get significant variables
    For_Modev = forward.sel(SNmodel3_sc, SN_soils, R2thresh=r$r.squared, adjR2thresh=r$adj.r.squared, alpha = 0.7)
    f_ev <- For_Modev
    
    For_Modsp = forward.sel(SNmodel3_sc,SNmem, R2thresh=r$r.squared, adjR2thresh=r$adj.r.squared, alpha = 0.7)
    f_sp <- For_Modsp
    
    
    
    ## calculate various fraction of variance explained using distance-based redundancy analysis (dbRDA)
    # dbRDA takes distance matrix, then makes an ordination (i.e., principal coordinates analysis), 
    # then a constrained ordination in which ordination axes are a linear function of environmental and spatial variables.
    # function "capscale" (using a formula with only a constant on the left hand side) 
    # takes the distance and creates a multivariate space that preserves those distances (i.e., performsvprincipal coordinates analysis);
    # the main message: the result of principal coordinate analysis (the object that below is named "PCO_model") 
    # is the response variable used for variance partitioning using function 'varpart'.
    PCO_model <- capscale(SN_dis ~ 1)
    
    # create two formulas, one for environment ("d") and  one for space ("e")
    # in c()below enter the columns corresponding to the env factor that are sig in the forward sel model
    d<-c(paste("unlist(SN_soils[,", as.character(f_ev$order), "])", sep=""))
    
    # for this model include significant spatial columns
    e<-c(paste("unlist(SNmem[,", f_sp$order, "])", sep=""))
    
    ## put it all together
    formula.environment <- as.formula(paste("SNmodel3_sc ~", paste(d, collapse= "+")))
    formula.environment
    formula.space <- as.formula(paste("SNmodel3_sc ~", paste(e, collapse= "+")))
    formula.space
    variance_explained_CAP_model_final <- varpart(SNmodel3_sc, formula.environment, formula.space) 
    variance_explained_CAP_model_final
    # [a] = pure enviro
    # [b] = E+S
    # [c] = pur space
    # {d} = Unexplained
    
    ##graph it
    showvarparts(2, bg = c("red", "blue"))
    plot(variance_explained_CAP_model_final, bg = c("red", "blue"), Xnames=c("env", "space"))
    
    
    out <- data.frame(trait = trait_s, site = sites, r2 =  variance_explained_CAP_model_final$part$indfract$Adj.R.squared, part = c("env", "both", "sp", "res"))
    
    vars_out <- bind_rows(f_ev, f_sp) %>% 
      mutate(site = sites,
             trait = trait_s)
    
    res <- bind_rows(res, out)
    vars <- bind_rows(vars, vars_out)
  }}

res %>% mutate(part = factor(part, levels = c("res",  "sp", "both",  "env"))) %>% 
  #mutate(trait = factor(trait, levels = c("BarkThickness", "LeafSize", "LWC", "SeedMassKEW", "SLA", "WoodDensity", "alkaloids", "benzenoids", "heterocyclics", "lipids", "organic_acids"))) %>% 
  ggplot(aes(x = trait, y = r2, fill = part)) +
  geom_bar(position = "fill", stat= "identity")+
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0)) +
  facet_wrap(~site)

write.csv(res, file = "output/varpart_res_comb.csv", row.names = F)
write.csv(vars, file = "output/vars_varpart_comb.csv", row.names = F)
