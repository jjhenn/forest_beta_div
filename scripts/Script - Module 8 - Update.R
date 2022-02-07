######   Module  8 Variation Partitioning and review    #####

## today we will review some of the concepts covered in class and conduct variation partitioning
library(vegan)
#library(tidyr)
library(FD)
library(adespatial)
library(tidyverse)

## start with the usual data cleaning

sites = "TRCP"

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




## let's set the row names
row.names(env) = env$code

##let's start with a simple analysis of comparing spatial and environmental drivers in two provinces
## we will start with NC and SN

## lets start with SN and isolate the spatial data
gx1 = env$gx
gy1 = env$gy
SPdata1 = data.frame(gx1,gy1)
row.names(SPdata1) = row.names(env)
head(SPdata1)
dim(SPdata1)

## and calculate our spatial eigenvalues and eigenvectors  

##calculate euclidean distances among plots
xyD1 = dist(SPdata1)

## construct the dbMEM eigenfunctions using the dbMEM function and make it a data.frame for further analyses
SNmem =  as.data.frame(dbmem(xyD1))
summary(SNmem)
##13 positive eigenvectors


################isolate soil variables
SN_soils = env[,6:25]
head(SN_soils)
######################  calculate dissimilarity ##############
## on taxonomic diversity
SN_dis <- vegdist(comp %>% select(-quadrat, -site, -code), method="bray")
SN_dis <- dist(cwm$lipids)

## merge predictor variables into a single dataframe
SNpred = cbind(SNmem, SN_soils)
head(SNpred)
dim(SNpred)


## run the analysis on all variables prior to forward selection
a <- names(SN_soils)
formula1 <- as.formula(paste("SN_dis ~ ", paste(a, collapse= "+")))
formula1

##  run a dbRDA on each data set
SNmodel1 <- dbrda(formula1, SNpred, add = T)

## check the model
SNmodel1
#examine model, see if there are negative eigenvalues that make up the non-Euclidean component of inertia, 
#given under the title "Imaginary". If so uses argument "add=T", there should be no imaginary component to the inertia.
## If still negative eigenvalues use sqrt transformation 

## plot the data and take a look
plot(SNmodel1)
# obtain adjusted R^2 for global model
RsquareAdj(SNmodel1)

###explore linear dependencies among variables by computing variance inflation factors (VIFs)
vif.cca(SNmodel1)
### VIFs above 20 indicate strong colinearity and suggests that variable reduction is justified

## looks like OM, ENR, MG, CA, K, and CEC and H

## so let's do a PCA on our soils
SN_pca = prcomp(SN_soils, scale. = TRUE)
biplot(SN_pca)

## lets looks at our PCA info
summary(SN_pca)
SN_pca

## first three axes explain 62% of the variation in the data
## axis 1 = OM, ENR, Bray P, Olsen P, Mg, NA
## axis 2 = K, CA, CEC,  
## axis 3 = pH, S, H,

## let's pull those out
SN_pc1 = scores(SN_pca, choices = 1)
SN_pc2 = scores(SN_pca, choices = 2)
SN_pc3 = scores(SN_pca, choices = 3)
SN_pc4 = scores(SN_pca, choices = 4)
SN_pc5 = scores(SN_pca, choices = 5)
SN_PCscores = cbind(SN_pc1,SN_pc2,SN_pc3, SN_pc4, SN_pc5)

## re-merge our predictors
SNpred1 = cbind(SNmem, SN_PCscores)
head(SNpred1)
dim(SNpred1)


## run the analysis on all varaibles prior to forward selection
a1 <- names(SNpred1)
formula2 <- as.formula(paste("SN_dis ~ ", paste(a1, collapse= "+")))
formula2

##  run a dbRDA on each data set
SNmodel2 <- dbrda(formula2, SNpred1)

## check the model
SNmodel2
#examine model, see if there are negative eigenvalues that make up the non-Euclidean component of inertia, 
#given under the title "Imaginary". If so uses argument "add=T", there should be no imaginary component to the inertia.
## If still negative eigenvalues use sqrt transformation 

## plot the data and take a look
plot(SNmodel2)
# obtain adjusted R^2 for global model
RsquareAdj(SNmodel2)

###explore linear dependencies among variables by computing variance inflation factors (VIFs)
vif.cca(SNmodel2)

## looks good!

# test significance of the whole model by permutation. 
#This is the "global test" in Blanchet et al. (2008) 
anova.cca(SNmodel2, alpha=0.05, beta=0.01, step=100, perm.max=99)


############## actual variation partitioning      ###################

## If the  model is significant, perform forward model choice 
## This model selection is based solely on adjusted R2 and P-value, following:
## Blanchet et al. (2008) Forward selection of explanatory variables. Ecology 89, 2623?2632.
## Once the final model is found, calculate various fractions of variance explained.

## first run the global model "by hand", using first function "capscale" 
#(using a formula with only a constant on the left hand side) for principal coordinates analysis,
# and then function "rda" on the site scores scaled by eigenvalues for redundancy analysis 
#(as described in the help page for function "capscale").
SNmodel3 <- capscale(SN_dis~ 1)
attributes(SNmodel3$CA)
SNmodel3_sc= scores(SNmodel3, choices = c(1), display = c("wa"), scaling =1)
formula1.rda2T <- as.formula(paste("SNmodel3_sc ~ ", paste(a1, collapse= "+")))
formula1.rda2T
RDA_mode2T <- rda(formula1.rda2T, SNpred1)
attributes(RDA_mode2T$CCA)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode2T)
plot(RDA_mode2T)

# test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(RDA_mode2T, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode2T)

## run forward selection model to get significant variables
For_Modev = forward.sel(SNmodel3_sc, SN_PCscores, R2thresh=0.8855448, adjR2thresh=0.7767001)
For_Modev

For_Modsp = forward.sel(SNmodel3_sc,SNmem, R2thresh=0.8855448, adjR2thresh=0.7767001)
For_Modsp


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
d<-c(paste("SN_PCscores[,", c(3,1,2), "]", sep=""))

# for this model include significant spatial columns
e<-c(paste("SNmem[,", c(17,6,3,14,44,20,22,16,1,7,19,2), "]", sep=""))

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

## lets run it again with a different set of variables
topo = SNenv[, 11:13]
head(topo)

toposoils = cbind(topo, SN_soils)
head (toposoils)

##  let's do a PCA on our soils
SN_pca2 = prcomp(toposoils, scale. = TRUE)
biplot(SN_pca2)

## lets looks at our PCA info
summary(SN_pca2)
SN_pca2

## first five axes explain 71% of the variation in the data
## axis 1 = Mg,CaMgRatio
## axis 2 = K, CA, CEC,  
## axis 3 = S, H
## axis 4 = Slope, OM, ENR, Bray P, pH,
## axis 5 = Aspect, Olsen P, NA

## let's pull those out
SN_pc11 = scores(SN_pca2, choices = 1)
SN_pc12 = scores(SN_pca2, choices = 2)
SN_pc13 = scores(SN_pca2, choices = 3)
SN_pc14 = scores(SN_pca2, choices = 4)
SN_pc15 = scores(SN_pca2, choices = 5)

SN_PCscores2 = cbind(SN_pc11,SN_pc12,SN_pc13, SN_pc14, SN_pc15)


## re-merge our predictors
SNpred2 = cbind(SNmem, SN_PCscores2)
head(SNpred2)
dim(SNpred2)


## run the analysis on all variables prior to forward selection
a2 <- names(SNpred2)
formula3 <- as.formula(paste("SN_dis ~ ", paste(a2, collapse= "+")))
formula3

##  run a dbRDA on each data set
SNmodel4 <- dbrda(formula3, SNpred2)

## check the model
SNmodel4
#examine model, see if there are negative eigenvalues that make up the non-Euclidean component of inertia, 
#given under the title "Imaginary". If so uses argument "add=T", there should be no imaginary component to the inertia.
## If still negative eigenvalues use sqrt transformation 

## plot the data and take a look
plot(SNmodel4)
# obtain adjusted R^2 for global model
RsquareAdj(SNmodel4)

###explore linear dependencies among variables by computing variance inflation factors (VIFs)
vif.cca(SNmodel4)

## looks good!

# test significance of the whole model by permutation. 
#This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(SNmodel4, alpha=0.05, beta=0.01, step=100, perm.max=99)


############## actual variation partitioning      ###################

## If the  model is significant, perform forward model choice 
## This model selection is based solely on adjusted R2 and P-value, following:
## Blanchet et al. (2008) Forward selection of explanatory variables. Ecology 89, 2623?2632.
## Once the final model is found, calculate various fractions of variance explained.

## first run the global model "by hand", using first function "capscale" 
#(using a formula with only a constant on the left hand side) for principal coordinates analysis,
# and then function "rda" on the site scores scaled by eigenvalues for redundancy analysis 
#(as described in the help page for function "capscale").
SNmodel4 <- capscale(SN_dis~ 1)
attributes(SNmodel4$CA)
SNmodel4_sc= scores(SNmodel4, choices = c(1:8), display = c("wa"), scaling =1)
formula1.rda2T1 <- as.formula(paste("SNmodel4_sc ~ ", paste(a2, collapse= "+")))
formula1.rda2T1
RDA_mode2T1 <- rda(formula1.rda2T1, SNpred2)
attributes(RDA_mode2T1$CCA)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode2T1)
plot(RDA_mode2T1)

# test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(RDA_mode2T1, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode2T1)

## run forward selection model to get significant variables
For_Modev = forward.sel(SNmodel4_sc, SN_PCscores2, R2thresh=0.5125844, adjR2thresh=0.4029159)
For_Modev

For_Modsp = forward.sel(SNmodel4_sc,SNmem, R2thresh=0.5125844, adjR2thresh=0.4029159)
For_Modsp


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
d2<-c(paste("SN_PCscores2[,", c(1,2,4), "]", sep=""))

# for this model include significant spatial columns
e2<-c(paste("SNmem[,", c(1,3,4,12,8,6,10,2,13,5), "]", sep=""))

## put it all together
formula.environment2 <- as.formula(paste("SNmodel4_sc ~", paste(d2, collapse= "+")))
formula.environment2
formula.space2 <- as.formula(paste("SNmodel4_sc ~", paste(e2, collapse= "+")))
formula.space2
variance_explained_CAP_model_final2 <- varpart(SNmodel4_sc, formula.environment2, formula.space2) #formula.environment,
variance_explained_CAP_model_final2
# [a] = pure enviro
# [b] = pure spatial
# [c] = E + S
# {d} = Unexplained

##graph it
showvarparts(2, bg = c("red", "blue"))
plot(variance_explained_CAP_model_final2, bg = c("red", "blue"), Xnames=c("env", "space"))

## compare the two
par(mfrow = c(1, 2))
plot(variance_explained_CAP_model_final, bg = c("red", "blue"), Xnames=c("env", "space"))
plot(variance_explained_CAP_model_final2, bg = c("red", "blue"), Xnames=c("env", "space"))


## let's run with a third circle for climate on the same dataset
clim = SNenv[,37:55]
head(clim)

## PCA on climate variables
C_PCA = prcomp(clim, scale.=TRUE)
biplot(C_PCA)
summary(C_PCA)
C_PCA
## first three axes explain 95% of the variation in the data
## axis 1 = 1,5,6,8,9,10,11,14,17,18,
## axis 2 = 4,12,13,15,16,19
## axis 3 = 2,3,7


## let's pull those out
C_pc1 = scores(C_PCA, choices = 1)
C_pc2 = scores(C_PCA, choices = 2)
C_pc3 = scores(C_PCA, choices = 3)

C1_PC = as.data.frame(cbind(C_pc1, C_pc2,C_pc3))
head(C1_PC)
names(C1_PC)[names(C1_PC) == "PC1"] <- "cPC1"
names(C1_PC)[names(C1_PC) == "PC2"] <- "cPC2"
names(C1_PC)[names(C1_PC) == "PC3"] <- "cPC3"
head(C1_PC)

## re-merge our predictors
SNpred3 = cbind(SNmem, SN_PCscores2, C1_PC)
head(SNpred3)
dim(SNpred3)


## run the analysis on all variables prior to forward selection
a3 <- names(SNpred3)
formula4 <- as.formula(paste("SN_dis ~ ", paste(a3, collapse= "+")))
formula4

##  run a dbRDA on each data set
SNmodel5 <- dbrda(formula4, SNpred3)

## check the model
SNmodel5
#examine model, see if there are negative eigenvalues that make up the non-Euclidean component of inertia, 
#given under the title "Imaginary". If so uses argument "add=T", there should be no imaginary component to the inertia.
## If still negative eigenvalues use sqrt transformation 

## plot the data and take a look
plot(SNmodel5)
# obtain adjusted R^2 for global model
RsquareAdj(SNmodel5)

###explore linear dependencies among variables by computing variance inflation factors (VIFs)
vif.cca(SNmodel5)

## looks good!

# test significance of the whole model by permutation. 
#This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(SNmodel5, alpha=0.05, beta=0.01, step=100, perm.max=99)


## first run the global model "by hand", using first function "capscale" 
#(using a formula with only a constant on the left hand side) for principal coordinates analysis,
# and then function "rda" on the site scores scaled by eigenvalues for redundancy analysis 
#(as described in the help page for function "capscale").
model5 <- capscale(SN_dis~ 1, sqrt.dist = T, add = T)
attributes(model5$CA)
model5_sc= scores(model5, choices = c(1:8), display = c("wa"), scaling =1)
formula1.rda5 <- as.formula(paste("model5_sc ~ ", paste(a3, collapse= "+")))
formula1.rda5
RDA_mode5 <- rda(formula1.rda5, SNpred3)
attributes(RDA_mode5$CCA)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode5)
plot(RDA_mode5)

# test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(RDA_mode5, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_mode5)


## run forward selection model to get significant variables
For_Modev3 = forward.sel(model5_sc, SN_PCscores2, R2thresh=0.5798347, adjR2thresh=0.4652441)
For_Modev3

For_Modsp3 = forward.sel(model5_sc,SNmem, R2thresh=0.5798347, adjR2thresh=0.4652441)
For_Modsp3

For_Modclim = forward.sel(model5_sc,C1_PC, R2thresh=0.5798347, adjR2thresh=0.4652441)
For_Modclim


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
d2<-c(paste("SN_PCscores2[,", c(1,2,4), "]", sep=""))

# for this model include significant spatial columns
e2<-c(paste("SNmem[,", c(1,3,4,12,8,6,10,2,13,5), "]", sep=""))

# for this model include significant microbial columns
c = c(paste("C1_PC[,", c(1,2,3), "]", sep=""))


formula.environment2 <- as.formula(paste("model5_sc ~", paste(d2, collapse= "+")))
formula.environment2
formula.space2 <- as.formula(paste("model5_sc ~", paste(e2, collapse= "+")))
formula.space2
formula.clim = as.formula(paste("model5_sc ~", paste(c, collapse= "+")))
formula.clim
variance_explained_CAP_model_final3 <- varpart(model5_sc, formula.environment2, formula.space2, formula.clim) 
variance_explained_CAP_model_final3
# [a] = pure enviro
# [b] = pure spatial
# [c] = pure climate
# [d] = E + S
# [e] = E + C
# [f] = S + C
# [g] = E + S + C
# {h} = Unexplained

showvarparts(3, bg = c("red", "blue", "yellow"))
plot(variance_explained_CAP_model_final3, bg = c("red", "blue", "yellow"), Xnames=c("env", "space", "climate"))


###########################  Exercise   ##################
## Run variation partitioning on another province
