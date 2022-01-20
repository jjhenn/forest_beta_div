######   Module  8 Variation Partitioning and review    #####

## today we will review some of the concepts covered in class and conduct variation partitioning
library(vegan)
library(tidyr)
library(FD)
library(adespatial)
library(tidyverse)

## start with the usual data cleaning

## first we will load all of our datasets
## 1. species composition data
dat = read.csv("CAstatewide_Local_Woody_Comdata_20190918.csv")

## now let's organize and calculate relative abundances
## make a site by species matrix
dat2 = spread(dat, Species, Abundance, fill = 0)

## fix is to make the row names the site names and remove the offending columns
row.names(dat2) = dat2$Loc.code

## and lets remove the Site info columns
abundances = dat2[,6:144]
head(abundances)

## Because we want to analyze community composition in a relative sense, 

## let’s relativize the data and focus on using Bray-Curtis dissimilarities:
comp <- decostand(abundances, "total")
head(comp)

### Again I am removing some plots here to make life easier, we will come back to this later
comp2=comp
comp2 = comp2[row.names(comp2) != c("ncwmsse"),]
comp2 = comp2[row.names(comp2) != c("kltmnns"),]
comp2 = comp2[row.names(comp2) != c("ncetnse"),]
comp2 = comp2[row.names(comp2) != c("sccisns"),]
comp2 = comp2[row.names(comp2) != c("ncrisns"),]
dim(comp2)

## 2. environmental data
env = read.csv("CAstatewide_Enviro_21090918.csv")
head(env)

## let's set the row names
row.names(env) = env$Loc.code

## We are also going to subset our data so that all the plots without soil data are excluded
env2 = na.omit(env)
dim(env2)
summary(env2)

##and subset our env data to match the comp data
##list of plots in comp data. Since we have our plots as row names, we can just use the row.names function
env.sites = row.names(comp2)
length(env.sites)

#sub-setting columns in ENV data frame to only have plots with comp data
selrow<-(is.element(row.names(env2), as.vector(env.sites)))
env3 = env2[selrow,]
dim(env3)


## and subset our comp data to match our env data
comp.sites = row.names(env3)
length(comp.sites)

comp3 = selrow<-(is.element(row.names(comp2), as.vector(comp.sites)))
comp3 = comp2[selrow,]
dim(comp3)

## now that everything is matched, we are going to further subset to start some "real analyses" and compare 
## drivers of variation among provinces

KLenv = env3[env3$Province == "kl",]
dim(KLenv)
NCenv = env3[env3$Province == "nc",]
dim(NCenv)
SNenv = env3[env3$Province == "sn",]
dim(SNenv)
SCenv = env3[env3$Province == "sc",]
dim(SCenv)


### now we match the composition data

## Kalamath
KL.sites = row.names(KLenv)
length(KL.sites)
KLcomp = selrow<-(is.element(row.names(comp3), as.vector(KL.sites)))
KLcomp =comp3[selrow,]
dim(KLcomp)
## remove species with no plots
x2<-colSums(KLcomp)
zero.cols=(!is.element(names(KLcomp), as.vector(names(x2[x2==0]))))
KLcomp2<-KLcomp[,zero.cols]
head(KLcomp2)
dim(KLcomp2)

## nor Cal
NC.sites = row.names(NCenv)
length(NC.sites)
NCcomp = selrow<-(is.element(row.names(comp3), as.vector(NC.sites)))
NCcomp =comp3[selrow,]
dim(NCcomp)
## remove species with no plots
x2<-colSums(NCcomp)
zero.cols=(!is.element(names(NCcomp), as.vector(names(x2[x2==0]))))
NCcomp2<-NCcomp[,zero.cols]
head(NCcomp2)
dim(NCcomp2)

## Sierra Nevada
SN.sites = row.names(SNenv)
length(SN.sites)
SNcomp = selrow<-(is.element(row.names(comp3), as.vector(SN.sites)))
SNcomp =comp3[selrow,]
dim(SNcomp)
## remove species with no plots
x2<-colSums(SNcomp)
zero.cols=(!is.element(names(SNcomp), as.vector(names(x2[x2==0]))))
SNcomp2<-SNcomp[,zero.cols]
head(SNcomp2)
dim(SNcomp2)

## SoCal
SC.sites = row.names(SCenv)
length(SC.sites)
SCcomp = selrow<-(is.element(row.names(comp3), as.vector(SC.sites)))
SCcomp =comp3[selrow,]
dim(SCcomp)
## remove species with no plots
x2<-colSums(SCcomp)
zero.cols=(!is.element(names(SCcomp), as.vector(names(x2[x2==0]))))
SCcomp2<-SCcomp[,zero.cols]
head(SCcomp2)
dim(SCcomp2)


##let's start with a simple analysis of comparing spatial and environmental drivers in two provinces
## we will start with NC and SN

##First we want to check our sample size
dim(SNcomp2)
## 99 plots in SN
dim(NCcomp2)
## 75 plots in NC

## these give us a fair amount of power to test things and probably don't warrant much data reductions
## but we will do some later

## lets start with SN and isolate the spatial data
gx1 = SNenv$x
gy1 = SNenv$y
SPdata1 = data.frame(gx1,gy1)
row.names(SPdata1) = row.names(SNenv)
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
SN_soils = SNenv[,24:36]
head(SN_soils)
######################  calculate dissimilarity ##############
## on taxonomic diversity
SN_dis <- vegdist(SNcomp2, method="bray")

## merge predictor variables into a single dataframe
SNpred = cbind(SNmem, SN_soils)
head(SNpred)
dim(SNpred)


## run the analysis on all variables prior to forward selection
a <- names(SNpred)
formula1 <- as.formula(paste("SN_dis ~ ", paste(a, collapse= "+")))
formula1

##  run a dbRDA on each data set
SNmodel1 <- dbrda(formula1, SNpred)

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
SN_PCscores = cbind(SN_pc1,SN_pc2,SN_pc3)

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
SNmodel3_sc= scores(SNmodel3, choices = c(1:8), display = c("wa"), scaling =1)
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
For_Modev = forward.sel(SNmodel3_sc, SN_PCscores, R2thresh=0.4798952, adjR2thresh=0.3784114)
For_Modev

For_Modsp = forward.sel(SNmodel3_sc,SNmem, R2thresh=0.4798952, adjR2thresh=0.3784114)
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
d<-c(paste("SN_PCscores[,", c(1,2), "]", sep=""))

# for this model include significant spatial columns
e<-c(paste("SNmem[,", c(1,3,4,12,8,6,10,2,13,5), "]", sep=""))

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
