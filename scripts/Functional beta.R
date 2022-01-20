## calcaulte Observed Functional and taxonomic beta among plots at each site
library(vegan)

## load CWM data
WR = read.table("WFDP_CWMV2.txt", header=TRUE)
head(WR)

TY = read.table("TRCP_CWM2s.txt", header=TRUE)
head(TY)

SC= read.table("SERC_CWMs.txt", header=TRUE)
head(SC) 

WR_DIS = vegdist(WR[,2:7], method="euclidean")
mean(WR_DIS)
seWR=sd(WR_DIS)/sqrt(640)

TY_DIS = vegdist(TY[,2:7], method="euclidean")
OBSB=mean(TY_DIS)
seTY=sd(TY_DIS)/sqrt(506)

SC_DIS = vegdist(SC, method="euclidean")
mean(SC_DIS)
seTY=sd(SC_DIS)/sqrt(379)


##Taxo beta 
WR_T = read.table("WFDP_comp.txt", header=TRUE)
head(WR_T)

TY_T = read.table("TRCP_comp.txt", header=TRUE)
head(TY_T)

SC_T = read.table("SERC_comp_2.txt", header=TRUE)
head(SC_T)


WR_DIS_T = vegdist(WR_T[,2:22], method="bray")
mean(WR_DIS_T)
seWR_T=sd(WR_DIS_T)/sqrt(640)

TY_DIS_T = vegdist(TY_T[,2:40], method="bray")
mean(TY_DIS_T)
seTY_T=sd(TY_DIS_T)/sqrt(506)

SC_DIS_T = vegdist(SC_T[,2:40], method="bray")
mean(SC_DIS_T)
seTY_T=sd(SC_DIS_T)/sqrt(379)

## variance partioning functional beta
############################################################################
##   Analysis of Tyson FDP
#########################################################################

##Examine ENV variables
FDPenv<-read.table("TRCP_Env.txt", header=T)
head(FDPenv)
dim(FDPenv)
plot(FDPenv[,2:3])

FDP_ENV = FDPenv[,c(4:21)]
head(FDP_ENV)
dim(FDP_ENV)

pcFDP<-prcomp(FDPenv[,c(4:21)], scale. = TRUE)
summary(pcFDP)
pcFDP$rotation #examine the contribution of each variable to each principal component
pcFDP$x #examine the score of each plot along each principal component; this is the matrix of environmental variables that will be used in the variance partition analysis
pc_scores <- data.frame(FDPenv$quadrat, pcFDP$x)
names(pc_scores)[1] <- "quadrat"
pc_scores <- pc_scores[order(pc_scores$quadrat),] #order plots according to their names, this is important to match plot characteristics and floristic disimilarity among plots in the variance partitioning analysis
write.table(pc_scores, file = "TFDP_PCA.txt", quote = TRUE, sep = ",")
biplot(pcFDP)

#mantel(dist(FDPenv[,2:3]), dist(pcFDP$x))

## plot environmental PCA against geographic distance
## note: great circle distances were calculated with function "distVincentyEllipsoid" of the "geosphere" package (above) in units of meters
## note: environmental distances are measured using all PCA axes

#plot(dist(FDPenv[,2:3]), dist(pcFDP$x), xlab="Geographic distance (meters)", ylab="Environmental dissimilarity (PCA)", cex.lab=1.5, cex.axis=1.5, main="Wind River", las = 1, ylim = c(), cex.main = 1.5)
#legend(x = "topleft", legend = c("P = 0.001", "r = 0.31"), cex = 1.3, bty = "n")
#lm_burned = lm(dist(pcFDP$x)~dist(FDPenv[,2:3]))
#abline(lm_burned, col="BLUE", cex = 5)

#library(ecodist)
#FDP_mgram = mgram(dist(pcFDP$x), dist(FDPenv[,2:3]))
#plot(FDP_mgram, col = "blue")




################################################################
# calculate spatial eigenvalues and eigenvectors #
################################################################

## extract plot coordinates and map them
coor_all <- cbind(FDPenv$gx, FDPenv$gy)
coor_all <- coor_all[order(FDPenv$quadrat),] #order plots by their respective names
plot(coor_all, col="red")

# Check that all plots have unique coordinates
## number of plots
n_coor_all<-length(coor_all[,1])
n_coor_all

## distance matrix
d<-as.matrix(dist(coor_all), method = "euclidean")

## number of pairwise comparisons in distance matrix
length(d[upper.tri(d)])

## check that all plots have unique coor_alldinates.
inspect_coor_all<-which(d==0, arr.ind=T)

## below is the number of elements in the distance matrix having zero distance
nrow(inspect_coor_all)

## below is the number of non-diagonal elements in the distance matrix having zero distance
(nrow(inspect_coor_all)-n_coor_all)/2

## if (and only if) the number above is higher than zero, then look up the columns and rows
## where the non-diagonal elements in the distance matrix having zero distance are. Make sure
## to correct the coordinates of the plots that have equal coor_alldinates (two plots cannot possibly
## have equal coordinates).
inspect_coor_all



##############################################################################################################################
# Construct spatial eigenvalues and eigenvectors; see Griffith & Peres-Neto (2006) Ecology pg. 2605.
# The procedure bellow follows description in pages 612-613 of Dormann et al. (2007, Ecography 30: 609-628).
##############################################################################################################################

## preliminaries
n_coor_all <- length(coor_all[,1])
n_coor_all

## Calculate euclidean distance between each pair of plots
Geodist = dist((FDPenv[,c(2:3)]), method = "euclidean", diag = FALSE)

length(Geodist) #this is the total number of distances, it should be equal to the number below
nrow(FDPenv)*(nrow(FDPenv)-1)/2 #this is the total number of plot pairs, it should be equal to the number above
Geodist.m = matrix(NA, nrow=nrow(FDPenv), ncol=nrow(FDPenv)) #create an empty square matrix that will hold geographic distances between plots
Geodist.m[lower.tri(Geodist.m)] <- Geodist #populate the lower triangle of the matrix with the respecitve geographic distances
Geodist.m[upper.tri(Geodist.m)] <- t(Geodist.m)[upper.tri(Geodist.m)] #populate the upper triangle of the matrix with the respecitve geographic distances
diag(Geodist.m) <- 0 #populate the diagonal of the matrix with zeros

#check the distance matrix, make sure it is square and symetric
Geodist.m[1:10, 1:10]
dim(Geodist.m)
isSymmetric(Geodist.m)
mst <- spantree(Geodist.m) #this claculates the minimum spanning tree


plot(coor_all, col="red")
lines(mst, coor_all)#this plots the minimum spaning tree on the map
ws <- (Geodist.m<=max(mst$dist))*(1-((Geodist.m/(4*(max(mst$dist))))^2))
w <- (Geodist.m>0)*ws
Id <- diag(n_coor_all)#construct an identity matrix with the length of the diagonal equal to the number of plots, i.e., n_coor_burned
ones <- rep(1, n_coor_all)
dim(ones) <- c(n_coor_all, 1)
l<-ones%*%t(ones)
l_n<-l/n_coor_all
res1<-Id-l_n
res2<-res1%*%w
res3<-res2%*%res1
ei_sim_all <- eigen(res3)#these are the spatial eigenvalues and eigenvectores
attributes(ei_sim_all)

## examine the spatial eigenvalues
ei_sim_all$"values"

## examine the first few rows of some of the spatial eigenvectors
ei_sim_all$"vectors"[1:10, 1:10]

## create text files with spatial eigenvalues and respective vectors for each quadrat
write.table(ei_sim_all$"values", file="SpatialEigenValues_TFDP.txt", quote=T, sep=",", row.names = F)
spatial.eigenvectors <- data.frame(FDPenv$quadrat[order(FDPenv$quadrat)], ei_sim_all$"vectors")
names(spatial.eigenvectors)[1] <- "quadrat" 
write.table(spatial.eigenvectors, file="SpatialEigenVectors_TFDP.txt", quote=T, sep=",", row.names = F)


##################################
## analysis based on raw env data
#################################
matrix1 = FDP_ENV
matrix1[1:2,]

## read the table with spatial eigenvalues and spatial eigenvectors, and examine it
# read file with spatial eigenvalues
spatial_eigen_values <- read.table("SpatialEigenValues_TFDP.txt", header=T, sep=",") 
spatial_eigen_values$x
# read file with spatial eigenvectors
spatial_eigen_vectors <- read.table("SpatialEigenVectors_TFDP.txt", header=T, sep=",")
spatial_eigen_vectors[1:10,1:10]
identical(labels(TY_DIS), as.character(spatial_eigen_vectors$Quadrat)) #test the match of tree plot names between the response variable and the environmental variables
spatial_eigen_vectors <- spatial_eigen_vectors[,-1] # exclude plot names if they matched the respective names in the response variable  
summary(spatial_eigen_vectors)
dim(spatial_eigen_vectors)

## read table with plot coordinates and examine it
plot_coor = coor_all
plot_coor[1:5,]
dim(plot_coor)

## build the spatial matrix concatenating the x y points of projected coordinates with the spatial eigen vectors that correspond to positive spatial eigenvalues.
spatial_matrix <- cbind(plot_coor, spatial_eigen_vectors[,spatial_eigen_values$x>0]) 
# examine at least a few rows of the resulting matrix
spatial_matrix[1:5,]
#names(spatial_matrix)[3:6] <- c("SEV1", "SEV2", "SEV3", "SEV4" )
#summary(spatial_matrix)
ncol(spatial_matrix)
nrow(spatial_matrix)
# examine corelations among spatial eigenvectors.
# spatial eigenvectors are orthogonal to each other
spatcor =cor(spatial_matrix)
spatcor[1:15, 1:15]

##select spatial egienvectors 3:255 or xy 1:2
spatial_matrix <- spatial_matrix[,c(3:255)]
dim(spatial_matrix)


## construct the formula for the global CAP model
ev <- cbind(matrix1, spatial_matrix)
head(ev)
a <- c(paste("ev[,", c(1:length(ev[1,])), "]", sep=""))
formula1 <- as.formula(paste("TY_DIS ~ ", paste(a, collapse= "+")))
formula1

# variable names above have to be short because code for model simplification might not parse otherwise; 
#so below is the match between short and long names for explanatory variables ev 11 and ev 12 are spatial variables
variable_names <- cbind(a, c(names(matrix1), names(spatial_matrix)))
colnames(variable_names) <- c("short_name", "long_name")



## run the global CAP model using function "capscale"
CAP_model0 <- capscale(formula1, sqrt.dist = T)
#examine CAP model, see if there are negative eigenvalues that make up the non-Euclidean component of inertia, 
#given under the title "Imaginary". Because the code line above uses argument "add=T", there should be no
#imaginary component to the inertia.
## If still negative eigenvalues use sqrt transforamtion 
##Taxodiversity use Add=T
CAP_model0
plot(CAP_model0)
# obtain adjusted R^2 for global model
RsquareAdj(CAP_model0)


# test significance of the whole model by permutation. 
#This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(CAP_model0, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
RsquareAdj(CAP_model0)


## run the global CAP model "by hand", using first function "capscale" 
#(using a formula with only a constant on the left hand side) for principal coordinates analysis,
# and then function "rda" on the site scores scaled by eigenvalues for redundancy analysis 
#(as described in the help page for function "capscale").
PCO_model <- capscale(TY_DIS ~ 1)
attributes(PCO_model$CA)
PCO_sc= scores(PCO_model, choices = c(1:88), display = c("wa"), scaling =1)
formula1.rda <- as.formula(paste("PCO_sc ~ ", paste(a, collapse= "+")))
formula1.rda
RDA_model0 <- rda(formula1.rda)
attributes(RDA_model0$CCA)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_model0)
plot(RDA_model0)

# test significance of the whole model by permutation. This is the "global test" in Blanchet et al. (2008) 
#Forward selection of explanatory variables. Ecology 89, 2623?2632.
anova.cca(RDA_model0, alpha=0.05, beta=0.01, step=100, perm.max=999)
# obtain adjusted R^2 for global model
RsquareAdj(RDA_model0)




##########################################################################################################################
##load packfor package
library(packfor)

For_Modev = forward.sel(PCO_sc,matrix1, R2thresh=0.921295, adjR2thresh=0.8293594)
For_Modev

For_Modsp = forward.sel(PCO_sc,spatial_matrix, R2thresh=0.921295, adjR2thresh=0.8293594)
For_Modsp

## calculate various fraction of variance explained using distance-based redundancy analysis (dbRDA)
# dbRDA takes distance matrix, then makes an ordination (i.e., principal coordinates analysis), 
# then a constrained ordination in which ordination axes are a linear function of environmental and spatial variables.
# function "capscale" (using a formula with only a constant on the left hand side) 
# takes the distance and creates a multivariate space that preserves those distances (i.e., performsvprincipal coordinates analysis);
# the main message: the result of principal coordinate analysis (the object that below is named "PCO_model") 
# is the response variable used for variance partitioning using function 'varpart'.
PCO_model <- capscale(TY_DIS ~ 1)
PCO_sc= scores(PCO_model, choices = c(1:88), display = c("wa"), scaling =1)

# create two formulas, one for environment ("d"), one for space ("e")
# in c()below enter the columns coresponding to the env factor that are sig in the forward sel model

d<-c(paste("matrix1[,", c(17,2,14,11,18,5,10,3,13,12,8,6,7,9), "]", sep=""))

# for this model include significant sptial columns
e<-c(paste("spatial_matrix[,", c(3,7,4,40,12,14,6,21,13,28,8,11,16,46,70,9,27,2,26,22,33,48,71,45,
                                 44,50,20,25,5,39,144,118,58,18,109,85,115,133,1,129,171,184,32,37,
                                 73,125,97,24,29,41,100,75,51,36,187,65,19,94,167,67,55,72,234,98,
                                 215,253,23,142,208,47,214,34,61,130,15,108,64,31,226,107,89,210,
                                 165,190), "]", sep=""))

formula.environment <- as.formula(paste("PCO_sc ~", paste(d, collapse= "+")))
formula.environment
formula.space <- as.formula(paste("PCO_sc ~", paste(e, collapse= "+")))
formula.space
variance_explained_CAP_model_final <- varpart(PCO_sc, formula.environment, formula.space)
variance_explained_CAP_model_final
# [a] = pure enviro
# [b] = envrio with space
# [c] = pure space
# [d] = unexplained




################################################################################################
##  Null model for functional beta diversity
###############################################################################################

#load required packages
library(vegan)
library(FD)

## load species compositon data
SC_Tn = read.table("TRCP_comp.txt", header=TRUE)
head(SC_Tn)

##reload existing env data indluding env and spatial data
evn = ev

## load in trait data
rsp = read.table("TRCP_Traits2s.txt", header=TRUE)

## set probablilites for sampling from Regional species pool
RegAbund = read.table("TRCP_RegAbund.txt", header=TRUE)
str(RegAbund)
RA = as.vector(RegAbund$Abund/100)

nreps<-999  

##create matrix to store loop results
res3= matrix(data=NA, nrow=1, ncol=nreps)
sim.floristic.dis.bc<-c()
sim.floristic.dis.TOC<-c()

##run loop
for(i in 1:nreps){
  
  ##Create a null community by randomizing the species matrix, while maintaining sample counts and species richness with each plot
  ##Only works if Relative abundance data is whole numbers##
  
  ##shuffles abundances to different species, but keeps row total and actual abundance numbers
  comm = permatfull(SC_Tn, fixedmar = "rows", shuffle = "samp", strata = NULL,
                    mtype = "count", times = 1)
  
  ##Put into usable form for FD package
  mat=data.frame(comm$perm)
  
  ## randomly draw species from RSP and name to match composition dataset
  Ntraits=rsp[sample(nrow(rsp), (dim(mat)[2]), replace=FALSE),]
  row.names(Ntraits) = names(mat)
  #paste("X", 1:nrow(Ntraits), sep="")
  mat2=as.matrix(mat)
  ###Calculate Null CWM trait values for each plot####
  res1 = functcomp(Ntraits, mat2)
  
  ## calcualte functional beta for each itteration and store it  
  #res2=vegdist(res1, method="euclidean")
  #res3[,i]=mean(res2)
  sim.floristic.dis.bc<-append(sim.floristic.dis.bc, vegdist(res1, method="euclidean"), after=length(sim.floristic.dis.bc))
  CWmS_dis=vegdist(res1, method="euclidean")
  name.SC = as.vector(rep("TRCPC", times=504))
  res2 = betadisper(CWmS_dis, name.SC, type = c("centroid"))  
  sim.floristic.dis.TOC<-append(sim.floristic.dis.TOC,  mean(res2$distances), after=length(sim.floristic.dis.TOC))   
  
  
  #res4=sd(res2)
  #ses[,i] = (OBSB-res3) 
}

write.table(as.vector(sim.floristic.dis.bc), "NULL_TRCP2.txt", quote=T, sep=",", row.names=F)
write.table(as.vector(sim.floristic.dis.TOC), "NULL_DTC_TRCP.txt", quote=T, sep=",", row.names=F)


#sim.floristic.dis.bc_MAT= as.matrix(sim.floristic.dis.bc)
#write.table(sim.floristic.dis.bc_MAT, "NULL_SERC_mat.txt", quote=T, sep=",", row.names=F)



