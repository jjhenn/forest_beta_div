
Borcard, D., and P. Legendre. 2002. All-scale spatial analysis of
ecological data by means of principal coordinates of
neighbour matrices. Ecological Modelling 153:51–68.
Borcard, D., P. Legendre, C. Avois-Jacquet, and H. Tuomisto.
2004. Dissecting the spatial structure of ecological data at
multiple scales. Ecology 85:1826–1832.
Borcard, D., P. Legendre, and P. Drapeau. 1992. Partialling out
the spatial component of ecological variation. Ecology 73:
  1045–1055.



################################################################
# calculate spatial eigenvalues and eigenvectors #
################################################################

library(tidyverse)
library(vegan)

#SERC, TRCP, WFDP
site_n = "WFDP"

coor_all <- read.csv("output/all_env.csv") %>% filter(variable %in% c("gx", "gy")) %>% 
  pivot_wider(names_from = variable, values_from = vals) %>% 
  filter(site == site_n) %>% 
  select(gx, gy)

coor_all1 <- read.csv("output/all_env.csv") %>% filter(variable %in% c("gx", "gy")) %>% 
  pivot_wider(names_from = variable, values_from = vals) %>% 
  filter(site == site_n) %>% 
  select(quadrat)

# Check that all plots have unique coordinates
## number of plots
n_coor_all<-length(as.matrix(coor_all)[,1])
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


mst <- spantree(d) #this claculates the minimum spanning tree


plot(coor_all, col="red")
lines(mst, coor_all)#this plots the minimum spaning tree on the map
ws <- (d<=max(mst$dist))*(1-((d/(4*(max(mst$dist))))^2))
w <- (d>0)*ws
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
write.table(ei_sim_all$"values", file=paste("output/SpatialEigenValues_", site_n, ".txt", sep = ""), quote=T, sep=",", row.names = F)
spatial.eigenvectors <- data.frame(coor_all1, ei_sim_all$"vectors")

write.table(spatial.eigenvectors, file=paste("output/SpatialEigenVectors_", site_n, ".txt", sep = ""), quote=T, sep=",", row.names = F)
