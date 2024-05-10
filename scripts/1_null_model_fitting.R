
#load required packages
library(tidyverse)
library(vegan)
library(FD)
library(readxl)

## load species compositon data
all_comp <- read.csv("output/all_comp.csv") 

all_traits <- read.csv("output/all_traits.csv") %>% 
  filter(trait != "Alcohols and polyols",
         trait != "Carbonyl compounds")

nreps<-999

#### test
for(i in c("morph", "meta")){
  #make the comp input file
  SC_Tn <- all_comp %>% 
    ungroup() %>% 
    #filter(site == j) %>% 
    select(-genus, -rel_abund, -code) %>% 
    pivot_wider(names_from = sp, values_from = abu, values_fill = 0) 
  
  site <- as.factor(SC_Tn$site)
  quadrat <- as.factor(SC_Tn$quadrat)
  
  ## make trait datafile
  rsp <- all_traits %>% 
    select(-code) %>% 
    na.omit() %>% 
    ungroup() %>% 
    filter(type == i) %>% 
    select(sp, trait, values) %>%
    distinct() %>% 
    filter(sp %in% colnames(SC_Tn)) %>% 
    pivot_wider(names_from = trait, values_from = values)
  
  
  SC_Tn <- SC_Tn %>% select(rsp$sp)
  
  ##create matrix to store loop results
  # res3= matrix(data=NA, nrow=1, ncol=nreps)
  res4 = matrix(data=NA, nrow=1545, ncol=nreps)
  # sim.floristic.dis.TOC<-c()
  
  
  ##run loop
  for(k in 1:nreps){
    
    comm_SERC = permatfull(SC_Tn %>% filter(site == "SERC") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_SERC = rsp %>% filter(sp %in% colnames(as.data.frame(comm_SERC$perm))) %>% select(-sp)
    
    comm_TRCP = permatfull(SC_Tn %>% filter(site == "TRCP") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_TRCP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_TRCP$perm))) %>% select(-sp)
    
    comm_WFDP = permatfull(SC_Tn %>% filter(site == "WFDP") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_WFDP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_WFDP$perm))) %>% select(-sp)
    
    ##Put into usable form for FD package
    mat_SERC=data.frame(comm_SERC$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    trt_SERC <- as.data.frame(trt_SERC)
    
    rownames(trt_SERC) = names(mat_SERC)
    
    mat2_SERC=as.matrix(mat_SERC)
    ###Calculate Null CWM trait values for each plot####
    res2_SERC = dbFD(trt_SERC, mat2_SERC, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_SERC = as.matrix(res2_SERC$FDis)
    
    
    ##Put into usable form for FD package
    mat_TRCP=data.frame(comm_TRCP$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    trt_TRCP <- as.data.frame(trt_TRCP)
    
    row.names(trt_TRCP) = names(mat_TRCP)
    
    mat2_TRCP=as.matrix(mat_TRCP)
    ###Calculate Null CWM trait values for each plot####
    res2_TRCP = dbFD(trt_TRCP, mat2_TRCP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_TRCP = as.matrix(res2_TRCP$FDis)
    
    
    ##Put into usable form for FD package
    mat_WFDP=data.frame(comm_WFDP$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    trt_WFDP <- as.data.frame(trt_WFDP)
    
    row.names(trt_WFDP) = names(mat_WFDP)
    
    mat2_WFDP=as.matrix(mat_WFDP)
    ###Calculate Null CWM trait values for each plot####
    res2_WFDP = dbFD(trt_WFDP, mat2_WFDP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_WFDP = as.matrix(res2_WFDP$FDis)
    
    res3 <- rbind(res3_WFDP, res3_TRCP, res3_SERC)
    
    res4[,k]=res3
  }
  
  out <- as.data.frame(res4) %>% mutate(type = i, site = site, quadrat = quadrat)
  
  write.csv(out, file = paste("output/NULL_fids1",i,".csv", sep = "_"), row.names=F)
  
}





out2 <- data.frame(NULL)
for(i in c("morph", "meta")){
  #make the comp input file
  SC_Tn <- all_comp %>% 
    ungroup() %>% 
    #filter(site == j) %>% 
    select(-genus, -rel_abund, -code) %>% 
    pivot_wider(names_from = sp, values_from = abu, values_fill = 0) 
  
  site <- as.factor(SC_Tn$site)
  quadrat <- as.factor(SC_Tn$quadrat)
  
  rsp <- all_traits %>% 
    select(-code) %>% 
    na.omit() %>% 
    ungroup() %>% 
    filter(type == i) %>% 
    select(sp, trait, values) %>%
    distinct() %>% 
    pivot_wider(names_from = trait, values_from = values) %>% 
    filter(sp %in% colnames(SC_Tn))
  
  SC_Tn <- SC_Tn %>% select(rsp$sp)

  ## randomly draw species from RSP and name to match composition dataset
  
  comm_SERC = SC_Tn %>% filter(site == "SERC") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_SERC = rsp %>% filter(sp %in% colnames(as.data.frame(comm_SERC))) %>% select(-sp) 
  
  comm_TRCP = SC_Tn %>% filter(site == "TRCP") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_TRCP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_TRCP))) %>% select(-sp)
  
  comm_WFDP = SC_Tn %>% filter(site == "WFDP") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_WFDP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_WFDP))) %>% select(-sp)
  
  ##Put into usable form for FD package
  mat_SERC=data.frame(comm_SERC) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  trt_SERC <- as.data.frame(trt_SERC)
  
  row.names(trt_SERC) = names(mat_SERC)
  
  mat2_SERC=as.matrix(mat_SERC)
  ###Calculate Null CWM trait values for each plot####
  res2_SERC = dbFD(trt_SERC, mat2_SERC, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_SERC = as.matrix(res2_SERC$FDis)
  
  
  ##Put into usable form for FD package
  mat_TRCP=data.frame(comm_TRCP) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  trt_TRCP <- as.data.frame(trt_TRCP)
  
  row.names(trt_TRCP) = names(mat_TRCP)
  
  mat2_TRCP=as.matrix(mat_TRCP)
  ###Calculate Null CWM trait values for each plot####
  res2_TRCP = dbFD(trt_TRCP, mat2_TRCP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_TRCP = as.matrix(res2_TRCP$FDis)
  
  
  ##Put into usable form for FD package
  mat_WFDP=data.frame(comm_WFDP) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  trt_WFDP <- as.data.frame(trt_WFDP)
  
  row.names(trt_WFDP) = names(mat_WFDP)
  
  mat2_WFDP=as.matrix(mat_WFDP)
  ###Calculate Null CWM trait values for each plot####
  res2_WFDP = dbFD(trt_WFDP, mat2_WFDP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_WFDP = as.matrix(res2_WFDP$FDis)
  
  out1 <- data.frame(dis = rbind(res3_WFDP, res3_TRCP, res3_SERC)) %>% 
    mutate(type = i, obs = "obs", site = site, quadrat = quadrat)
  out2 <- bind_rows(out1, out2)
}

out1 %>% 
  ggplot(aes(x = site, y = dis)) +
  geom_boxplot()



#### visualize results ####
file.list <- list.files("output", pattern = "NULL_fids1_m", full.names = T)

# null_dat <- file.list %>% 
#   map_dfr(~read.csv(.)) %>% 
#   pivot_longer(-c(site, quadrat, type), names_to = "rep", values_to = "dis_exp") %>% 
#   #mutate(obs = "exp") %>% 
#   select(-rep)

all_null <- file.list %>% 
  map_dfr(~read.csv(.)) %>% 
  pivot_longer(cols = -c(type, site, quadrat), values_to = "fdis", names_to = "run") %>% 
  left_join(out2 %>% select(-obs)) %>%
  group_by(type, site, quadrat, dis) %>% 
  summarize(mean_fdis = mean(fdis),
            sd_fdis = sd(fdis),
            high_ci = quantile(fdis, 0.975),
            low_ci = quantile(fdis, 0.025),
            p = ifelse(dis > mean_fdis, sum(fdis > dis), sum(fdis < dis))/999) %>% 
  mutate(fdis_dev = (dis - mean_fdis)/sd_fdis,
         high_ci_dev = (high_ci - mean_fdis)/sd_fdis,
         low_ci_dev = (low_ci - mean_fdis)/sd_fdis,
         sig = ifelse(fdis_dev > high_ci_dev | fdis_dev < low_ci_dev, "sig", "not"),
         alt_ses = VGAM::probitlink(1 - p),
         sig2 = ifelse((p*2)< 0.05, "sig", "not")) %>% 
  mutate(alt_ses = ifelse(dis > mean_fdis, alt_ses, -alt_ses)) %>% 
  filter(!is.na(sig)) %>% 
  distinct()

all_null %>% 
  ggplot(aes(x = fdis_dev, y = alt_ses, color = sig2)) +
  geom_point() +
  geom_abline(color = "red") +
  geom_hline(aes(yintercept = -1.96)) +
  geom_hline(aes(yintercept = 1.96)) +
  facet_grid(type~site)

cor.test(all_null$fdis_dev, all_null$alt_ses)

write.csv(all_null, file = "output/all_null_fdis.csv", row.names=F)

#look at the distribution of the null model results
dist_dat <- file.list %>% 
  map_dfr(~read.csv(.)) %>%   
  pivot_longer(cols = -c(type, site, quadrat), values_to = "fdis", names_to = "run")

dist_test <- dist_dat %>% 
  group_by(site, type, quadrat) %>% 
  summarize(skew = asbio::skew(fdis))

norm_test <- dist_dat %>% 
  group_by(site, type, quadrat) %>%
  filter(sd(fdis) > 0) %>% 
  do(test = shapiro.test(.$fdis)) %>% 
  mutate(p = test$p.value)

dist_test %>% 
  ggplot(aes(x = fdis)) +
  geom_histogram() +
  facet_grid(site~type, scales = "free")

all_null <-  read.csv("output/all_null_fdis.csv")

all_null %>% 
  ggplot(aes(x = site, y = fdis_dev)) +
  geom_boxplot() +
  facet_wrap(~type)

pc_points <- read.csv("output/pca_points.csv")
# env <- read.csv("output/all_env.csv") %>% 
#   filter(variable == "TEB")

fdis_pc_dat <- all_null %>% 
  left_join(pc_points) 

fdis_pc_dat %>% filter(fdis_dev < 10) %>% 
  ggplot(aes(x = PC2, y = fdis_dev)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(color = site)) +
  stat_smooth(method = "lm", color = "black") +
  facet_grid(type ~ site, scales = "free") +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text= element_text(size = 14))

#### test for spatial autocorrelation in models ####

library(spdep)
points <- read.csv("output/all_env.csv") %>% 
  filter(variable %in% c("gx", "gy")) %>% 
  pivot_wider(names_from = variable, values_from = vals)

fdis_pc_dat <- fdis_pc_dat %>% left_join(points)

rs <- data.frame(rs = rstandard(lm(fdis_dev~PC1 + PC2, data = 
                                     fdis_pc_dat %>% 
                     filter(site == "SERC") %>% 
                     filter(type == "meta"))))

pt_data <- fdis_pc_dat %>% filter(site == "SERC", type == "meta") %>% 
  select(gx, gy)
pt_data <- cbind(pt_data, rs)
coordinates(pt_data) <- c("gx", "gy")
bubble(pt_data, "rs", col = c("blue", "orange"), main = "Residuals", 
       xlab = "X-coordinates", ylab = "Y-coordinates")

library(ncf)
fit1 <- correlog(x = pt_data$gx, y = pt_data$gy, z = pt_data$rs, increment = 5, 
                 resamp = 100, quiet = F)

par(mfrow = c(1,1))
plot(fit1, ylab = "Moran's I", xlab = "Mean distance (coordinate degrees)", 
     main = "Local stability", ylim = c(-.7, .6))
abline(h = 0, lty = 2)
#there is spatial autocorrelation

library(nlme)
m2 <- gls(fdis_dev~PC1 + PC2, data = fdis_pc_dat %>% 
            filter(site == "SERC") %>% 
            filter(type == "meta") %>% 
            na.omit())
vario2 <- Variogram(m2, form = ~gx + gy, resType = "pearson")
plot(vario2, smooth = TRUE)

# Exponential correlation structure
m3 <- gls(fdis_dev~PC1 + PC2, 
          correlation = corExp(form = ~gx + gy, nugget = TRUE), 
          data = fdis_pc_dat %>% 
            filter(site == "SERC") %>% 
            filter(type == "meta") %>% 
            na.omit())

# Gaussian correlation structure
m4 <- gls(fdis_dev~PC1 + PC2, 
          correlation = corGaus(form = ~gx + gy, nugget = TRUE), 
          data = fdis_pc_dat %>% 
            filter(site == "SERC") %>% 
            filter(type == "meta") %>% 
            na.omit())

# spherical correlation strcucture
m5 <- gls(fdis_dev~PC1 + PC2, 
          correlation = corSpher(form = ~gx + gy, nugget = TRUE), 
          data = fdis_pc_dat %>% 
            filter(site == "SERC") %>% 
            filter(type == "meta") %>% 
            na.omit())

# compare models with AIC, lower value has more support
AIC(m2, m3, m4, m5)

#exponential model has best fit. I will use that for all of the models
rs <- data.frame(rs = residuals(m3, type = "normalized"))

pt_data <- fdis_pc_dat %>% filter(site == "SERC", type == "meta") %>% 
  select(gx, gy)
pt_data <- cbind(pt_data, rs)
coordinates(pt_data) <- c("gx", "gy")
bubble(pt_data, "rs", col = c("blue", "orange"), main = "Residuals", 
       xlab = "X-coordinates", ylab = "Y-coordinates")

library(ncf)
fit1 <- correlog(x = pt_data$gx, y = pt_data$gy, z = pt_data$rs, increment = 5, 
                 resamp = 100, quiet = F)

par(mfrow = c(1,1))
plot(fit1, ylab = "Moran's I", xlab = "Mean distance (coordinate degrees)", 
     main = "Local stability", ylim = c(-.7, .6))
abline(h = 0, lty = 2)
#autocorrelation disappears

#### overal model output ####
library(broom.mixed)
points <- read.csv("output/all_env.csv") %>% 
  filter(variable %in% c("gx", "gy")) %>% 
  pivot_wider(names_from = variable, values_from = vals)

fdis_pc_dat <- fdis_pc_dat %>% left_join(points)

gls_mod <- fdis_pc_dat %>% 
  group_by(type, site) %>% 
  filter(fdis_dev < 40) %>% 
  do(tidy(gls(fdis_dev ~ PC1 + PC2, data = ., correlation = corExp(form = ~gx + gy, nugget = TRUE)))) 

slopes_PC <- gls_mod %>% 
  filter(term %in% c("PC1", "PC2")) %>% 
  select(type, site, estimate, std.error, p.value, term)

write.csv(slopes_PC, file = "output/overall_slopes.csv", row.names = F)


#### create trait-specific FDis results ####
all_traits <- read.csv("output/all_traits.csv")

out2 <- data.frame(NULL)
trait <- unique(all_traits$trait)
for(i in trait){
  #make the comp input file
  SC_Tn <- all_comp %>% 
    ungroup() %>% 
    #filter(site == j) %>% 
    select(-genus, -rel_abund, -code) %>% 
    pivot_wider(names_from = sp, values_from = abu, values_fill = 0) 
  
  rsp <- all_traits %>% 
    select(-code) %>% 
    na.omit() %>% 
    ungroup() %>% 
    filter(trait == i) %>% 
    select(sp, trait, values) %>%
    distinct() %>% 
    pivot_wider(names_from = trait, values_from = values) %>% 
    filter(sp %in% colnames(SC_Tn))
  
  SC_Tn <- SC_Tn %>% select(site, quadrat, rsp$sp) %>% filter(rowSums(.[,-c(1:2)]) > 0) 
  
  site <- as.factor(SC_Tn$site)
  quadrat <- as.factor(SC_Tn$quadrat)
  
  SC_Tn <- SC_Tn %>% select(-site, -quadrat)
  comm_SERC = SC_Tn %>% filter(site == "SERC") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_SERC = rsp %>% filter(sp %in% colnames(as.data.frame(comm_SERC))) %>% select(-sp) 
  
  comm_TRCP = SC_Tn %>% filter(site == "TRCP") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_TRCP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_TRCP))) %>% select(-sp)
  
  comm_WFDP = SC_Tn %>% filter(site == "WFDP") %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  trt_WFDP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_WFDP))) %>% select(-sp)
  
  ##Put into usable form for FD package
  mat_SERC=data.frame(comm_SERC) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  row.names(trt_SERC) = names(mat_SERC)
  
  mat2_SERC=as.matrix(mat_SERC)
  ###Calculate Null CWM trait values for each plot####
  res2_SERC = dbFD(trt_SERC, mat2_SERC, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_SERC = as.matrix(cbind(res2_SERC$FDis, res2_SERC$CWM))
  
  
  ##Put into usable form for FD package
  mat_TRCP=data.frame(comm_TRCP) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  row.names(trt_TRCP) = names(mat_TRCP)
  
  mat2_TRCP=as.matrix(mat_TRCP)
  ###Calculate Null CWM trait values for each plot####
  res2_TRCP = dbFD(trt_TRCP, mat2_TRCP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_TRCP = as.matrix(cbind(res2_TRCP$FDis, res2_TRCP$CWM))
  
  
  ##Put into usable form for FD package
  mat_WFDP=data.frame(comm_WFDP) %>% 
    replace(is.na(.), 0) %>% 
    select(order(colnames(.)))
  
  row.names(trt_WFDP) = names(mat_WFDP)
  
  mat2_WFDP=as.matrix(mat_WFDP)
  ###Calculate Null CWM trait values for each plot####
  res2_WFDP = dbFD(trt_WFDP, mat2_WFDP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
  res3_WFDP = as.matrix(cbind(res2_WFDP$FDis, res2_WFDP$CWM))
  
  out1 <- data.frame(rbind(res3_WFDP, res3_TRCP, res3_SERC)) %>% 
    mutate(trait = i, obs = "obs", site = site, quadrat = quadrat)
  colnames(out1) <- c("FDis", "CWM", "trait", "obs", "site", "quadrat")
  out2 <- bind_rows(out1, out2)
}

trait_fdis <- out2 

write.csv(trait_fdis, "output/trait_fdis.csv",  row.names = F)


#### create null results for each trait ####
nreps<-999
trait <- unique(all_traits$trait)[c(13)]
for(i in trait){
  #make the comp input file
  SC_Tn <- all_comp %>% 
    ungroup() %>% 
    #filter(site == j) %>% 
    select(-genus, -rel_abund, -code) %>% 
    pivot_wider(names_from = sp, values_from = abu, values_fill = 0) 
  

  ## make trait datafile
  rsp <- all_traits %>% 
    select(-code) %>% 
    na.omit() %>% 
    ungroup() %>% 
    filter(trait == i) %>% 
    select(sp, trait, values) %>%
    distinct() %>% 
    pivot_wider(names_from = trait, values_from = values) %>% 
    filter(sp %in% colnames(SC_Tn))
  
  SC_Tn <- SC_Tn %>% select(site, quadrat, rsp$sp) %>% filter(rowSums(.[,-c(1:2)]) > 0) 
  
  site <- as.factor(SC_Tn$site)
  quadrat <- as.factor(SC_Tn$quadrat)
  
  SC_Tn <- SC_Tn %>% select(-site, -quadrat)
  ##create matrix to store loop results
  # res3= matrix(data=NA, nrow=1, ncol=nreps)
  res4 = matrix(data=NA, nrow=nrow(SC_Tn), ncol=nreps)
  # sim.floristic.dis.TOC<-c()
  
  
  ##run loop
  for(k in 1:nreps){
    
    comm_SERC = permatfull(SC_Tn %>% filter(site == "SERC") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_SERC = rsp %>% filter(sp %in% colnames(as.data.frame(comm_SERC$perm))) %>% select(-sp)
    
    comm_TRCP = permatfull(SC_Tn %>% filter(site == "TRCP") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_TRCP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_TRCP$perm))) %>% select(-sp)
    
    comm_WFDP = permatfull(SC_Tn %>% filter(site == "WFDP") %>% select_if(~ !is.numeric(.) || sum(.) != 0), fixedmar = "rows", shuffle = "samp", strata = NULL, mtype = "count", times = 1)
    
    trt_WFDP = rsp %>% filter(sp %in% colnames(as.data.frame(comm_WFDP$perm))) %>% select(-sp)
    
    ##Put into usable form for FD package
    mat_SERC=data.frame(comm_SERC$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    row.names(trt_SERC) = names(mat_SERC)
    
    mat2_SERC=as.matrix(mat_SERC)
    ###Calculate Null CWM trait values for each plot####
    res2_SERC = dbFD(trt_SERC, mat2_SERC, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_SERC = as.matrix(res2_SERC$FDis)
    
    
    ##Put into usable form for FD package
    mat_TRCP=data.frame(comm_TRCP$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    row.names(trt_TRCP) = names(mat_TRCP)
    
    mat2_TRCP=as.matrix(mat_TRCP)
    ###Calculate Null CWM trait values for each plot####
    res2_TRCP = dbFD(trt_TRCP, mat2_TRCP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_TRCP = as.matrix(res2_TRCP$FDis)
    
    
    ##Put into usable form for FD package
    mat_WFDP=data.frame(comm_WFDP$perm) %>% 
      replace(is.na(.), 0) %>% 
      select(order(colnames(.)))
    
    row.names(trt_WFDP) = names(mat_WFDP)
    
    mat2_WFDP=as.matrix(mat_WFDP)
    ###Calculate Null CWM trait values for each plot####
    res2_WFDP = dbFD(trt_WFDP, mat2_WFDP, calc.FRic = FALSE, calc.FDiv = FALSE, print.pco = FALSE, messages = FALSE, corr = "lingoes")
    res3_WFDP = as.matrix(res2_WFDP$FDis)
    
    res3 <- rbind(res3_WFDP, res3_TRCP, res3_SERC)
    
    res4[,k]=res3
  }
  
  out <- as.data.frame(res4) %>% mutate(trait = i, site = site, quadrat = quadrat)
  
  write.csv(out, file = paste("output/NULL_fdis",i,".csv", sep = "_"), row.names=F)
  
}

trait_fdis <- read.csv("output/trait_fdis.csv")

files_in <- list.files(path = "output", pattern = "NULL_fdis", full.names = T)

ses_calc <- files_in %>% 
  map_dfr(~read.csv(.)) %>% 
  pivot_longer(cols = -c(trait, site, quadrat), values_to = "fdis", names_to = "run") %>% 
  group_by(trait, site, quadrat) %>% 
  summarize(mean_fdis = mean(fdis),
            sd_fdis = sd(fdis),
            high_ci = quantile(fdis, 0.975),
            low_ci = quantile(fdis, 0.025)) %>% 
  left_join(trait_fdis) %>% 
  mutate(fdis_dev = (FDis - mean_fdis)/sd_fdis,
         high_ci_dev = (high_ci - mean_fdis)/sd_fdis,
         low_ci_dev = (low_ci - mean_fdis)/sd_fdis,
         sig = ifelse(fdis_dev > high_ci_dev | fdis_dev < low_ci_dev, "sig", "not"))

write.csv(ses_calc, "output/trait_fdis_ses.csv", row.names = F)

ses_calc %>% 
  ggplot(aes(x = site, y = fdis_dev)) +
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = "free")

ses_calc <- read.csv("output/trait_fdis_ses.csv") %>% 
  left_join(pc_points)

ses_calc %>% 
  ggplot(aes(x = PC1, y = CWM)) +
  geom_point()





lm_mod_trait <- ses_calc %>% 
  group_by(trait, site) %>% 
  filter(fdis_dev < 40) %>% 
  do(tidy(lm(fdis_dev ~ PC1, data = .))) 

slopes_trait <- lm_mod_trait %>% 
  filter(term == "PC1") %>% 
  select(trait, site, estimate, p.value)
