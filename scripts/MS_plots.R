library(tidyverse)
library(broom)

sp_traits <- read.csv("output/sp_traits.csv") %>% select(-Carbonyl.compounds, -ChlorophyllContent, -Alcohols.and.polyols) %>%  na.omit()

sp_pca <- prcomp(sp_traits[,c(3:15)], center = T, scale. = T)
# sp_pcoa <- ape::pcoa(dist(sp_traits[,c(3:16)]))

sp_pca_points <- bind_cols(as.data.frame(sp_pca$x)[c(1:3)], sp_traits[,c(1:2,16)])%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))

sp_pca_vecs <- as.data.frame(sp_pca$rotation)[c(1:3)] %>% 
  mutate(trait = rownames(.)) %>% 
  mutate(trait = plyr::mapvalues(trait, from = c("Alcohols.and.polyols", "Benzenoids","Organic.acids.and.derivatives", "Organoheterocyclic.compounds", "Other.glycosides", "Phenolic.glycosides", "Phenylpropanoids.and.polyketides", "Prenol.lipids..terpenoids", "Steroids.and.steroid.derivatives", "LeafSize", "SLA", "LWC", "BarkThickness", "WoodDensity", "SeedMassKEW"), to = c("Alcohols", "Benzenoids","Organic acids", "Organoheterocyclics", "Other glycosides", "Phenolic glycosides", "Phenylpropanoids", "Prenol lipids", "Steroids", "Leaf Mass", "Specific Leaf Area", "Leaf Water Content", "Bark Thickness", "Wood Density", "Seed Mass"))) %>% 
  mutate(type = ifelse(trait %in% c("Alcohols", "Benzenoids","Organic acids", "Organoheterocyclics", "Other glycosides", "Phenolic glycosides", "Phenylpropanoids", "Prenol lipids", "Steroids", "Glycosides"), "Chemical Trait", "Morphological Trait"))

pca_plot <- ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(x = PC1, y = PC2, fill = site), data = sp_pca_points, size = 2, pch = 21) +
  geom_segment(aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6), data = sp_pca_vecs, arrow=arrow(length=unit(0.25,"cm"))) +
  ggrepel::geom_text_repel(data=sp_pca_vecs, aes(x=PC1*6, y=PC2*6, label=trait, color = type), size = 3, show.legend = F) +
  theme_classic() +
  scale_color_manual(values = c("black", "navy")) +
  scale_fill_manual(values = c("darkslategray3", "brown2", "gold2")) +
  theme(legend.title = element_blank(),
        text= element_text(size = 12),
        legend.position = c(0.85,0.2)) +
  xlab(paste("Trait PC1 (", round(summary(sp_pca)$importance[2]*100, 1), "%", ")", sep = "")) +
  ylab(paste("Trait PC2 (", round(summary(sp_pca)$importance[5]*100, 1), "%", ")", sep = ""))

jpeg(file = "output/trait_pca.jpeg",
    width = 8,
    height = 6,
    unit= "in",
    res = 400)
pca_plot
dev.off()



#### make trait group ses plots ####
fdis_ses <- read.csv("output/all_null_fdis.csv")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP"))) %>% 
  mutate(type = plyr::mapvalues(type, from = c("morph", "meta"), to = c("Morphological", "Metabolomic")))

pc_points <- read.csv("output/pca_points.csv")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))

cwm <- read.csv("output/cwm_all.csv")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))%>% 
  mutate(type = plyr::mapvalues(type, from = c("func", "meta"), to = c("Morphological", "Metabolomic")))

trait_ses <- fdis_ses %>% 
  left_join(pc_points) #%>% 
# left_join(cwm) 

# lm_mod_trait <- trait_ses %>% 
#   group_by(type, site) %>% 
#   filter(fdis_dev < 40) %>% 
#   do(tidy(lm(fdis_dev ~ PC1, data = .))) 
# 
# slopes_trait <- lm_mod_trait %>% 
#   filter(term == "PC1") %>% 
#   select(type, site, estimate, p.value) %>% 
#   filter(p.value < 0.05) %>% 
#   select(-estimate)

slopes_trait_PC1 <- read.csv("output/overall_slopes.csv") %>% 
  filter(p.value < 0.05) %>% 
  filter(term == "PC1")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))%>% 
  mutate(type = plyr::mapvalues(type, from = c("morph", "meta"), to = c("Morphological", "Metabolomic")))

#calculate proportion of points that are significantly greater and lower than expected
sig_prop <- trait_ses %>% 
  group_by(site, type) %>% 
  filter(fdis_dev < 10) %>% 
  mutate(sig_type = ifelse(sig == "sig" & fdis_dev < 0, "sig_neg", sig)) %>% 
  summarize(prop_pos = sum(fdis_dev > 0, na.rm = T)/n(),
            prop_neg = sum(fdis_dev < 0, na.rm = T)/n(),
            average = mean(fdis_dev, na.rm = T))%>% 
  mutate(type = plyr::mapvalues(type, from = c("morph", "meta"), to = c("Morphological", "Metabolomic")))


sig_text <- sig_prop %>% 
  mutate(text = paste("Overdispersion:", round(prop_pos, 2), "\nUnderdispersion:", round(prop_neg, 2), "\nAverage Value:", round(average, 1)))

all_fdis_pc1 <- trait_ses %>% filter(fdis_dev < 10) %>% 
  #mutate(sig = ifelse(sig == "not", 0.75, 1)) %>% 
  ggplot(aes(x = PC1, y = fdis_dev)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm", data = slopes_trait_PC1 %>% left_join(trait_ses) %>% filter(fdis_dev < 10) ) +
  geom_text(aes(x = -Inf, y = Inf, label = text, hjust = 0, vjust = 1.1), data = sig_text, size = 3) +
  geom_hline(aes(yintercept = 1.96), color = "lightgray", linetype = "dashed") +
  geom_hline(aes(yintercept = -1.96), color = "lightgray", linetype = "dashed") +
  facet_grid(type~site, scales = "free") +
  #scale_fill_manual(values = c("white", NULL)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Functional Dispersion Deviation")


jpeg(file = "output/all_fdis_pc1.jpeg",
    width = 12,
    height = 7,
    unit = "in",
    res = 400)
all_fdis_pc1
dev.off()



slopes_trait_PC2 <- read.csv("output/overall_slopes.csv") %>% 
  filter(p.value < 0.05) %>% 
  filter(term == "PC2")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))%>% 
  mutate(type = plyr::mapvalues(type, from = c("morph", "meta"), to = c("Morphological", "Metabolomic")))

all_fdis_pc2 <- trait_ses %>% filter(fdis_dev < 10) %>% 
  #mutate(sig = ifelse(sig == "not", 0.75, 1)) %>% 
  ggplot(aes(x = PC2, y = fdis_dev)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm", data = slopes_trait_PC2 %>% left_join(trait_ses) %>% filter(fdis_dev < 10) ) +
  #geom_text(aes(x = -Inf, y = Inf, label = text, hjust = 0, vjust = 1.1), data = sig_text, size = 3) +
  facet_grid(type~site, scales = "free") +
  #scale_fill_manual(values = c("white", NULL)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Functional Dispersion Deviation")


jpeg(file = "output/all_fdis_pc2.jpeg",
     width = 12,
     height = 7,
     unit = "in",
     res = 400)
all_fdis_pc2
dev.off()

#### trait figures ####
pc_points <- read.csv("output/pca_points.csv")%>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP")))

ses_calc_morph <- read.csv("output/trait_fdis_ses.csv") %>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP"))) %>% 
  left_join(pc_points) %>% 
  filter(!trait %in% c("Alcohols and polyols", "Carbonyl compounds", "Other glycosides", "Phenolic glycosides"  )) %>% 
  filter(trait %in% c("LWC", "SLA", "SeedMassKEW", "LeafSize", "WoodDensity", "BarkThickness")) %>% 
  mutate(trait = plyr::mapvalues(trait, from = c("LWC", "SLA", "SeedMassKEW", "LeafSize", "WoodDensity", "BarkThickness"), to = c("Leaf Water\nContent", "Specific Leaf\nArea", "Seed Mass", "Leaf Area", "Wood Density", "Bark Thickness")))

ses_calc_meta <- read.csv("output/trait_fdis_ses.csv") %>% 
  mutate(site = plyr::mapvalues(site, from = c("SERC", "TRCP", "WFDP"), to = c("SERC", "TRCP", "WFDP"))) %>% 
  left_join(pc_points) %>% 
  filter(!trait %in% c("Alcohols and polyols", "Carbonyl compounds", "Other glycosides", "Phenolic glycosides"  )) %>% 
  filter(!trait %in% c("LWC", "SLA", "SeedMassKEW", "LeafSize", "WoodDensity", "BarkThickness")) %>% 
  mutate(trait = plyr::mapvalues(trait, from = c("Benzenoids", "Glycosides", "Organic acids and derivatives", "Organoheterocyclic compounds", "Phenylpropanoids and polyketides", "Prenol lipids, terpenoids", "Steroids and steroid derivatives"), to = c("Benzenoids", "Glycosides", "Organic acids\nand derivatives", "Organoheterocyclic\ncompounds", "Phenylpropanoids\nand polyketides", "Prenol lipids,\nterpenoids", "Steroids and\nsteroid derivatives")))

#### CWM ####
morph_cwm_pc1 <- ses_calc_morph %>% 
  ggplot(aes(x = PC1, y = CWM)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Community Weighted Mean")

jpeg(file = "output/morph_cwm_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_cwm_pc1
dev.off()


morph_cwm_pc2 <- ses_calc_morph %>% 
  ggplot(aes(x = PC2, y = CWM)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Community Weighted Mean")

jpeg(file = "output/morph_cwm_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_cwm_pc2
dev.off()

meta_cwm_pc1 <- ses_calc_meta %>% 
  ggplot(aes(x = PC1, y = CWM)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Community Weighted Mean")

jpeg(file = "output/meta_cwm_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_cwm_pc1
dev.off()

meta_cwm_pc2 <- ses_calc_meta %>% 
  ggplot(aes(x = PC2, y = CWM)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Community Weighted Mean")

jpeg(file = "output/meta_cwm_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_cwm_pc2
dev.off()


#### FDis ####
morph_fdis_pc1 <- ses_calc_morph %>% 
  ggplot(aes(x = PC1, y = FDis)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Functional Dispersion")

jpeg(file = "output/morph_fdis_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_fdis_pc1
dev.off()


morph_fdis_pc2 <- ses_calc_morph %>% 
  ggplot(aes(x = PC2, y = FDis)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Functional Dispersion")

jpeg(file = "output/morph_fdis_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_fdis_pc2
dev.off()

meta_fdis_pc1 <- ses_calc_meta %>% 
  ggplot(aes(x = PC1, y = FDis)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Functional Dispersion")

jpeg(file = "output/meta_fdis_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_fdis_pc1
dev.off()

meta_fdis_pc2 <- ses_calc_meta %>% 
  ggplot(aes(x = PC2, y = FDis)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Functional Dispersion")

jpeg(file = "output/meta_fdis_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_fdis_pc2
dev.off()

#### FDis Deviation ####
morph_fdis_dev_pc1 <- ses_calc_morph %>% 
  ggplot(aes(x = PC1, y = fdis_dev)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Functional Dispersion Deviation")

jpeg(file = "output/morph_fdis_dev_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_fdis_dev_pc1
dev.off()


morph_fdis_dev_pc2 <- ses_calc_morph %>% 
  ggplot(aes(x = PC2, y = fdis_dev)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Functional Dispersion Deviation")

jpeg(file = "output/morph_fdis_dev_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
morph_fdis_dev_pc2
dev.off()

meta_fdis_dev_pc1 <- ses_calc_meta %>% 
  ggplot(aes(x = PC1, y = fdis_dev)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC1") +
  ylab("Functional Dispersion Deviation")

jpeg(file = "output/meta_fdis_dev_pc1.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_fdis_dev_pc1
dev.off()

meta_fdis_dev_pc2 <- ses_calc_meta %>% 
  ggplot(aes(x = PC2, y = fdis_dev)) +
  geom_point(size = 2, aes(color = site)) +
  stat_smooth(color = "black", method = "lm") +
  facet_grid(trait~site, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  scale_color_manual(values = c("darkslategray3", "brown2", "gold2")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Topo-edaphic PC2") +
  ylab("Functional Dispersion Deviation")

jpeg(file = "output/meta_fdis_dev_pc2.jpeg",
     width = 5,
     height = 8,
     unit = "in",
     res = 400)
meta_fdis_dev_pc2
dev.off()
