library(tidyverse)

#### create PCAs of env. characteristics for each plot ####

env <- read.csv("output/all_env.csv") 

#SERC
serc_pc_in <- env %>% 
  filter(site == "SERC") %>% 
  pivot_wider(names_from = variable, values_from = vals)

serc_pca <- prcomp(serc_pc_in %>% select(-site, -quadrat, -gx, -gy, -aspect), center = T, scale. = T)

serc_points <- bind_cols(serc_pc_in %>% select(quadrat, site), PC1 =  serc_pca$x[,c(1)], PC2 = serc_pca$x[,2], PC3 = serc_pca$x[,3])

serc_loadings <- as.data.frame(serc_pca$rotation) %>% 
  mutate(env = row.names(.)) %>% 
  mutate(site = "SERC")

#WFDP
WFDP_pc_in <- env %>% 
  filter(site == "WFDP") %>% 
  pivot_wider(names_from = variable, values_from = vals)

WFDP_pca <- prcomp(WFDP_pc_in %>% select(-site, -quadrat, -gx, -gy, -aspect), center = T, scale. = T)

WFDP_points <- bind_cols(WFDP_pc_in %>% select(quadrat, site), PC1 =  WFDP_pca$x[,c(1)], PC2 = WFDP_pca$x[,2], PC3 = WFDP_pca$x[,3])

WFDP_loadings <- as.data.frame(WFDP_pca$rotation) %>% 
  mutate(env = row.names(.)) %>% 
  mutate(site = "WFDP")

#TRCP
TRCP_pc_in <- env %>% 
  filter(site == "TRCP") %>% 
  pivot_wider(names_from = variable, values_from = vals)

TRCP_pca <- prcomp(TRCP_pc_in %>% select(-site, -quadrat, -gx, -gy, -aspect), center = T, scale. = T)

TRCP_points <- bind_cols(TRCP_pc_in %>% select(quadrat, site), PC1 =  TRCP_pca$x[,c(1)], PC2 = TRCP_pca$x[,2], PC3 = TRCP_pca$x[,3])

TRCP_loadings <- as.data.frame(TRCP_pca$rotation) %>% 
  mutate(env = row.names(.)) %>% 
  mutate(site = "TRCP")

all_points <- bind_rows(serc_points, WFDP_points, TRCP_points)
all_loadings <- bind_rows(serc_loadings, WFDP_loadings, TRCP_loadings)

all_loadings_table <- all_loadings %>% 
  select(site, env, PC1, PC2) %>% 
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "vals") %>% 
  mutate(name = paste(PC, site)) %>% 
  select(-site, -PC) %>% 
  pivot_wider(names_from = name, values_from = vals)
  

write.csv(all_loadings_table, file = "output/all_env_rotation.csv",
          row.names = F)
  
  
summary(serc_pca)
summary(TRCP_pca)
summary(WFDP_pca)

write.csv(all_points, file = "output/pca_points.csv")
write.csv(all_loadings, file = "output/pca_loadings.csv")

serc <- ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(data = serc_points %>% mutate(site = "SERC"), aes(x = PC1, y = PC2, color = site), alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), data = serc_loadings, arrow=arrow(length=unit(0.2,"cm"))) +
  ggrepel::geom_text_repel(data=serc_loadings, aes(x=PC1*10, y=PC2*10, label=env), size = 3) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_color_manual(values = "darkslategray3") +
  xlab(paste("PC1 (", round(summary(serc_pca)$importance[2]*100, 1), "%", ")", sep = "")) +
  ylab(paste("PC2 (", round(summary(serc_pca)$importance[5]*100, 1), "%", ")", sep = ""))

trcp <- ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(data = TRCP_points %>% mutate(site = "TRCP"), aes(x = PC1, y = PC2, color = site), alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), data = TRCP_loadings, arrow=arrow(length=unit(0.2,"cm"))) +
  ggrepel::geom_text_repel(data=TRCP_loadings, aes(x=PC1*10, y=PC2*10, label=env), size =3) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_color_manual(values = "brown2") +
  xlab(paste("PC1 (", round(summary(TRCP_pca)$importance[2]*100, 1), "%", ")", sep = "")) +
  ylab(paste("PC2 (", round(summary(TRCP_pca)$importance[5]*100, 1), "%", ")", sep = ""))

wfdp <- ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_point(data = WFDP_points %>% mutate(site = "WFDP"), aes(x = PC1, y = PC2, color = site), alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), data = WFDP_loadings, arrow=arrow(length=unit(0.2,"cm"))) +
  ggrepel::geom_text_repel(data=WFDP_loadings, aes(x=PC1*10, y=PC2*10, label=env), size = 3) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_color_manual(values = "gold2") +
  xlab(paste("PC1 (", round(summary(WFDP_pca)$importance[2]*100, 1), "%", ")", sep = "")) +
  ylab(paste("PC2 (", round(summary(WFDP_pca)$importance[5]*100, 1), "%", ")", sep = "")) 

env_pca_plot <- ggpubr::ggarrange(serc, trcp, wfdp, nrow = 3, common.legend = F)

jpeg(file = "output/env_pca_plot.jpeg",
     width = 3.5,
     height = 9,
     unit = "in",
     res = 600)
env_pca_plot
dev.off()
