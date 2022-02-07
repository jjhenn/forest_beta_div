#explore trait data from forests

library(tidyverse)

cwm <- read.csv("output/cwm_all.csv") %>% 
  group_by(trait) %>% 
  mutate(cwm_scale = scale(cwm))

cwm %>% filter(type == "func") %>% 
  ggplot(aes(x = cwm_scale, fill = site, color = site)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~trait, scales = "free")

cwm %>% filter(type == "meta") %>% 
  ggplot(aes(x = cwm_scale, fill = site, color = site)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~trait, scales = "free")

env <- read.csv("output/all_env.csv")

dat <- cwm %>% 
  left_join(env)


dat %>% 
  ggplot(aes(x = vals, y = cwm_scale, color = site, group = site)) +
  stat_smooth() +
  facet_grid(trait~variable, scales = "free")
