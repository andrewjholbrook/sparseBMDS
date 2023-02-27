# Import libraries
library(tidyverse)
library(invgamma)
library(mvtnorm)
library(truncnorm)
library(ggpubr)
library(wesanderson)

# Set working direction
#setwd("/Users/amisheth/sparseBMDS/code")

# Import data
# grid_data <- read_csv("txt_simdata/error_grid_results3.txt",
#                        col_names = c("SD", "BAND", "ESS", "X", "Y", "color"))
grid_data <- read_csv("txt_simdata/error_grid_cmdsD_results.txt",
                       col_names = c("SD", "BAND", "ESS", "X", "Y", "color"))

# RESULTS -> N(0, 1) PRIOR
# RESULTS2 -> N(0, 10) PRIOR
# RESULTS3 -> N(0, 100) PRIOR

# cmds -> initial position using cmds on dist_test2 (distance matrix)
# cmdsD -> initial position using cmds on D.noise (dist matrix with noise)

grid_data %>% filter(ESS < 100) #%>% print(n = 25)

# Graphs for Metropolis

#notes: iter = 300000, burnin = 1, stepsize = .1, N(0, 1) prior for latent variables 
#gives all min ess > 100, but issues with 5 bands 

#### cmds on D.noise as initial position

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = rel(1)),
             legend.position = "none",
             #aspect.ratio = 1,
             axis.text.x = element_text(size = 7),
             axis.text.y = element_text(size = 7))
             #plot.margin = unit(c(-0.75, 0, -.8, 0, "mm")))

pal <- wes_palette(name = "Darjeeling1", n = 5, type = "discrete")

# 25 bands
g1.25.add <- grid_data %>% filter(SD == .05 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) +
  labs(x = NULL,
       y = "25 BANDS",
       title = "SD = .05") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g2.25.add <- grid_data %>% filter(SD == .1 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "SD = .1") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g3.25.add <- grid_data %>% filter(SD == .2 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "", 
       title = "SD = .2") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g4.25.add <- grid_data %>% filter(SD == .4 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "SD = .4") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

# 15 bands
g1.15.add <- grid_data %>% filter(SD == .05 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "15 BANDS",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g2.15.add <- grid_data %>% filter(SD == .1 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "", 
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g3.15.add <- grid_data %>% filter(SD == .2 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g4.15.add <- grid_data %>% filter(SD == .4 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

# 10 bands
g1.10.add <- grid_data %>% filter(SD == .05 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "10 BANDS",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g2.10.add <- grid_data %>% filter(SD == .1 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g3.10.add <- grid_data %>% filter(SD == .2 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g4.10.add <- grid_data %>% filter(SD == .4 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

# 5 bands
g1.5.add <- grid_data %>% filter(SD == .05 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "5 BANDS",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g2.5.add <- grid_data %>% filter(SD == .1 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g3.5.add <- grid_data %>% filter(SD == .2 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

g4.5.add <- grid_data %>% filter(SD == .4 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  scale_color_manual(values = pal) + 
  labs(x = NULL,
       y = "",
       title = "") +
  xlim(c(-2, 2)) + 
  ylim(c(-2, 2))

gg_all <- ggarrange(g1.25.add, g2.25.add, g3.25.add, g4.25.add, 
                    g1.15.add, g2.15.add, g3.15.add, g4.15.add, #ncol = 4, nrow = 2) 
                    g1.10.add, g2.10.add, g3.10.add, g4.10.add, 
                    g1.5.add, g2.5.add, g3.5.add, g4.5.add,
                    nrow = 4, ncol = 4)

ggsave("error_grid_graphs.png", plot = gg_all, path = "results_graphs/", bg = "white",
       height = 8, width = 8)

# ggsave("error_grid_graphs.png", 
#        plot = ggarrange(g1.25.add, g2.25.add, g3.25.add, g4.25.add, 
#                         g1.15.add, g2.15.add, g3.15.add, g4.15.add, 
#                         g1.10.add, g2.10.add, g3.10.add, g4.10.add, 
#                         g1.5.add, g2.5.add, g3.5.add, g4.5.add, 
#                         nrow = 4, ncol = 4), 
#        path = "results_graphs/", bg = "white")
