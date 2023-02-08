# Import libraries
library(tidyverse)
library(invgamma)
library(mvtnorm)
library(truncnorm)
library(gridExtra)

# Import data
grid_data <- read_csv("txt_simdata/error_grid_results.txt", 
                      col_names = c("SD", "BAND", "ESS", "X", "Y", "color"))

grid_data %>% filter(ESS < 100)

# Graphs for Metropolis

theme_update(plot.title = element_text(hjust = 0.5, size = rel(.75)),
             legend.position = "none")

# 25 bands
g1.25.add <- grid_data %>% filter(SD == .1 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "25 bands \n sd = .1",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g2.25.add <- grid_data %>% filter(SD == .25 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "25 bands \n sd = .25",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g3.25.add <- grid_data %>% filter(SD == .4 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "25 bands \n sd = .4",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g4.25.add <- grid_data %>% filter(SD == .55 & BAND == 25) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "25 bands \n sd = .55",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

# 15 bands
g1.15.add <- grid_data %>% filter(SD == .1 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "15 bands \n sd = .1",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g2.15.add <- grid_data %>% filter(SD == .25 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "15 bands \n sd = .25",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g3.15.add <- grid_data %>% filter(SD == .4 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "15 bands \n sd = .4",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g4.15.add <- grid_data %>% filter(SD == .55 & BAND == 15) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "15 bands \n sd = .55",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

# 10 bands
g1.10.add <- grid_data %>% filter(SD == .1 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "10 bands \n sd = .1",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g2.10.add <- grid_data %>% filter(SD == .25 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "10 bands \n sd = .25",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g3.10.add <- grid_data %>% filter(SD == .4 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "10 bands \n sd = .4",
       x = "",
       y = "") + 
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g4.10.add <- grid_data %>% filter(SD == .55 & BAND == 10) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "10 bands \n sd = .55",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

# 5 bands
g1.5.add <- grid_data %>% filter(SD == .1 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "5 bands \n sd = .1",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g2.5.add <- grid_data %>% filter(SD == .25 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "5 bands \n sd = .25",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g3.5.add <- grid_data %>% filter(SD == .4 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "5 bands \n sd = .4",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

g4.5.add <- grid_data %>% filter(SD == .55 & BAND == 5) %>% 
  ggplot(data = .) + 
  geom_point(aes(X, Y, color = color)) +
  labs(title = "5 bands \n sd = .55",
       x = "",
       y = "") +
  coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))

grid.arrange(g1.25.add, g2.25.add, g3.25.add, g4.25.add,
             g1.15.add, g2.15.add, g3.15.add, g4.15.add,
             ncol = 4)
grid.arrange(g1.10.add, g2.10.add, g3.10.add, g4.10.add, 
             g1.5.add, g2.5.add, g3.5.add, g4.5.add, 
             ncol = 4)
