# ===================== #
# ---- General setup ----
# ===================== #
library(tidyverse)
library(viridis)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ======================== #
# ---- Define functions ----
# ======================== #
basetheme <- function(p)
{
  p + theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2)),
            strip.text = element_text(color = "white"),
            panel.spacing.x = unit(1, "lines")) 
}


fitness <- function(R, r, c, k, z_0, z)
{
  w <- k*(1 - exp(-c*(z - z_0))) + r*R*k*((1 - exp(-c*(z - z_0))))/z
  w
}


# =================================== #
# ---- Preliminary data processing ----
# =================================== #

d <- read.table(str_c("data.txt"), T) %>% as_tibble()

# Calculate inclusive fitness of organism when specific entities control investment 
dp <- d %>% mutate(w_pat2 = fitness(R, 0, c, k, z_0, z_m),
                   w_org2 = fitness(R, 0.25, c, k, z_0, z_m),
                   w_mat2 = fitness(R, 0.5, c, k, z_0, z_m),
                   r = factor(r),
                   c = factor(c))

# Tidy up data frame for easier figure plotting etc.
dp <- dp %>% mutate(fc = str_c("italic(c)==", c)) %>% 
  gather(w_pat2, w_org2, w_mat2, key = "Element", value = "Fitness") %>%  
  mutate(Element = factor(Element))


levels(dp$r) <- c("Paternal", "Organismal", "Maternal")
levels(dp$Element) <- c("Maternal", "Organismal", "Paternal")

# Calculate relative fitness (used for Supplementary Figure 2)
dp2 <- dp %>% group_by(R, c, k, z_0, Element) %>% 
  summarize(Fitness = Fitness, Optimality = Fitness / max(Fitness), r = r)


# Filtered dataframe with optimality under maternal/paternal genome control
dp3 <- dp2 %>% filter(r == 'Organismal', Element != 'Organismal') %>% 
  mutate(fc = str_c("italic(c)==", c))

# Fitness unity data frame
dp4 <- dp2 %>% filter(Element == "Organismal") %>% 
  select(-Optimality) %>% 
  spread(key = r, value = Fitness)

dp4 <- dp4 %>% mutate(P = (Organismal - Paternal)/Organismal,
                      M = (Organismal - Maternal)/Organismal,
                      Uw = 1 - abs(P)/2- abs(M) / 2,
                      fc = str_c("italic(c)==", c, ""))

# Trait unity data frame
dp5 <- dp %>% filter(Element == "Organismal") %>% 
  select(R, r, z_m, c, z_0, fc) %>% 
  spread(key = r, value = z_m)

dp5 <- dp5 %>% mutate(devP = (Paternal - Organismal)/Organismal,
                      devM = (Maternal - Organismal)/Organismal,
                      
                      
                      P = ifelse(devP > 0, abs(devP / (R - Organismal)), 
                                 ifelse(devP < 0, abs(devP / (Organismal - z_0)), 0)),
                      M = ifelse(devM > 0, abs(devM / (R - Organismal)), 
                                 ifelse(devM < 0, abs(devM / (Organismal - z_0)), 0)),
                      
                      Uz = exp(-(P+M)/2),
                      Uz2 = exp(-(P+M)/2) / (1 + exp(-(P+M)/2)))


# ============================== #
# ---- Supplementary Figure 1 ----
# ============================== #

p1 <- dp %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = z_m), breaks = c(1:20), color = "black", linewidth = 0.25) +
  facet_grid(fc~r, labeller = label_parsed) + 
  scale_fill_viridis(discrete = T) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4, 20, 4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote("Optimal allocation ("*{italic(z[i]^'*')}*")"))

p1p <- basetheme(p1)





pdf(file = "Supplementary_Figure_1.pdf", height = 8, width = 9)
p1p
dev.off()


png(file = "Supplementary_Figure_1.png", height = 8, width = 9, res = 1000, units = "in")
p1p
dev.off()

# ============================== #
# ---- Supplementary Figure 2 ----
# ============================== #

# Start with alternate version (continuous values, not binned)
p2 <- dp3 %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_tile(aes(fill = Optimality)) +
  facet_grid(Element ~ fc, labeller = label_parsed) + 
  scale_fill_viridis() + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote("Organismal optimality"))

p2p <- basetheme(p2) + 
  guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "top")) + 
  theme(legend.position = "bottom")

pdf(file = "Supplementary_Figure_2_continuous.pdf", height = 5, width = 9)
p2p
dev.off()


png(file = "Supplementary_Figure_2_continuous.png", height = 5, width = 9, res = 1000, units = "in")
p2p
dev.off()

# Version actually used in manuscript
p3 <- dp3 %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = Optimality), breaks = seq(0,1.05, 0.05), color = "black", linewidth = 0.25) +
  facet_grid(Element ~ fc, labeller = label_parsed) + 
  scale_fill_viridis(discrete = T) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote("Organismal optimality"))

p3p <- basetheme(p3) + 
  guides(fill = guide_legend(nrow = 2, title.hjust = 0.5, title.position = "top")) + 
  theme(legend.position = "bottom")


pdf(file = "Supplementary_Figure_2.pdf", height = 5, width = 9)
p3p
dev.off()

png(file = "Supplementary_Figure_2.png", height = 5, width = 9, res = 1000, units = "in")
p3p
dev.off()

# ============================== #
# ---- Supplementary Figure 3 ----
# ============================== #

# Part 3A
pz <- dp5 %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = 1 - (P+M)/2), color = "black", breaks = seq(-0.025, 1.025, 0.05),  linewidth = 0.25) +
  facet_wrap(~fc, ncol = length(unique(dp5$c)), labeller = label_parsed) + 
  scale_fill_viridis(discrete = T, drop = F) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote(italic(U[z]))) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1, "lines"))

# Part 3B
pw <- dp4 %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = Uw), color = "black", breaks = seq(-0.025, 1.025, 0.05), linewidth = 0.25) +
  facet_wrap(~fc, ncol = length(unique(dp4$c)), labeller = label_parsed) + 
  scale_fill_viridis(discrete = T, drop = F, labels= seq(0,1,0.05)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = "Unity") +
  guides(fill = guide_legend(ncol = 3, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1, "lines"))

# Extract legend
legend <- get_legend(pw)

# Construct figure
p1 <- plot_grid(pz + theme(legend.position = "none"),
                pw + theme(legend.position = "none"),
                ncol = 1, labels = LETTERS[1:2])
p2 <- plot_grid(p1, legend, ncol = 2, rel_widths = c(5,1.5))

# Generate output files
pdf(file = "Supplementary_Figure_3.pdf", width = 9, height = 4.5)
p2
dev.off()

png(file = "Supplementary_Figure_3.png", width = 9, height = 4.5, units = "in", res = 1000)
p2
dev.off()

# ================ #
# ---- Figure 2 ----
# ================ #

# Disclaimer: We are producing multiple versions of this figure in this script.

dp6a <- dp5 %>% filter(c %in% c(0.5, 1.5)) %>% 
  mutate(Return = ifelse(c == 0.5, "Low efficiency", "High efficiency")) %>%
  mutate(Return = factor(Return, levels = c("Low efficiency", "High efficiency")))

dp6b <- dp4 %>% filter(c %in% c(0.5, 1.5)) %>% 
  mutate(Return = ifelse(c == 0.5, "Low efficiency", "High efficiency")) %>%
  mutate(Return = factor(Return, levels = c("Low efficiency", "High efficiency")))


# Part 2 (top row)
pz <- dp6a %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = 1 - (P+M)/2), color = "black", breaks = seq(-0.05, 1.05, 0.1),  linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6a$c))) + 
  scale_fill_viridis(discrete = T, drop = F) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote(italic(U[z])),
       title = bquote("Trait unity ("*italic(U[z])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        text = element_text(size = 18),
        axis.title = element_text(hjust = 0.5, size = 17),
        plot.title = element_text(hjust = 0.5, size = 17),
        legend.title = element_text(hjust = 0.5, size = 17))

# Part 2 (bottom row)
pw <- dp6b %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = Uw), color = "black", breaks = seq(-0.05, 1.05, 0.1), linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6b$c))) + 
  scale_fill_viridis(discrete = T, drop = F, labels= seq(0,1,0.1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = "Unity\n(±0.025)",
       title = bquote("Fitness unity ("*italic(U[w])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        text = element_text(size = 18),
        axis.title = element_text(hjust = 0.5, size = 17),
        plot.title = element_text(hjust = 0.5, size = 17),
        legend.title = element_text(hjust = 0.5, size = 17))

# Extract legend
legend <- get_legend(pw)

# Construct figure
p1 <- plot_grid(pz + theme(legend.position = "none"),
                NULL,
                pw + theme(legend.position = "none"),
                ncol = 1, rel_heights = c(10,0.5,10))
p2 <- plot_grid(p1, legend, ncol = 2, rel_widths = c(5,1))

p2

# Generate output files
pdf(file = "Figure_2.pdf", width = 7, height = 7)
p2
dev.off()

png(file = "Figure_2.png", width = 7, height = 7, units = "in", res = 1000)
p2
dev.off()


# =============================== #
# ---- Figure 2 - fine-grained ----
# =============================== #

# Part 2 (top row)
pz <- dp6a %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = 1 - (P+M)/2), color = "black", breaks = seq(-0.025, 1.025, 0.05),  linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6a$c))) + 
  scale_fill_viridis(discrete = T, drop = F) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote(italic(U[z])),
       title = bquote("Trait unity ("*italic(U[z])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        text = element_text(size = 18),
        axis.title = element_text(hjust = 0.5, size = 17),
        plot.title = element_text(hjust = 0.5, size = 17),
        legend.title = element_text(hjust = 0.5, size = 17))

# Part 3  (bottom row)
pw <- dp6b %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = Uw), color = "black", breaks = seq(-0.025, 1.025, 0.05), linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6b$c))) + 
  scale_fill_viridis(discrete = T, drop = F, labels= seq(0,1,0.05)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = "Unity\n(±0.025)",
       title = bquote("Fitness unity ("*italic(U[w])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        text = element_text(size = 18),
        axis.title = element_text(hjust = 0.5, size = 17),
        plot.title = element_text(hjust = 0.5, size = 17),
        legend.title = element_text(hjust = 0.5, size = 17))

# Extract legend
legend <- get_legend(pw)

# Construct figure
p1 <- plot_grid(pz + theme(legend.position = "none"),
                NULL,
                pw + theme(legend.position = "none"),
                ncol = 1, rel_heights = c(10,0.5,10))
p2 <- plot_grid(p1, legend, ncol = 2, rel_widths = c(5,1))

p2

# Generate output files
pdf(file = "Figure_2_fine.pdf", width = 7, height = 7)
p2
dev.off()

png(file = "Figure_2_fine.png", width = 7, height = 7, units = "in", res = 1000)
p2
dev.off()

# =============================== #
# ---- Figure 2 - fine-grained ----
# =============================== #

# Part 2 (top row)
pz <- dp6a %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = 1 - (P+M)/2), color = "black", breaks = seq(-0.025, 1.025, 0.05),  linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6a$c))) + 
  scale_fill_viridis(discrete = T, drop = F) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = bquote(italic(U[z])),
       title = bquote("Trait unity ("*italic(U[z])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Part 3  (bottom row)
pw <- dp6b %>% filter(k == 1) %>% ggplot(aes(R, z_0)) + 
  geom_contour_filled(aes(z = Uw), color = "black", breaks = seq(-0.025, 1.025, 0.05), linewidth = 0.25) +
  facet_wrap(~Return, ncol = length(unique(dp6b$c))) + 
  scale_fill_viridis(discrete = T, drop = F, labels= seq(0,1,0.05)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(4,20), breaks = seq(4,20,4)) + 
  labs(x = bquote("Resources ("*italic(R)*")"),
       y = bquote("Minimal cost ("*italic(z[0]*")")),
       fill = "Unity\n(±0.025)",
       title = bquote("Fitness unity ("*italic(U[w])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Extract legend
legend <- get_legend(pw)

# Construct figure
p1 <- plot_grid(pz + theme(legend.position = "none"),
                NULL,
                pw + theme(legend.position = "none"),
                ncol = 1, rel_heights = c(10,0.5,10))
p2 <- plot_grid(p1, legend, ncol = 2, rel_widths = c(5,1))

p2

# Generate output files
pdf(file = "Figure_2_fine_small.pdf", width = 7, height = 7)
p2
dev.off()

png(file = "Figure_2_fine_small.png", width = 7, height = 7, units = "in", res = 1000)
p2
dev.off()
