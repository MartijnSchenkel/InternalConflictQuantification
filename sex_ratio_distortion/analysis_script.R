# ============================ #
# ---- Administrative setup ----
# ============================ #
library(tidyverse)
library(viridis)
library(cowplot)
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ======================= #
# ---- Data processing ----
# ======================= #

# Load raw data
d <- read.table("data.txt", T) %>% as_tibble()

# Split out entity-specific optima
dz <- d %>% gather(3:8, key = "Element", value = "zG") %>% 
  mutate(Element = str_sub(Element, 4, 7))

# Split out individual-level fitness scores at optima of different entities
dw <- d %>% gather(10:15, key = "Element", value = "wG") %>% 
  mutate(Element = str_sub(Element, 4, 7),
         wI = wI_star)

# Drop variables that are not used in trait- vs. fitness-based analyses
dw2 <- dw %>% select(-(3:8))
dz2 <- dz %>% select(-(3:9))

# Combine into 1 data frame, calculate deviations from fitness and trait optima
d2 <- dw2 %>% full_join(dz2)
d2 <- d2 %>% mutate(dw = (wI - wG) / wI,
              dz = ifelse(zG < zI, abs(zI - zG) / zI, abs(zI - zG) / (1 - zI)),
              Gene = str_sub(Element, 1, 1),
              Genotype = str_sub(Element, 3, 4),
              Count = ifelse(Genotype %in% c("XX", "ZZ"), 2, 1)) %>%
  select(-Element)


# Collect data for XX and XY individuals
d3z_XY <- d2 %>% filter(Genotype %in% c("XX", "XY")) %>% 
  spread(Gene, dz) %>% rename(zX = X, zY = Y) %>% 
  mutate(zY = replace_na(zY, 0),
         zX = replace_na(zX, 0)) %>% 
  group_by(n, Genotype) %>% 
  summarize(zX = sum(zX), zY = sum(zY), Count = Count)


d3w_XY <- d2 %>% filter(Genotype %in% c("XX", "XY")) %>% 
  spread(Gene, dw) %>% rename(wX = X, wY = Y) %>% 
  mutate(wY = replace_na(wY, 0),
         wX = replace_na(wX, 0)) %>% 
  group_by(n, Genotype) %>% 
  summarize(wX = sum(wX), wY = sum(wY), wI = wI, Count = Count)

# Collect data for ZZ and ZW individuals
d3z_ZW <- d2 %>% filter(Genotype %in% c("ZZ", "ZW")) %>% 
  spread(Gene, dz) %>% rename(zZ = Z, zW = W) %>% 
  mutate(zW = replace_na(zW, 0),
         zZ = replace_na(zZ, 0)) %>% 
  group_by(n, Genotype) %>% 
  summarize(zZ = sum(zZ), zW = sum(zW), Count = Count)

d3w_ZW <- d2 %>% filter(Genotype %in% c("ZZ", "ZW")) %>% 
  spread(Gene, dw) %>% rename(wZ = Z, wW = W) %>% 
  mutate(wW = replace_na(wW, 0),
         wZ = replace_na(wZ, 0)) %>% 
  group_by(n, Genotype) %>% 
  summarize(wZ = sum(wZ), wW = sum(wW), wI = wI, Count = Count)

# Calculate relative optimality (fitness at faction optimum z_G* / fitness at 
# individual optimum z_I^*)
d3o <- d2 %>% mutate(O = wG / wI) %>% 
  spread(Gene, O) %>% rename(oX = X, oY = Y, oZ = Z, oW = W) %>% 
  mutate(oY = replace_na(oY, 0),
         oX = replace_na(oX, 0),
         oZ = replace_na(oZ, 0),
         oW = replace_na(oW, 0)) %>% 
  group_by(n, Genotype) %>% 
  summarize(oX = sum(oX), oY = sum(oY), oZ = sum(oZ), oW = sum(oW),
            wI, Count = Count) %>% select(-wI)


# Recombine into XX/XY and ZZ/ZW data frames
d3_XY <- d3w_XY %>% full_join(d3z_XY) %>% ungroup()
d3_ZW <- d3w_ZW %>% full_join(d3z_ZW) %>% ungroup()

d4z_XY <- d4z_ZW <- tibble(Genotype = character(),
             n = numeric(),
             Uz = numeric(),
             p = numeric())
d4w_XY <- d4w_ZW <- tibble(Genotype = character(),
              n = numeric(),
              Uw = numeric(),
              p = numeric())

# Calculate resulting unity scores when
for(p in seq(0, 1, 0.005))
{
  if(p %in% seq(0,1, 0.1)) {print(p)}
  dt <- d3z_XY %>% mutate(Uz = 1 - p * Count * (zX + zY) / 2, p = p) %>%
    select(Genotype, n, p, Uz)
  d4z_XY <- d4z_XY %>% full_join(dt, by = names(dt))

  dt <- d3w_XY %>% mutate(Uw = 1 - p * Count * (wX + wY) / 2, p = p) %>%
    select(Genotype, n, p, Uw)
  d4w_XY <- d4w_XY %>% full_join(dt, by = names(dt))

  dt <- d3z_ZW %>% mutate(Uz = 1 - p * Count * (zZ + zW) / 2, p = p) %>%
    select(Genotype, n, p, Uz)
  d4z_ZW <- d4z_ZW %>% full_join(dt, by = names(dt))

  dt <- d3w_ZW %>% mutate(Uw = 1 - p * Count * (wZ + wW) / 2, p = p) %>%
    select(Genotype, n, p, Uw)
  d4w_ZW <- d4w_ZW %>% full_join(dt, by = names(dt))
}

# ========================== #
# ---- Data visualization ----
# ========================== #

# Plot trait unity values in XX/XY individuals
pz_XY <- d4z_XY %>% ggplot(aes(n, p)) + 
  geom_contour_filled(aes(z = Uz), breaks = seq(-0.025, 1.025, 0.05), color = rgb(0.2,0.2,0.2)) + 
  scale_fill_viridis(discrete = T, drop = F, labels = seq(0,1,0.05)) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  labs(x = bquote("Number of foundresses per deme ("*italic(n)*")"), 
       y = bquote("Fraction of genome ("*italic(p)*")"), 
       fill = bquote(italic(U[z]))) + #,
       #title = bquote("Trait unity ("*italic(U[z])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5))+ 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Plot fitness unity values in XX/XY individuals
pw_XY <- d4w_XY %>% ggplot(aes(n, p)) + 
  geom_contour_filled(aes(z = Uw), breaks = seq(-0.025, 1.025, 0.05), color = rgb(0.2,0.2,0.2)) + 
  scale_fill_viridis(discrete = T, drop = F, labels = seq(0,1,0.05)) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  labs(x = bquote("Number of foundresses per deme ("*italic(n)*")"), 
       y = bquote("Fraction of genome ("*italic(p)*")"), 
       fill = bquote(italic(U[w]))) + #,
       # title = bquote("Fitness unity ("*italic(U[w])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Plot trait unity values in ZZ/ZW individuals
pz_ZW <- d4z_ZW %>% ggplot(aes(n, p)) + 
  geom_contour_filled(aes(z = Uz), breaks = seq(-0.025, 1.025, 0.05), color = rgb(0.2,0.2,0.2)) + 
  scale_fill_viridis(discrete = T, drop = F, labels = seq(0,1,0.05)) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  labs(x = bquote("Number of foundresses per deme ("*italic(n)*")"), 
       y = bquote("Fraction of genome ("*italic(p)*")"), 
       fill = bquote(italic(U[z]))) + #,
       # title = bquote("Trait unity ("*italic(U[z])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5))+ 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Plot fitness unity values in ZZ/ZW individuals
pw_ZW <- d4w_ZW %>% ggplot(aes(n, p)) + 
  geom_contour_filled(aes(z = Uw), breaks = seq(-0.025, 1.025, 0.05), color = rgb(0.2,0.2,0.2)) + 
  scale_fill_viridis(discrete = T, drop = F, labels = seq(0,1,0.05)) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  labs(x = bquote("Number of foundresses per deme ("*italic(n)*")"), 
       y = bquote("Fraction of genome ("*italic(p)*")"), 
       fill = bquote(italic(U[w]))) + #,
       # title = bquote("Fitness unity ("*italic(U[w])*")")) +
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2, 0.2, 0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Generate plot for XX/XY, i.e. Figure 3.
p1_XY <- plot_grid(pz_XY + theme(legend.position = "none"),
                   NULL,
                   pw_XY + theme(legend.position = "none"), ncol = 1,
                   rel_heights = c(10,0.5,10), labels = c("A", "", "B", ""))

# Generate faux legend that uses entire range of values between 0 and 1, rather than just observed values.
xf <- seq(0.025, 1.025, 0.05)
y <- rep(1, length(xf))

dt <- tibble(xf, y)
pt <- ggplot(data = dt, aes(xf, y, fill = factor(xf))) + 
  geom_col(aes(z = xf), color = "black") +
  scale_fill_viridis(discrete = T, drop = F, labels= seq(0,1,0.05)) + 
  labs(fill = "Unity\n(±0.025)") + 
  guides(fill = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5))


legend <- get_legend(pt)

p2_XY <- plot_grid(p1_XY, legend, ncol = 2, rel_widths = c(5,1))

pdf(file = "Figure_4.pdf", width = 7, height = 7)
p2_XY
dev.off()

png(file = "Figure_4.png",
    width = 7, height = 7, units = "in", res = 1000)
p2_XY
dev.off()

# Generate plot for ZZ/ZW, i.e. Supplementary Figure 4.
p1_ZW <- plot_grid(pz_ZW + theme(legend.position = "none"), 
                   NULL,
                   pw_ZW + theme(legend.position = "none"), ncol = 1,
                   rel_heights = c(10,0.5,10), labels = c("A", "", "B", ""))

p2_ZW <- plot_grid(p1_ZW, legend, ncol = 2, rel_widths = c(5,1))

pdf(file = "Supplementary_Figure_4.pdf", width = 7, height = 7)
p2_ZW
dev.off()

png(file = "Supplementary_Figure_4.png", width = 7, height = 7, units = "in", res = 1000)
p2_ZW
dev.off()

# ========================== #
# ---- Optimal strategies ----
# ========================== #
p1 <- d2 %>% filter(Genotype == "XX" & Gene == "X" | 
                      Genotype == "XY" & Gene %in% c("X", "Y") | 
                      Genotype == "ZZ" & Gene == "Z" | 
                      Genotype == "ZW" & Gene %in% c("Z", "W")) %>% 
  ggplot(aes(n, zG, color = Gene)) +  
  geom_line(aes(n, zI), linewidth = 1, color = "black", linetype = 2) + 
  geom_line(linewidth = 1) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.minor = element_blank()) +
  labs(x = bquote("Foundresses ("*italic(n)*")"), y = bquote("Optimal strategy ("*italic(z)[P]^"*"*")"), fill = bquote(italic(U[w]))) +
  guides(color = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  scale_color_viridis(option = "A", begin = 0.1, end = 0.8, discrete = T)

# Relative optimality under full control
p2 <- d3o %>% rename(X = oX, Y = oY, Z = oZ, W = oW) %>% 
  gather(X, Y, Z, W, key = "Element", value = "Optimality") %>% 
  filter(Genotype == "XX" & Element == "X" | 
           Genotype == "XY" & Element %in% c("X", "Y") | 
           Genotype == "ZZ" & Element == "Z" | 
           Genotype == "ZW" & Element %in% c("Z", "W")) %>% 
  ggplot(aes(n, Optimality, color = Element)) + 
  geom_line(linewidth = 1) +
  facet_wrap(~Genotype) + 
  scale_x_continuous(limits = c(1,20), expand = c(0,0), breaks = c(1,2,3,5,10,15,20)) +
  scale_y_continuous(limits = c(0,1)) + 
  theme(strip.background = element_rect(fill = rgb(0.2, 0.2, 0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.minor = element_blank()) +
  labs(x = bquote("Foundresses ("*italic(n)*")"), y = "Relative fitness", fill = bquote(italic(U[w]))) +
  guides(color = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  scale_color_viridis(option = "A", begin = 0.1, end = 0.8, discrete = T)

legend <- get_legend(p1)
pt <- plot_grid(p1 + theme(legend.position = "none"),
                p2 + theme(legend.position = "none"),
        
                ncol = 1, labels = LETTERS[1:2])
pt2 <- plot_grid(pt, legend, ncol = 2, rel_widths = c(5,1))

pdf(file = "Supplementary_Figure_3.pdf", width = 9, height = 9)
pt2
dev.off()

png(file = "Supplementary_Figure_3.png", width = 9, height = 9, units = "in", res = 1000)
pt2
dev.off()
