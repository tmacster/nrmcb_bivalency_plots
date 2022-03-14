# code to create promoter profile plots for chromatin marks
# imports output data from deeptools, applies "upstream tail normalization" to correct for global enrichment differences

library(tidyverse)
library(viridis)
library(fs)
library(jsonlite)

setwd("~/profile_plots")
#https://www.biostars.org/p/230599/


# H3K27me3 ------------------------------------------------------------
# read in, clean up plotProfile data
col_k27 <- c("bins", "pre_2cell", "pre_8cell", "epiblast", "ICM", "zygote", "gamete")

k27_maternal_bivalent <- t(read.delim("H3K27me3/k27_maternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_maternal_bivalent) <- col_k27

k27_maternal_random <- t(read.delim("H3K27me3/k27_maternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_maternal_random) <- col_k27

k27_mat <- as_tibble(k27_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k27_maternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_8cell, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

k27_paternal_bivalent <- t(read.delim("H3K27me3/k27_paternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_paternal_bivalent) <- col_k27

k27_paternal_random <- t(read.delim("H3K27me3/k27_paternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_paternal_random) <- col_k27

k27_pat <- as_tibble(k27_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k27_paternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_8cell, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

k27_all <- bind_rows("maternal" = k27_mat, "paternal" = k27_pat, .id = "parent")

k27_all %>%
  ggplot(aes(x = bins, y = counts, color = stage), fill = "white") +
  geom_smooth() +
  facet_grid(parent ~ geneset) +
  coord_fixed(ratio = 50) +
  theme_classic() +
  theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"), 
        axis.ticks.length=unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), name = "Normalized counts") + xlab("Bins (vs TSS)") + ggtitle("Maternal H3K27me3")

# bivalent/random -----------------------------------------------------------
k27_wide <- k27_all %>%
  pivot_wider(names_from = geneset, values_from = counts) %>%
  mutate(norm_counts = bivalent/random)

k27_wide %>%
  ggplot(aes(bins, norm_counts, color = stage)) +
  geom_smooth() +
  facet_grid(~parent) +
  theme_classic()

# "upstream tail" normalization (see: https://stackoverflow.com/questions/27117429/scale-relative-to-a-value-in-each-group-via-dplyr)
k27_wide_upstream <- k27_wide %>%
  select(-norm_counts) %>%
  group_by(parent, stage) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), random_norm = random/first(random, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "pre_2cell", "pre_8cell", "ICM", "epiblast"))
  
# plot maternal/paternal - bivalent/random
k27_wide_upstream %>%
  ggplot(aes(bins, norm_counts, color = stage, fill = stage)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(se = FALSE) +
  facet_grid(geneset ~ parent) +
  theme_classic() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(size=12, color="black"), axis.text.y = element_text(size=12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(0.75, 1.5), expand = expansion(mult = c(0, .1)), name = "Normalized counts") + 
  xlab("Bins (vs TSS)") + ggtitle(" H3K27me3")

# plot by stage
k27_wide_upstream %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts, color = geneset)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage) +
  theme_classic() +
  scale_color_manual(values = c("lightgrey", "goldenrod")) +
  theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"), 
        axis.ticks.length=unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(0.75, 1.5), expand = expansion(mult = c(0, .1)), name = "Normalized counts") + 
  xlab("Distance from TSS (kb)") + ggtitle(" H3K27me3") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0, 100, 200), labels = c("-10", "TSS", "+10"))



# Mei H2Aub read-in -------------------------------------------------------
mat_col_h2a <- c("bins", "zygote", "oocyte_7d", "ICM", "early_2cell", "FGO", "late_2cell", "MII", "morula")
pat_col_h2a <- c("bins", "zygote", "ICM", "early_2cell", "late_2cell", "morula")

h2a_maternal_bivalent <- t(read.delim("H2Aub/h2a_mei_maternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
h2a_maternal_random <- t(read.delim("H2Aub/h2a_mei_maternal_random_data.txt", skip = 1, header = F))[-1:-2,]
h2a_paternal_bivalent <- t(read.delim("H2Aub/h2a_mei_paternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
h2a_paternal_random <- t(read.delim("H2Aub/h2a_mei_paternal_random_data.txt", skip = 1, header = F))[-1:-2,]

# H2Aub  ------------------------------------------------------------------
col_h2a <- c("bins", "pre_2cell", "pre_4cell", "epiblast", "ICM", "morula", "zygote", "gamete")

h2a_maternal_bivalent <- t(read.delim("H2Aub/h2a_maternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
#read.delim("H2Aub/h2a_maternal_bivalent_data.txt", skip = 1, header = F)[,1] #CHECK COLUMN ORDER
colnames(h2a_maternal_bivalent) <- col_h2a

h2a_maternal_random <- t(read.delim("H2Aub/h2a_maternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_maternal_random) <- col_h2a

h2a_mat <- as_tibble(h2a_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(h2a_maternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_4cell, morula, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

h2a_paternal_bivalent <- t(read.delim("H2Aub/h2a_paternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_paternal_bivalent) <- col_h2a

h2a_paternal_random <- t(read.delim("H2Aub/h2a_paternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_paternal_random) <- col_h2a

h2a_pat <- as_tibble(h2a_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(h2a_paternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_4cell, morula, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

h2a_all <- bind_rows("maternal" = h2a_mat, "paternal" = h2a_pat, .id = "parent")

h2a_all %>%
  ggplot(aes(x = bins, y = counts, color = stage), fill = "white") +
  geom_smooth() +
  facet_grid(parent ~ geneset) +
  theme_classic() +
  theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"), 
        axis.ticks.length=unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"),
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), name = "Normalized counts") + xlab("Bins (vs TSS)") + ggtitle("H2Aub") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)))

# "upstream tail" normalization (see: https://stackoverflow.com/questions/27117429/scale-relative-to-a-value-in-each-group-via-dplyr)
h2a_wide <- h2a_all %>%
  pivot_wider(names_from = geneset, values_from = counts)

h2a_wide_upstream <- h2a_wide %>%
  group_by(parent, stage) %>%
  #mutate(bivalent_norm = bivalent/nth(bivalent, 5, order_by = bins), random_norm = random/nth(random, 5, order_by = bins)) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), random_norm = random/first(random, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "pre_2cell", "pre_4cell", "morula", "ICM", "epiblast"))

# plot maternal/paternal - bivalent/random
h2a_wide_upstream %>%
  ggplot(aes(bins, norm_counts, color = stage)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(size = 0.8, span = 0.3, se = FALSE) +
  facet_grid(parent ~ geneset) +
  theme_classic() +
  #scale_color_viridis_d() +
  scale_color_manual(values = c("#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#91003f")) +
  #scale_color_manual(values = colorRampPalette(brewer.pal(7, "YlGnBu"))(12)[6:12]) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(0.8, 1.75), expand = expansion(mult = c(0, .1)), name = "Normalized counts") + 
  xlab("Bins (vs TSS)") + ggtitle(" H2Aub") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0)), breaks = c(0, 100, 200), labels = c("-10", "TSS", "+10"))

# plot by stage
h2a_wide_upstream %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts)) +
  #geom_rect(aes(xmin = 75, xmax = 125, ymin = 0.75, ymax = 1.8), fill = "grey92") +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(aes(color = geneset), span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage) +
  theme_classic() +
  scale_color_manual(values = c("lightgrey", "maroon")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0.75, 1.8), name = "Normalized counts") + 
  xlab("Distance from TSS (kb)") + ggtitle("H2Aub") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0, 100, 200), labels = c("-10", "TSS", "+10"))


# H3K4me3 -----------------------------------------------------------------
col_k4 <- c("bins", "early_2cell", "late_2cell", "pre_4cell", "pre_8cell", "ICM", "zygote", "gamete")
k4_maternal_bivalent <- t(read.delim("H3K4me3/k4_maternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
#read.delim("H3K4me3/k4_maternal_bivalent_data.txt", skip = 1, header = F)[,1] #CHECK COLUMN ORDER
colnames(k4_maternal_bivalent) <- col_k4

k4_maternal_random <- t(read.delim("H3K4me3/k4_maternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_maternal_random) <- col_k4

k4_mat <- as_tibble(k4_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k4_maternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, early_2cell, late_2cell, pre_4cell, pre_8cell, ICM, geneset) %>%
  mutate_at(vars(bins:ICM), as.numeric) %>%
  pivot_longer(cols = gamete:ICM, names_to = "stage", values_to = "counts")

k4_paternal_bivalent <- t(read.delim("H3K4me3/k4_paternal_bivalent_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_paternal_bivalent) <- col_k4

k4_paternal_random <- t(read.delim("H3K4me3/k4_paternal_random_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_paternal_random) <- col_k4

k4_pat <- as_tibble(k4_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k4_paternal_random), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, early_2cell, late_2cell, pre_4cell, pre_8cell, ICM, geneset) %>%
  mutate_at(vars(bins:ICM), as.numeric) %>%
  pivot_longer(cols = gamete:ICM, names_to = "stage", values_to = "counts")

k4_all <- bind_rows("maternal" = k4_mat, "paternal" = k4_pat, .id = "parent")

k4_all %>%
  ggplot(aes(x = bins, y = counts, color = stage), fill = "white") +
  geom_smooth() +
  facet_grid(parent ~ geneset) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), name = "Normalized counts") + xlab("Bins (vs TSS)") + ggtitle("H3K4me3")

# "upstream tail" normalization
k4_wide <- k4_all %>%
  pivot_wider(names_from = geneset, values_from = counts)

k4_wide_upstream <- k4_wide %>%
  group_by(parent, stage) %>%
  #mutate(bivalent_norm = bivalent/nth(bivalent, 5, order_by = bins), random_norm = random/nth(random, 5, order_by = bins)) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), random_norm = random/first(random, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "early_2cell", "late_2cell", "pre_4cell", "pre_8cell", "ICM"))

# plot maternal/paternal - bivalent/random
k4_wide_upstream %>%
  ggplot(aes(bins, norm_counts, color = stage, fill = stage)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "darkgrey") +
  geom_smooth(se = FALSE) +
  facet_grid(geneset ~ parent) +
  theme_classic() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(0.75, 1.5), expand = expansion(mult = c(0, .1)), name = "Normalized counts") + 
  xlab("Bins (vs TSS)") + ggtitle("H3K4me3")

# plot by stage
k4_wide_upstream %>%
  filter(stage != "late_2cell") %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts, color = geneset)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage) +
  theme_classic() +
  scale_color_manual(values = c("lightgrey", "#5F9595")) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0.5, 3), name = "Normalized counts") + 
  xlab("Distance from TSS (kb)") + ggtitle("H3K4me3") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0, 100, 200), labels = c("-10", "TSS", "+10"))



# Pointrange quantification -----------------------------------------------
k27_quants <- k27_wide_upstream %>%
  group_by(parent, geneset, stage) %>%
  slice(75:125) %>%
  summarize(median = mean(norm_counts), sd = sd(norm_counts))

k27_quants %>%
  ggplot(aes(stage, median, color = geneset)) +
  geom_point() + geom_line(aes(group = geneset)) +
  #geom_pointrange(aes(ymin = median - sd, ymax = median + sd), position = position_dodge2(width = 0.5), fatten = 2) +
  facet_grid(parent ~ .) +
  theme_classic() +
  scale_color_manual(values = c("goldenrod", "lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.1), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  ylab("Norm. counts (median)") + xlab("Stage") + ggtitle("H3K27me3 median counts") + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(labels = c("gamete", "zygote", "2C", "8C", "ICM", "E6.5"))

# H2Aub 
h2a_quants <- h2a_wide_upstream %>%
  group_by(parent, geneset, stage) %>%
  slice(75:125) %>%
  summarize(median = mean(norm_counts), sd = sd(norm_counts))

h2a_quants %>%
  ggplot(aes(stage, median, color = geneset)) +
  geom_point() + geom_line(aes(group = geneset)) +
  #geom_pointrange(aes(ymin = median - sd, ymax = median + sd), position = position_dodge2(width = 0.5), fatten = 2) +
  facet_grid(parent ~ .) +
  theme_classic() +
  scale_color_manual(values = c("maroon", "lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.1), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  ylab("Norm. counts (median)") + xlab("Stage") + ggtitle("H2Aub median counts") + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(labels = c("gamete", "zygote", "2C", "4C", "morula", "ICM", "E6.5"))

# H3K4me3 
k4_quants <- k4_wide_upstream %>%
  group_by(parent, geneset, stage) %>%
  slice(75:125) %>%
  summarize(median = mean(norm_counts), sd = sd(norm_counts))

k4_quants %>%
  ggplot(aes(stage, median, color = geneset)) +
  geom_point() + geom_line(aes(group = geneset)) +
  #geom_pointrange(aes(ymin = median - sd, ymax = median + sd), position = position_dodge2(width = 0.5), fatten = 2) +
  facet_grid(parent ~ .) +
  theme_classic() +
  scale_color_manual(values = c("#5F9595", "lightgrey")) +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.1), axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside") +
  ylab("Norm. counts (median)") + xlab("Stage") + ggtitle("H3K4me3 median counts") + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(labels = c("gamete", "zygote", "early 2C", "late 2C", "4C", "8C", "morula", "ICM"))

##
