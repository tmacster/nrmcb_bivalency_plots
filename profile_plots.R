# code to create promoter profile plots for chromatin marks
# imports output data from deeptools, applies "upstream tail normalization" to correct for global enrichment differences

library(tidyverse)
library(viridis)
library(fs)
library(jsonlite)
library(lemon)

setwd("profileplots/")
#https://www.biostars.org/p/230599/

labels <- as_labeller(c('gamete' = "gamete", 
                        'zygote' = "zygote", 
                        'early_2cell' = "2 cell", 'pre_2cell' = "2 cell", 
                        'pre_4cell' = "4 cell", 
                        'pre_8cell' = "8 cell",
                        'morula' = "morula",
                        'ICM' = "ICM", 
                        'epiblast' = "epiblast"))

# H3K27me3 ------------------------------------------------------------
col_k27 <- c("bins", "pre_2cell", "pre_8cell", "epiblast", "ICM", "zygote", "gamete")

# read in plotProfile data
k27_maternal_bivalent <- t(read.delim("H3K27me3/k27_maternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_maternal_bivalent) <- col_k27

k27_maternal_random <- t(read.delim("H3K27me3/k27_maternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_maternal_random) <- col_k27

k27_maternal_matched <- t(read.delim("H3K27me3/k27_maternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_maternal_matched) <- col_k27

# paternal
k27_paternal_bivalent <- t(read.delim("H3K27me3/k27_paternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_paternal_bivalent) <- col_k27

k27_paternal_random <- t(read.delim("H3K27me3/k27_paternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_paternal_random) <- col_k27

k27_paternal_matched <- t(read.delim("H3K27me3/k27_paternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k27_paternal_matched) <- col_k27

# clean up
k27_mat <- as_tibble(k27_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k27_maternal_random), 
            "matched" = as_tibble(k27_maternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_8cell, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

k27_pat <- as_tibble(k27_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k27_paternal_random), 
            "matched" = as_tibble(k27_paternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_8cell, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

k27_all <- bind_rows("maternal" = k27_mat, "paternal" = k27_pat, .id = "parent")

# "upstream tail" normalization (see: https://stackoverflow.com/questions/27117429/scale-relative-to-a-value-in-each-group-via-dplyr)
k27_wide <- k27_all %>%
  pivot_wider(names_from = geneset, values_from = counts)

k27_wide_upstream <- k27_wide %>%
  group_by(parent, stage) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), 
         random_norm = random/first(random, order_by = bins),
         matched_norm = matched/first(matched, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "pre_2cell", "pre_8cell", "ICM", "epiblast"))

# plot by stage

k27_wide_upstream %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "matched_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts, color = geneset)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage, 
             labeller = labeller(stage = labels)) +
  scale_color_manual(values = c("lightgrey", "darkgrey", "goldenrod")) +
  theme(axis.text.x = element_text(size = 9, color="black"), 
        axis.text.y = element_text(size = 10, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.x = element_text(size = 12, color = "goldenrod"),
        strip.text.y = element_text(size = 12, color = "black"),
        aspect.ratio = 1, 
        panel.spacing = unit(1, "lines")) +
  scale_y_continuous(limits = c(0.5, 2), expand = expansion(mult = c(0, 0.1)), 
                     breaks = c(0.6, 1.2, 1.8), name = "Normalized counts") + 
  ggtitle("H3K27me3") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(1, 100, 200), 
                     labels = c("-10", "TSS", "+10"), name = "Distance from TSS (kb)")
#ggsave("figures/k27_profile.pdf", device = "pdf", dpi = 300)

# H2Aub  ------------------------------------------------------------------
col_h2a <- c("bins", "pre_2cell", "pre_4cell", "epiblast", "ICM", "morula", "zygote", "gamete")

# read in mat data
h2a_maternal_bivalent <- t(read.delim("H2Aub/h2a_maternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
#read.delim("H2Aub/h2a_maternal_bivalent_data.txt", skip = 1, header = F)[,1] #CHECK COLUMN ORDER
colnames(h2a_maternal_bivalent) <- col_h2a

h2a_maternal_random <- t(read.delim("H2Aub/h2a_maternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_maternal_random) <- col_h2a

h2a_maternal_matched <- t(read.delim("H2Aub/h2a_maternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_maternal_matched) <- col_h2a

# read in pat data
h2a_paternal_bivalent <- t(read.delim("H2Aub/h2a_paternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_paternal_bivalent) <- col_h2a

h2a_paternal_random <- t(read.delim("H2Aub/h2a_paternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_paternal_random) <- col_h2a

h2a_paternal_matched <- t(read.delim("H2Aub/h2a_paternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(h2a_paternal_matched) <- col_h2a

# clean up
h2a_mat <- as_tibble(h2a_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(h2a_maternal_random), 
            "matched" = as_tibble(h2a_maternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_4cell, morula, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

h2a_pat <- as_tibble(h2a_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(h2a_paternal_random), 
            "matched" = as_tibble(h2a_paternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, pre_2cell, pre_4cell, morula, ICM, epiblast, geneset) %>%
  mutate_at(vars(bins:epiblast), as.numeric) %>%
  pivot_longer(cols = gamete:epiblast, names_to = "stage", values_to = "counts")

h2a_all <- bind_rows("maternal" = h2a_mat, "paternal" = h2a_pat, .id = "parent")

# "upstream tail" normalization (see: https://stackoverflow.com/questions/27117429/scale-relative-to-a-value-in-each-group-via-dplyr)
h2a_wide <- h2a_all %>%
  pivot_wider(names_from = geneset, values_from = counts)

h2a_wide_upstream <- h2a_wide %>%
  group_by(parent, stage) %>%
  #mutate(bivalent_norm = bivalent/nth(bivalent, 5, order_by = bins), random_norm = random/nth(random, 5, order_by = bins)) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), random_norm = random/first(random, order_by = bins), 
         matched_norm = matched/first(matched, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "pre_2cell", "pre_4cell", "morula", "ICM", "epiblast"))

# plot by stage

h2a_wide_upstream %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "matched_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(aes(color = geneset), span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage, 
             labeller = labeller(stage = labels)) +
  scale_color_manual(values = c("lightgrey", "darkgrey", "maroon")) +
  theme(axis.text.x = element_text(size = 9, color = "black"), 
        axis.text.y = element_text(size = 10, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.x = element_text(size = 12, color = "maroon"),
        strip.text.y = element_text(size = 12, color = "black"),
        aspect.ratio = 1, 
        panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0.8, 2), 
                     breaks = c(1, 1.5, 2), name = "Normalized counts") + 
  xlab("Distance from TSS (kb)") + ggtitle("H2Aub") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(1, 100, 200), 
                     labels = c("-10", "TSS", "+10"))
#ggsave("figures/h2a_profile.pdf", device = "pdf", dpi = 300)


# H3K4me3 -----------------------------------------------------------------
col_k4 <- c("bins", "early_2cell", "late_2cell", "pre_4cell", "pre_8cell", "ICM", "zygote", "gamete")

# read in mat
k4_maternal_bivalent <- t(read.delim("H3K4me3/k4_maternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_maternal_bivalent) <- col_k4

k4_maternal_random <- t(read.delim("H3K4me3/k4_maternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_maternal_random) <- col_k4

k4_maternal_matched <- t(read.delim("H3K4me3/k4_maternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_maternal_matched) <- col_k4

# read in pat data
k4_paternal_bivalent <- t(read.delim("H3K4me3/k4_paternal_bivalent_3868_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_paternal_bivalent) <- col_k4

k4_paternal_random <- t(read.delim("H3K4me3/k4_paternal_random_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_paternal_random) <- col_k4

k4_paternal_matched <- t(read.delim("H3K4me3/k4_paternal_matched_cpg_sample_data.txt", skip = 1, header = F))[-1:-2,]
colnames(k4_paternal_matched) <- col_k4

# clean up
k4_mat <- as_tibble(k4_maternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k4_maternal_random), 
            "matched" = as_tibble(k4_maternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, early_2cell, pre_4cell, pre_8cell, ICM, geneset) %>%
  mutate_at(vars(bins:ICM), as.numeric) %>%
  pivot_longer(cols = gamete:ICM, names_to = "stage", values_to = "counts")

k4_pat <- as_tibble(k4_paternal_bivalent) %>%
  bind_rows("bivalent" = ., "random" = as_tibble(k4_paternal_random), 
            "matched" = as_tibble(k4_paternal_matched), .id = "geneset") %>%
  na.omit() %>%
  select(bins, gamete, zygote, early_2cell, pre_4cell, pre_8cell, ICM, geneset) %>%
  mutate_at(vars(bins:ICM), as.numeric) %>%
  pivot_longer(cols = gamete:ICM, names_to = "stage", values_to = "counts")

k4_all <- bind_rows("maternal" = k4_mat, "paternal" = k4_pat, .id = "parent")

# "upstream tail" normalization
k4_wide <- k4_all %>%
  pivot_wider(names_from = geneset, values_from = counts)

k4_wide_upstream <- k4_wide %>%
  group_by(parent, stage) %>%
  mutate(bivalent_norm = bivalent/first(bivalent, order_by = bins), random_norm = random/first(random, order_by = bins), 
         matched_norm = matched/first(matched, order_by = bins)) %>%
  ungroup() %>%
  pivot_longer(cols = matches("_norm"), names_to = "geneset", values_to = "norm_counts") %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "early_2cell", "pre_4cell", "pre_8cell", "ICM"))

# plot by stage

k4_wide_upstream %>%
  mutate(geneset = fct_relevel(geneset, "random_norm", "matched_norm", "bivalent_norm")) %>%
  ggplot(aes(bins, norm_counts)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(aes(color = geneset), span = 0.3, size = 0.8, se = FALSE) +
  facet_grid(parent ~ stage, 
             labeller = labeller(stage = labels)) +
  scale_color_manual(values = c("lightgrey", "darkgrey", "#5F9595")) +
  theme(axis.text.x = element_text(size = 9, color = "black"), 
        axis.text.y = element_text(size = 10, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.x = element_text(size = 12, color = "#5F9595"),
        strip.text.y = element_text(size = 12, color = "black"),
        aspect.ratio = 1, 
        panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0.5, 4.25),
                     breaks = c(1, 2, 3, 4), name = "Normalized counts") + 
  xlab("Distance from TSS (kb)") + ggtitle("H3K4me3") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(1, 100, 200), 
                     labels = c("-10", "TSS", "+10"))

#ggsave("figures/k4_profile.pdf", device = "pdf", dpi = 300)

