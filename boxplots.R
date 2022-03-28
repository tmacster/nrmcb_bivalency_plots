# R code for quantifying histone mark enrichment over promoters
# takes output from deeptools, filters uncovered regions, plots, performs stats

library(tidyverse)
library(pheatmap)
library(viridis)
library(ggpubr)

setwd("dir/")


genesets <- c("bivalent_3868", "matched_cpg_sample", "random_sample")

# H2Aub -------------------------------------------------------------------

h2a_files <- str_glue("counts/chen_H2Aub_filt_{genesets}.txt")
h2a_all <- h2a_files %>%
  map_dfr(read_tsv, 
          .id = "geneset", skip = 1, 
          col_names = c("chr", "start", "end", "mat_zygote", "pat_zygote", 
                        "mat_2C", "pat_2C", "mat_4C", "pat_4C", "mat_epiblast", "pat_epiblast", 
                        "mat_ICM", "pat_ICM", "mat_oocyte", "mat_morula", "pat_morula", "pat_sperm")) %>%
  select(-chr:-end) %>%
  mutate(geneset = case_when(geneset == '1' ~ "bivalent", geneset == '2' ~ "matched", geneset == '3' ~ "random"))

h2a_long <- h2a_all %>%
  pivot_longer(cols = mat_zygote:pat_sperm, names_to = "stage1", values_to = "coverage", values_drop_na = TRUE) %>%
  separate(stage1, c("parent", "stage")) %>%
  mutate(stage = case_when(stage == "oocyte" ~ "gamete", stage == "sperm" ~ "gamete", TRUE ~ stage)) %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "2C", "4C", "morula", "ICM", "epiblast"))
  #mutate(parent = case_when(grepl("mat", stage1) ~ "maternal", grepl("pat", stage1) ~ "paternal")) %>%
  

# K27me3 ------------------------------------------------------------------

k27_files <- str_glue("counts/zheng_H3K27me3_filt_{genesets}.txt")
k27_all <- k27_files %>%
  map_dfr(read_tsv, 
          .id = "geneset", skip = 1, 
          col_names = c("chr", "start", "end", 
                        "mat_2C", "pat_2C", "mat_8C", "pat_8C", "mat_epiblast", "pat_epiblast", "mat_ICM", "pat_ICM", "mat_oocyte", 
                        "mat_morula", "pat_morula", "mat_zygote", "pat_zygote", "pat_sperm", "pat_sperm_chen")) %>%
  select(-chr:-end, -pat_sperm_chen, -matches("morula")) %>% #remove morula data (from a different study)
  mutate(geneset = case_when(geneset == '1' ~ "bivalent", geneset == '2' ~ "matched", geneset == '3' ~ "random"))

k27_long <- k27_all %>%
  pivot_longer(cols = mat_2C:pat_sperm, names_to = "stage1", values_to = "coverage", values_drop_na = TRUE) %>%
  separate(stage1, c("parent", "stage")) %>%
  mutate(stage = case_when(stage == "oocyte" ~ "gamete", stage == "sperm" ~ "gamete", TRUE ~ stage)) %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "2C", "8C", "ICM", "epiblast"))


# K4me3 -------------------------------------------------------------------

k4_files <- str_glue("counts/zhang_H3K4me3_filt_{genesets}.txt")
k4_all <- k4_files %>%
  map_dfr(read_tsv, 
          .id = "geneset", skip = 1, 
          col_names = c("chr", "start", "end", 
                        "mat_early_2C", "pat_early_2C", "mat_late_2C", "pat_late_2C", "mat_4C", "pat_4C", 
                        "mat_8C", "pat_8C", "mat_ICM", "pat_ICM", "mat_oocyte", "pat_sperm", "mat_zygote", "pat_zygote")) %>%
  select(-chr:-end) %>%
  mutate(geneset = case_when(geneset == '1' ~ "bivalent", geneset == '2' ~ "matched", geneset == '3' ~ "random"))

k4_long <- k4_all %>%
  pivot_longer(cols = mat_early_2C:pat_zygote, names_to = "stage1", values_to = "coverage", values_drop_na = TRUE) %>%
  separate(stage1, c("parent", "stage"), extra = "merge") %>%
  mutate(stage = case_when(stage == "oocyte" ~ "gamete", stage == "sperm" ~ "gamete", TRUE ~ stage))


# combine & export  ------------------------------------------------------------
all_data_long <- bind_rows("H3K4me3" = k4_long, "H3K27me3" = k27_long, "H2AK119ub1" = h2a_long,
                           .id = "histone") %>%
  unite("chip", histone:geneset, remove = FALSE) %>%
  mutate(stage = fct_relevel(stage, "gamete", "zygote", "early_2C", "2C", "late_2C", "4C", "8C", "morula", "ICM", "epiblast"))



# print plots -------------------------------------------------------------

# set themes
parent <- c(mat = "maternal", pat = "paternal")
h2a_colors <- c("maroon", "darkgrey", "lightgrey")
k27_colors <- c("#EFBC68", "darkgrey", "lightgrey")
k4_colors <- c("#5F9595", "darkgrey", "lightgrey")

theme_var <- theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color="black"), 
        axis.text.y = element_text(size = 12, color="black"), 
        axis.ticks.length = unit(.25, "cm"), 
        axis.line = element_line(size = 0.5, color = "black"), 
        panel.border = element_rect(fill = "transparent", color = "transparent"), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        panel.spacing = unit(1, "lines"), 
        legend.position = "top")

# will plot individually for export
all_data_long %>% 
  filter(histone == "H2AK119ub1") %>%
  ggplot(aes(stage, log2(coverage + 0.1), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~parent, scales = "free", 
             labeller = labeller(parent = parent)) +
  scale_fill_manual(values = h2a_colors) +
  scale_x_discrete(name = "stage") +
  ylim(c(-3,6)) +
  theme_var
ggsave("/figures/h2a_boxplots.pdf", device = "pdf", dpi = 300, width = 6, height = 3, units = "in")

all_data_long %>% 
  filter(histone == "H3K27me3") %>%
  ggplot(aes(stage, log2(coverage + 0.1), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~parent, scales = "free", 
             labeller = labeller(parent = parent)) +
  scale_fill_manual(values = k27_colors) +
  scale_x_discrete(name = "stage") +
  ylim(c(-3,6)) +
  theme_var
ggsave("/figures/k27_boxplots.pdf", device = "pdf", dpi = 300, width = 6, height = 3, units = "in")


all_data_long %>% 
  filter(histone == "H3K4me3") %>%
  ggplot(aes(stage, log2(coverage + 0.1), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~parent, scales = "free", 
             labeller = labeller(parent = parent)) +
  scale_fill_manual(values = k4_colors) +
  scale_x_discrete(name = "stage") +
  ylim(c(-2,5)) +
  theme_var
ggsave("/figures/k4_boxplots.pdf", device = "pdf", dpi = 300, width = 6, height = 3, units = "in")


# can also plot together
all_data_long %>% 
  ggplot(aes(stage, log2(coverage + 0.1), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(parent ~ histone, scales = "free", 
             labeller = labeller(parent = parent)) +
  scale_fill_manual(values = c("maroon", "darkgrey", "lightgrey", "#EFBC68", "darkgrey", "lightgrey", "#5F9595", "darkgrey", "lightgrey")) +
  scale_x_discrete(name = "stage") +
  ylim(c(-3,6)) +
  theme_var


# statistics --------------------------------------------------------------
library(rstatix)

stats_aov <- all_data_long %>% 
  group_by(histone, parent, stage) %>%
  anova_test(coverage ~ chip)
#write_tsv(stats_aov, "stats_FigS1_anova.txt")


stats_hsd <- all_data_long %>% 
  group_by(histone, parent, stage) %>%
  tukey_hsd(coverage ~ chip) 
#write_tsv(stats_hsd, "stats_FigS1_posthoc.txt")


stats_hsd_subset <- stats_hsd %>%
  filter(if_any(contains("group"), ~ str_detect(., "_bivalent")))  #filter only bivalent comparisons

stats_hsd_subset %>%
  filter(grepl("matched", group2)) %>% #filter only matched cpg control
  filter(p.adj.signif == "ns")

stats_hsd_subset_sig <- stats_hsd_subset %>%
  filter(grepl("matched", group2)) %>% #filter only matched cpg control
  filter(p.adj.signif != "ns" & p.adj.signif < 2.2e-16)
