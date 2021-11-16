# R code for quantifying histone mark enrichment over promoters
# takes output from deeptools, filters uncovered regions, plots, performs stats

library(tidyverse)
library(pheatmap)
library(viridis)
library(ggpubr)

setwd("~/")

####################
# H2Aub - Chen ------------------------------------------------------------
##########

# read in data with counts over TSS of ~3500 bivalent genes (defined in Sachs et al)
chen_H2Aub <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/chen_H2Aub_filt.txt", skip = 1, col_names = 
                         c("chr", "start", "end", "zygote_mat", "zygote_pat", "mat_2C", "pat_2C", "mat_4C", "pat_4C", "mat_epiblast", "pat_epiblast", 
                           "mat_ICM", "pat_ICM", "mat_MII_oocyte", "mat_morula", "pat_morula", "sperm")) %>%
  rowid_to_column(var = "ID") #%>% filter_at(vars(-chr:-end), any_vars(. > 0))

# sort by max enrichment in epiblast...
chen_H2Aub_mat <- chen_H2Aub %>% select(-chr:-end) %>%
  select(ID, A_oocyte = mat_MII_oocyte, B_zygote = zygote_mat, C_2C = mat_2C, D_4C= mat_4C, E_morula = mat_morula, F_ICM = mat_ICM, G_epiblast = mat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_epiblast))

chen_H2Aub_pat <- chen_H2Aub %>% select(-chr:-end) %>%
  select(ID, A_sperm = sperm, B_zygote = zygote_pat, C_2C = pat_2C, D_4C= pat_4C, E_morula = pat_morula, F_ICM = pat_ICM, G_epiblast = pat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_epiblast))

# and plot
#pheatmap(log2(chen_H2Aub_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H2Aub, Chen et al")
#pheatmap(log2(chen_H2Aub_pat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H2Aub, Chen et al")


##########
# H3K27me3 ----------------------------------------------------------------
##########

h3k27me3 <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/zheng_H3K27me3_filt.txt", skip = 1, col_names = 
                       c("chr", "start", "end", "mat_2C", "pat_2C", "mat_8C", "pat_8C", "mat_epiblast", "pat_epiblast", "mat_ICM", "pat_ICM", "oocyte", 
                         "mat_morula", "pat_morula", "mat_zygote", "pat_zygote", "sperm", "sperm_chen"))

# trying to sort by max enrichment in epiblast...
h3k27me3_mat <- h3k27me3 %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_oocyte = oocyte, B_zygote = mat_zygote, C_2C = mat_2C, D1_8C = mat_8C, D2_morula = mat_morula, E_ICM = mat_ICM, F_epiblast = mat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(F_epiblast))

h3k27me3_pat <- h3k27me3 %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A1_sperm = sperm, A2_sperm_chen = sperm_chen, B_zygote = pat_zygote, C_2C = pat_2C, D1_8C = pat_8C, D2_morula = pat_morula, E_ICM = pat_ICM, F_epiblast = pat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(F_epiblast))

# and plot
#pheatmap(log2(h3k27me3_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 5, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H3K27me3")

#pheatmap(log2(h3k27me3_pat[,-c(1, 3)]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 5, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H3K27me3")

##########
# H3K4me3 -----------------------------------------------------------------
##########

h3k4me3 <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/zhang_H3K4me3_filt.txt", skip = 1, col_names = 
                       c("chr", "start", "end", "mat_early_2C", "pat_early_2C", "mat_late_2C", "pat_late_2C", "mat_4C", "pat_4C", 
                         "mat_8C", "pat_8C", "mat_ICM", "pat_ICM", "oocyte", "sperm", "mat_zygote", "pat_zygote"))

h3k4me3_mat <- h3k4me3 %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_oocyte = oocyte, B_zygote = mat_zygote, C_early_2C = mat_early_2C, D_late_2C = mat_late_2C, E_4C = mat_4C, F_8C = mat_8C, G_ICM = mat_ICM) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_ICM))

h3k4me3_pat <- h3k4me3 %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_sperm = sperm, B_zygote = pat_zygote, C_early_2C = pat_early_2C, D_late_2C = pat_late_2C, E_4C = pat_4C, F_8C = pat_8C, G_ICM = pat_ICM) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_ICM))

# and plot
#pheatmap(log2(h3k4me3_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H3K4me3")

#pheatmap(log2(h3k4me3_pat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H3K4me3")

####################
# H2Aub random ------------------------------------------------------------
####################

chen_H2Aub_random <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/chen_H2Aub_random_filt.txt", skip = 1, col_names = 
                         c("chr", "start", "end", "zygote_mat", "zygote_pat", "mat_2C", "pat_2C", "mat_4C", "pat_4C", "mat_epiblast", "pat_epiblast", 
                           "mat_ICM", "pat_ICM", "mat_MII_oocyte", "mat_morula", "pat_morula", "sperm")) %>%
  rowid_to_column(var = "ID") #%>%  filter_at(vars(-chr:-end), any_vars(. > 0))

chen_H2Aub_random_mat <- chen_H2Aub_random %>% select(-chr:-end) %>%
  select(ID, A_oocyte = mat_MII_oocyte, B_zygote = zygote_mat, C_2C = mat_2C, D_4C= mat_4C, E_morula = mat_morula, F_ICM = mat_ICM, G_epiblast = mat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_epiblast))

chen_H2Aub_random_pat <- chen_H2Aub_random %>% select(-chr:-end) %>%
  select(ID, A_sperm = sperm, B_zygote = zygote_pat, C_2C = pat_2C, D_4C= pat_4C, E_morula = pat_morula, F_ICM = pat_ICM, G_epiblast = pat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_epiblast))

#pheatmap(log2(chen_H2Aub_random_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H2Aub, Chen et al")
#pheatmap(log2(chen_H2Aub_random_pat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H2Aub, Chen et al")

# H3K27me3 random ---------------------------------------------------------
h3k27me3_random <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/zheng_H3K27me3_random_filt.txt", skip = 1, col_names = 
                       c("chr", "start", "end", "mat_2C", "pat_2C", "mat_8C", "pat_8C", "mat_epiblast", "pat_epiblast", "mat_ICM", "pat_ICM", "oocyte", 
                         "mat_morula", "pat_morula", "mat_zygote", "pat_zygote", "sperm", "sperm_chen"))

h3k27me3_random_mat <- h3k27me3_random %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_oocyte = oocyte, B_zygote = mat_zygote, C_2C = mat_2C, D1_8C = mat_8C, D2_morula = mat_morula, E_ICM = mat_ICM, F_epiblast = mat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(F_epiblast))

h3k27me3_random_pat <- h3k27me3_random %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A1_sperm = sperm, A2_sperm_chen = sperm_chen, B_zygote = pat_zygote, C_2C = pat_2C, D1_8C = pat_8C, D2_morula = pat_morula, E_ICM = pat_ICM, F_epiblast = pat_epiblast) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(F_epiblast))

#pheatmap(log2(h3k27me3_random_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H3K27me3")
#pheatmap(log2(h3k27me3_random_pat[,-c(1, 3)]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H3K27me3")

# H3K4me3 random ----------------------------------------------------------
h3k4me3_random <- read_tsv("~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/zhang_H3K4me3_random_filt.txt", skip = 1, col_names = 
                      c("chr", "start", "end", "mat_early_2C", "pat_early_2C", "mat_late_2C", "pat_late_2C", "mat_4C", "pat_4C", 
                        "mat_8C", "pat_8C", "mat_ICM", "pat_ICM", "oocyte", "sperm", "mat_zygote", "pat_zygote"))

h3k4me3_random_mat <- h3k4me3_random %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_oocyte = oocyte, B_zygote = mat_zygote, C_early_2C = mat_early_2C, D_late_2C = mat_late_2C, E_4C = mat_4C, F_8C = mat_8C, G_ICM = mat_ICM) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_ICM))

h3k4me3_random_pat <- h3k4me3_random %>% 
  rowid_to_column(var = "ID") %>%
  select(ID, A_sperm = sperm, B_zygote = pat_zygote, C_early_2C = pat_early_2C, D_late_2C = pat_late_2C, E_4C = pat_4C, F_8C = pat_8C, G_ICM = pat_ICM) %>%
  #filter_at(vars(-ID), any_vars(. > 1)) %>%
  arrange(desc(G_ICM))

#pheatmap(log2(h3k4me3_random_mat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "maternal H3K4me3")
#pheatmap(log2(h3k4me3_random_pat[,-1]), color = viridis(20), scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
 #        fontsize_row = 1, cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, main = "paternal H3K4me3")

########################################
# Boxplots ----------------------------------------------------------------

# bind maternal
maternal <- bind_rows("H2Aub" = chen_H2Aub_mat, "k27" = h3k27me3_mat, "k4" = h3k4me3_mat, "H2Aub_random" = chen_H2Aub_random_mat, "k27_random" = h3k27me3_random_mat, "k4_random" = h3k4me3_random_mat, .id = "chip")
paternal <- bind_rows("H2Aub" = chen_H2Aub_pat, "k27" = h3k27me3_pat[,-3], "k4" = h3k4me3_pat, "H2Aub_random" = chen_H2Aub_random_pat, "k27_random" = h3k27me3_random_pat, "k4_random" = h3k4me3_random_pat, .id = "chip")

maternal_long <- maternal %>% select(-ID) %>%
  pivot_longer(cols = A_oocyte:G_ICM, names_to = "stage", values_to = "coverage", values_drop_na = TRUE) %>%
  mutate(histone = case_when(grepl("H2Aub", chip) ~ "H2Aub", grepl("k27", chip) ~ "k27", grepl("k4", chip) ~ "k4")) %>%
  mutate(geneset = case_when(grepl("random", chip) ~ "random", TRUE ~ "bivalent"))

paternal_long <- paternal %>% select(-ID) %>%
  pivot_longer(cols = A_sperm:G_ICM, names_to = "stage", values_to = "coverage", values_drop_na = TRUE) %>%
  mutate(histone = case_when(grepl("H2Aub", chip) ~ "H2Aub", grepl("k27", chip) ~ "k27", grepl("k4", chip) ~ "k4"))

# plot by parent
maternal_long %>% ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~histone, scales = "free") +
  scale_fill_manual(values = c("maroon", "lightgrey", "#EFBC68", "lightgrey", "#5F9595", "lightgrey")) +
  #ylim(c(-2.5, 5)) +
  theme_classic() +
  ggtitle("maternal coverage") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

paternal_long %>% ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~histone, scales = "free") +
  scale_fill_manual(values = c("maroon", "lightgrey", "#EFBC68", "lightgrey", "#5F9595", "lightgrey")) +
  theme_classic() + ggtitle("paternal coverage") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# plot by histone
all_data_long <- bind_rows("maternal" = maternal_long, "paternal" = paternal_long, .id = "parent")
plot <- all_data_long %>%
  filter(histone == "H2Aub") %>%
  ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ parent, scales = "free") +
  scale_fill_manual(values = c("maroon", "lightgrey")) +
  scale_x_discrete(name = "stage", labels = c("gamete", "zygote", "2C", "8C", "morula", "ICM", "epiblast")) +
  theme_classic() +
  ylim(c(-3,6)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")

ggsave("~/Dropbox/NRMCB Bivalency Review 2021/figures and tables/H2A_boxplots.pdf", plot, device = "pdf", dpi = 300, width = 8, height = 3, units = "in")

plot <- all_data_long %>%
  filter(histone == "k27") %>%
  ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ parent, scales = "free") +
  scale_fill_manual(values = c("#EFBC68", "lightgrey")) +
  scale_x_discrete(name = "stage", labels = c("gamete", "zygote", "2C", "8C", "morula", "ICM", "epiblast")) +
  theme_classic() +
  ylim(c(-6,4.5)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")
ggsave("~/Dropbox/NRMCB Bivalency Review 2021/figures and tables/K27_boxplots.pdf", plot, device = "pdf", dpi = 300, width = 8, height = 3, units = "in")

plot <- all_data_long %>%
  filter(histone == "k4") %>%
  ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ parent, scales = "free") +
  scale_fill_manual(values = c("#5F9595", "lightgrey")) +
  theme_classic() + ggtitle("H3K4me3") +
  scale_x_discrete(name = "stage", labels = c("gamete", "zygote", "early 2C", "late 2C", "4C", "8C", "ICM")) +
  ylim(c(-2,4)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")
ggsave("~/Dropbox/NRMCB Bivalency Review 2021/figures and tables/K4_boxplots.pdf", plot, device = "pdf", dpi = 300, width = 8, height = 3, units = "in")


# plot together
all_data_long %>%
  ggplot(aes(stage, log2(coverage), fill = chip)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(parent ~ histone, scales = "free") +
  scale_fill_manual(values = c("maroon", "lightgrey", "#EFBC68", "lightgrey", "#5F9595", "lightgrey")) +
  theme_classic() +
  #ylim(c(-5,6)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(color = "black", size = "12"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.ticks.length = unit(.25, "cm"),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside") #+
  #ggsave("boxplots.pdf", device = "pdf", dpi = 300, scale = 1)

########################################
# Boxplot stats -------------------------------------------------------------------
stats <- compare_means(coverage ~ chip, group.by = c("parent", "histone", "stage"), data = all_data_long, method = "t.test")
write_tsv(stats, "~/Dropbox/NRMCB Bivalency Review 2021/TM_analysis/ChIP_stats_ttest.txt")
