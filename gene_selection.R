# script to select mm10 TSS controls

library(tidyverse)
library(readxl)

setwd("~/working_dir/")
set.seed(0)

# read in files -----------------------------------------------------------

bivalent <- read_tsv("mas_bivalent_genes.txt", col_names = c("GeneID")) %>%
  unique()
nrow(bivalent)

mm10 <- read_tsv("~/Documents/bioinf/mm10_tss_merge.bed", col_names = c("chr", "start", "end", "GeneID")) %>%
  group_by(GeneID) %>%
  slice_sample(n = 1) %>%
  ungroup()

mm10_bivalent <- mm10 %>%
  filter(GeneID %in% bivalent$GeneID)

mm10_non_bivalent <- mm10 %>%
  filter(!GeneID %in% bivalent$GeneID)

# select a random gene sample
#random_sample <- mm10 %>%
#  filter(!GeneID %in% bivalent$GeneID) %>%
#  slice_sample(n = 3868)

#random_sample %>% 
#  select(GeneID) %>%
#  write_tsv("mm10_non_bivalent_sample_3868.txt", col_names = FALSE)

random_sample_genes <- read_tsv("mm10_non_bivalent_sample_3868.txt", col_names = "GeneID")
random_sample <- mm10 %>%
  filter(GeneID %in% random_sample_genes$GeneID)

random_sample %>% count(GeneID) %>% filter(n > 1)

# compare by CpG content --------------------------------------------------
# read in data from Mikkelsen et al 2007
genes_cpg <- read_excel("Mikkelsen_CpG_content.xlsx", 
                             col_names = c("TSS", "strand", "GeneID", "RefSeq", "class", "ESC", "NPC", "MEF")) %>%
  select(GeneID, class)

#293 genes have >1 promoter listed in this table
nrow(genes_cpg %>% count(GeneID) %>% filter(n > 1))

high_cpg_genes <- genes_cpg %>%
  filter(class == "HCP") %>%
  select(GeneID)

# annotate above subsets by CpG content
bivalent_cpg <- mm10_bivalent %>%
  left_join(genes_cpg)
bivalent_cpg %>% group_by(class) %>% tally()
#HCP 2535 / ICP 274 / LCP 128 / NA 954

random_sample_cpg <- random_sample %>%
  left_join(genes_cpg)
random_sample_cpg %>% group_by(class) %>% tally()
#HCP 876 / ICP 333 / LCP 300 / NA 2383

# define control based on matched CpG content
 high_cpg_2535 <- mm10_non_bivalent %>%
  filter(GeneID %in% high_cpg_genes$GeneID) %>%
  slice_sample(n = 2535) 

matched_cpg_sample <- mm10_non_bivalent %>%
  filter(!GeneID %in% high_cpg_2535$GeneID) %>%
  slice_sample(n = 3868-2535) %>%
  bind_rows(high_cpg_2535) %>%
  left_join(genes_cpg)


# write files -------------------------------------------------------------

#write_tsv(mm10_bivalent, "~/mas/bivalent_3868.bed", col_names = FALSE)
#matched_cpg_sample %>% select(-class) %>% write_tsv("mas/matched_cpg_sample.bed", col_names = FALSE)
#write_tsv(random_sample, "mas/random_sample.bed", col_names = FALSE)

# plots of cpg content ----------------------------------------------------
all_subsets <- bind_rows("matched_ctrl" = matched_cpg_sample, 
                         "random_ctrl" = random_sample_cpg, 
                         "bivalent" = bivalent_cpg, 
                         .id = "geneset") %>%
  mutate(class = replace_na(class, "other")) 

all_subsets %>%
  group_by(geneset, class) %>%
  tally() %>%
  mutate(class = factor(class, levels = c("other", "LCP", "ICP", "HCP"))) %>%
  ggplot(aes(x = geneset, y = n, fill = class)) + 
  geom_bar(stat = "identity", width = 0.8, color = "darkgrey", size = 0.2) + 
  labs("Mikkelsen class") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4000)) +
  coord_flip() +
  ggtitle("CpG island density") + ylab("Number of genes") +
  scale_x_discrete(name = "Gene set", 
                   labels = c("Bivalent", "CpG-matched", "Random")) +
  scale_fill_manual(values = c("#eaeaea", "#fae2a2","#ffc214","#d95f0e")) +
  theme(axis.line = element_line(size = 0.5, colour = "black", linetype = 1),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black", size = 14), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.ticks.length = unit(.25, "cm"), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "transparent", size = 0)) 
#ggsave("~/figures/cpg_density.pdf", device = "pdf", dpi = 300)


# venn diagram ------------------------------------------------------------
library(VennDiagram)
venn.diagram(x = list(mm10_bivalent$GeneID, random_sample$GeneID, matched_cpg_sample$GeneID), 
             category.names = c("bivalent", "random", "matched"), 
             filename = "~/Dropbox/NRMCB Bivalency Review 2021/THIRD SUBMISSION/figures/venn_diagram.png", 
             output = TRUE)
