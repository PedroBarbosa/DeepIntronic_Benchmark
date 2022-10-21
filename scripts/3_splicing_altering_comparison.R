library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggbump)

setwd(getwd())
list_files <- list(branchpoint = '../out/out_branchpoint/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   pe_new_donor = '../out/out_pe_new_donor/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   pe_donor_associated = '../out/out_pe_donor_associated/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   pe_acceptor = '../out/out_pe_acceptor_associated/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   pe_sre  = '../out/out_pe_change_in_sre/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   elong_new_donor = '../out/out_elong_new_donor/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   elong_donor_associated = '../out/out_elong_donor_associated/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   elong_acceptor = '../out/out_elong_acceptor_associated/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv',
                   elong_sre = '../out/out_elong_change_in_sre/tools_benchmark/all_types/results_tsv/statistics_all_types_all.tsv')

# Concatenate dzata
df <- map2_df(list_files, 
              names(list_files),
              ~read_tsv(.x) %>% 
                mutate(analysis = .y) %>% 
                select(c('tool', 'norm_mcc', 'weighted_norm_mcc', 'weighted_F1', 'total_p', 'total_n',
                         'fraction_nan', 'weighted_accuracy', 'auROC', 'pr_auROC', 'analysis'))) 

df <- df %>% drop_na()

df$analysis <- factor(df$analysis, levels = c("branchpoint", "pe_acceptor", "pe_sre", "pe_new_donor", "pe_donor_associated",
                                              "elong_acceptor", "elong_sre", "elong_new_donor", "elong_donor_associated"))



pe <- df %>% filter(str_detect(analysis, "pe_"))
elong <- df %>% filter(str_detect(analysis, "elong_"))

# Keep only tools that predict all categories
pe_all <- pe %>% 
  group_by(tool) %>% 
  filter(n() >= 4)

elong_all <- elong %>% 
  group_by(tool) %>% 
  filter(n() >= 4)


pal <- c("#B15928", "#377EB8", "#E5D8BD", "#FBB4AE", "#984EA3", "#CCEBC5", "#999999", "#E7298A")

# "#A6CEE3", "#A6761D", "#D95F02", "#FDCDAC", "#A6D854", "#984EA3", "#D9D9D9",
# "#CCEBC5", "#FFFFB3", "#E7298A", "#FFED6F", "#66A61E", "#666666", "#BEAED4", "#E6F5C9"

pe_all$tool = factor(pe_all$tool, levels = df %>% filter(analysis == "pe_acceptor") %>% arrange(-pr_auROC) %>% pull(tool))
labels_pe <- c(paste0("Acceptor associated\nN pos=", 
                      pe_all %>% filter(analysis == "pe_acceptor") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_acceptor") %>% pull(total_n) %>% first()),  
               paste0("SRE\nN pos=", 
                      pe_all %>% filter(analysis == "pe_sre") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_sre") %>% pull(total_n) %>% first()),
               paste0("New splice donor\nN pos=", 
                      pe_all %>% filter(analysis == "pe_new_donor") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_new_donor") %>% pull(total_n) %>% first()),
               paste0("Donor associated\nN pos=", 
                      pe_all %>% filter(analysis == "pe_donor_associated") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_donor_associated") %>% pull(total_n) %>% first()))

ggplot(pe_all, aes(x=analysis, y=pr_auROC, group=tool, color =tool)) +
  geom_bump(smooth = 15, size = 2, alpha = 1) + 
  geom_point(size=5) +
  scale_color_manual(values = pal) +
  labs(title = "Pseudoexon activation", x = '', y = "Average precision") +
  scale_x_discrete(labels= labels_pe) +
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=15, hjust = 0.5),
        legend.text = element_text(size=12)) +
  ylim(0.2, 1) 



elong_all$tool = factor(elong_all$tool, levels = df %>% filter(analysis == "pe_acceptor") %>% arrange(-pr_auROC) %>% pull(tool))
labels_elong <- c(paste0("Acceptor associated\nN pos=", 
                      elong_all %>% filter(analysis == "elong_acceptor") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      elong_all %>% filter(analysis == "elong_acceptor") %>% pull(total_n) %>% first()),  
               paste0("SRE\nN pos=", 
                      elong_all %>% filter(analysis == "elong_sre") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      elong_all %>% filter(analysis == "elong_sre") %>% pull(total_n) %>% first()),
               paste0("New splice donor\nN pos=", 
                      elong_all %>% filter(analysis == "elong_new_donor") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      elong_all %>% filter(analysis == "elong_new_donor") %>% pull(total_n) %>% first()),
               paste0("Donor associated\nN pos=", 
                      elong_all %>% filter(analysis == "elong_donor_associated") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      elong_all %>% filter(analysis == "elong_donor_associated") %>% pull(total_n) %>% first()))

ggplot(elong_all, aes(x=analysis, y=pr_auROC, group=tool, color =tool)) +
  geom_bump(smooth = 15, size = 2, alpha = 1) + 
  geom_point(size=5) +
  scale_color_manual(values = pal) +
  labs(title = "Exon elongation", x = '', y = "Average precision") +
  scale_x_discrete(labels= labels_elong) +
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=15, hjust = 0.5),
        legend.text = element_text(size=12)) +
  ylim(0.2, 1) 

