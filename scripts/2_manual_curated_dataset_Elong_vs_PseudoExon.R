library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggpubr)

setwd("~/git_repos/paper_intronic_benchmark/data/manual_curation/")
data_pathogenic <- read_tsv('manual_curated.tsv')
data_benign <- read_tsv('benign.tsv')

###############################
###### Calculate offsets ######
###############################
data_pathogenic$offset <- unlist(lapply(data_pathogenic$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
data_pathogenic$offset[which(data_pathogenic$HGVSc == "ENST00000231188.5:c.1355-587dup")] <- 587
data_pathogenic$offset[which(data_pathogenic$HGVSc == "ENST00000396625.3:c.72-26_72-23del")] <- 23
data_pathogenic$offset[which(data_pathogenic$HGVSc == "ENST00000278616.4:c.2839-579_2839-576del")] <- 576

data_benign$offset <- unlist(lapply(data_benign$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
data_benign$offset[which(data_benign$HGVSc == "ENST00000370225.3:c.2919-575del")] <- 575
data_benign$offset[which(data_benign$HGVSc == "ENST00000280362.3:c.84-715_84-713del")] <- 713 
data_benign$offset[which(data_benign$HGVSc == "ENST00000240651.9:c.415-852del")] <- 852
data_benign$offset[which(data_benign$HGVSc == "ENST00000240651.9:c.415-625dup")] <- 625
data_benign$offset[which(data_benign$HGVSc == "ENST00000552810.1:c.2991+1179del ")] <- 1179
data_benign$offset[which(data_benign$HGVSc == "ENST00000359596.3:c.14647-1630dup")] <- 1630
data_benign$offset[which(data_benign$HGVSc == "ENST00000359596.3:c.14647-1452_14647-1451del")] <- 1451
data_benign$offset[which(data_benign$HGVSc == "ENST00000317799.5:c.443-394dup")] <- 394
data_benign$offset[which(data_benign$HGVSc == "ENST00000270142.6:c.357+89del")] <- 89
data_benign$offset[which(data_benign$HGVSc == "ENST00000510224.1:c.2077-163_2077-162dup ")] <- 162
data_benign$offset[which(data_benign$HGVSc == "ENST00000282516.8:c.4560+2107del")] <- 2107
data_benign$offset[which(data_benign$HGVSc == "ENST00000282516.8:c.4561-2202dup")] <- 2202
data_benign$offset[which(data_benign$HGVSc == "ENST00000231188.5:c.1355-928del")] <- 928
data_benign$offset[which(data_benign$HGVSc == "ENST00000262186.5:c.3330+648_3330+653del")] <- 648
data_benign$offset[which(data_benign$HGVSc == "ENST00000342228.3:c.289+1128del")] <- 1128
data_benign$offset[which(data_benign$HGVSc == "ENST00000375708.3:c.98-1019_98-1018insAT")] <- 1018
data_benign$offset[which(data_benign$HGVSc == "ENST00000467482.1:c.885+584dup")] <- 584
data_benign$offset[which(data_benign$HGVSc == "ENST00000467482.1:c.885+385del")] <- 385
data_benign$offset[which(data_benign$HGVSc == "ENST00000357033.4:c.3603+2253del")] <- 2253
data_benign$offset[which(data_benign$HGVSc == "ENST00000298542.4:c.206-60_206-59del")] <- 59
data_benign$offset[which(data_benign$HGVSc == "ENST00000298556.7:c.485+1626dup")] <- 1626
data_benign$offset[which(data_benign$HGVSc == "ENST00000552810.1:c.2991+1179del")] <- 471
data_benign$offset[which(data_benign$HGVSc == "ENST00000360256.4:c.2113+471_2113+473dup")] <- 471
data_benign$offset[which(data_benign$HGVSc == "ENST00000510224.1:c.2077-163_2077-162dup")] <- 162
data_benign$offset[which(data_benign$HGVSc == "ENST00000342228.3:c.289+1128dup")] <- 1128

data_benign <- data_benign %>% arrange(offset)
##########################################################
#### Split variants: Pseudoexon vs Partial retention #####
##########################################################
pseudoexon <- data_pathogenic %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon$group <- "Pseudoexon activation"

exon_extension <- data_pathogenic %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                           str_detect(Functional_consequence, 'intron_retention')))
exon_extension$group <- "Partial intron retention"

df_pathogenic <- bind_rows(pseudoexon, exon_extension) %>% select(c(Effect_category, Functional_consequence, offset, group))


# Selet benign variants according to the offset
benign_to_exon_elong <- head(data_benign, nrow(exon_extension))
benign_to_exon_elong$group <- "Partial intron retention"
#write.table(file = "pos_benign_exon_elong.txt", x = benign_to_exon_elong$POS, quote=F, row.names = F, col.names = F)

#benign_to_pseudoexon <- tail(data_benign, nrow(data_benign) - nrow(benign_to_exon_elong)) %>%  filter(offset < 10000) %>% sample_n(size=nrow(pseudoexon))
#benign_to_pseudoexon$group <- "Pseudoexon activation"
#write.table(file = "pos_benign_pseudoexon.txt", x =benign_to_pseudoexon$POS, quote=F, row.names = F, col.names = F)

#df_benign <- bind_rows(benign_to_exon_elong, benign_to_pseudoexon) 


# Density plot
df_pathogenic %>% drop_na(offset) %>% 
  ggplot() + 
  geom_density(aes(x=log10(offset), y=..scaled.., fill=group), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  xlab('Log10(Distance to splice site)') +
  ylab('Scaled density') + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.87, 0.25))

# Histogram
df_pathogenic %>% drop_na(offset) %>% 
  ggplot() + 
  geom_histogram(aes(x=offset, fill=group), bins=20, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  xlim(0, 1000) + 
  xlab('Distance to splice site') +
  ylab('Counts') + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.87, 0.75))

# Overall dist 
df_pathogenic %>% drop_na(offset) %>% 
  ggplot() + 
  geom_histogram(aes(x=offset), bins=70, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  xlab('Distance to splice site') +
  ylab('Counts') + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.87, 0.75)) +
  xlim(0,5000)


#######################################
### Compare performance PE vs Elong ###
#######################################
pe <- read_tsv("pseudoexon_activation/statistics_all_types_all.tsv")
pe$analysis <- 'Pseudoexon activation'
retent <- read_tsv("partial_intron_retention/statistics_all_types_all.tsv")
retent$analysis <- 'Partial intron retention'
df <- bind_rows(pe, retent)
df$tool = factor(df$tool, levels = df %>% filter(analysis != "Pseudoexon activation") %>% arrange(-weighted_norm_mcc) %>% pull(tool))

pal <- c("#A6CEE3", "#F1E2CC", "#A6761D", "#FC8D62", "#7FC97F", "#0430B5", "#984EA3", "#999999", "#FFED6F")
# library(hues)
# hues::swatch(pal)

# Keep only tools with metrics in both analysis
paired_tools <- df %>% group_by(tool) %>% tally() %>% filter(n==2) %>% pull(tool)
df <- df[df$tool %in% paired_tools,]

# Remove tools with high fraction of missing data
df <- df %>% filter(fraction_nan < 0.5)

# Sort df by levels so that a paired statistical comparison can be performed
df <- df %>% arrange(factor(tool, levels = levels(df$tool)))
compare_means(weighted_norm_mcc ~ analysis, data = df, paired = TRUE) #, alternative = 'greater')

ggplot(df, aes(x = analysis, y=fraction_nan)) + 
  geom_boxplot(fill = 'antiquewhite2', alpha=0.6) +
  geom_point(aes(fill=tool, group=tool), alpha=0.8, size=3.5, shape=21, position = position_dodge(0.2)) +
  geom_line(aes(group=tool), linetype = "dashed", size = 0.25, position = position_dodge(0.2)) + 
  scale_fill_manual(values =pal) +
  theme_bw() +
  ylim(0, 0.5) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12)) +
  ylab('Fraction of missing predictions') +
  xlab('')


ggplot(df, aes(x = analysis, y=weighted_norm_mcc)) + 
  geom_boxplot(fill = 'antiquewhite2', outlier.shape = NA, alpha=0.6) +
  geom_point(aes(fill=tool, group=tool), alpha=0.8, size=3.5, shape=21, position = position_dodge(0.2)) +
  geom_line(aes(group=tool), linetype = "dashed", size = 0.25, position = position_dodge(0.2)) + 
  stat_compare_means(paired = T, method = 'wilcox.test', label.x = 1.4, label.y = 0.95) +
  scale_fill_manual(values =pal) +
  theme_bw() +
  ylim(0.6,1) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12)) +
  ylab('Weighted normalized MCC') +
  xlab('')

