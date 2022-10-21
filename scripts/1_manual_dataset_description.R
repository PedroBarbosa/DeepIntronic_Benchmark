library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

setwd(getwd())
data <- read_tsv('../data/manual_curation/manual_curated.tsv')

data <- data %>% filter(Source == "pbarbosa_2022")

###################
# Disease barplot #
###################
data <- data %>% separate(., col=`Phenotype (OMIM)`,into = "disease", sep = "\\(", remove = F)
counts <- data %>% dplyr::select(disease) %>% group_by(disease) %>% tally()
counts <- counts %>% mutate(disease_filt = case_when(n >= 3 ~ disease,
                                                     n < 3 ~ "Other"))

# Diseases with more than 3 variants
no_other <- counts %>% filter(disease_filt != "Other")
ggplot(no_other) +
  geom_bar(mapping = aes(x = n_occur, y = reorder(disease, n_occur)), stat = 'identity',fill = 'antiquewhite2', colour='black') +
  labs(y='', x="Variant counts") + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
 xlim(0, 20)

# All diseases
aux <- filter(counts, disease_filt == "Other")
other_df <- tibble_row(disease = c("Other"),
                       n_occur = sum(aux$n_occur),
                       disease_filt = c("Other"))
all <- bind_rows(no_other, other_df)

ggplot(all) +
  geom_bar(mapping = aes(x = n_occur, y = reorder(disease, n_occur)), stat = 'identity',fill = 'lightblue', colour='black') +
  labs(y='', x="Variant counts") + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlim(0, 100)

###################
## AF histogram ###
###################
data <- data %>% mutate(gnomADg_AF = as.numeric(gnomADg_AF)) %>% 
  mutate(gnomADg_AF = ifelse(is.na(gnomADg_AF), 0, gnomADg_AF))

ggplot(data, aes(gnomADg_AF)) +
  geom_histogram(bins = 200, fill='coral4', alpha=0.8) +
  #geom_density(fill='coral4', alpha=0.7) +
  xlab('gnomAD v2.1 frequency') +
  ylab('Count')
  
##################
###  Overlaps  ###
##################
data <- read_tsv("../data/manual_curation/overlaps.tsv")
counts <- data %>% group_by(source) %>% summarize(n=n())

our_curation = counts %>% filter(source == "pbarbosa") %>% pull(n)
in_clinvar = counts %>% filter(source == "in_clinvar") %>% pull(n)
in_clinvar_and_patho_or_likely = in_clinvar - counts %>% filter(source == "in_clinvar_wrong_interpt") %>% pull(n)
in_gnomad = counts %>% filter(source == "in_gnomad") %>% pull(n)


v <- c(our_curation, in_clinvar, in_clinvar_and_patho_or_likely, in_gnomad)
cols <- brewer.pal(n = 4, name = "PuBuGn")
names(v) <- c("Our curation",
             "In ClinVar", 
             "In Clinvar & Pathogenic/Likely_pathogenic",
             "In gnomAD v2.1")
          
barplot(v,
        ylim = c(0, 180),
        col = cols,
        xaxt='n')
legend(locator(1), legend = names(v),fill = cols)

# library(VennDiagram)
# flog.threshold(ERROR)
#Make the plot
# venn.diagram(
#   x = list(data %>% filter(source=="pbarbosa") %>% unlist(),
#            data %>% filter(source=="in_clinvar") %>% unlist(),
#            data %>% filter(source=="in_clinvar_wrong_interpt") %>% unlist(),
#            data %>% filter(source=="in_gnomad") %>% unlist()),
#   category.names = c("Our curation" , "In ClinVar" , 
#                      "In ClinVar with wrong interpretations",
#                      "In gnomAD v2.1"),
#   filename='~/Desktop/venn.png')
 

###############################
###### Bin distribution #######
###############################
data$offset <- unlist(lapply(data$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
data$offset[which(data$HGVSc == "ENST00000231188.5:c.1355-587dup")] <- 587
data$offset[which(data$HGVSc == "ENST00000396625.3:c.72-26_72-23del")] <- 23
data$offset[which(data$HGVSc == "ENST00000278616.4:c.2839-579_2839-576del")] <- 576

pseudoexon <- data %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon$group <- "pseudoexon"
exon_extension <- data %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                         str_detect(Functional_consequence, 'intron_retention')))
exon_extension$group <- "exon_extension"
write.table(pseudoexon$HGVSg, file ='../data/manual_curation/pseudoexon/pos.txt', 
            quote = F,
            row.names = F,
            col.names = F)

df <- bind_rows(pseudoexon, exon_extension) %>% select(c(Effect_category, Functional_consequence, offset, group))
# Density plot
df %>% drop_na(offset) %>% 
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
df %>% drop_na(offset) %>% 
  ggplot() + 
  geom_histogram(aes(x=offset, fill=group), bins=30, alpha = 0.8, position = "identity") +
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
df %>% drop_na(offset) %>% 
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
pe <- read_tsv("../data/manual_curation/pe_activation/statistics_all_types_all.tsv")
pe$analysis <- 'Pseudoexon activation'
elong <- read_tsv("../data/manual_curation/exon_elongation/statistics_all_types_all.tsv")
elong$analysis <- 'Exon elongation'
df <- bind_rows(pe, elong)

library(RColorBrewer)
n <- 25
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
palette=sample(col_vector, n)

pal <- c("#A6CEE3", "#F1E2CC", "#A6761D", "#FC8D62", "#7FC97F", "#D95F02", "#FDCDAC", "#FFFF33", "#A6D854", "#984EA3", "#999999", "#D9D9D9",
         "#B15928", "#377EB8", "#E5D8BD", "#FBB4AE", "#CCEBC5", "#FFFFB3", "#E7298A", "#FFED6F", "#66A61E", "#666666", "#BEAED4", "#E6F5C9",
         "#FFFF99")

df$tool = factor(df$tool, levels = df %>% filter(analysis != "Pseudoexon activation") %>% arrange(-weighted_norm_mcc) %>% pull(tool))
ggplot(df, aes(x = analysis, y=weighted_norm_mcc)) +
  geom_boxplot(fill = 'antiquewhite2', alpha=0.8) +
#  geom_point(aes(color = tool, size =2) , size = 2) +
  geom_point(aes(fill=tool, group=tool),size=5, shape=21, position = position_dodge(0.2)) +
  geom_line(aes(group=tool), position = position_dodge(0.2)) + 
  scale_fill_manual(values =palette) +
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
        legend.text = element_text(size=12)) +
  ylab('Weighted normalized MCC') +
  xlab('')
