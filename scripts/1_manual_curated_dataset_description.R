library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggpubr)

setwd("~/git_repos/DeepIntronic_Benchmark/data/splicing_pathogenic_manual_curation/1_dataset_description")
data <- read_tsv('manual_curated.tsv')
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
  geom_bar(mapping = aes(x = n, y = reorder(disease, n)), stat = 'identity',fill = 'antiquewhite2', colour='black') +
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
                       n_occur = sum(aux$n),
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
# data_af <- data %>% mutate(gnomADg_AF = as.numeric(gnomADg_AF)) %>% 
#   mutate(gnomADg_AF = ifelse(is.na(gnomADg_AF), 0, gnomADg_AF))
data_af <- data %>% dplyr::select(gnomADg_AF) %>% mutate(gnomADg_AF = as.numeric(gnomADg_AF)) %>% 
  mutate(gnomADg_AF = ifelse(is.na(gnomADg_AF), 0.00001, gnomADg_AF)) %>% 
  mutate(gnomADg_AF = -log10(gnomADg_AF))
  
ggplot(data_af, aes(gnomADg_AF)) +
  geom_histogram(bins = 100, fill='coral4', alpha=0.8) +
  #geom_density(fill='coral4', alpha=0.7) +
  xlab('-log10 (gnomAD v2.1 frequency)') +
  ylab('Count') +
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
  ylim(0, 150)
  
##################
###  Overlaps  ###
##################
data <- read_tsv("overlaps.tsv")
counts <- data %>% group_by(source) %>% dplyr::summarize(n=n())

our_curation = counts %>% filter(source == "pbarbosa") %>% pull(n)
in_clinvar = counts %>%  filter(source == "in_clinvar" | source == "in_clinvar_wrong_interpt") %>% summarize(sum=sum(n)) %>% pull(sum)
in_clinvar_and_patho_or_likely = in_clinvar - counts %>% filter(source == "in_clinvar_wrong_interpt") %>% pull(n)
in_gnomad = counts %>% filter(source == "in_gnomad") %>% pull(n)


v <- c(our_curation, in_clinvar, in_clinvar_and_patho_or_likely, in_gnomad)
cols <- brewer.pal(n = 4, name = "PuBuGn")
names(v) <- c("Our curation",
             "In ClinVar", 
             "In Clinvar & Pathogenic/Likely_pathogenic",
             "In gnomAD v2.1")
          
barplot(rev(v),
        xlim = c(0, 180),
        col = rev(cols),
        yaxt='n',
        horiz = T)

legend(locator(1), legend = names(rev(v)),fill = rev(cols))

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
