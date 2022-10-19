library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(plyr)
library(ggplot2)
library("RColorBrewer")

# Read and subset data
setwd("~/git_repos/paper_intronic_benchmark/scripts")
data <- read_csv('../data/clinvar/clinvar_to_VETA_evaluation.vcf.gz.tsv')

target_cols <- c('index', 'ref', 'alt', 'SYMBOL', 'CLNREVSTAT', 'location', 'intron_bin', 'intron_offset', 'which_ss', 'outcome')
data <- data %>% dplyr::select(all_of(target_cols))
data$intron_offset_log <- log(data$intron_offset, 10)
data$intron_bin <- factor(data$intron_bin, levels=c("1-2","3-10","11-40", "41-200", "201-500", "501-1000", "1000+"))

# Remove 0 star variants
data <- data %>%filter(CLNREVSTAT != "no_assertion_criteria_provided" &
               CLNREVSTAT != "no_assertion_provided" &
               CLNREVSTAT != "no_interpretation_for_the_single_variant") %>% dplyr::select(-CLNREVSTAT)

data %>% group_by(outcome, intron_bin) %>% tally()
#########################
# All variants per bin ##
#########################
data %>% drop_na(intron_bin) %>%
  ggplot() + 
  geom_density(aes(x=intron_offset, y=..scaled.., fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

#################
# All but 1-2 ###
#################
data %>% drop_na(intron_bin) %>%
  filter(intron_bin != "1-2") %>%
  ggplot() + 
  geom_density(aes(x=intron_offset_log, y=..scaled.., fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  xlab('Log10(Intron offset)') 

########################
#### + 1000bp bin ######
########################
d <- data %>% filter(intron_bin == "1000+") %>% arrange(intron_offset)
d %>% ggplot() + 
  geom_density(aes(x=intron_offset, y=..scaled.., fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
 # xlim(3, 7) + 
  scale_x_continuous(labels = scales::comma) +
  xlab('Log10(Intron offset)')

########################
#### Per ss per bin ####
########################
data %>% drop_na(intron_bin) %>%
  filter(which_ss == "acceptor") %>%
  ggplot() + 
  geom_density(aes(x=intron_offset, y=..scaled.., fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

data %>% drop_na(intron_bin) %>%
  filter(which_ss == "donor") %>%
  ggplot() + 
  geom_density(aes(x=intron_offset, y=..scaled.., fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

##############################
### Nucleotide composition ###
##############################
data$allele <- paste0(data$ref, ">", data$alt)

data %>% drop_na(intron_bin) %>%
  ggplot(aes(x = outcome, fill = allele)) + 
  geom_bar(stat = "count", position = "fill", color = "black",size=0.5) +
  facet_wrap(~intron_bin) +
  scale_fill_brewer(palette = "Set3") +
  ylab('Fraction') +
  xlab('')
  
