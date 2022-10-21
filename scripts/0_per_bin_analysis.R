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

target_cols <- c('index', 'ref', 'alt', 'SYMBOL', 'HGVSc', 'CLNREVSTAT', 'location', 'intron_bin', 'intron_offset', 'which_ss', 'outcome')
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
  
#################
### + 1000bp ####
#################
# Check overlaps with other transcripts of the same gene
plus_1000 <- read_tsv("../data/clinvar/plus_1000/plus_1000_all_consequences.tsv")

check_overlaps <- function(df, data_evaluated){
  conseq_used <- df %>% filter(HGVSc %in% data_evaluated$HGVSc)
  df <- df %>% filter(Gene %in% conseq_used$Gene)
  if (nrow(df) == 1) {
    return(conseq_used %>% add_column(MoreTx = NA))
  }
  else{
    exon <- df %>% filter(!is.na(EXON))
    if (nrow(exon) >= 1){
      return(conseq_used %>% add_column(MoreTx = -1))  
    }
    else{
      df$offset <- unlist(lapply(df$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
      df <-  df %>% slice(which.min(offset))
      return(conseq_used %>% add_column(MoreTx = df$offset))
    }
  }
}

out <- plus_1000 %>% group_by(`#CHROM`, POS, REF, ALT) %>% group_modify(~check_overlaps(.x, data))
out$offset <- unlist(lapply(out$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))

exonic <- out %>% filter(MoreTx == -1)
n_exonic <- nrow(out %>% filter(MoreTx == -1))

intronic <- out %>% filter(is.na(MoreTx) | MoreTx >= 1)

one_conseq <- intronic %>% filter(is.na(MoreTx))
n_one_conseq <- nrow(intronic %>% filter(is.na(MoreTx)))
more_than_1_conseq <- intronic %>% filter(!is.na(MoreTx))
more_than_1_conseq$new_offset <- more_than_1_conseq$offset - more_than_1_conseq$MoreTx

more_than_1_conseq_at_original_offset <- more_than_1_conseq %>% filter(new_offset == 0)
n_more_than_1_conseq_at_original_offset <- nrow(more_than_1_conseq_at_original_offset
                                                )
more_than_1_conseq_at_shorter_offset <- more_than_1_conseq %>% filter(new_offset != 0)
more_than_1_conseq_at_shorter_offset_higher_1000 <- more_than_1_conseq_at_shorter_offset %>% filter(MoreTx >= 1000)
n_more_than_1_conseq_at_shorter_offset_higher_1000 <- nrow(more_than_1_conseq_at_shorter_offset %>% filter(MoreTx >= 1000))

more_than_1_conseq_at_shorter_offset_shorter_1000 <- more_than_1_conseq_at_shorter_offset %>% filter(MoreTx < 1000)
n_more_than_1_conseq_at_shorter_offset_shorter_1000 <- nrow(more_than_1_conseq_at_shorter_offset_shorter_1000)

n_exonic + n_one_conseq + n_more_than_1_conseq_at_original_offset + n_more_than_1_conseq_at_shorter_offset_higher_1000 + n_more_than_1_conseq_at_shorter_offset_shorter_1000
counts <- c(n_exonic, n_one_conseq, n_more_than_1_conseq_at_original_offset, n_more_than_1_conseq_at_shorter_offset_higher_1000, n_more_than_1_conseq_at_shorter_offset_shorter_1000)
names <- c("Exonic", "No other transcript",  "> 1 Tx but same offset", "> 1 Tx smaller offset > 1000bp", "> 1 Tx smaller offset < 1000bp")

to_plot <- tibble(names, counts)
to_plot$names <- factor(to_plot$names, levels = rev(c("No other transcript",
                                                  "> 1 Tx but same offset",
                                                  "> 1 Tx smaller offset > 1000bp",
                                                  "> 1 Tx smaller offset < 1000bp", 
                                                  "Exonic")))
ggplot(to_plot) +
  geom_bar(mapping = aes(x = counts, y = names), stat = 'identity',fill = 'antiquewhite2', alpha=0.7, colour='black') +
  labs(y='', x="Variant counts") + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(more_than_1_conseq_at_shorter_offset_shorter_1000) + 
  geom_histogram(aes(MoreTx), alpha=0.7, fill = 'antiquewhite2', colour = 'black', bins = 50) + 
  theme_bw() + 
  xlab('Distance to splice site') +
  ylab('Counts') +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

to_exclude <- bind_rows(exonic, more_than_1_conseq_at_shorter_offset_shorter_1000)
write.table(to_exclude$HGVSc, file="../data/clinvar/plus_1000/hgvs_to_exclude_exonic_and_shorter_1000.txt", row.names = F, col.names = F, quote = F)
#left_join(to_exclude, data %>% select(HGVSc,outcome), on="HGVSc")

write.table(one_conseq$HGVSc, file="../data/clinvar/plus_1000/hgvs_one_conseq.txt", row.names = F, col.names = F, quote = F)
write.table(more_than_1_conseq_at_original_offset$HGVSc, file="../data/clinvar/plus_1000/hgvs_more_than_1_conseq_same_offset.txt", row.names = F, col.names = F, quote = F)
write.table(exonic$HGVSc, file="../data/clinvar/plus_1000/hgvs_exonic.txt", row.names = F, col.names = F, quote = F)
write.table(more_than_1_conseq_at_shorter_offset_higher_1000$HGVSc, file="../data/clinvar/plus_1000/hgvs_more_than_1_conseq_shorter_offset_higher_1000.txt", row.names = F, col.names = F, quote = F)
write.table(more_than_1_conseq_at_shorter_offset_shorter_1000$HGVSc, file="../data/clinvar/plus_1000/hgvs_more_than_1_conseq_shorter_offset_lower_1000.txt", row.names = F, col.names = F, quote = F)

out_exonic <- read_tsv("../out/TEST_exonic/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_exonic$analysis <- "Exonic"
out_1_consq <- read_tsv("../out/TEST_1_conseq/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_1_consq$analysis <- "No other transcript"
out_more_1_same_offset <- read_tsv("../out/TEST_more_1_conseq_same_offset/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_more_1_same_offset$analysis <- "> 1 Tx but same offset"
out_more_1_shorter_higher_1000 <- read_tsv("../out/TEST_more_1_conseq_shorter_offset_higher_1000/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_more_1_shorter_higher_1000$analysis <- "> 1 Tx smaller offset > 1000bp"
out_more_1_shorter_smaller_1000 <- read_tsv("../out/TEST_more_1_conseq_shorter_offset_lower_1000/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_more_1_shorter_smaller_1000$analysis <- "> 1 Tx smaller offset < 1000bp"

outt <- bind_rows(list(out_exonic, out_1_consq, out_more_1_same_offset, out_more_1_shorter_higher_1000, out_more_1_shorter_smaller_1000)) %>% select(tool, weighted_F1, analysis)
outt$analysis <- factor(outt$analysis, levels = rev(c("No other transcript",
                                                      "> 1 Tx but same offset",
                                                      "> 1 Tx smaller offset > 1000bp",
                                                      "> 1 Tx smaller offset < 1000bp", 
                                                      "Exonic")))

ggplot(outt %>% arrange(-weighted_F1)) +  
  geom_boxplot(aes(y=analysis, x=weighted_F1), alpha = 0.7, fill = 'antiquewhite2') +
  geom_jitter(data= outt %>% filter(weighted_F1 > 0.8), aes(y=analysis, x=weighted_F1, size=0.1, colour=tool), height = 0.3) +
  geom_vline(xintercept = 0.8, linetype='dashed') +
  scale_color_manual(values = brewer.pal(12, "Set3")) +
  theme_bw() + 
  xlab('Weighted F1') +
  ylab('') +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12)) +
  guides(colour = guide_legend(override.aes = list(size=7)))

