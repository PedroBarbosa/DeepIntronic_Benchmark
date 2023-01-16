library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(plyr)
library(ggplot2)
library("RColorBrewer")

# Read and subset data
setwd("~/git_repos/paper_intronic_benchmark/")
data <- read_csv('data/clinvar/clinvar_to_VETA_evaluation.vcf.gz.tsv')

target_cols <- c('index','chr', 'pos', 'ref', 'alt', 'SYMBOL', 'HGVSc', 'CLNREVSTAT', 'location', 'intron_bin', 'intron_offset', 'which_ss', 'outcome')
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
  geom_density(aes(x=intron_offset, y=after_stat(scaled), fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

#################
# All but 1-2 ###
#################
data %>% drop_na(intron_bin) %>%
  filter(intron_bin != "1-2") %>%
  ggplot() + 
  geom_density(aes(x=intron_offset_log, y=after_stat(scaled), fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  xlab('Log10(Intron offset)') 

########################
#### + 1000bp bin ######
########################
d <- data %>% filter(intron_bin == "1000+") %>% arrange(intron_offset)
d %>% ggplot() + 
  geom_density(aes(x=intron_offset, y=after_stat(scaled), fill=outcome), alpha = 0.8, position = "identity") +
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
  geom_density(aes(x=intron_offset, y=after_stat(scaled), fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

data %>% drop_na(intron_bin) %>%
  filter(which_ss == "donor") %>%
  ggplot() + 
  geom_density(aes(x=intron_offset, y=after_stat(scaled), fill=outcome), alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(vars(intron_bin), scales='free_x')

##############################
### Nucleotide composition ###
##############################
# Need to fix alleles based on strand info
# data$allele <- paste0(data$ref, ">", data$alt)
# 
# data %>% drop_na(intron_bin) %>%
#   ggplot(aes(x = outcome, fill = allele)) + 
#   geom_bar(stat = "count", position = "fill", color = "black",size=0.5) +
#   facet_wrap(~intron_bin) +
#   scale_fill_brewer(palette = "Set3") +
#   ylab('Fraction') +
#   xlab('')
  
#############################
### 501-1000 and +1000bp ####
##############################
# Check overlaps with other transcripts of the same gene
plus_500 <- read_tsv("data/clinvar/2_plus_500/plus_500_all_consequences.tsv")

check_overlaps <- function(df, data_evaluated){
  conseq_used <- df %>% filter(HGVSc %in% data_evaluated$HGVSc)
  df <- df %>% filter(Gene %in% conseq_used$Gene)

  # No other tx, same conseq returned
  if (nrow(df) == 1) {
    return(conseq_used %>% add_column(MoreTx = NA))
  }
  else{
    # Other tx, at least 1 exonic
    exon <- df %>% filter(!is.na(EXON))
    if (nrow(exon) >= 1){
      return(conseq_used %>% add_column(MoreTx = -1))  
    }
    
    else{
      
      # Other intronic, check for smaller offset
      intron <- df %>% filter(!is.na(INTRON))
      if (nrow(intron) >= 1){
        intron$offset <- unlist(lapply(intron$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
        intron <- intron %>% slice(which.min(offset))
        return(conseq_used %>% add_column(MoreTx = intron$offset))
      }
      else{
        show("Problem here. Only upstream/downstream consequences retrieved.")
      }
    }
  }
}

out <- plus_500 %>% group_by(`#CHROM`, POS, REF, ALT) %>% group_modify(~check_overlaps(.x, data))
out <- left_join(out, data %>% select(chr, pos, ref, alt, intron_bin, outcome), by=c("#CHROM" = "chr", "POS" = "pos", "REF" = "ref", "ALT" = "alt"))
out$offset <- unlist(lapply(out$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))

exonic <- out %>% filter(MoreTx == -1)
n_exonic <- nrow(out %>% filter(MoreTx == -1))

intronic <- out %>% filter(is.na(MoreTx) | MoreTx >= 1)

one_conseq <- intronic %>% filter(is.na(MoreTx))
n_one_conseq <- nrow(one_conseq)

more_than_1_conseq <- intronic %>% filter(!is.na(MoreTx))
more_than_1_conseq$new_offset_diff <- more_than_1_conseq$offset - more_than_1_conseq$MoreTx

more_than_1_conseq_at_original_offset <- more_than_1_conseq %>% filter(new_offset_diff == 0)
n_more_than_1_conseq_at_original_offset <- nrow(more_than_1_conseq_at_original_offset)
                                                
more_than_1_conseq_at_shorter_offset <- more_than_1_conseq %>% filter(new_offset_diff != 0)
n_more_than_1_conseq_at_shorter_offset <- nrow(more_than_1_conseq_at_shorter_offset)

#more_than_1_conseq_at_shorter_offset_higher_500 <- more_than_1_conseq_at_shorter_offset %>% filter(MoreTx >= 500)
#n_more_than_1_conseq_at_shorter_offset_higher_500 <- nrow(more_than_1_conseq_at_shorter_offset %>% filter(MoreTx >= 500))

#more_than_1_conseq_at_shorter_offset_shorter_500 <- more_than_1_conseq_at_shorter_offset %>% filter(MoreTx < 500)
#n_more_than_1_conseq_at_shorter_offset_shorter_500 <- nrow(more_than_1_conseq_at_shorter_offset_shorter_500)

n_exonic + n_one_conseq + n_more_than_1_conseq_at_original_offset + n_more_than_1_conseq_at_shorter_offset
counts <- c(n_exonic, n_one_conseq, n_more_than_1_conseq_at_original_offset, n_more_than_1_conseq_at_shorter_offset)
names <- c("Exonic", "No other transcript",  "> 1 transcript (same offset)", "> 1 transcript (smaller offset)")

to_plot <- tibble(names, counts)
to_plot$names <- factor(to_plot$names, levels = rev(c("No other transcript",
                                                  "> 1 transcript (same offset)",
                                                  "> 1 transcript (smaller offset)",
                                                  "Exonic")))


ggplot(to_plot, aes(x=counts, y=names)) +
  geom_bar(stat = 'identity',fill = 'antiquewhite2', alpha=0.7, colour='black') +
  geom_text(aes(label=counts), size = 4, position=position_dodge(width=0.9), hjust=-0.25) + 
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

# Limit the range to 1000bp for clear depiction
# 41 1000+ variants are still > 1000+; 5 501-1000 are still > 500
gg1 <- ggplot() + 
  geom_histogram(data = more_than_1_conseq_at_shorter_offset %>% filter(intron_bin == "1000+"), aes(x = MoreTx), bins = 100, color = 'black', fill="antiquewhite2") +
  geom_label( aes(x =700, y = 20, label="Variants originally assigned to '1000+' bin"), size = 6, color="black") + 
  xlim(0, 1000) + 
  theme_bw() +
  xlab('Distance to splice site') +
  ylab('Counts') +
  xlim(0,1000) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

gg2 <- ggplot() +  geom_histogram(data = more_than_1_conseq_at_shorter_offset %>% filter(intron_bin == "501-1000"), aes(x = MoreTx, fill=''), bins =100, color = 'black', fill= "#404080") +
  geom_label( aes(x =700, y = 4, label="Variants originally assigned to '501-1000' bin"), size = 6, color="black") + 
  xlim(0, 1000) +
  scale_y_reverse() + 
  theme_bw() +
  xlab('Distance to splice site') +
  ylab('Counts') +
  xlim(0,1000) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggpubr::ggarrange(gg1, gg2,
          ncol = 1, nrow = 2)


####################################################
######### Write variants per category to file ######
####################################################
# To run VETA on each one separately

write.table(one_conseq$HGVSc, file="data/clinvar/plus_500/hgvs_one_conseq.txt", row.names = F, col.names = F, quote = F)
write.table(more_than_1_conseq_at_original_offset$HGVSc, file="data/clinvar/plus_500/hgvs_more_than_1_conseq_same_offset.txt", row.names = F, col.names = F, quote = F)
write.table(exonic$HGVSc, file="data/clinvar/plus_500/hgvs_exonic.txt", row.names = F, col.names = F, quote = F)
write.table(more_than_1_conseq_at_shorter_offset$HGVSc, file="data/clinvar/plus_500/hgvs_more_than_1_conseq_smaller_offset.txt", row.names = F, col.names = F, quote = F)


############################################
########### Upload VETA results ############
############################################

out_exonic <- read_tsv("out/out_clinvar_plus_500_exonic/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_exonic$analysis <- "Exonic"

out_1_conseq <- read_tsv("out/out_clinvar_plus_500_intronic_1_conseq/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_1_conseq$analysis <- "No other transcript"

out_more_1_same_offset <- read_tsv("out/out_clinvar_plus_500_intronic_moreThan1Tx_same_offset/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_more_1_same_offset$analysis <- "> 1 transcript (same offset)"

out_more_1_smaller_offset <- read_tsv("out/out_clinvar_plus_500_intronic_moreThan1Tx_smaller_offset/1s_l/tools_benchmark/snps/results_tsv/statistics_snps_all.tsv")
out_more_1_smaller_offset$analysis <- "> 1 transcript (smaller offset)"

outt <- bind_rows(list(out_exonic, out_1_conseq, out_more_1_same_offset, out_more_1_smaller_offset)) %>% select(tool, weighted_F1, average_precision, fraction_nan, analysis)
outt$analysis <- factor(outt$analysis, levels = rev(c("No other transcript",
                                                      "> 1 transcript (same offset)",
                                                      "> 1 transcript (smaller offset)", 
                                                      "Exonic")))

# Remove splicing tools with range limits 
outt <- outt %>% filter(!tool %in% c("HAL", "SPIDEX", "dbscSNV", "S-CAP", "kipoiSplice4", "MLCsplice", "MMSplice", "MaxEntScan"))

ggplot() +  
  geom_boxplot(data = outt, aes(y=analysis, x=weighted_F1), outlier.shape = NA, alpha = 0.7, fill = 'antiquewhite2') +
  geom_point(data= outt %>% filter(weighted_F1 < 0.65), 
             position = position_dodge(width=1), 
             size = 1,
             aes(y=analysis, x=weighted_F1)) +
  
  geom_point(data= outt %>% filter(weighted_F1 >= 0.65), 
             size = 4,
             position = position_jitter(width = 0, height = 0.2), 
             aes(y=analysis, x=weighted_F1, colour=tool)) +
  geom_vline(xintercept = 0.65, linetype='dashed') +
  scale_color_manual(values = brewer.pal(12, "Set3")) +
  theme_bw() + 
  xlim(0, 1) +
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
  guides(colour = guide_legend(override.aes = list(size=6)))


###################################################
### Compare last two bins with previous results ###
###################################################
out_500_old <- read_tsv('out/out_clinvar_all/1s_l/intron_analysis/results_tsv/statistics_all_501-1000.tsv') %>% select(tool, weighted_F1) %>% dplyr::rename(old = weighted_F1)
out_1000_old <- read_tsv('out/out_clinvar_all/1s_l/intron_analysis/results_tsv/statistics_all_1000+.tsv') %>% select(tool, weighted_F1) %>% dplyr::rename(old = weighted_F1)
out_500_new <- read_tsv('out/out_clinvar_excluding_exonic_and_closer_to_ss_all_tools/1s_l/intron_analysis/results_tsv/statistics_all_501-1000.tsv') %>% select(tool, weighted_F1) %>% dplyr::rename(new = weighted_F1)
out_1000_new <- read_tsv('out/out_clinvar_excluding_exonic_and_closer_to_ss_all_tools/1s_l/intron_analysis/results_tsv/statistics_all_1000+.tsv') %>% select(tool, weighted_F1) %>% dplyr::rename(new = weighted_F1)

out_500 <- left_join(out_500_old, out_500_new)
out_500$bin <- "500-1000"

out_1000 <- left_join(out_1000_old, out_1000_new)
out_1000$bin <- "1000+"


splicing_tools_to_remove <- c("HAL", "SPIDEX", "dbscSNV", "S-CAP", "kipoiSplice4", "MLCsplice", "MMSplice", "MaxEntScan")
to_plot <- bind_rows(out_500, out_1000) %>% filter(!tool %in% splicing_tools_to_remove)
to_plot$difference <- to_plot$new - to_plot$old
to_plot <- to_plot %>% arrange(-difference)

to_plot$bin = factor(to_plot$bin, levels=c("500-1000","1000+")) 

ggplot(to_plot) +
  geom_point(aes(x = difference, y = reorder(tool, difference), fill = "antiquewhite2"), size = 3, pch=21, color = 'black') +
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.75) +
  scale_fill_manual(values=c("antiquewhite2")) +
  facet_wrap(~bin) +
  geom_text(data=subset(to_plot, difference > 0), aes(x=difference, y=reorder(tool, difference), label=new), hjust=-0.5, cex=4, col='black') +
  xlim(c(-0.6,0.4)) +
  ylab('') +
  xlab('Performance difference (measured by Weighted F1)') +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13),
        strip.text.x = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

