library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(dutchmasters)

setwd("~/git_repos/giga_science_reviews/scripts")

# Manual curation files
curation <- read_tsv('../data/splicing_pathogenic_manual_curation/1_dataset_description/manual_curated.tsv') %>% mutate(`#CHROM` = as.character(`#CHROM`),
                                                                                                                          POS = as.integer(POS))

#####################
###### SPiP #########
#####################
# Correct SPiP predictions
spip <- read_tsv('../data/splicing_pathogenic_manual_curation/2_interpretability/spip/SPiP_data.tsv')


# Join data
spip <- left_join(spip, curation %>% select(HGVSg, 
                                            Effect_category,
                                            Functional_consequence), 
                  by = 'HGVSg')

# Get interpretation and score columns
spip <- spip %>% separate(SPiP, into=c(NA, NA, "Interpretation", "CI", "Score", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "Region", "Extra"), 
                          extra = 'merge', 
                          sep = "[|]") %>% select(-Extra)

spip <- spip %>% filter(!is.na(Effect_category))
spip <- spip %>% mutate(Score=as.numeric(Score))

pseudoexon <- spip %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon$group <- "Pseudoexon activation"
exon_extension <- spip %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                         str_detect(Functional_consequence, 'intron_retention')))
exon_extension$group <- "Partial intron retention"

exon_skipping <- spip %>% filter(Functional_consequence == "exon_skipping")
exon_skipping$group <- "Exon skipping"

spip <- bind_rows(pseudoexon, exon_extension, exon_skipping)
counts_spip <- spip %>% group_by(Interpretation, Effect_category, group) %>% tally()

assign_spip_category <- function(row){
  spip_ <- row[['Interpretation']]
  true_label <- row[['Effect_category']]
  group <- row[['group']]

  if(spip_ == "NTR"){
    return('No interpretation')
  }
  else if(spip_ == "Alter BP" && true_label == "branchpoint_associated"){
    return('Correct')
  }
  else if(spip_ == "Alter BP"){
    return('Incorrect')
  }
  else if(spip_ == "Alter by complex event"){
    return('Not informative')
  }
  else if(spip_ == "Alter by create New Exon" & group == "Pseudoexon activation"){
    return('Correct')
  }
  else if(spip_ == "Alter by create New Exon"){
    return('Incorrect')
  }
  else if(grepl('Alter by create New splice site', spip_) & true_label %in% c("new_splice_acceptor", 
                                                                                  "new_splice_donor",
                                                                                  "strengthening_donor",
                                                                                  "strengthening_acceptor")){
    return('Correct')
  }
  else if(grepl('Alter by create New splice site', spip_)){
    return('Incorrect')
  }
  else if(spip_ == "Alter ESR" & true_label == "change_sre"){
    return('Correct')
  }
  else if(spip_ == "Alter ESR"){
    return('Incorrect')
  }
  else if(grepl('Alter by MES', spip_) & true_label %in% c("strengthening_acceptor", 
                                                          "weakening_acceptor")){
    return('Correct')
  }
  else if(grepl('Alter by MES', spip_)){
    return('Incorrect')
  }
  else{
   return('Never seen') 
  }
  
}

spip$final_interpretation <- apply(spip, 1, assign_spip_category)
# counts_spip$final_interpretation <- apply(counts_spip, 1, assign_spip_category)
# counts_spip <- counts_spip[rep(seq(nrow(counts_spip)), counts_spip$n),]
# counts_spip <- counts_spip %>% select(-n)

#####################
###### SQUIRLS#######
#####################
# Correct SQUIRLS predictions
squirls <- read_tsv('../data/splicing_pathogenic_manual_curation/2_interpretability/squirls/SQUIRLS_data.tsv')


# Join data
squirls <- left_join(squirls, curation %>% select(HGVSg, 
                                            Effect_category,
                                            Functional_consequence), 
                  by = 'HGVSg')

#squirls <- squirls %>% filter(!is.na(Effect_category))

pseudoexon <- squirls %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon$group <- "Pseudoexon activation"
exon_extension <- squirls %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                           str_detect(Functional_consequence, 'intron_retention')))
exon_extension$group <- "Partial intron retention"

exon_skipping <- squirls %>% filter(Functional_consequence == "exon_skipping")
exon_skipping$group <- "Exon skipping"

squirls <- bind_rows(pseudoexon, exon_extension, exon_skipping)

counts_squirls <- squirls %>% group_by(SQUIRLS_interpretation, SQUIRLS_n_bases_info, Effect_category, group) %>% tally()

assign_squirls_category <- function(row){
  squirls_ <- row[['SQUIRLS_interpretation']]
  true_label <- row[['Effect_category']]
  group <- row[['group']]
  
  if(squirls_ == "No explanation"){
    return('No interpretation')
  }
  else if(squirls_ == "Not informative"){
    return('Not informative')
  }
  else if(squirls_ == "New splice donor" && true_label == "new_splice_donor"){
    return('Correct')
  }
  else if(squirls_ == "New splice donor"){
    return('Incorrect')
  }
  else if(squirls_ == "New splice acceptor" & true_label == "new_splice_acceptor"){
    return('Correct')
  }
  else if(squirls_ == "New splice acceptor"){
    return('Incorrect')
  }
  else if(squirls_ == "Activate cryptic donor" & true_label == "strengthening_donor"){
    return('Correct')
  }
  else if(squirls_ == "Activate cryptic donor"){
    return('Incorrect')
  }
  else if(squirls_ == "Activate cryptic acceptor" & true_label == "strengthening_acceptor"){
    return('Correct')
  }
  else if(squirls_ == "Activate cryptic acceptor"){
    return('Incorrect')
  }
  else{
    return('Never seen') 
  }
}

squirls$final_interpretation <- apply(squirls, 1, assign_squirls_category)
# counts_squirls$final_interpretation <- apply(counts_squirls, 1, assign_squirls_category)
# counts_squirls <- counts_squirls[rep(seq(nrow(counts_squirls)), counts_squirls$n),]
# counts_squirls <- counts_squirls %>% select(-n)

##########################
###### SpliceVault #######
##########################
# Possible variants to analyse (those that do not refer to pseudoexon activation outcomes)
splicevault <- read_tsv('../data/splicing_pathogenic_manual_curation/2_interpretability/splicevault/SpliceVault_data.tsv')
splicevault <- splicevault %>% filter(!str_detect(Functional_consequence, 'pseudoexon'))
counts_splicevault <- splicevault %>% group_by(Analysis,Seen_in_SpliceVault, Rank) %>% tally()

assign_splicevault_category <- function(row){
  in_splicevault <- row[['Seen_in_SpliceVault']]
  
  if(str_detect(in_splicevault, "Yes")){
    rank <- row[['Rank']]
    if(rank <= 4){
      return('Correct')
    }
    else if(rank > 4){
      return('Incorrect')
    }
    else{
      return('Never seen')
    }
  }
  else{
    return('No interpretation')
  }
  
}

splicevault$final_interpretation <- apply(splicevault, 1, assign_splicevault_category)
sv_no_explan <- splicevault %>% filter(final_interpretation == "No interpretation") 

###############################
### Plotting tools together ###
##############################
spip_simple <- spip %>% select(c(Score, group, final_interpretation))
spip_simple$tool <- "SPiP"
colnames(spip_simple)[1] <- 'score'
squirls_simple <- squirls %>% select(c(`SQUIRLS (>0.016)`, group, final_interpretation))
squirls_simple$tool <- "SQUIRLS"
colnames(squirls_simple)[1] <- 'score'
splicevault_simple <- splicevault %>% select(final_interpretation)
splicevault_simple$tool <- "SpliceVault"

#f <- df %>%mutate_all(funs(str_replace(., "No interpretation", "Absent")))
df <- bind_rows(spip_simple, squirls_simple, splicevault_simple)

df$tool <- factor(df$tool, levels=c("SQUIRLS", "SPiP", "SpliceVault"))
df %>% ggplot(aes(tool, fill=final_interpretation)) + 
  geom_bar(position="fill", stat="count", colour="black") + 
  scale_fill_dutchmasters(palette = "pearl_earring") +
  ylab('Fraction of interpretations') +
  xlab('') +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        legend.text=element_text(size=10))

df %>% filter(tool != "SpliceVault") %>% ggplot(aes(tool, score, fill=final_interpretation)) + 
  geom_boxplot(alpha=0.9, outlier.shape=NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
             aes(fill = final_interpretation), pch = 21) +
  scale_fill_dutchmasters(palette = "pearl_earring") +
  ylab('Prediction score') +
  xlab('') +
  theme_classic() +
  guides(fill=guide_legend(title="Interpretation")) +
  theme(axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        legend.text=element_text(size=10))

ggplot(sv_no_explan, 
       aes(`Distance_to_cryptic_splice_site(or canonical, for branchpoint variants)`)) +
  geom_histogram(fill = 'coral4', color = 'black', alpha = 0.9) +
  xlab('Distance to cryptic splicing event') +
  ylab('Count') +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 11)) +
  scale_y_continuous(breaks = seq(0, 12, 3), labels = seq(0, 12, 3)) +
  scale_x_continuous(breaks = seq(0, 30, 1), labels = seq(0, 30, 1))
sv_no_explan
