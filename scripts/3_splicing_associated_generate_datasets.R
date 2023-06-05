library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
setwd("~/git_repos/giga_science_reviews/scripts")

#########################
### Splicing altering ###
#########################
df <- read_tsv('../data/splicing_altering/splicing_variants.tsv') %>% 
  mutate(Chromosome = as.character(Chromosome), 
         Position_hg19 = as.integer(Position_hg19),
         Position_hg38 = as.integer(Position_hg38))

df_positive <- df %>% filter(Affects_splicing)
df_neutral <- df %>% filter(!Affects_splicing)

## PSEUDOEXON-INDUCING ##
pseudoexon <- df_positive %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon_sre <- pseudoexon %>% filter(Analysis == 'Exonic-like')
pseudoexon_new_donor <- pseudoexon %>% filter(Analysis == 'New splice donor' &
                                                Distance_to_cryptic_splice_site %in% c(0, 1))

pseudoexon_new_acceptor <- pseudoexon %>% filter(Analysis == 'Acceptor associated' &
                                                   Distance_to_cryptic_splice_site %in% c(0, 1))

pseudoexon_donor_associated <- pseudoexon %>% filter(Analysis == 'Donor downstream' &
                                                       !Distance_to_cryptic_splice_site %in% c(0, 1))

pseudoexon_acceptor_associated <-  pseudoexon %>% filter(Analysis == 'Acceptor associated' &
                                                           !Distance_to_cryptic_splice_site %in% c(0, 1))

pseudoexon_branchpoint <- pseudoexon %>% filter(Analysis == 'Branchpoint associated')
pseudoexon$offset <- unlist(lapply(pseudoexon$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))

write.table(pseudoexon_sre$Position_hg19, file ='../data/splicing_altering/per_category/pseudoexon_activation/exonic_like/position_positive.txt',
          quote = F,
          row.names = F,
          col.names = F)
write.table(pseudoexon_new_donor$Position_hg19, file ='../data/splicing_altering/per_category/pseudoexon_activation/new_donor/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(pseudoexon_donor_associated$Position_hg19, file ='../data/splicing_altering/per_category/pseudoexon_activation/donor_downstream/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(pseudoexon_new_acceptor$Position_hg19, file ='../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_positive_new_acceptor.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(pseudoexon_acceptor_associated$Position_hg19, file ='../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_positive_acceptor_upstream.txt',
            quote = F,
            row.names = F,
            col.names = F)

## EXON ELONGATION / PARTIAL INTRON RETENTION
# Intron retention Jung == partial intron retention == exon_extension 
exon_extension <- df_positive %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                         str_detect(Functional_consequence, 'intron_retention')))

exon_extension_sre <- exon_extension %>% filter(Analysis == 'Exonic-like')
exon_extension_new_donor <- exon_extension %>% filter(Analysis == 'New splice donor' &
                                                        Distance_to_cryptic_splice_site %in% c(0, 1))
exon_extension_donor_associated <- exon_extension  %>% filter(Analysis == 'Donor downstream' &
                                                                !Distance_to_cryptic_splice_site %in% c(0, 1))

exon_extension_new_acceptor <- exon_extension %>% filter(Analysis == 'Acceptor associated' &
                                                            Distance_to_cryptic_splice_site %in% c(0, 1))

exon_extension_acceptor_associated <- exon_extension %>% filter(Analysis == 'Acceptor associated' &
                                                                  !Distance_to_cryptic_splice_site %in% c(0, 1))

write.table(exon_extension_sre$Position_hg19, file ='../data/splicing_altering/per_category/partial_intron_retention/exonic_like/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(exon_extension_new_donor$Position_hg19, file ='../data/splicing_altering/per_category/partial_intron_retention/new_donor/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(exon_extension_donor_associated$Position_hg19, file ='../data/splicing_altering/per_category/partial_intron_retention/donor_downstream/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(exon_extension_new_acceptor$Position_hg19, file ='../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_positive_new_acceptor.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(exon_extension_acceptor_associated$Position_hg19, file ='../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_positive_acceptor_upstream.txt',
            quote = F,
            row.names = F,
            col.names = F)

## Others (most exon skipping) ##
other <- setdiff(df_positive, bind_rows(pseudoexon %>% select(-offset), exon_extension))
branchpoint_disruption <- other %>% filter(Major_group == "Branchpoint associated")


# Branchpoint-associated (inclues both pseudoexon-activating and exon-skipping|partial intron retention)
branchpoint <- df_positive %>% filter(Major_group == "Branchpoint associated")

write.table(branchpoint$Position_hg19, file ='../data/splicing_altering/per_category/branchpoint_associated/position_positive.txt',
            quote = F,
            row.names = F,
            col.names = F)


skip <- df_positive %>% filter(str_detect(Functional_consequence, "skipping"))
skip %>% group_by(Major_group, Analysis) %>% count()

###################
##### Neutral #####
###################
df_positive %>% group_by(Major_group, Analysis) %>% count()

df_neutral$offset <- unlist(lapply(df_neutral$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))
with_offset <- df_neutral %>% filter(!is.na(offset))
no_offset <- df_neutral %>% filter(is.na(offset))
no_offset$offset <- ifelse(no_offset$Analysis == "Exonic-like", 0, NA)
no_offset[no_offset$HGVSc == "ENST00000358273.4:c.4110+681_4110+693del",]$offset <- 681
no_offset[no_offset$HGVSc == "ENST00000544455.1:c.-39-137_-39-135del",]$offset <- 135
no_offset[no_offset$HGVSc == "ENST00000544455.1:c.632-1136dup",]$offset <- 1136
no_offset[no_offset$HGVSc == "ENST00000544455.1:c.9256+2960_9256+2991dup",]$offset <- 2960 
no_offset[no_offset$HGVSc == "ENST00000471181.2:c.4548-156_4548-155insAGTGAACATT",]$offset <- 155
no_offset[no_offset$HGVSc == "ENST00000471181.2:c.81-1360del",]$offset <- 1360
no_offset[no_offset$HGVSc == "ENST00000355703.3:c.2099+10del",]$offset <- 10
no_offset[no_offset$HGVSc == "ENST00000371953.3:c.209+2257del",]$offset <- 2257
no_offset[no_offset$HGVSc == "ENST00000370225.3:c.768+6457del",]$offset <- 6457
no_offset[no_offset$HGVSc == "ENST00000278616.4:c.8419-54_8419-49del",]$offset <- 49
no_offset[no_offset$HGVSc == "ENST00000544455.1:c.632-943dup" ,]$offset <- 943
no_offset[no_offset$HGVSc == "ENST00000471181.2:c.5257-1687_5257-1676del",]$offset <- 1676
no_offset[no_offset$HGVSc == "ENST00000471181.2:c.4186-2448_4186-2447insGGA",]$offset <- 2447
no_offset[no_offset$HGVSc == "ENST00000233146.2:c.2005+42_2005+44del",]$offset <- 42
no_offset[no_offset$HGVSc == "ENST00000334583.6:c.4037-22_4037-21del",]$offset <- 21
no_offset[no_offset$HGVSc == "ENST00000371953.3:c.802-51_802-14del",]$offset <- 14
no_offset[no_offset$HGVSc == "ENST00000278616.4:c.6453-418del",]$offset <- 418
no_offset[no_offset$HGVSc == "ENST00000408978.4:c.3370+22_3370+23del",]$offset <- 22

df_neutral <- bind_rows(with_offset, no_offset)
  
#################
## Branchpoint ##
#################
branchpoint_neutral <- df_neutral %>% filter(Analysis == "Branchpoint associated")
branchpoint_neutral$Major_group_new <- "Branchpoint associated"
write.table(branchpoint_neutral$Position_hg19, file ='../data/splicing_altering/per_category/branchpoint_associated/position_negative.txt',
           quote = F,
           row.names = F,
           col.names = F)

#################
## Exonic-like ##
#################
# Split neutral datasets to Pseudoexon | Exon elongation
exonic_like_neutral <- df_neutral %>% filter(Analysis == "Exonic-like") %>% dplyr::sample_frac(size=1)
exonic_like_neutral_pe <- head(exonic_like_neutral, 84) 
exonic_like_neutral_pe$Major_group_new <- "Pseudoexon activation"
exonic_like_neutral_partial_ir <- tail(exonic_like_neutral, 35) 
exonic_like_neutral_partial_ir$Major_group_new <- "Partial intron retention"
write.table(exonic_like_neutral_pe %>% pull(Position_hg19), file ='../data/splicing_altering/per_category/pseudoexon_activation/exonic_like/position_negative.txt',
            quote = F,
            row.names = F,
            col.names = F)
write.table(exonic_like_neutral_partial_ir %>% pull(Position_hg19), file ='../data/splicing_altering/per_category/partial_intron_retention/exonic_like/position_negative.txt',
            quote = F,
            row.names = F,
            col.names = F)

other_neutral <- setdiff(df_neutral, bind_rows(branchpoint_neutral %>% select(-Major_group_new), exonic_like_neutral))

###############
## Acceptors ##
###############
# Acceptor associated merged, given small datasets size 
acceptor_neutral <- other_neutral %>% filter(Analysis == "Acceptor associated")
acceptor_neutral_pe <- acceptor_neutral %>% filter(offset >= 100)
acceptor_neutral_pe$Major_group_new <- "Pseudoexon activation"
acceptor_neutral_partial_ir <- acceptor_neutral %>% filter(offset < 100)
acceptor_neutral_partial_ir$Major_group_new <- "Partial intron retention"
write.table(acceptor_neutral_pe %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_negative.txt', quote=F, row.names = F, col.names = F)
write.table(acceptor_neutral_partial_ir %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_negative.txt', quote=F, row.names = F, col.names = F)

############
## Donors ##
############
new_donor_neutral <- other_neutral %>% filter(Analysis == "New splice donor") %>% dplyr::sample_frac(size=1)
new_donor_neutral_pe <- head(new_donor_neutral, 161)
new_donor_neutral_pe$Major_group_new <- "Pseudoexon activation"
new_donor_neutral_partial_ir <- tail(new_donor_neutral, 36)
new_donor_neutral_partial_ir$Major_group_new <- "Partial intron retention"
write.table(new_donor_neutral_pe %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/pseudoexon_activation/new_donor/position_negative.txt', quote=F, row.names = F, col.names = F)
write.table(new_donor_neutral_partial_ir %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/partial_intron_retention/new_donor/position_negative.txt', quote=F, row.names = F, col.names = F)

########################
### Donor downstream ###
########################
donor_downstream_neutral <- other_neutral %>% filter(Analysis == "Donor downstream")
donor_downstream_neutral_pe <- donor_downstream_neutral %>% filter(offset >= 20)
donor_downstream_neutral_pe$Major_group_new <- "Pseudoexon activation"
donor_downstream_neutral_partial_ir <- donor_downstream_neutral %>% filter(offset < 20)
donor_downstream_neutral_partial_ir$Major_group_new <- "Partial intron retention"
write.table(donor_downstream_neutral_pe %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/pseudoexon_activation/donor_downstream/position_negative.txt', quote=F, row.names = F, col.names = F)
write.table(donor_downstream_neutral_partial_ir %>% pull(Position_hg19), file = '../data/splicing_altering/per_category/partial_intron_retention/donor_downstream/position_negative.txt', quote=F, row.names = F, col.names = F)


aux <- bind_rows(branchpoint_neutral, 
                  exonic_like_neutral_pe,
                  exonic_like_neutral_partial_ir,
                  acceptor_neutral_pe,
                  acceptor_neutral_partial_ir,
                  new_donor_neutral_pe,
                  new_donor_neutral_partial_ir,
                  donor_downstream_neutral_pe,
                  donor_downstream_neutral_partial_ir)

df %>% left_join(aux) %>% 
  select(-offset) %>%
  writexl::write_xlsx('~/Desktop/splicing_variants_with_neutral_major_group.xlsx')
