library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

setwd("~/git_repos/paper_intronic_benchmark/scripts")

#########################
### Splicing altering ###
#########################
list_files <- list(pbarbosa = '../data/splicing_altering/per_study/pbarbosa_and_vazDrago2017/pbarbosa_tabular.tsv',
                   vazDrago = '../data/splicing_altering/per_study/pbarbosa_and_vazDrago2017/vazDrago_tabular.tsv',
                   petersen = '../data/splicing_altering/per_study/petersen_et_al_2022/petersen_2021_tabular.tsv',
                   keegan = '../data/splicing_altering/per_study/keegan_et_al_2022/keegan_2022_tabular.tsv',
                   jung  = '../data/splicing_altering/per_study/jung_et_all_2021/jung_2021_tabular.tsv',
                   moles_fernandez = '../data/splicing_altering/per_study/moles_fernandez_2021/moles_fernandez_2021_splicing_altering_tabular.tsv',
                   tubeuf = '../data/splicing_altering/per_study/tubeuf_et_al_2020/tubeuf_2020_tabular.tsv')

# Concatenate data
df <- map2_df(list_files, 
        names(list_files),
        ~read_tsv(.x) %>% 
          select(-any_of(c('Phenotype (OMIM)', 'Intron_number', 'gnomADg_AF'))) %>%
          mutate(`#CHROM` = as.character(`#CHROM`),
                 POS = as.integer(POS)))

## PSEUDOEXON-INDUCING ##
pseudoexon <- df %>% filter(str_detect(Functional_consequence, 'pseudoexon'))
pseudoexon_sre <- pseudoexon %>% filter(Effect_category == 'change_sre')
pseudoexon_new_donor <- pseudoexon %>% filter(Effect_category == 'new_splice_donor')
pseudoexon_new_acceptor <- pseudoexon %>% filter(Effect_category == 'new_splice_acceptor')
pseudoexon_donor_associated <- pseudoexon %>% filter(Effect_category == 'strengthening_donor')
pseudoexon_acceptor_associated <- pseudoexon %>% filter(Effect_category == 'strengthening_acceptor')
pseudoexon_branchpoint <- pseudoexon %>% filter(Effect_category == 'branchpoint_associated')
pseudoexon$offset <- unlist(lapply(pseudoexon$HGVSc, function(x) as.integer(sub("[A-Z>A-Z ]+", "", sub(".*[+-]", "", x)))))

# write.table(pseudoexon_sre$POS, file ='../data/splicing_altering/per_category/pseudoexon_activation/sre_associated/position_positive.txt', 
#           quote = F,
#           row.names = F,
#           col.names = F)
# write.table(pseudoexon_new_donor$POS, file ='../data/splicing_altering/per_category/pseudoexon_activation/new_donor/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(pseudoexon_donor_associated$POS, file ='../data/splicing_altering/per_category/pseudoexon_activation/donor_associated/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(pseudoexon_new_acceptor$POS, file ='../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_positive_new_acceptor.txt',
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(pseudoexon_acceptor_associated$POS, file ='../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_positive_acceptor_associated.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)

## EXON ELONGATION / PARTIAL INTRON RETENTION
# Intron retention Jung == partial intron retention == exon_extension 
exon_extension <- df %>% filter((!str_detect(Functional_consequence, 'pseudoexon')) & (str_detect(Functional_consequence, 'exon_extension') | 
                                                                                         str_detect(Functional_consequence, 'intron_retention')))

exon_extension_sre <- exon_extension %>% filter(Effect_category == 'change_sre')
exon_extension_new_donor <- exon_extension %>% filter(Effect_category == 'new_splice_donor')
exon_extension_donor_associated <- exon_extension %>% filter(Effect_category == 'strengthening_donor')
exon_extension_new_acceptor <- exon_extension %>% filter(Effect_category == 'new_splice_acceptor')
exon_extension_acceptor_associated <- exon_extension %>% filter(Effect_category == 'strengthening_acceptor')

# write.table(exon_extension_sre$POS, file ='../data/splicing_altering/per_category/partial_intron_retention/sre_associated/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(exon_extension_new_donor$POS, file ='../data/splicing_altering/per_category/partial_intron_retention/new_donor/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(exon_extension_donor_associated$POS, file ='../data/splicing_altering/per_category/partial_intron_retention/donor_associated/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(exon_extension_new_acceptor$POS, file ='../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_positive_new_acceptor.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)
# write.table(exon_extension_acceptor_associated$POS, file ='../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_positive_acceptor_associated.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)

## Others (most exon skipping) ##
other <- setdiff(df, bind_rows(pseudoexon %>% select(-offset), exon_extension))
other <- other %>%
  filter(!if_any(c(Effect_category, Functional_consequence), is.na) &
         ! Effect_category == "-")
branchpoint_disruption <- other %>% filter(Effect_category == "branchpoint_associated")


# Branchpoint-associated (inclues both pseudoexon-activating and exon-skipping|partial intron retention)
branchpoint <- df %>% filter(Effect_category == "branchpoint_associated")

# write.table(branchpoint$POS, file ='../data/splicing_altering/per_category/branchpoint_associated/position_positive.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)


skip <-df %>% filter(str_detect(Functional_consequence, "skipping"))
skip %>% group_by(Effect_category) %>% count()

###################
##### Neutral #####
###################
list_files_neutral <- list(moles_fernandez_neutral = '../data/splicing_altering/per_study/moles_fernandez_2021/moles_fernandez_2021_no_splicing_effect_tabular.tsv',
                           adamson_neutral_exonic = '../data/splicing_altering/per_study/adamson_2018_Vex_seq/adamson_2018_exonic_no_splicing_effect_tabular.tsv',
                           adamson_neutral_intronic = '../data/splicing_altering/per_study/adamson_2018_Vex_seq/adamson_2018_intronic_no_splicing_effect_tabular.tsv',
                           cheung_neutral = '../data/splicing_altering/per_study/cheung_2019_MFASS/cheung_2019_no_splicing_altering_tabular.tsv',
                           new_donor_neutral = '/Users/pbarbosa/git_repos/paper_intronic_benchmark/data/splicing_altering/per_study/neutral_gnomAD/new_donor_neutral_tabular.tsv',
                           new_acceptor_neutral = '/Users/pbarbosa/git_repos/paper_intronic_benchmark/data/splicing_altering/per_study/neutral_gnomAD/new_acceptor_neutral_tabular.tsv')

# Concatenate data
df_neutral <- map2_df(list_files_neutral, 
              names(list_files_neutral),
              ~read_tsv(.x) %>% 
                select(-any_of(c('Intron_number', 'Exon_number', 'gnomADg_AF'))) %>%
                mutate(`#CHROM` = as.character(`#CHROM`),
                       POS = as.integer(POS)))

branchpoint_neutral <- df_neutral %>% filter(Effect_category == "branchpoint_associated")
# write.table(branchpoint_neutral$POS, file ='../data/splicing_altering/per_category/branchpoint_associated/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F)

# Split neutral datasets to Pseudoexon | Exon elongation
sre_associated <- df_neutral %>% filter(Effect_category == "change_sre")  
sre_associated <- sre_associated %>% dplyr::sample_frac(size=1)

# write.table(head(sre_associated, 90) %>% pull(POS), file ='../data/splicing_altering/per_category/pseudoexon_activation/sre_associated/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 
# write.table(tail(sre_associated, 33) %>% pull(POS), file ='../data/splicing_altering/per_category/partial_intron_retention/sre_associated/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 

other_neutral <- setdiff(df_neutral, bind_rows(branchpoint_neutral, sre_associated))
other_neutral$offset <- unlist(lapply(other_neutral$HGVSc, function(x) as.integer(str_extract(strsplit(x, split="[+-]") %>% sapply( "[", 2 ), '\\d+'))))

# Acceptor associated merged, given small datasets size 
other_neutral_pseudoexon <- other_neutral %>% filter(Effect_category == "new_splice_acceptor" | (Effect_category == "closer_to_acceptor" & offset >= 150))
other_neutral_partial_IR <- other_neutral %>% filter(Effect_category == "strengthening_acceptor" | (Effect_category == "closer_to_acceptor" & offset < 150)) 
# write.table(other_neutral_pseudoexon %>% pull(POS), file = '../data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/position_negative.txt', quote=F, row.names = F, col.names = F)
# write.table(other_neutral_partial_IR %>% pull(POS), file = '../data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/position_negative.txt', quote=F, row.names = F, col.names = F)

############
## Donors ##
############
new_donor <- other_neutral %>% filter(Effect_category == "new_splice_donor") %>% dplyr::sample_frac(size=1)
# PE
pseudoexon_donor_associated <- other_neutral %>% filter((Effect_category == "strengthening_donor" & offset >= 20) | (Effect_category == "closer_to_donor" & offset >= 20))
pseudoexon_new_donor <- head(new_donor, 171)
# write.table(pseudoexon_donor_associated %>% pull(POS), file ='../data/splicing_altering/per_category/pseudoexon_activation/donor_associated/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 
# 
# write.table(pseudoexon_new_donor %>% pull(POS), file ='../data/splicing_altering/per_category/pseudoexon_activation/new_donor/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 

# Exon elongation
elongation_donor_associated <- other_neutral %>% filter((Effect_category == "strengthening_donor" & offset < 20) | (Effect_category == "closer_to_donor" & offset < 20))
elongation_new_donor <- tail(new_donor, 40)
# write.table(elongation_donor_associated %>% pull(POS), file ='../data/splicing_altering/per_category/partial_intron_retention/donor_associated/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 
# 
# write.table(elongation_new_donor %>% pull(POS), file ='../data/splicing_altering/per_category/partial_intron_retention/new_donor/position_negative.txt', 
#             quote = F,
#             row.names = F,
#             col.names = F) 

