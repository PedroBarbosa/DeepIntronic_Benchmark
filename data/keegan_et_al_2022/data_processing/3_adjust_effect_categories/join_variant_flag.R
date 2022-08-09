library(readr)
library(dplyr)
setwd("~/git_repos/paper_intronic_benchmark/data/keegan_et_al_2022/data_processing/3_compare_keegan_with_our_flag/")

latest_processed_data <- read_tsv('1_processed_data.tsv') %>% select(hgvs,mutation_flag)

################
#### keegan  ###
################
keeganTab <- read_tsv('keegan_2022_tabular_original.tsv')
joined <- left_join(keeganTab, latest_processed_data, by =c("ID" = "hgvs")) %>% distinct() 
joined$source <- "keegan_2022"
write.table(joined, file='keegan_2022_tabular.tsv', sep = "\t", quote =  F, row.names = F, col.names = T)

################
### vazDrago ###
################
vazDragoTab <- read_tsv('vazDrago_tabular_original.tsv')
keegan_in_vazDrago <- read_tsv('keegan_in_vazDrago.tsv')

joined <- left_join(vazDragoTab, keegan_in_vazDrago, by="HGVSc")
joined2 <- left_join(joined, latest_processed_data, by =c("ID" = "hgvs")) %>% distinct() %>% select(-c(ID, mutation_flag))
# Dont write now because I manually curated the variants
#joined2 %>% write_tsv(file='vazDrago_tabular.tsv', col_names =  T)

################
### pbarbosa ###
################
pbarbosaTab_fromVCF <- read_tsv('pbarbosa_2022_tabular_original.tsv')
pbarbosa_all_tab <- read_tsv('pbarbosa_deep_intronic.tsv')
merged <- left_join(pbarbosa_all_tab, pbarbosaTab_fromVCF, by=c("HGVS_to_VEP"="ID"))
merged <- merged %>% select(c(colnames(pbarbosaTab_fromVCF)[! colnames(pbarbosaTab_fromVCF) %in% c('ID')],
                              "HGVS_to_VEP", "Effect_category")) %>% write_tsv(file='pbarbosa_2022_tabular.tsv', col_names =  T)

# Check mutation flag of keegan mutations
keegan_in_pbarbosa <- read_tsv('keegan_in_pbarbosa.tsv')
keegan_in_pbarbosa_with_mutation_class <- left_join(keegan_in_pbarbosa, 
                                                    latest_processed_data, by=c('ID'='hgvs')) %>% distinct()
keegan_in_pbarbosa_with_mutation_class %>% count(ID) %>% filter(n>1)

merged_with_keegan <- left_join(merged, 
                                keegan_in_pbarbosa_with_mutation_class, 
                                by="HGVSc") %>% select(c(HGVS_to_VEP, HGVSc, 
                                                         ID,
                                                         Effect_category,
                                                         mutation_flag)) %>% drop_na()
diff <- merged_with_keegan[merged_with_keegan$Effect_category != merged_with_keegan$mutation_flag,]


##############
## Petersen ##
##############
petersenTab <- read_tsv('petersen_2021_tabular_original.tsv')
keegan_in_peterson <- read_tsv('keegan_in_petersen.tsv')
joined <- left_join(petersenTab, keegan_in_peterson, by="HGVSc")
