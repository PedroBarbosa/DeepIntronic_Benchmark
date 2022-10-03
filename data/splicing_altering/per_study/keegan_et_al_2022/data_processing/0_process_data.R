library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
setwd("~/git_repos/paper_intronic_benchmark/data/deep_intronic/keegan_et_al_2022/data_processing/")

data <- read_xlsx('0_supp_tables.xlsx', sheet = 1, skip= 1) %>% select(c(ID,
                                                                         Gene,
                                                                         Intron,
                                                                         Size,
                                                                         `Transcript ID`, 
                                                                         `Instigating mutation(s)`,
                                                                         `Mutation type`,
                                                                         `Mutation type updated`,
                                                                         `RNA source`))
names(data) <- tolower(names(data))
data <- data %>% rename(tx_id = `transcript id`, 
                        mutation = `instigating mutation(s)`,
                        mutation_type = `mutation type`,
                        mutation_type_updated = `mutation type updated`,
                        rna_source = `rna source`)

# Melt variants triggering the same pseudoexons
data <- data %>% mutate(mutation = strsplit(mutation, "; ")) %>% unnest(mutation)

# Remove variants with rsIDs with low condifence
rsIDs_to_remove <- c("rs2014886", "rs2545162", "rs2075356")
data <- filter(data, !grepl(paste(rsIDs_to_remove, collapse="|"), mutation))

# Remove spaces from mutation column
data$mutation <- sapply(strsplit(data$mutation, "\\s+"), "[", 1)

# Select rows where mutation starts with c.*
data <- data %>% filter(str_detect(mutation, "^c."))

# Generate HGVSc
data$hgvs <- paste0(data$tx_id, ":", data$mutation)

# PROCESS MUTATION TYPES

# Mutation types: 
#A-SNV/A-SNP = Acceptor motif, SNV or SNP -> new_acceptor
#A-OTH = Acceptor motif, non-substitution mutation; -> new_acceptor
#D-SNV/D-SNP = Donor motif, SNV or SNP; -> new_donor
#D-OTH = Donor motif, non-susbtitution mutation;  -> new_donor
#BPD = Mutation (any size) affecting definition of the pseudoexon BP -> branchpoint
#DIST3 = Mutation 3´ of pseudoexon donor site; -> sre_intronic -> strenghtning_donor
#DIST5 = Mutation 5´ of pseudoexon acceptor and branch point motifs; -> strenghtning_acceptor
#ENCO = Insertion or inversion mutation where the inserted/inverted sequence completely encompasses the pseudoexon and its splice motifs; 
#M-SNV/M-SNP = Mid-pseudoexon, SNV or SNP;-> sre
#M-OTH = Mid-pseudoexon, non-substitution mutation; -> sre
#UNKN = Mutation not determined.

# Remove ENCO and DIST3/5 and tricky mutations
data <- data %>% filter(mutation_type != "ENCO" & mutation_type != "DIST3/5")
data <- data %>% filter(id != "PMS2-7-1") 
data <- data %>% filter(mutation != "c.1390-12_1390-13insL1(6044)")

# Change type of composite mutation cases
# (HADH-5-1) A-SNP + D-SNV
data[data$mutation == "c.636+385A>G",]$mutation_type <- 'A-SNP'
data[data$mutation == "c.636+385A>G",]$mutation_type_updated <- 'A-SNP'
data[data$mutation == "c.636+471G>T",]$mutation_type <- 'D-SNV'
data[data$mutation == "c.636+471G>T",]$mutation_type_updated <- 'D-SNV'

# (NTRK1-8-1a) D-SNV + M-SNV
data[data$mutation == "c.851-794C>G",]$mutation_type <- 'D-SNV'
data[data$mutation == "c.851-794C>G",]$mutation_type_updated <- 'D-SNV'
data[data$mutation == "c.851-798C>T",]$mutation_type <- 'M-SNV'
data[data$mutation == "c.851-798C>T",]$mutation_type_updated <- 'M-SNV'

# (DMD-2-3) M-OTH
data[data$mutation == "c.94-78858C>G",]$mutation_type <- 'M-SNV'
data[data$mutation == "c.94-78858C>G",]$mutation_type_updated <- 'M-SNV'
data[data$mutation == "c.94-78836T>G",]$mutation_type <- 'M-SNV'
data[data$mutation == "c.94-78836T>G",]$mutation_type_updated <- 'M-SNV'

# (COL4A3-48-1) D-SNV
data[data$mutation == "c.4463-523C>G",]$mutation_type <- 'D-SNV'
data[data$mutation == "c.4463-523C>G",]$mutation_type_updated <- 'D-SNV'
data[data$mutation == "c.4463-537A>G",]$mutation_type <- 'M-SNV'
data[data$mutation == "c.4463-537A>G",]$mutation_type_updated <- 'M-SNV'

# (COL4A5-29-1) BPD 
# One of the composite mutations appear to influence more 
#the polypirimide tract close to the pseudoexon cryptic acceptor
data[data$mutation == "c.2395+1292G>T",]$mutation_type_updated <- 'A-SNV_outDinuc'

# Delete one OTC-9-1: same variants originates 2 slightly different pseudoexons.
# The one that the mutation corresponds to the donor site is going to be kept
data <- data %>% filter(id != "OTC-9-1a")

# Delete MYBPC3-12-1a and MYBPC3-20-2: same variant originates two pseudoexons.
# The one that the mutation corresponds to the donor site is going to be kept
data <- data %>% filter(!id %in% c("MYBPC3-12-1a", "MYBPC3-20-2"))

# Delete DMD-32-1b same variant originates two slightly different pseudoexons
# The one that the mutation corresponds to the acceptor site is going to be kept
data <- data %>% filter(id != c("DMD-32-1b"))

# Delete mutation that had a wrong allele 
data <- data %>% filter(id != c("PSMC3-10-1"))

#############################
## Simplify mutation types ##
#############################
# Intronic mutations outside the dinucleotide sequences that compose the pseudoexon
#were treated as strenghtening variants

# First and last pseudoexon positions that were mutated (and treated as A-SNV and D-SNV in the paper)
# were kept as they are.
map_dict <- c('A-OTH'='new_splice_acceptor',
              'A-SNP'='new_splice_acceptor',
              'A-SNV'='new_splice_acceptor',
              'A-SNV_outDinuc' = 'strengthening_acceptor',
              'A-SNP_outDinuc' = 'strengthening_acceptor',
              'A-OTH_outDinuc' = 'strengthening_acceptor',
              'D-OTH'='new_splice_donor',
              'D-SNP'='new_splice_donor',
              'D-SNV'='new_splice_donor',
              'D-SNV_outDinuc' = 'strengthening_donor',
              'D-SNP_outDinuc' = 'strengthening_donor',
              'D-OTH_outDinuc' = 'strengthening_donor',
              'DIST3'='strengthening_donor',
              'DIST5'='strengthening_acceptor',
              'BPD'='branchpoint_associated',
              'M-OTH'='change_sre',
              'M-SNV'='change_sre',
              'M-SNP'='change_sre')
data <- data %>% 
  mutate(mutation_flag = recode(mutation_type_updated, !!!map_dict, .default = NA_character_))

data %>% count(mutation_type_updated)
data %>% count(mutation_flag)
# Write output
data %>% write_tsv(file='1_processed_data.tsv')
write.csv(unique(data$hgvs), file="1_hgvs.txt", row.names=F, quote = F)

