library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
setwd("~/git_repos/paper_intronic_benchmark/data/keegan_et_al_2022")

data <- read_xlsx('supp_tables.xlsx', sheet = 1, skip= 1) %>% select(c(Gene,
                                                                       Intron,
                                                                       Size,
                                                                       `Transcript ID`, 
                                                                       `Instigating mutation(s)`,
                                                                       `Mutation type`,
                                                                       `RNA source`))
names(data) <- tolower(names(data))
data <- data %>% rename(tx_id = `transcript id`, 
                        mutation = `instigating mutation(s)`,
                        mutation_type = `mutation type`,
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

# Mutation types: 
#A-SNV/A-SNP = Acceptor motif, SNV or SNP;
#A-OTH = Acceptor motif, non-substitution mutation; 
#D-SNV/D-SNP = Donor motif, SNV or SNP;
#D-OTH = Donor motif, non-susbtitution mutation; 
#BPD = Mutation (any size) affecting definition of the pseudoexon BP
#DIST3 = Mutation 3´ of pseudoexon donor site;
#DIST5 = Mutation 5´ of pseudoexon acceptor and branch point motifs;
#ENCO = Insertion or inversion mutation where the inserted/inverted sequence completely encompasses the pseudoexon and its splice motifs; 
#M-SNV/M-SNP = Mid-pseudoexon, SNV or SNP; 
#M-OTH = Mid-pseudoexon, non-substitution mutation; 
#UNKN = Mutation not determined.

colnames(data)
data %>% filte