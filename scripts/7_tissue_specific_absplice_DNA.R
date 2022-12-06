library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)
library(dutchmasters)

setwd("~/git_repos/paper_intronic_benchmark/scripts")
tissues <- c("Adipose_Subcutaneous",
                 "Adipose_Visceral_Omentum",
                 "Adrenal_Gland",
                 "Artery_Aorta",
                 "Artery_Coronary",
                 "Artery_Tibial",
                 "Brain_Amygdala",
                 "Brain_Anterior_cingulate_cortex_BA24",
                 "Brain_Caudate_basal_ganglia",
                 "Brain_Cerebellar_Hemisphere",
                 "Brain_Cerebellum",
                 "Brain_Cortex",
                 "Brain_Frontal_Cortex_BA9",
                 "Brain_Hippocampus",
                 "Brain_Hypothalamus",
                 "Brain_Nucleus_accumbens_basal_ganglia",
                 "Brain_Putamen_basal_ganglia",
                 "Brain_Spinal_cord_cervical_c_1",
                 "Brain_Substantia_nigra",
                 "Breast_Mammary_Tissue",
                 "Cells_Cultured_fibroblasts",
                 "Cells_EBV_transformed_lymphocytes",
                 "Colon_Sigmoid",
                 "Colon_Transverse",
                 "Esophagus_Gastroesophageal_Junction",
                 "Esophagus_Mucosa",
                 "Esophagus_Muscularis",
                 "Heart_Atrial_Appendage",
                 "Heart_Left_Ventricle",
                 "Kidney_Cortex",
                 "Liver",
                 "Lung",
                 "Minor_Salivary_Gland",
                 "Muscle_Skeletal",
                 "Nerve_Tibial",
                 "Ovary",
                 "Pancreas",
                 "Pituitary",
                 "Prostate",
                 "Skin_Not_Sun_Exposed_Suprapubic",
                 "Skin_Sun_Exposed_Lower_leg",
                 "Small_Intestine_Terminal_Ileum",
                 "Spleen",
                 "Stomach",
                 "Testis",
                 "Thyroid",
                 "Uterus",
                 "Vagina",
                 "Whole_Blood")

variants <- read_tsv('../data/manual_curation/3_tissue_specificity/manual_curated_with_closest_GTEx_tissue.tsv')
absplice <- read_tsv('../data/manual_curation/3_tissue_specificity/AbSplice_DNA_data.tsv')
absplice_names <- c("Gene", tissues)

select_gene <- function(row){
  
  if(is.na(row[['gene2']])){
    return(row[['gene1']])
  }
  else{
    preds <- unlist(row[["gene1"]] %>% str_split(., pattern = "[|]"))
    preds <- preds[-1]
    preds <- as.numeric(preds)
    max1 = max(preds, na.rm = T)
    
    preds <- unlist(row[["gene2"]] %>% str_split(., pattern = "[|]"))
    preds <- preds[-1]
    preds <- as.numeric(preds)
    max2 = max(preds, na.rm = T)
    
    if(max1 >= max2){
      return(row[['gene1']])
    }
    else{
      return(row[["gene2"]])

  }
  }
}

# Select max prediction per variant
absplice <-absplice %>% separate(AbSplice_DNA, into = c("gene1", "gene2"), sep = "[,]")
absplice$AbSplice_DNA <- apply(absplice, 1, select_gene)
absplice <- absplice %>% select(-c(gene1, gene2))

# Separate preds per tissue into different columns
absplice_wide <- absplice %>% separate(AbSplice_DNA, into = absplice_names, extra = 'merge', 
                                  sep = "[|]")

# Merge variants info
absplice_wide <- left_join(absplice_wide, variants %>% select(`#CHROM`, POS, REF, ALT, `Phenotype (OMIM)`, `Closest GTEx Major Tissue affected`))
absplice_wide <- absplice_wide %>% 
  rename(Disease = `Phenotype (OMIM)`, Disease_major_tissue = `Closest GTEx Major Tissue affected`) 
absplice_wide <- absplice_wide %>% 
  mutate(Disease=str_trim(sub("\\(.*", "", Disease)))

# Remove Diseases with no tissue  
absplice_wide <- absplice_wide %>% filter(Disease_major_tissue != "-")

# Select eye-associated diseases
# eye_retina_diseases <- c("Achromatopsia", "Choroideremia", "Cone-rod dystrophy",
#                           "FRMD7-related infantile nystagmus", "Gyrate atrophy of choroid and retina",
#                           "Leber congenital amaurosis", "Ocular albinism", "Oguchi disease", "Optic atrophy", "Photoreceptor Dystrophy",
#                           "Retinal dystrophy", "Retinitis pigmentosa", "Stargardt disease",
#                           "Retinoblastoma", "Stargardt disease", "Usher syndrome")
#  
#  absplice_wide <- absplice_wide %>% filter(Disease %in% eye_retina_diseases)


##################
#### Heatmap #####
##################
library(ComplexHeatmap)

heat <- absplice_wide %>% select(-c(`#CHROM`, POS, REF, ALT, Consequence, SYMBOL, HGVSc, HGVSg, Gene))

# Select diseases associated with 1 tissue
heat <- heat %>% filter(!str_detect(Disease_major_tissue, ";"))
# Select diseases associated with more than 1 tissue
#heat <- heat %>% filter(str_detect(Disease_major_tissue, ";"))

heat <- heat %>% 
  group_by(Disease) %>% mutate(n=1:n()) %>% 
  group_by(Disease_major_tissue) %>% 
  mutate(Display_group = if(n()>3) Disease_major_tissue else "Other") %>% 
  unite("Disease", c("Disease", "n"), remove = T) %>% 
  column_to_rownames(var = "Disease") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  arrange(Display_group) 

disease_group_map <- heat %>% rownames_to_column(var='Disease') %>% select(c(Disease, Display_group))

heat <- heat %>% select(-c(Display_group, Disease_major_tissue)) %>% 
  mutate_if(is.character, as.numeric)

heat <- t(scale(t(heat)))

# Remove Diseases with no tissue differences (SD is 0)
heat <- as.data.frame(heat) %>% rownames_to_column %>% as_tibble() %>% column_to_rownames() %>%  filter(if_any(everything(), ~ !is.na(.)))

# Join Tissue group
heat <- left_join(heat %>% rownames_to_column(var='Disease'), disease_group_map)

labels <- unique(heat$Display_group)
fill <- c(dutchmasters$pearl_earring[[11]], 
          dutchmasters$pearl_earring[[10]],
          dutchmasters$pearl_earring[[9]],
          dutchmasters$pearl_earring[[1]],
          dutchmasters$pearl_earring[[2]], 
          dutchmasters$milkmaid[[9]])

# fill <- c(dutchmasters$pearl_earring[[11]], 
#           dutchmasters$pearl_earring[[10]],
#           dutchmasters$pearl_earring[[9]])

ha = rowAnnotation(Muscle = anno_block(gp = gpar(fill = fill),
                                       labels = labels))

aux <- heat %>% group_by(Display_group) %>% summarise(n=n())
split = rep(aux %>% pull(Display_group), times = aux %>% pull(n))

final_heat <- heat %>% column_to_rownames('Disease') %>% select(-Display_group)

Heatmap(final_heat,      
        name="z-score", 
        na_col = "white",
        row_split = split,
        cluster_row_slices = F,
        cluster_columns = F,
        row_gap = unit(5, "mm"),
        row_title = NULL,
        rect_gp = gpar(col = "white", lwd = 2),
        border = T,
        column_names_rot = 75,
        left_annotation = ha)

Heatmap(final_heat,      
        name="z-score", 
        na_col = "white",
        cluster_columns = F,
        cluster_rows =F ,
        rect_gp = gpar(col = "white", lwd = 2),
        border = T,
        column_names_rot = 75)


###############################################
############### Non- heatmap ##################
###############################################
# Tidy format
# absplice <- absplice_wide %>% pivot_longer(cols = tissues, names_to = 'tissue', values_to = 'pred')
# 
# # Get major tissue
# absplice$major_tissue <- gsub("_.*$","",absplice$tissue)
# absplice <- absplice %>% mutate(major_tissue = recode(major_tissue, 'Whole' = 'Blood', 
#                                                       'Cells' = 'Cultured_cells', 
#                                                       'Small' = 'Small_intestine',
#                                                       'Minor' = 'Minor_salivary_gland'))
# 
# # Merge variants info
# #absplice <- left_join(absplice, variants %>% select(`#CHROM`, POS, REF, ALT, `Phenotype (OMIM)`, `Closest GTEx Major Tissue affected`))
# #absplice <- absplice %>% rename(Disease = `Phenotype (OMIM)`, Disease_major_tissue = `Closest GTEx Major Tissue affected`)
# 
# # Drop variants with no GTEx tissue associated with disease
# absplice <- absplice %>% filter(Disease_major_tissue != "-")
# absplice <- absplice %>% mutate_at('pred', as.numeric)
# 
# #Remove blood temporarily
# absplice <- absplice %>% filter(Disease_major_tissue != "Blood")
# #absplice <- absplice %>% filter(!str_detect(Disease_major_tissue, ";"))
# 
# # Compare scores in disease tissues vs non-diseases
# get_average_per_group <- function(group){
#   disease_tissues <- unlist(str_split(group %>% slice(1) %>% pull(Disease_major_tissue), pattern = "[;]"))
#   median_disease <- subset(group, major_tissue %in% disease_tissues) %>% summarise(median = median(pred, na.rm = T)) %>% pull(median)
#   median_non_disease <- subset(group, !major_tissue %in% disease_tissues) %>% summarise(median = median(pred, na.rm = T)) %>% pull(median)
#   out <- group %>% slice(1) %>% select(c(HGVSc, Gene, Disease, Disease_major_tissue))
#   
#   return(out %>% add_column(Median_disease = median_disease, Median_non_disease=median_non_disease))
# }
# 
# 
# to_scatter <- absplice %>% group_by(`#CHROM`, POS, REF, ALT) %>% group_modify(~get_average_per_group(.x)) %>% drop_na()
# to_scatter %>% ggplot(aes(x=Median_non_disease, y=Median_disease)) +
#   geom_point(fill="bisque2", colour='black', size = 2, pch = 21, alpha=0.8) +
#   #geom_smooth(method='lm') +
#   geom_abline(linetype='dashed') +
#   theme_classic() +
#   xlim(0,0.2) +
#   ylab('Median AbSplice-DNA (Disease tissues)') + 
#   xlab('Median AbSplice-DNA (Non-Disease tissues)') +
#   theme(legend.title = element_blank(), 
#         axis.text.y = element_text(size=11),
#         axis.text.x = element_text(size=11),
#         axis.title.y = element_text(size=12),
#         axis.title.x = element_text(size=12))
#   

