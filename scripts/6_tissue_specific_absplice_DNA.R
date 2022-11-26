library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(dutchmasters)

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

major_tissues <- gsub("_.*$","",tissues)
major_tissues <- replace(major_tissues, major_tissues=="Whole", "Blood")
major_tissues <- replace(major_tissues, major_tissues=="Cells", "Cultured_cells")
major_tissues <- replace(major_tissues, major_tissues=="Small", "Small_Intestine")
major_tissues <- replace(major_tissues, major_tissues=="Minor", "Minor_Salivary_Gland")