library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(ggbump)
library(forcats)

setwd("~/git_repos/paper_intronic_benchmark/")

list_files <- list(branchpoint = 'data/splicing_altering/per_category/branchpoint_associated/statistics_all_types_all.tsv',
                   pe_new_donor = 'data/splicing_altering/per_category/pseudoexon_activation/new_donor/statistics_all_types_all.tsv',
                   pe_donor_downstream = 'data/splicing_altering/per_category/pseudoexon_activation/donor_downstream/statistics_all_types_all.tsv',
                   pe_acceptor_associated = 'data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/statistics_all_types_all.tsv',
                   pe_sre  = 'data/splicing_altering/per_category/pseudoexon_activation/sre_associated/statistics_all_types_all.tsv',
                   partial_ir_new_donor = 'data/splicing_altering/per_category/partial_intron_retention/new_donor/statistics_all_types_all.tsv',
                   partial_ir_donor_downstream = 'data/splicing_altering/per_category/partial_intron_retention/donor_downstream/statistics_all_types_all.tsv',
                   partial_ir_acceptor_associated = 'data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/statistics_all_types_all.tsv',
                   partial_ir_sre = 'data/splicing_altering/per_category/partial_intron_retention/sre_associated/statistics_all_types_all.tsv')

# Concatenate dzata
df <- map2_df(list_files, 
              names(list_files),
              ~read_tsv(.x) %>% 
                mutate(analysis = .y) %>% 
                select(c('tool', 'norm_mcc', 'weighted_norm_mcc', 'weighted_F1', 'total_p', 'total_n',
                         'fraction_nan', 'weighted_accuracy', 'auROC', 'pr_auROC', 'analysis'))) 
to_heatmap <- df
df <- df %>% drop_na()

df$analysis <- factor(df$analysis, levels = c("branchpoint", "pe_acceptor_associated", "pe_sre", "pe_new_donor", "pe_donor_downstream",
                                              "partial_ir_acceptor_associated", "partial_ir_sre", "partial_ir_new_donor", "partial_ir_donor_downstream"))


####################################
########## Per class boxplot #######
####################################
major_label_f <- function(x){

  if (grepl("acceptor", x)){
    return("Acceptor associated")
  }
  else if (grepl("sre", x)){
    return("Splicing regulatory elements")
  }
  else if (grepl("new_donor", x)){
    return("New splice donor")
  }
  else if (grepl("donor_downstream",x )){
    return("Donor downstream")
  }
}

pe <- df %>% filter(str_detect(analysis, "pe_"))
pe$group <- 'Pseudoexon activation'
pe$major_group <- sapply(pe$analysis, major_label_f)

partial_ir <- df %>% filter(str_detect(analysis, "partial_ir_"))
partial_ir$group <- 'Partial intron retention'
partial_ir$major_group <- sapply(partial_ir$analysis, major_label_f)
df <- bind_rows(pe, partial_ir)
   
ggplot(df %>% filter(analysis != "branchpoint"), aes(x=group, y=pr_auROC,fill=group)) +
  geom_boxplot(alpha=0.5) + 
  geom_point() +
  facet_wrap(~major_group) + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=11),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=11),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 10, colour = "black"),
        plot.title = element_text(size=15, hjust = 0.5),
        legend.text = element_text(size=12)) +
  ylim(0.2, 1) +
  labs(x = '', y = 'Average Precision')

ggplot(df, aes(x=reorder(tool, -pr_auROC, median), y=pr_auROC, fill = group)) +
  geom_boxplot(alpha=0.5, outlier.colour = NULL) + 
  #geom_point(shape = 21) +
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 75),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.text = element_text(size=12)) +
  ylim(0.2, 1) +
  labs(x = '', y = 'Average Precision')

#####################################
########   All categories  ##########
#####################################
# Keep only tools that predict all categories
pe_all <- pe %>% 
  group_by(tool) %>% 
  filter(n() >= 4)

partial_ir_all <- partial_ir %>% 
  group_by(tool) %>% 
  filter(n() >= 4 & tool != "MMSplice")


# "#A6CEE3", "#A6761D", "#D95F02", "#FDCDAC", "#A6D854", "#984EA3", "#D9D9D9",
# "#CCEBC5", "#FFFFB3", "#E7298A", "#FFED6F", "#66A61E", "#666666", "#BEAED4", "#E6F5C9"

pe_all$tool = factor(pe_all$tool, levels = df %>% filter(analysis == "pe_acceptor_associated") %>% arrange(-pr_auROC) %>% pull(tool))
labels_pe <- c(paste0("Acceptor associated\nN pos=", 
                      pe_all %>% filter(analysis == "pe_acceptor_associated") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_acceptor_associated") %>% pull(total_n) %>% first()),  
               paste0("Exonic-like\nN pos=", 
                      pe_all %>% filter(analysis == "pe_sre") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_sre") %>% pull(total_n) %>% first()),
               paste0("New splice donor\nN pos=", 
                      pe_all %>% filter(analysis == "pe_new_donor") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_new_donor") %>% pull(total_n) %>% first()),
               paste0("Donor downstream\nN pos=", 
                      pe_all %>% filter(analysis == "pe_donor_downstream") %>% pull(total_p) %>% first(), 
                      ";\nN neg=", 
                      pe_all %>% filter(analysis == "pe_donor_downstream") %>% pull(total_n) %>% first()))

labels_partial_ir <- c(paste0("Acceptor associated\nN pos=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_acceptor_associated") %>% pull(total_p) %>% first(), 
                              ";\nN neg=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_acceptor_associated") %>% pull(total_n) %>% first()),  
                       paste0("Exonic-like\nN pos=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_sre") %>% pull(total_p) %>% first(), 
                              ";\nN neg=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_sre") %>% pull(total_n) %>% first()),
                       paste0("New splice donor\nN pos=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_new_donor") %>% pull(total_p) %>% first(), 
                              ";\nN neg=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_new_donor") %>% pull(total_n) %>% first()),
                       paste0("Donor associated\nN pos=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_donor_downstream") %>% pull(total_p) %>% first(), 
                              ";\nN neg=", 
                              partial_ir_all %>% filter(analysis == "partial_ir_donor_downstream") %>% pull(total_n) %>% first()))

label_map <- tibble(analysis=c("pe_acceptor_associated",
                                  "pe_sre", 
                                  "pe_new_donor", 
                                  "pe_donor_downstream",
                                  "partial_ir_acceptor_associated",
                                  "partial_ir_sre", 
                                  "partial_ir_new_donor", 
                                  "partial_ir_donor_downstream"),
                       label=c(labels_pe, labels_partial_ir))

pal <- c("#B15928", "#377EB8", "#E5D8BD", "#FBB4AE", "#E7298A", "#984EA3", "#7dc99e", "#999999")

all <- bind_rows(pe_all, partial_ir_all)
all <- left_join(all, label_map)
all$tool <- factor(all$tool, levels = all %>% filter(analysis == "pe_acceptor_associated") %>% arrange(-pr_auROC) %>% pull(tool))
all$label <- factor(all$label, levels = c(labels_pe, labels_partial_ir))
all$group <- factor(all$group, levels=c('Pseudoexon activation','Partial intron retention'))
ggplot(all, aes(x=label, y=pr_auROC, group=tool, color =tool)) +
  geom_bump(smooth = 15, size = 1.2, alpha = 1) + 
  geom_point(size=5, position = position_dodge(0.15)) +
  scale_color_manual(values = pal) +
  facet_grid(~group, scales = "free_x") +
  labs(x = '', y = "Average precision") +
  #scale_x_discrete(labels = c(labels_pe, labels_partial_ir)) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = c(0.5, 0.6),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 17),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(5, "lines"),
        panel.background = element_blank(),
        plot.title = element_text(size=15, hjust = 0.5),
        legend.text = element_text(size=12)) +
  ylim(0.2, 1) 



#################################
########## Heatmap ##############
#################################
library(ComplexHeatmap)
to_heatmap <- to_heatmap %>% filter(tool != "S-CAP")
unique_tools <- unique(to_heatmap$tool)

fill_missing_analysis <- function(i){
  group_name <- group_names[[i]] 
  df_ <- list_dfs[[i]] %>% select(tool,auROC, pr_auROC)
  missing_tools <- setdiff(unique_tools, df_$tool)
  
  not_perfomed <- tibble(tool = missing_tools, 
                           auROC = rep(c(-1), times = length(missing_tools)),
                           pr_auROC = rep(c(-1), times = length(missing_tools)))
  df_ <- bind_rows(df_, not_perfomed)
  df_$analysis <- group_name
  return(df_)
}

list_dfs <- to_heatmap %>% 
  group_split(analysis) %>%
  setNames(unique(to_heatmap$analysis))

group_names <- sort(unique(to_heatmap$analysis))
out <- lapply(seq_along(list_dfs), fill_missing_analysis) %>% bind_rows() %>%
  mutate(auROC = replace_na(auROC, -2),
         pr_auROC = replace_na(pr_auROC, -2))

# Change MLCsplice and dbscSNV. They were run, but had all predictions missed, b
# eing wrongly assigned to the "no run for variant region" class
out <-out %>% mutate(pr_auROC=case_when(tool=="MLCsplice" & analysis == "pe_new_donor" ~ -2,
                                   tool=="dbscSNV" & analysis %in% c("pe_new_donor", "partial_ir_new_donor") ~ -2,
                                   TRUE ~ pr_auROC))

out <-
  out %>% mutate(pr_auROC_bins = cut(
    pr_auROC,
    breaks = c(-2.1,-1.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    labels = c(
      "Too many missing predictions",
      "Not run for variant region",
      "> 0",
      "> 0.1",
      "> 0.2",
      "> 0.3",
      "> 0.4",
      "> 0.5",
      "> 0.6",
      "> 0.7",
      "> 0.8",
      "> 0.9"
    )
  ))

final_heat <- out %>% select(-auROC, -pr_auROC) %>% pivot_wider(names_from = analysis, values_from = pr_auROC_bins) %>% column_to_rownames(var = 'tool')
tool_order <- c("SpliceAI", "Pangolin", "CI-SpliceAI", "ConSpliceML", "AbSplice-DNA", "SPiP", "SQUIRLS", "TraP", "regSNP-intron", "MMSplice", "SPIDEX", "dbscSNV",
                "kipoiSplice4", "MLCsplice", "HAL", "LaBranchoR", "BPP", "SVM_BP_finder", "BPHunter", "IntSplice2", "SpliceRover", "Spliceator", "DSSP", "HEXplorer", "ESRseq")
col_order <- c("branchpoint", "pe_acceptor_associated", "partial_ir_acceptor_associated", "pe_sre",
               "partial_ir_sre", "pe_new_donor", "partial_ir_new_donor", "pe_donor_downstream", 
               "partial_ir_donor_downstream")
final_heat <- final_heat[tool_order,col_order]

##################################
##### Colors of Heatmap values ###
##################################
value_labels <- rev(c("Too many missing predictions",
                  "Not run for variant region",
                  "> 0",
                  "> 0.1",
                  "> 0.2",
                  "> 0.3",
                  "> 0.4",
                  "> 0.5",
                  "> 0.6",
                  "> 0.7",
                  "> 0.8",
                  "> 0.9"))
fill_labels <- c(rev(brewer.pal(n = 11, name = 'BrBG')[c(1,2,3,5,7,8,9,10,11)]),"#03045E",  "#F5F5F5", "darkgrey")
#fill_labels <- c("#03045E", "#6891C3", "#A4C3D2", "#4f6d7a", "#e0fbfc", "#03045E", 1:3, "#e0fbfc", "#ece2d0", "#FDAAAA")
#"#ece2d0"
#0096c7
colors <- structure(fill_labels, names=value_labels)

##################################
### Variant regions annotation ###
##################################
labels <- c("BP", "Acceptor associated", "Exonic-like", "New splice donor", "Donor downstream")
fill <- rep("lightgrey", times=5) #c("#D8E2DC", "#FFE5D9", "#FFCAD4", "#F4ACB7", "#9D8189")
regions_block <- anno_block(gp = gpar(fill = fill),
                                       labels = labels,
                            labels_gp = gpar(col = 'black', fontsize = 8))

split <- rep(factor(c("bp", "acc", "sre", "new_don", "d_down"), 
                    levels = c("bp", "acc", "sre", "new_don", "d_down")), 
             times=c(1, 2, 2, 2, 2))

#################################
##### Major group annotation ####
#################################
major_group_annot <- c("Branchpoint associated", rep(c("Pseudoexon activation", "Partial intron retention"), times=4))

ht_opt$COLUMN_ANNO_PADDING = unit(0.2, "cm")

htmp <- Heatmap(
  final_heat,
  cluster_column_slices = F,
  cluster_row_slices = F,
  cluster_rows = F,
  rect_gp = gpar(col = "black", lwd = 2),
  name = "Average precision",
  border = T,
  col = colors,
  row_names_side = 'left',
  column_names_side = 'top',
  column_split = split,
  column_gap = unit(4, "mm"),
  cluster_columns = F,
  show_column_names = F,
  column_title = NULL,
  top_annotation = HeatmapAnnotation(
    regions = regions_block,
    "Major variant group" = major_group_annot,
    show_annotation_name = F,
    gap = unit(0.1, "cm"),
    border = T,
    simple_anno_size = unit(3, "mm"),
    col = list(
      "Major variant group" = c(
        "Branchpoint associated" = "#E8F6EF",
        "Pseudoexon activation" = "#874C62",
        "Partial intron retention" = "#3F3B6C"
      )
    )
  )
)
draw(htmp,  merge_legend=T)

