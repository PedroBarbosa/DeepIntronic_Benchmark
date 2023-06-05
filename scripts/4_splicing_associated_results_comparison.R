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
library(broom)
setwd("~/git_repos/giga_science_reviews/")

list_files <- list(branchpoint = 'data/splicing_altering/per_category/branchpoint_associated/statistics_all_types_all.tsv',
                   pe_new_donor = 'data/splicing_altering/per_category/pseudoexon_activation/new_donor/statistics_all_types_all.tsv',
                   pe_donor_downstream = 'data/splicing_altering/per_category/pseudoexon_activation/donor_downstream/statistics_all_types_all.tsv',
                   pe_acceptor_associated = 'data/splicing_altering/per_category/pseudoexon_activation/acceptor_associated/statistics_all_types_all.tsv',
                   pe_sre  = 'data/splicing_altering/per_category/pseudoexon_activation/exonic_like/statistics_all_types_all.tsv',
                   partial_ir_new_donor = 'data/splicing_altering/per_category/partial_intron_retention/new_donor/statistics_all_types_all.tsv',
                   partial_ir_donor_downstream = 'data/splicing_altering/per_category/partial_intron_retention/donor_downstream/statistics_all_types_all.tsv',
                   partial_ir_acceptor_associated = 'data/splicing_altering/per_category/partial_intron_retention/acceptor_associated/statistics_all_types_all.tsv',
                   partial_ir_sre = 'data/splicing_altering/per_category/partial_intron_retention/exonic_like/statistics_all_types_all.tsv')

# Concatenate dzata
df <- map2_df(list_files, 
              names(list_files),
              ~read_tsv(.x) %>% 
                mutate(analysis = .y) %>% 
                select(c('tool', 'norm_mcc', 'weighted_norm_mcc', 'weighted_F1', 'total_p', 'total_n',
                         'tp', 'fp', 'tn', 'fn',
                         'fraction_nan', 'weighted_accuracy', 'auROC', 'average_precision', 'analysis'))) 

# Remove AP score from tools with more than 50% of missing data since ROC curves are not drawn in such cases
df$average_precision <- ifelse(df$tool == "MLCsplice" & df$analysis == "branchpoint", NA, df$average_precision)
df$average_precision <- ifelse(df$tool == "IntSplice2" & df$analysis == "branchpoint", NA, df$average_precision)

to_heatmap <- df
df <- df %>% drop_na()

df$analysis <- factor(df$analysis, levels = c("branchpoint", "pe_acceptor_associated", "pe_sre", "pe_new_donor", "pe_donor_downstream",
                                              "partial_ir_acceptor_associated", "partial_ir_sre", "partial_ir_new_donor", "partial_ir_donor_downstream"))


#######################
###### Boxplots #######
#######################
major_label_f <- function(x){

  if (grepl("acceptor", x)){
    return("Acceptor associated")
  }
  else if (grepl("sre", x)){
    return("Exonic-like")
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
pe$region <- sapply(pe$analysis, major_label_f)

partial_ir <- df %>% filter(str_detect(analysis, "partial_ir_"))
partial_ir$group <- 'Partial intron retention'
partial_ir$region <- sapply(partial_ir$analysis, major_label_f)
df <- bind_rows(pe, partial_ir)

###########################
###### Per region #########
###########################
ggplot(df %>% filter(analysis != "branchpoint"), aes(x=group, y=average_precision,fill=group)) +
  geom_boxplot(alpha=0.5) + 
  scale_fill_manual(values=c("#619CFF", "#F8766D")) +
  geom_point() +
  facet_wrap(~region) + 
  theme_bw() +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size=16, hjust = 0.5),
        legend.text = element_text(size=13)) +
  ylim(0.2, 1) +
  labs(x = '', y = 'Average Precision')

########################
###### Per tool ########
########################
fisher_test <- function(data){
  if (nrow(data) == 1){
    col_names <- c("estimate", "p.value", "conf.low", "conf.high", "method", "alternative", "tool", "region")
    values <- c(NA, NA, NA, NA, NA, NA, data %>% pull(tool), data %>% pull(region))
    fisher <- tibble(!!!setNames(values, nm = col_names)) %>% 
      mutate_at(vars(estimate, p.value, conf.low, conf.high), as.double)
    return(fisher)
  }
  else{
    pe <- data %>% filter(group == "Pseudoexon activation")
    pir <- data %>% filter(group == "Partial intron retention")
    pe_success <- pe %>% pull(tp) + pe %>% pull(tn)
    pe_failure <- pe %>% pull(fp) + pe %>% pull(fn)
    pir_success <- pir %>% pull(tp) + pir %>% pull(tn)
    pir_failure <- pir %>% pull(fp) + pir %>% pull(fn)
    matrix <- matrix(c(pe_success, pe_failure, pir_success, pir_failure), nrow = 2)
    fisher = tidy(fisher.test(matrix))
    fisher$region = unique(data$region)
    fisher$tool = unique(data$tool)
    return(fisher)
  }
}

color_df <- data.frame(colors = c("coral4", "aquamarine4", "darkslateblue", "darkgrey"),
                       region = c("Acceptor associated", "Donor downstream", "Exonic-like", "New splice donor"))

fisher_df <- bind_rows(df %>% group_by(tool, region) %>% 
  group_map(~fisher_test(.x), .keep = T))
pval_corr <- fisher_df %>% filter(!is.na(p.value))
pval_corr$p.adj <- p.adjust(pval_corr$p.value, method = "holm")
fisher_df <- left_join(fisher_df, pval_corr)%>% left_join(color_df)

table(pval_corr$p.value < 0.05)
table(pval_corr$p.adj < 0.05)


df <- left_join(df, fisher_df %>% select(tool, region, p.value, p.adj, colors))

df_ <- df %>% mutate(p.value = round(p.value, digits = 5))
df_ <- df %>% mutate(p.adj = round(p.adj, digits = 5))
df_ <- df_ %>% mutate(p.adj_lower_1 = ifelse(p.adj < 1, p.adj, ""))

df_$p.value <- as.character(df_$p.value) 
df_$p.adj <- as.character(df_$p.adj) 
df_$p.adj_lower_1 <- as.character(df_$p.adj_lower_1) 

df_ <- df_ %>% mutate(p.value = sub("^", "pval=", p.value))
df_ <- df_ %>% mutate(p.adj = sub("^", "pval=", p.adj))
df_ <- df_ %>% mutate(p.adj_lower_1 = ifelse(p.adj_lower_1 == "", "", sub("^", "p=", p.adj_lower_1)))

df_$region <- factor(df_$region, levels = c("Acceptor associated", "Exonic-like", "New splice donor", "Donor downstream"))
df_ <- df_ %>% filter(!is.na(p.value))

ggplot(df_, aes(x=reorder(tool, -average_precision, median))) +
  geom_boxplot(aes(x=group, y=average_precision, fill=group, alpha=0.5), outlier.shape=NA) + 
  geom_point(aes(x=group, y=average_precision, fill = region), position = position_dodge(width=0.2), size=3, shape = 21)  +
  scale_fill_manual(values=c("coral4", "aquamarine4", "darkslateblue", "darkgrey", "#619CFF", "#F8766D")) +
  geom_text(data = subset(df_, region == "Acceptor associated"), aes(x = 1.5, y = 1.21, label = p.adj, color=colors), size = 2.5) +
  geom_text(data = subset(df_, region == "Exonic-like"), aes(x = 1.5, y = 1.09, label = p.adj, color=colors), size = 2.5) +
  geom_text(data = subset(df_, region == "New splice donor"), aes(x = 1.5, y = 1.03, label = p.adj, color=colors), size = 2.5) +
  geom_text(data = subset(df_, region == "Donor downstream"), aes(x = 1.5, y = 1.15, label = p.adj, color=colors), size = 2.5) +
  scale_colour_manual(values = c("aquamarine4", "coral4", "darkgrey", "darkslateblue")) +
  facet_wrap(reorder(tool, -average_precision, median) ~ ., nrow=2) +
  theme_bw() +
  ylim(0.2, 1.21) +
  ylab('Average precision') +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=7),
        panel.background = element_blank(),
        legend.text = element_text(size=13)) 


#####################################
######## All regions together #######
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

pe_all$tool = factor(pe_all$tool, levels = df %>% filter(analysis == "pe_acceptor_associated") %>% arrange(-average_precision) %>% pull(tool))
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
all$tool <- factor(all$tool, levels = all %>% filter(analysis == "pe_acceptor_associated") %>% arrange(-average_precision) %>% pull(tool))
all$label <- factor(all$label, levels = c(labels_pe, labels_partial_ir))
all$group <- factor(all$group, levels=c('Pseudoexon activation','Partial intron retention'))
ggplot(all, aes(x=label, y=average_precision, group=tool, color =tool)) +
  geom_bump(smooth = 15, size = 0.5, alpha = 1) + 
  geom_point(size=5, position = position_dodge(0.15)) +
  scale_color_manual(values = pal) +
  facet_grid(~group, scales = "free_x") +
  labs(x = '', y = "Average precision") +
  #scale_x_discrete(labels = c(labels_pe, labels_partial_ir)) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = c(0.5, 0.6),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
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
  df_ <- list_dfs[[i]] %>% select(tool,auROC, average_precision)
  missing_tools <- setdiff(unique_tools, df_$tool)
  
  not_perfomed <- tibble(tool = missing_tools, 
                           auROC = rep(c(0.18), times = length(missing_tools)),
                           average_precision = rep(c(0.18), times = length(missing_tools)))
  df_ <- bind_rows(df_, not_perfomed)
  df_$analysis <- group_name
  return(df_)
}

list_dfs <- to_heatmap %>% 
  group_split(analysis) %>%
  setNames(unique(to_heatmap$analysis))

group_names <- sort(unique(to_heatmap$analysis))
out <- lapply(seq_along(list_dfs), fill_missing_analysis) %>% bind_rows() %>%
  mutate(auROC = replace_na(auROC, 0.19),
         average_precision = replace_na(average_precision, 0.19))

# Change MLCsplice and dbscSNV. They were run, but had all predictions missed, b
# eing wrongly assigned to the "no run for variant region" class
# out <-out %>% mutate(average_precision=case_when(tool=="MLCsplice" & analysis == "pe_new_donor" ~ 0.19,
#                                    tool=="dbscSNV" & analysis %in% c("pe_new_donor", "partial_ir_new_donor") ~ 0.19,
#                                    TRUE ~ average_precision))
out <-out %>% mutate(average_precision=case_when(tool=="IntSplice2" & analysis %in% c("partial_ir_acceptor_associated", "pe_acceptor_associated") ~ 0.19,
                                                 TRUE ~ average_precision))

out <-
  out %>% mutate(average_precision_bins = cut(
    average_precision,
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

final_heat <- out %>% select(-auROC, -average_precision) %>% pivot_wider(names_from = analysis, values_from = average_precision_bins) %>% column_to_rownames(var = 'tool')
final_heat <- out %>% select(-auROC, -average_precision_bins) %>% pivot_wider(names_from = analysis, values_from = average_precision) %>% column_to_rownames(var = 'tool')
tool_order <- c("SpliceAI", "Pangolin", "CI-SpliceAI", "ConSpliceML", "AbSplice-DNA", "PDIVAS", "SPiP", "SQUIRLS", "TraP", "MMSplice", "SPIDEX", 
                "kipoiSplice4", "MLCsplice","LaBranchoR", "BPHunter", "BPP", "SVM_BP_finder", "IntSplice2", "MaxEntScan", "SpliceRover", "Spliceator", "DSSP", "HAL", "HEXplorer", "ESRseq", "ESEfinder")
col_order <- c("branchpoint", "pe_acceptor_associated", "partial_ir_acceptor_associated", "pe_sre",
               "partial_ir_sre", "pe_new_donor", "partial_ir_new_donor", "pe_donor_downstream", 
               "partial_ir_donor_downstream")
final_heat <- final_heat[tool_order,col_order]

##################################
### Variant regions annotation ###
##################################
labels <- c("BP", "Acceptor associated", "Exonic-like", "New splice donor", "Donor downstream")
fill <- rep("grey90", times=5) #c("#D8E2DC", "#FFE5D9", "#FFCAD4", "#F4ACB7", "#9D8189")
regions_block <- anno_block(gp = gpar(fill = fill),
                            labels = labels,
                            labels_gp = gpar(col = 'black', fontsize = 8.75))

split <- rep(factor(c("bp", "acc", "sre", "new_don", "d_down"), 
                    levels = c("bp", "acc", "sre", "new_don", "d_down")), 
             times=c(1, 2, 2, 2, 2))

################
### Palette ####
################
# mypal <- carto.pal(pal1 = "pastel.pal", n1 = 10, transparency = F)
# k <- length(mypal)
# image(1:k, 1, as.matrix(1:k), col =mypal, xlab = paste(k," classes",sep=""),
#       ylab = "", xaxt = "n", yaxt = "n",bty = "n")
library(circlize)
ccc <- colorRampPalette(c("#0b253f", "#e2e8eb"))
fill_labels <- c(ccc(35))[c(1,4,8,12,16,20,24,28,32)]
col_fun = colorRamp2(c(0.18, 0.19, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
                     c(c("white", "#FFE8EA"), fill_labels))
#library(cartography)
#ccc <- carto.pal(pal1 = "blue.pal" ,n1 = 20, middle = TRUE, transparency = TRUE)[c(1,3,5,7,9,11,13,15,18)] 
# col_fun = colorRamp2(c(0.18, 0.19, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
#                      c(c("white", "#FFE8EA"), rev(ccc)))

region_annot <- c("Branchpoint associated", rep(c("Pseudoexon activation", "Partial intron retention"), times=4))

htmp <- Heatmap(
  final_heat,
  cluster_column_slices = F,
  cluster_row_slices = F,
  cluster_rows = F,
  rect_gp = gpar(col = "black", lwd = 1),
  name = "Average precision",
  border = T,
  col = col_fun,
  row_names_side = 'left',
  column_names_side = 'top',
  column_split = split,
  column_gap = unit(4, "mm"),
  cluster_columns = F,
  show_column_names = F,
  column_title = NULL,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v = pindex(final_heat, i, j)
    l<-which(lapply(v,length)>0)
    if (final_heat[i, j] > 0.2){
      if(final_heat[i, j] > 0.6){
        grid.text(sprintf("%.2f", final_heat[i, j]), x, y, gp = gpar(col="white", fontsize = 8))
      }
      else{
        grid.text(sprintf("%.2f", final_heat[i, j]), x, y, gp = gpar(fontsize = 8))}
    }
  },
  top_annotation = HeatmapAnnotation(
    regions = regions_block,
    "Major variant group" = region_annot,
    show_annotation_name = F,
    gap = unit(0.1, "cm"),
    border = T,
    simple_anno_size = unit(3, "mm"),
    col = list(
      "Major variant group" = c(
        "Branchpoint associated" = "#C99FC6",
        "Pseudoexon activation" = "#8AAF96",
        "Partial intron retention" = "#8A5946"
      )
    )
  ),
)

lgd = Legend(
  labels = c("Not run for variant region", "Too many missing predictions"), 
  legend_gp = gpar(fill = c("white", "#FFE8EA")),
  title = "", 
  border = "black")

draw(htmp,  merge_legend=T, annotation_legend_list = list(lgd))
  
