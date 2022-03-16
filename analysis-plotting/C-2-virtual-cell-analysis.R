# setup -------------------------------------------------------------------

source("function-library.R")

args <- commandArgs(trailingOnly=TRUE)
output_str <- ifelse(length(args)>0,args[1],"subset")
dir.create(paste('figures',output_str,sep = "/"),showWarnings = F)

library(ggrepel)
library(umap)
library(cowplot)
library(dendextend)

multi_objective_matching <- c("anemia_binary","binary","mds_binary","disease_binary")

conversion_list <- list(
  `Disease detection` = c(
    Normal = "Normal",`Non-SF3B1-mutant` = "Disease",
    `SF3B1-mutant` = "Disease",`SRSF2-mutant` = "Disease",
    `RUNX1-mutant` = "Disease",`U2AF1-mutant` = "Disease",
    `Iron deficiency` = "Disease",
    `Megaloblastic` = "Disease"),
  `Disease classification` = c(
    `Non-SF3B1-mutant` = "MDS",`SF3B1-mutant` = "MDS",
    `SRSF2-mutant` = "MDS",`RUNX1-mutant` = "MDS",`U2AF1-mutant` = "MDS",
    `Iron deficiency` = "Anaemia",`Megaloblastic` = "Anaemia"),
  `SF3B1mut detection` = c(
    `SF3B1-mutant` = "SF3B1-mutant",`SRSF2-mutant` = "SF3B1-wildtype",
    `RUNX1-mutant` = "SF3B1-wildtype",`U2AF1-mutant` = "SF3B1-wildtype",
    `Non-SF3B1-mutant` = "SF3B1-wildtype"),
  `Anaemia classification` = c(
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"),
  `Multi-objective` = c(
    Normal = "Normal",`SF3B1-mutant` = "SF3B1-mutant",`SRSF2-mutant` = "SF3B1-wildtype",
    `RUNX1-mutant` = "SF3B1-wildtype",`U2AF1-mutant` = "SF3B1-wildtype",
    `Non-SF3B1-mutant` = "SF3B1-wildtype",
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"))

all_conditions <- rbind(
  read_csv(
    "data/all_classes.csv",progress = F,
    col_types = c(col_character(),col_character(),col_character())) %>%
    select(slide_id = id,fine_class,coarse_class) %>% 
    mutate(fine_class = factor(fine_class_conversion[fine_class],fine_class_conversion),
           coarse_class = factor(class_conversion[coarse_class],class_conversion)),
  read_csv(
    "data/all_classes_adden_1.csv",progress = F,
    col_types = c(col_character(),col_character(),col_character())) %>%
    select(slide_id,fine_class,coarse_class) %>% 
    mutate(fine_class = factor(fine_class_conversion[fine_class],fine_class_conversion),
           coarse_class = factor(class_conversion[coarse_class],class_conversion)) %>%
    mutate(slide_id = gsub('\\.','',str_match(slide_id,'[0-9_-]+.'))),
  read_csv(
    "data/all_classes_adden_2.csv",progress = F,
    col_types = c(col_character(),col_character(),col_character())) %>%
    select(slide_id,fine_class,coarse_class) %>% 
    mutate(fine_class = factor(fine_class_conversion[fine_class],fine_class_conversion),
           coarse_class = factor(class_conversion[coarse_class],class_conversion))
)

read_final_layer <- function(path) {
  df <- read_csv(path,col_names=F)
  colnames(df) <- c("model_name","fold","class_id","task_idx","virtual_cell_type","coefficient")
  df$max_vc <- as.numeric(str_match(df$model_name,'[0-9]+'))
  
  df <- df %>%
    na.omit %>%
    mutate(model_name = gsub("models/","",model_name),
           virtual_cell_type = virtual_cell_type + 1) %>% 
    group_by(model_name) %>%
    filter(virtual_cell_type <= (max_vc*2)) %>%
    mutate(cell_type = ifelse(virtual_cell_type <= max_vc,"WBC","RBC")) %>% 
    mutate(
      virtual_cell_type = ifelse(
        virtual_cell_type > max_vc,virtual_cell_type - max_vc,virtual_cell_type)) %>% 
    mutate(data_type = ifelse(grepl("bc$",model_name),"Morphology + B.C.","Morphology")) %>%
    ungroup %>% 
    spread(key = "class_id",value = "coefficient") %>% 
    mutate(class_difference = as.numeric(`1`) - as.numeric(`0`))
  
  return(df)} 

read_vcq_layer <- function(path,wbc_subset_path=NULL,rbc_subset_path=NULL) {
  df <- read_csv(path,col_names=F)
  colnames(df) <- c("model_name","fold","feature","cell_type","virtual_cell_type","coefficient")
  df$max_vc <- as.numeric(str_match(df$model_name,'[0-9]+'))
  if (!is.null(wbc_subset_path)) {
    wbc_feature_subset <- unlist(read_csv(wbc_subset_path,col_names = F))
  } else {
    wbc_feature_subset <- sort(unique(df$feature[df$cell_type=="dataset_0"])) + 1
  }
  if (!is.null(rbc_subset_path)) {
    rbc_feature_subset <- unlist(read_csv(rbc_subset_path,col_names = F))
  } else {
    rbc_feature_subset <- sort(unique(df$feature[df$cell_type=="dataset_1"])) + 1
  }
  
  df <- df %>%
    na.omit %>%
    mutate(model_name = gsub("models/","",model_name),
           virtual_cell_type = 1+virtual_cell_type) %>% 
    group_by(model_name) %>%
    mutate(cell_type = ifelse(cell_type == 'dataset_0',"WBC","RBC")) %>% 
    mutate(feature = ifelse(
      cell_type == 'WBC',
      c(features_all,features_nuclear)[wbc_feature_subset][feature+1],
      features_all[rbc_feature_subset][feature+1]
    )) %>%
    mutate(
      virtual_cell_type = ifelse(
        virtual_cell_type > max_vc,virtual_cell_type - max_vc,virtual_cell_type)) %>% 
    mutate(data_type = ifelse(grepl("bc$",model_name),"Morphology + B.C.","Morphology")) %>%
    ungroup 
  
  return(df)}

read_proportions <- function(path) {
  df <- read_csv(path,col_names = F)
  colnames(df)[c(1,2)] <- c("slide_id","cell_type") 
  
  df <- df %>%
    mutate(cell_type = ifelse(cell_type == "cell_type_0","WBC","RBC")) %>%
    gather(key = "virtual_cell_type",value = "proportion",-slide_id,-cell_type) %>%
    mutate(virtual_cell_type = gsub("X","",virtual_cell_type)) %>%
    mutate(virtual_cell_type = as.numeric(virtual_cell_type)) %>%
    mutate(virtual_cell_type = virtual_cell_type - min(virtual_cell_type)+1) %>%
    mutate(max_vc = max(virtual_cell_type))
  
  return(df)}

coefficients_plot <- function(x) {
  x %>%
    ggplot(aes(x = virtual_cell_type_fctr,
               y = feature,
               fill = coefficient)) + 
    geom_tile() + 
    facet_wrap( ~ sprintf('%s (%s)',data_type,cell_type),scales = "free") + 
    theme_pretty(base_size = 6) + 
    theme(legend.position = "bottom") +
    xlab("Virtual cell type") + 
    ylab("Feature") +
    scale_fill_gradient2(low = 'blue4',mid = "white",high = "red4",
                         name = "Coefficient") + 
    theme(legend.key.height = unit(0.1,"cm"),
          strip.text = element_text(margin = ggplot2::margin(t = 1,b = 2)),
          axis.text = element_text(colour = "black")) +
    scale_x_discrete(expand = c(0,0),labels = function(x) str_match(x,'[0-9]+$')) +
    scale_y_discrete(expand = c(0,0)) %>%
    return
}

subset_vcq <- function(best_vcq_subset,S) {
  best_vcq_subset %>%
    subset(virtual_cell_type_fctr %in% S) %>%
    mutate(virtual_cell_type_fctr = factor(
      virtual_cell_type_fctr,levels=S,ordered = T)) %>%
    return
}

# data loading and processing ---------------------------------------------

so_layers <- paste(Filter(
  nchar,c("../mile-vice/best_models/best_layers",output_str)),collapse="_")
mo_layers <- paste(Filter(
  nchar,c("../mile-vice/best_models/best_layers_mo",output_str)),collapse="_")
vcq_layers <- paste(Filter(
  nchar,c("../mile-vice/best_models/best_vcq_layers",output_str)),collapse="_")
vcq_layers_mo <- paste(Filter(
  nchar,c("../mile-vice/best_models/best_vcq_layers_mo",output_str)),collapse="_")
wbc_subset_path <- NULL
rbc_subset_path <- NULL
if (grepl("subset",output_str)){
  wbc_subset_path <- "../mile-vice/scripts/wbc_feature_subset"
  rbc_subset_path <- "../mile-vice/scripts/rbc_feature_subset"
} 
if (grepl("unbiased",output_str)) {
  wbc_subset_path <- "../mile-vice/scripts/wbc_feature_subset_unbiased"
  rbc_subset_path <- "../mile-vice/scripts/rbc_feature_subset_unbiased"
}

feat_conv <- rev(gsub('\n',' ',features_conversion))

best_layers_subset <- read_final_layer(path = so_layers)
best_layers_subset_mo <- read_final_layer(path = mo_layers)
best_vcq_subset <- read_vcq_layer(
  path = vcq_layers,
  wbc_subset_path = wbc_subset_path,
  rbc_subset_path = rbc_subset_path) %>%
  mutate(virtual_cell_type_fctr = paste(
    decode_model_name(model_name),data_type,cell_type,virtual_cell_type)) %>%
  mutate(feature = factor(feat_conv[feature],levels = feat_conv))
best_vcq_subset_mo <- read_vcq_layer(
  path = vcq_layers_mo,
  wbc_subset_path = wbc_subset_path,
  rbc_subset_path = rbc_subset_path) %>%
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>%
  mutate(feature = factor(feat_conv[feature],levels = feat_conv))

all_cell_proportions <- list()
all_cell_proportions_mo <- list()

for (model in unique(best_layers_subset$model_name)) {
  path <- sprintf("../mile-vice/cell-proportions/mll/%s.csv",model)
  proportions <- read_proportions(path)
  proportions$model_name <- model
  all_cell_proportions[[model]] <- proportions
}

for (model in unique(best_layers_subset_mo$model_name)) {
  path <- sprintf("../mile-vice/cell-proportions/mll/%s.csv",model)
  proportions <- read_proportions(path) %>%
    distinct
  proportions$model_name <- model
  all_cell_proportions_mo[[model]] <- proportions
}

cell_proportions_df <- do.call(rbind,all_cell_proportions) %>%
  merge(all_conditions,by = "slide_id")
cell_proportions_mo_df <- do.call(rbind,all_cell_proportions_mo) %>%
  merge(all_conditions,by = "slide_id")

full_proportions_cell_type <- merge(
  cell_proportions_df,best_layers_subset,
  by=c("cell_type","virtual_cell_type","max_vc","model_name")) %>%
  mutate(effect = class_difference * proportion) %>%
  mutate(direction = ifelse(effect < 0,-1,1),ae = abs(effect)) %>%
  group_by(direction,cell_type,model_name,task_idx) %>%
  mutate(normalized_effect = direction*(ae-min(ae))/diff(range(ae)))
full_proportions_cell_type_mo <- merge(
  cell_proportions_mo_df,best_layers_subset_mo,
  by=c("cell_type","virtual_cell_type","max_vc","model_name")) %>%
  mutate(effect = class_difference * proportion) %>%
  mutate(direction = ifelse(effect < 0,-1,1),ae = abs(effect)) %>%
  group_by(direction,cell_type,model_name,task_idx) %>%
  mutate(normalized_effect = direction*(ae-min(ae))/diff(range(ae)))
full_proportions_cell_type_mo$task_idx <- multi_objective_matching[full_proportions_cell_type_mo$task_idx+1]

# single objective plots --------------------------------------------------

N_MIN <- 12
H <- 2
W <- 5
M <- 1
mm <- 0.5

generic_objects <- list( 
  theme = theme_pretty(base_size = 6) + theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.key.size = unit(0,"cm"),legend.position = "top",
    strip.text = element_text(margin = ggplot2::margin()),
    strip.background = element_rect(fill = NA,color = NA),
    legend.box.spacing = unit(0.1,"cm"),panel.spacing = unit(0.6,"cm"),
    legend.box.background = element_rect(fill=NA,colour="grey90")),
  facet = facet_wrap(~ sprintf('%s\n(%s)',data_type, cell_type),scales = "free_x",nrow = 1),
  average = geom_point(position = position_dodge(width = 0.5),size = 0.75,alpha = 0.5),
  interval_1 = geom_linerange(aes(ymin = y25,ymax = y75),position = position_dodge(width = 0.5),
                              size = 0.6,alpha = 0.5),
  interval_2 = geom_errorbar(position = position_dodge(width = 0.5),size = 0.25,alpha = 0.5,width=0.3),
  proportion = geom_point(colour = "grey10",
                          mapping = aes(y = log_ratio,shape = "Proportion ratio"),size = 1.0),
  scale_y = scale_y_continuous(
    sec.axis = sec_axis(~ 1*., name = "Average proportion ratio (\u25C7)",
                        breaks = mm*c(-2,-1,0,1,2),
                        labels = function(x) 2^(x/mm))),
  scale_x = scale_x_discrete(labels = function(x) str_match(x,'[0-9]+$')),
  scale_shape = scale_shape_manual(values = c(`Proportion ratio` = 5),name = NULL),
  label_y = ylab("Normalized effect size (\u25CF)"),
  label_x = xlab("Virtual cell type"),
  axis_line = geom_hline(yintercept = 0,size = 0.25,alpha = 1,colour = "grey90"))

# Disease detection

tmp_data <- full_proportions_cell_type %>%
  mutate(model_name = decode_model_name(model_name)) %>% 
  subset(model_name == "Disease detection") %>% 
  mutate(coarse_class = ifelse(coarse_class == "Normal","Normal","Disease")) %>% 
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,coarse_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(model_name,data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(coarse_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[2] - y[1])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>%
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio))

tmp_data %>% 
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = factor(coarse_class,c("Normal","Disease")))) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = fine_colours,name = NULL,breaks = c("Normal","Disease")) +
  ggsave(
    sprintf("figures/%s/mile-vice-so-disease-detection.pdf",output_str),device=cairo_pdf,height=H,width=W)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-disease-detection-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-disease-detection-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# Disease classification 

tmp_data <- full_proportions_cell_type %>%
  mutate(model_name = decode_model_name(model_name)) %>% 
  subset(model_name == "Disease classification" & coarse_class != "Normal") %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,coarse_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(model_name,data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(coarse_class) %>% 
  mutate(ratio = median_prop[1]/median_prop[2]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[1] - y[2])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = coarse_class)) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = fine_colours,name = NULL,breaks = c("MDS","Anaemia")) +
  ggsave(sprintf("figures/%s/mile-vice-so-disease-classification.pdf",output_str),device=cairo_pdf,height=H,width=W)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-disease-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-disease-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# SF3B1mut detection

tmp_data <- full_proportions_cell_type %>%
  mutate(model_name = decode_model_name(model_name)) %>% 
  subset(model_name == "SF3B1mut detection" & coarse_class == "MDS") %>% 
  mutate(fine_class = ifelse(fine_class == "SF3B1-mutant","SF3B1-mutant","SF3B1-wildtype")) %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,fine_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(model_name,data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(fine_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[1] - y[2])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  mutate(AD = abs(y[1] - y[2])) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = factor(fine_class,c("SF3B1-mutant","SF3B1-wildtype")))) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = c(fine_colours,`log(ratio)` = "grey10"),name = NULL,
                      breaks = c("SF3B1-mutant","SF3B1-wildtype")) + 
  ggsave(sprintf("figures/%s/mile-vice-so-mds-classification.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-mds-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-mds-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# Anaemia detection

tmp_data <- full_proportions_cell_type %>%
  mutate(model_name = decode_model_name(model_name)) %>% 
  subset(model_name == "Anaemia classification" & coarse_class == "Anaemia") %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,fine_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(model_name,data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(fine_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[2] - y[1])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = fine_class)) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = c(fine_colours,`log(ratio)` = "grey10"),name = NULL,
                      breaks = c("Iron deficiency","Megaloblastic")) + 
  ggsave(sprintf("figures/%s/mile-vice-so-anaemia-classification.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-anaemia-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-so-anaemia-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# multiple objective plots ------------------------------------------------

groups_of_cells <- list()

# Disease detection

tmp_data <- full_proportions_cell_type_mo %>%
  mutate(model_name = decode_model_name(model_name)) %>% 
  subset(task_idx == "binary") %>% 
  mutate(coarse_class = ifelse(coarse_class == "Normal","Normal","Disease")) %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,coarse_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(coarse_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[2] - y[1])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

groups_of_cells[[length(groups_of_cells)+1]] <- tmp_data %>% 
  select(cell_type,virtual_cell_type,y) %>%
  distinct %>%
  mutate(task_idx = "Disease detection")

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = factor(coarse_class,c("Normal","Disease")))) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = fine_colours,name = NULL,breaks = c("Normal","Disease")) +
  ggsave(sprintf("figures/%s/mile-vice-mo-disease-detection.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset_mo,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-disease-detection-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-disease-detection-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# Disease classification 

tmp_data <- full_proportions_cell_type_mo %>%
  subset(task_idx == "disease_binary" & coarse_class != "Normal") %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,coarse_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(coarse_class) %>% 
  mutate(ratio = median_prop[1]/median_prop[2]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[1] - y[2])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

groups_of_cells[[length(groups_of_cells)+1]] <- tmp_data %>% 
  select(cell_type,virtual_cell_type,y) %>%
  distinct %>%
  mutate(task_idx = "Disease classification")

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = coarse_class)) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = fine_colours,name = NULL,breaks = c("MDS","Anaemia")) + 
  ggsave(sprintf("figures/%s/mile-vice-mo-disease-classification.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset_mo,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-disease-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-disease-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# SF3B1mut detection

tmp_data <- full_proportions_cell_type_mo %>%
  subset(task_idx == "mds_binary" & coarse_class == "MDS") %>% 
  mutate(fine_class = ifelse(fine_class == "SF3B1-mutant","SF3B1-mutant","SF3B1-wildtype")) %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,fine_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(fine_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[1] - y[2])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

groups_of_cells[[length(groups_of_cells)+1]] <- tmp_data %>% 
  select(cell_type,virtual_cell_type,y) %>%
  distinct %>%
  mutate(task_idx = "SF3B1mut detection")

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = factor(fine_class,c("SF3B1-mutant","SF3B1-wildtype")))) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = c(fine_colours,`log(ratio)` = "grey10"),name = NULL,
                      breaks = c("SF3B1-mutant","SF3B1-wildtype")) + 
  ggsave(sprintf("figures/%s/mile-vice-mo-mds-classification.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset_mo,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-mds-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-mds-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

# Anaemia detection

tmp_data <- full_proportions_cell_type_mo %>%
  subset(task_idx == "anemia_binary" & coarse_class == "Anaemia") %>% 
  mutate(effect = class_difference * proportion) %>%
  group_by(model_name,data_type,virtual_cell_type,
           cell_type,fine_class) %>%
  summarise(ymin = min(normalized_effect),ymax = max(normalized_effect),
            y = median(normalized_effect),
            y25 = quantile(normalized_effect,0.25),y75 = quantile(normalized_effect,0.75),
            median_prop = median(proportion)) %>%
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>% 
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  arrange(fine_class) %>% 
  mutate(ratio = median_prop[2]/median_prop[1]) %>%
  mutate(log_ratio = mm*log(ratio,2),d = abs(y[2] - y[1])) %>% 
  mutate(abs_log_ratio = abs(log_ratio)) %>%
  group_by(data_type,cell_type) %>% 
  filter(d > rev(sort(d))[N_MIN]) %>% 
  mutate(virtual_cell_type_fctr = reorder(virtual_cell_type_fctr,log_ratio)) 

groups_of_cells[[length(groups_of_cells)+1]] <- tmp_data %>% 
  select(cell_type,virtual_cell_type,y) %>%
  distinct %>%
  mutate(task_idx = "Anaemia classification")

tmp_data %>%
  ggplot(aes(x = virtual_cell_type_fctr,
             y = y,ymin = ymin,ymax = ymax,
             colour = fine_class)) + 
  generic_objects$axis_line +
  generic_objects$proportion +
  generic_objects$average + 
  generic_objects$interval_1 +
  generic_objects$interval_2 +
  generic_objects$theme +
  generic_objects$facet +
  generic_objects$label_x +
  generic_objects$label_y +
  generic_objects$scale_x +
  generic_objects$scale_y + 
  generic_objects$scale_shape +
  scale_colour_manual(values = c(fine_colours,`log(ratio)` = "grey10"),name = NULL,
                      breaks = c("Iron deficiency","Megaloblastic")) +
  ggsave(sprintf("figures/%s/mile-vice-mo-anaemia-classification.pdf",output_str),height=H,width=W,device=cairo_pdf)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset_mo,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-anaemia-classification-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-anaemia-classification-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)

groups_of_cells %>%
  do.call(what = rbind) %>%
  group_by(virtual_cell_type,data_type,cell_type) %>%
  mutate(total_relevance = length(task_idx) + mean(abs(y))) %>% 
  mutate(virtual_cell_type = paste(data_type,cell_type,virtual_cell_type)) %>%
  mutate(task_idx = factor(task_idx,levels = c("Disease detection","Disease classification",
                                               "SF3B1mut detection","Anaemia classification") %>% rev)) %>% 
  subset(data_type == "Morphology + B.C.") %>%
  ggplot(aes(x = reorder(virtual_cell_type,-total_relevance),fill = y,y = task_idx)) +
  geom_tile() + 
  facet_grid( ~ data_type + cell_type,scales = "free_x") + 
  theme_pretty(base_size = 6) + 
  theme(strip.text = element_text(margin = ggplot2::margin()))  +
  scale_x_discrete(labels = function(x) str_match(x,"[0-9]+")) + 
  scale_fill_gradient2(low = "blue4",high = "red4",mid = "white",midpoint = 0,name = NULL) + 
  theme(axis.title.y = element_blank(),legend.key.width = unit(0.2,"cm")) + 
  xlab("Virtual cell type") + 
  ggsave(
    sprintf("figures/%s/mile-vice-mo-map.pdf",output_str),height=1.5,width=4.5)

# multiple objective correlations -----------------------------------------

tmp_data <- full_proportions_cell_type_mo %>%
  mutate(model_name = decode_model_name(task_idx)) %>%
  mutate(fine_class = ifelse(
    coarse_class == "MDS",
    ifelse(grepl('SF3B1',fine_class),"SF3B1-mutant MDS","non-SF3B1-mutant MDS"),
    as.character(fine_class))) %>%
  select(cell_type,virtual_cell_type,task_idx,slide_id,model_name,proportion,fine_class,data_type,class_difference) %>%
  distinct %>% 
  mutate(effect = class_difference * proportion) %>%
  select(-class_difference,-proportion) %>%
  group_by(virtual_cell_type,data_type,cell_type) %>% 
  summarise(average_effect_size = abs(median(effect))) %>%
  subset(abs(average_effect_size) > 0.05) %>% 
  mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>% 
  subset(data_type == "Morphology + B.C.")

tmp_data %>% 
  ggplot(aes(x = reorder(virtual_cell_type_fctr,average_effect_size),y = average_effect_size)) + 
  geom_point(size = 0.5) + 
  facet_grid(~ cell_type,scales = "free_x",space = "free") + 
  theme_pretty(base_size = 6) +
  xlab("Virtual cell type") + 
  ylab("Absolute average effect size") + 
  scale_x_discrete(labels = function(x) str_match(x,"[0-9]+$")) +
  scale_y_continuous(limits = c(0,1),expand = c(0,0,0.05,0)) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.key.size = unit(0,"cm"),
    strip.text = element_text(margin = ggplot2::margin()),
    strip.background = element_rect(fill = NA,color = NA)) + 
  ggsave(sprintf("figures/%s/mile-vice-multiple-objective-median-effect.pdf",output_str),height=1.5,width=3)

S <- sort(unique(tmp_data$virtual_cell_type_fctr))

curr_vcq <- subset_vcq(best_vcq_subset_mo,S)

subset(curr_vcq,cell_type == 'RBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-coefficients-rbc.pdf",output_str),height=3.3,width=4.5)

subset(curr_vcq,cell_type == 'WBC') %>% 
  coefficients_plot +
  ggsave(
    sprintf("figures/%s/mile-vice-mo-coefficients-wbc.pdf",output_str),height=4.2,width=4.5)


# slide-similarity heatmap ------------------------------------------------

all_heatmaps <- list()
for (M in grep("\\.bc",unique(full_proportions_cell_type$model_name),value = T)) {
  S <- full_proportions_cell_type %>%
    ungroup %>% 
    select(-direction) %>%
    subset(cell_type == "RBC") %>% 
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class) %>%
    subset(model_name == M) %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    select(-virtual_cell_type,-cell_type) %>% 
    spread(key = "VC",value = "proportion") %>%
    mutate(fine_class = conversion_list[[decode_model_name(M)]][as.character(fine_class)])
  # cos_sim <- as.matrix(S[,-c(1:5)] / sqrt(rowSums(S[,-c(1:5)] * S[,-c(1:5)])))
  # cos_sim <- cos_sim %*% t(cos_sim)
  # dm <- as.dist(1-cos_sim)
  # dm <- dist(S[,-c(1:5)])
  dm <- as.dist(1 - cor(t(S[,-c(1:5)])))
  hc <- hclust(dm, method = "ward.D")
  hc$labels <- S$slide_id
  hc <- dendextend::rotate(hc,arrange(S,fine_class)$slide_id)
  O <-  hc$order
  heatmap_data <- data.frame(
    slide_id_1 = S$slide_id,
    model_name = decode_model_name(S$model_name),
    fine_class = S$fine_class,
    data_type = S$data_type,
    as.data.frame(as.matrix(dm)))
  colnames(heatmap_data)[5:ncol(heatmap_data)] <- S$slide_id
  heatmap_data_long <- heatmap_data %>%
    gather(key = "slide_id_2",value = "distance",-slide_id_1,-model_name,-fine_class,-data_type) %>%
    mutate(slide_id_1 = factor(slide_id_1,levels = S$slide_id[O]),
           slide_id_2 = factor(slide_id_2,levels = S$slide_id[O]))
  CM <- heatmap_data_long %>%
    select(slide_id_1,fine_class) %>%
    distinct %>%
    ggplot(aes(x = slide_id_1,y = 1,fill = fine_class)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm"),
          axis.title = element_blank(),axis.line = element_blank(),
          legend.position = "top",legend.key.height = unit(0.2,"cm"),
          legend.margin = ggplot2::margin()) + 
    scale_fill_manual(values = fine_colours,name = NULL,
                      breaks = names(fine_colours)[names(fine_colours) %in% unique(S$fine_class)])
  HM <- heatmap_data_long %>% 
    ggplot(aes(x = slide_id_1,slide_id_2,fill = distance)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm")) + 
    xlab("Slide A") + 
    ylab("Slide B") + 
    scale_fill_distiller(direction = 1)
  all_heatmaps[[decode_model_name(M)]] <- plot_grid(CM,HM,rel_heights = c(0.2,1),align = "hv",axis = "lr",ncol = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-distance-rbc.pdf",output_str),height = 4.5,width = 5)

all_heatmaps <- list()
for (M in grep("\\.bc",unique(full_proportions_cell_type$model_name),value = T)) {
  S <- full_proportions_cell_type %>%
    ungroup %>% 
    select(-direction) %>%
    subset(cell_type == "WBC") %>% 
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class) %>%
    subset(model_name == M) %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    select(-virtual_cell_type,-cell_type) %>% 
    spread(key = "VC",value = "proportion") %>%
    mutate(fine_class = conversion_list[[decode_model_name(M)]][as.character(fine_class)])
  dm <- as.dist(1 - cor(t(S[,-c(1:5)])))
  hc <- hclust(dm, method = "ward.D")
  hc$labels <- S$slide_id
  hc <- dendextend::rotate(hc,arrange(S,fine_class)$slide_id)
  O <-  hc$order
  heatmap_data <- data.frame(
    slide_id_1 = S$slide_id,
    model_name = decode_model_name(S$model_name),
    fine_class = S$fine_class,
    data_type = S$data_type,
    as.data.frame(as.matrix(dm)))
  colnames(heatmap_data)[5:ncol(heatmap_data)] <- S$slide_id
  heatmap_data_long <- heatmap_data %>%
    gather(key = "slide_id_2",value = "distance",-slide_id_1,-model_name,-fine_class,-data_type) %>%
    mutate(slide_id_1 = factor(slide_id_1,levels = S$slide_id[O]),
           slide_id_2 = factor(slide_id_2,levels = S$slide_id[O]))
  CM <- heatmap_data_long %>%
    select(slide_id_1,fine_class) %>%
    distinct %>%
    ggplot(aes(x = slide_id_1,y = 1,fill = fine_class)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm"),
          axis.title = element_blank(),axis.line = element_blank(),
          legend.position = "top",legend.key.height = unit(0.2,"cm"),
          legend.margin = ggplot2::margin()) + 
    scale_fill_manual(values = fine_colours,name = NULL,
                      breaks = names(fine_colours)[names(fine_colours) %in% unique(S$fine_class)])
  HM <- heatmap_data_long %>% 
    ggplot(aes(x = slide_id_1,slide_id_2,fill = distance)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm")) + 
    xlab("Slide A") + 
    ylab("Slide B") + 
    scale_fill_distiller(direction = 1)
  all_heatmaps[[decode_model_name(M)]] <- plot_grid(CM,HM,rel_heights = c(0.2,1),align = "hv",axis = "lr",ncol = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-distance-wbc.pdf",output_str),height = 4.5,width = 5)

all_heatmaps <- list()
for (M in grep("\\.bc",unique(full_proportions_cell_type$model_name),value = T)) {
  S <- full_proportions_cell_type %>%
    ungroup %>% 
    select(-direction) %>%
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class) %>%
    subset(model_name == M) %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    select(-virtual_cell_type,-cell_type) %>% 
    spread(key = "VC",value = "proportion") %>%
    mutate(fine_class = conversion_list[[decode_model_name(M)]][as.character(fine_class)])
  dm <- as.dist(1 - cor(t(S[,-c(1:5)])))
  hc <- hclust(dm, method = "ward.D")
  hc$labels <- S$slide_id
  hc <- dendextend::rotate(hc,arrange(S,fine_class)$slide_id)
  O <-  hc$order
  heatmap_data <- data.frame(
    slide_id_1 = S$slide_id,
    model_name = decode_model_name(S$model_name),
    fine_class = S$fine_class,
    data_type = S$data_type,
    as.data.frame(as.matrix(dm)))
  colnames(heatmap_data)[5:ncol(heatmap_data)] <- S$slide_id
  heatmap_data_long <- heatmap_data %>%
    gather(key = "slide_id_2",value = "distance",-slide_id_1,-model_name,-fine_class,-data_type) %>%
    mutate(slide_id_1 = factor(slide_id_1,levels = S$slide_id[O]),
           slide_id_2 = factor(slide_id_2,levels = S$slide_id[O]))
  CM <- heatmap_data_long %>%
    select(slide_id_1,fine_class) %>%
    distinct %>%
    ggplot(aes(x = slide_id_1,y = 1,fill = fine_class)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm"),
          axis.title = element_blank(),axis.line = element_blank(),
          legend.position = "top",legend.key.height = unit(0.2,"cm"),
          legend.margin = ggplot2::margin()) + 
    scale_fill_manual(values = fine_colours,name = NULL,
                      breaks = names(fine_colours)[names(fine_colours) %in% unique(S$fine_class)])
  HM <- heatmap_data_long %>% 
    ggplot(aes(x = slide_id_1,slide_id_2,fill = distance)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm")) + 
    xlab("Slide A") + 
    ylab("Slide B") + 
    scale_fill_distiller(direction = 1)
  all_heatmaps[[decode_model_name(M)]] <- plot_grid(CM,HM,rel_heights = c(0.2,1),align = "hv",axis = "lr",ncol = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-distance.pdf",output_str),height = 4.5,width = 5)

all_heatmaps <- list()
for (M in grep("\\.bc",unique(full_proportions_cell_type$model_name),value = T)) {
  S <- full_proportions_cell_type_mo %>%
    ungroup %>% 
    select(-direction) %>%
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class) %>%
    subset(grepl("\\.bc",model_name)) %>%
    distinct %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    select(-virtual_cell_type,-cell_type) %>% 
    spread(key = "VC",value = "proportion") %>%
    mutate(fine_class = conversion_list[[decode_model_name(M)]][as.character(fine_class)]) %>%
    subset(!is.na(fine_class))
  dm <- as.dist(1 - cor(t(S[,-c(1:5)])))
  hc <- hclust(dm, method = "ward.D")
  hc$labels <- S$slide_id
  hc <- dendextend::rotate(hc,arrange(S,fine_class)$slide_id)
  O <-  hc$order
  heatmap_data <- data.frame(
    slide_id_1 = S$slide_id,
    model_name = decode_model_name(S$model_name),
    fine_class = S$fine_class,
    data_type = S$data_type,
    as.data.frame(as.matrix(dm)))
  colnames(heatmap_data)[5:ncol(heatmap_data)] <- S$slide_id
  heatmap_data_long <- heatmap_data %>%
    gather(key = "slide_id_2",value = "distance",-slide_id_1,-model_name,-fine_class,-data_type) %>%
    mutate(slide_id_1 = factor(slide_id_1,levels = S$slide_id[O]),
           slide_id_2 = factor(slide_id_2,levels = S$slide_id[O]))
  CM <- heatmap_data_long %>%
    select(slide_id_1,fine_class) %>%
    distinct %>%
    ggplot(aes(x = slide_id_1,y = 1,fill = fine_class)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm"),
          axis.title = element_blank(),axis.line = element_blank(),
          legend.position = "top",legend.key.height = unit(0.2,"cm"),
          legend.margin = ggplot2::margin()) + 
    scale_fill_manual(values = fine_colours,name = NULL,
                      breaks = names(fine_colours)[names(fine_colours) %in% unique(S$fine_class)])
  HM <- heatmap_data_long %>% 
    ggplot(aes(x = slide_id_1,slide_id_2,fill = distance)) + 
    geom_tile() + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),axis.ticks = element_blank(),
          legend.key.width = unit(0.2,"cm")) + 
    xlab("Slide A") + 
    ylab("Slide B") + 
    scale_fill_distiller(direction = 1)
  all_heatmaps[[decode_model_name(M)]] <- plot_grid(CM,HM,rel_heights = c(0.2,1),align = "hv",axis = "lr",ncol = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-distance-mo.pdf",output_str),height = 4.5,width = 5)

# virtual cell composition heatmap ----------------------------------------

label_str <- function(x) {
  as.vector(str_match(x,"[WR]BC[0-9]+"))
}

all_heatmaps <- list()
for (M in grep("\\.bc",unique(full_proportions_cell_type$model_name),value = T)) {
  S <- full_proportions_cell_type %>%
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class) %>%
    subset(model_name == M) %>%
    mutate(cell_type = factor(cell_type,levels = c("RBC","WBC"))) %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    mutate(VC = reorder(VC,as.numeric(cell_type)*100 + as.numeric(virtual_cell_type))) %>%
    mutate(fine_class = conversion_list[[decode_model_name(M)]][as.character(fine_class)]) %>%
    group_by(VC,cell_type) %>%
    filter(any(proportion > 0.05)) %>% 
    mutate(VC_size = sum(proportion)) %>%
    ungroup %>%
    mutate(VC = reorder(VC,VC_size)) %>%
    group_by(slide_id) %>%
    mutate(slide_order = as.numeric(VC[which.max(proportion)]) + max(proportion))
  l <- S %>%
    select(slide_id,fine_class,slide_order) %>%
    distinct
  VC_order <- c("Class",as.character(unique(S$VC[order(S$VC_size)])))
  PH <- data.frame(NA,NA,NA,NA,NA,M,NA,NA,NA,NA,"Class",NA,NA)
  colnames(PH) <- colnames(S)
  CP <- S %>%
    select(slide_id,proportion,fine_class,VC,VC_size,model_name) %>%
    distinct %>%
    group_by(VC,VC_size,fine_class,model_name) %>% 
    summarise(M = median(proportion)) %>%
    ggplot(aes(x = M,y = reorder(VC,VC_size),colour = fine_class)) +
    geom_point(size = 1,alpha = 0.5) + 
    scale_colour_manual(values = fine_colours,name = NULL,breaks = unique(S$fine_class)) + 
    theme_pretty(base_size = 6) + 
    theme(axis.text.y = element_blank(),axis.title.y = element_blank(),
          legend.position = "bottom",legend.key.size = unit(0,"cm"),
          panel.grid.major.y = element_line(colour = "grey90"),
          strip.text = element_blank()) + 
    xlab("Median proportion") + 
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8)) + 
    guides(colour = guide_legend(nrow = 2)) + 
    facet_wrap(~ decode_model_name(model_name))
  HM <- S %>% 
    rbind.data.frame(PH) %>% 
    ggplot(aes(x = reorder(slide_id,-slide_order),y = factor(VC,VC_order),
               colour = fine_class,
               fill = proportion)) + 
    geom_tile(colour=NA) + 
    geom_point(aes(x = reorder(slide_id,-slide_order),y = factor("Class",VC_order),
                   colour = fine_class),fill=NA,
               data = l,alpha = 0.5,size = 0.25) +
    theme_pretty(base_size = 6) + 
    theme(axis.text.x = element_blank(),legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.3,"cm"),axis.ticks.x = element_blank(),
          legend.position = "bottom") + 
    xlab("Slide") + 
    ylab("Virtual cell type") + 
    scale_fill_distiller(direction = 1,name = "Proportion",breaks = c(0,0.2,0.4,0.6,0.8)) + 
    facet_wrap(~ decode_model_name(model_name)) + 
    scale_colour_manual(values = fine_colours,guide=F)
  all_heatmaps[[decode_model_name(M)]] <- plot_grid(HM,CP,rel_widths = c(1,0.4),align = "h",axis = "lr",nrow = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-vc.pdf",output_str),height = 7,width = 5.5)

all_heatmaps <- list()
for (curr_task in unique(full_proportions_cell_type_mo$task_idx)) {
  S <- full_proportions_cell_type_mo %>%
    select(slide_id,virtual_cell_type,cell_type,model_name,proportion,max_vc,data_type,fine_class,task_idx) %>%
    subset(grepl("\\.bc",model_name) & task_idx == curr_task) %>%
    mutate(cell_type = factor(cell_type,levels = c("RBC","WBC"))) %>%
    mutate(fine_class = conversion_list[[decode_model_name(curr_task)]][as.character(fine_class)]) %>%
    subset(!is.na(fine_class)) %>%
    mutate(VC = paste0(cell_type,virtual_cell_type)) %>%
    mutate(VC = reorder(VC,as.numeric(cell_type)*100 + as.numeric(virtual_cell_type))) %>%
    group_by(VC,cell_type) %>%
    filter(any(proportion > 0.1)) %>% 
    mutate(VC_size = sum(proportion)) %>%
    ungroup %>%
    mutate(VC = reorder(VC,VC_size)) %>%
    group_by(slide_id) %>%
    mutate(slide_order = as.numeric(VC[which.max(proportion)]) + max(proportion)) 
  l <- S %>%
    select(slide_id,fine_class,slide_order) %>%
    distinct
  VC_order <- c("Class",as.character(unique(S$VC[order(S$VC_size)])))
  PH <- data.frame(NA,NA,NA,NA,curr_task,NA,NA,NA,NA,NA,"Class",NA,NA)
  colnames(PH) <- colnames(S)
  CP <- S %>%
    select(slide_id,proportion,fine_class,VC,VC_size,model_name) %>%
    distinct %>%
    group_by(VC,VC_size,fine_class,model_name) %>% 
    summarise(M = median(proportion)) %>%
    ggplot(aes(x = M,y = reorder(VC,VC_size),colour = fine_class)) +
    geom_point(size = 1,alpha = 0.5) + 
    scale_colour_manual(values = fine_colours,name = NULL,breaks = unique(S$fine_class)) + 
    theme_pretty(base_size = 6) + 
    theme(axis.text.y = element_blank(),axis.title.y = element_blank(),
          legend.position = "bottom",legend.key.size = unit(0,"cm"),
          panel.grid.major.y = element_line(colour = "grey90"),
          strip.text = element_blank()) + 
    xlab("Median proportion") + 
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8)) + 
    guides(colour = guide_legend(nrow = 2))
  HM <- S %>% 
    rbind.data.frame(PH) %>% 
    ggplot(aes(x = reorder(slide_id,-slide_order),y = factor(VC,VC_order),
               colour = fine_class,
               fill = proportion)) + 
    geom_tile(colour=NA) + 
    geom_point(aes(x = reorder(slide_id,-slide_order),y = factor("Class",VC_order),
                   colour = fine_class),fill=NA,
               data = l,alpha = 0.3,size = 0.25) +
    theme_pretty(base_size = 6) + 
    theme(axis.text.x = element_blank(),legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.3,"cm"),axis.ticks.x = element_blank(),
          legend.position = "bottom") + 
    xlab("Slide") + 
    ylab("Virtual cell type") + 
    scale_fill_distiller(direction = 1,name = "Proportion",breaks = c(0,0.2,0.4,0.6,0.8)) + 
    scale_colour_manual(values = fine_colours,guide=F)
  all_heatmaps[[curr_task]] <- plot_grid(HM,CP,rel_widths = c(1,0.4),align = "h",axis = "lr",nrow = 1)
}

plot_grid(plotlist = all_heatmaps,ncol = 2) + 
  ggsave(sprintf("figures/%s/mile-vice-heatmap-vc-mo.pdf",output_str),height = 6.5,width = 5.5)

# proportion correlations (validation) ------------------------------------

conversion_list <- list(
  `Disease detection` = c(
    Normal = "Normal",`U2AF1-mutant` = "Disease",`RUNX1-mutant` = "Disease",
    `SRSF2-mutant` = "Disease",`SF3B1 mutant` = "Disease",`Iron deficiency` = "Disease",
    `Megaloblastic` = "Disease"),
  `Disease classification` = c(
    `U2AF1-mutant` = "MDS",`RUNX1-mutant` = "MDS",
    `SRSF2-mutant` = "MDS",`SF3B1 mutant` = "MDS",`Iron deficiency` = "Anaemia",
    `Megaloblastic` = "Anaemia"),
  `SF3B1mut detection` = c(
    `SRSF2-mutant` = "SRSF2-mutant",`SF3B1 mutant` = "SF3B1-mutant",
    `SRSF2 mutant` = "SRSF2-mutant",`SF3B1-mutant` = "SF3B1-mutant"),
  `Anaemia classification` = c(
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"))

conversion_list_order <- c(
  "Normal","Disease","MDS","Anaemia","SF3B1-mutant","SRSF2-mutant",
  "Iron deficiency","Megaloblastic")

all_cell_proportions_val <- list()

for (model in unique(best_layers_subset$model_name)) {
  path <- sprintf("../mile-vice/cell-proportions/adden_2/%s.csv",model)
  proportions_val <- read_proportions(path)
  proportions_val$model_name <- model
  all_cell_proportions_val[[model]] <- proportions_val
}

cell_proportions_val_df <- do.call(rbind,all_cell_proportions_val) %>%
  merge(all_conditions,by = "slide_id")

cell_prop_summary <- cell_proportions_df %>%
  rowwise() %>%
  mutate(fine_class = conversion_list[[
    decode_model_name(model_name)]][as.character(fine_class)]) %>% 
  na.omit() %>% 
  group_by(model_name,fine_class,cell_type,virtual_cell_type,max_vc) %>%
  summarise(q05 = quantile(proportion,0.05),q50 = quantile(proportion,0.5),
            q95 = quantile(proportion,0.95))

cell_prop_summary_val <- cell_proportions_val_df %>%
  mutate(fine_class = ifelse(fine_class=="SF3B1-mutant","SF3B1 mutant",
                             as.character(fine_class))) %>%
  rowwise() %>%
  mutate(fine_class = conversion_list[[decode_model_name(model_name)]][fine_class]) %>% 
  na.omit() %>%
  group_by(model_name,fine_class,cell_type,virtual_cell_type,max_vc) %>%
  summarise(q05_val = quantile(proportion,0.05),q50_val = quantile(proportion,0.5),
            q95_val = quantile(proportion,0.95))

comparison_df <- merge(
  cell_prop_summary,cell_prop_summary_val,
  by = c("model_name","fine_class","virtual_cell_type","cell_type")) %>%
  group_by(cell_type,model_name) %>%
  mutate(label = ifelse(virtual_cell_type %in% virtual_cell_type[abs(q50-q50_val)>0.1],
                        virtual_cell_type,NA)) %>%
  mutate(fine_class = factor(fine_class,conversion_list_order))

comparison_df %>%  
  mutate(data_type = ifelse(grepl('\\.bc',model_name),'Morphology + B.C.',"Morphology"),
         model_name = decode_model_name(model_name)) %>%
  group_by(cell_type,data_type) %>%
  mutate(p.val.global = cor.test(q50,q50_val)$p.val,
         R.global = cor.test(q50,q50_val)$estimate) %>%
  group_by(cell_type,model_name,data_type,R.global,p.val.global) %>%
  summarise(p.val = cor.test(q50,q50_val)$p.val,
            R = cor.test(q50,q50_val)$estimate) %>%
  ggplot(aes(x = R,y = -log(p.val),colour = model_name)) + 
  geom_point(aes(shape = data_type,alpha = p.val < (0.05 / length(R))),
             size = 0.5) + 
  geom_hline(aes(yintercept = -log(0.05 / length(R))),size = 0.25,linetype = 2) + 
  geom_point(aes(x = R.global,y = -log(p.val.global),shape = data_type),
             colour = "black",size = 1) +
  facet_wrap(~ cell_type) + 
  scale_alpha_manual(values = c(0.3,1),name = "Stat. significant",
                     labels = c("No","Yes")) + 
  theme_pretty(base_size = 6) + 
  ylab("-log(p-value)") + 
  theme(legend.key.size = unit(0.1,"cm"),legend.margin = ggplot2::margin()) + 
  scale_shape(name = "Dataset") + 
  scale_colour_aaas(name = NULL) + 
  ggsave(sprintf("figures/%s/mile-vice-proportion-correlations.pdf",output_str),height = 1.5,width = 4)

comparison_df %>%
  mutate(data_type = ifelse(grepl('\\.bc',model_name),'Morphology + B.C.',"Morphology"),
         model_name = decode_model_name(model_name)) %>%
  subset(data_type == "Morphology") %>% 
  ggplot(aes(x = q50,y = q50_val, colour = cell_type)) + 
  geom_abline(slope = 1,size = 0.25,linetype = 2) +
  geom_point(aes(shape = fine_class),size = 0.5,alpha = 0.5) + 
  geom_linerange(aes(ymin = q05_val,ymax = q95_val),size = 0.25,
                 alpha = 0.5) +
  geom_linerange(aes(xmin = q05,xmax = q95),size = 0.25,
                 alpha = 0.5) +
  geom_text_repel(aes(label = label),colour = "black",size = 2,
                  segment.size = 0.25,min.segment.length = 0.1) + 
  facet_grid(cell_type ~ model_name) + 
  theme_pretty(base_size = 6) +
  scale_shape_manual(name = NULL,values = c(1,2,6,7,10,11,15,17)) + 
  xlab("Proportion distribution in training data") + 
  ylab("Proportion distribution in validation data") + 
  scale_colour_manual(values = c("red4","darkorchid"),guide = F) + 
  theme(legend.key.size = unit(0,"cm")) +
  ggsave(sprintf("figures/%s/mile-vice-proportion-correlations-morph.pdf",output_str),
         width = 5,height = 2.5)

comparison_df %>%
  mutate(data_type = ifelse(grepl('\\.bc',model_name),'Morphology + B.C.',"Morphology"),
         model_name = decode_model_name(model_name)) %>%
  subset(data_type == "Morphology" & model_name == "Disease detection") %>% 
  ggplot(aes(x = q50,y = q50_val, colour = cell_type)) + 
  geom_abline(slope = 1,size = 0.25,linetype = 2) +
  geom_point(aes(shape = fine_class),size = 1,alpha = 0.5) + 
  geom_linerange(aes(ymin = q05_val,ymax = q95_val),size = 0.25,
                 alpha = 0.5) +
  geom_linerange(aes(xmin = q05,xmax = q95),size = 0.25,
                 alpha = 0.5) +
  geom_text_repel(aes(label = label),colour = "black",size = 2,
                  segment.size = 0.25,min.segment.length = 0.1) + 
  facet_grid(. ~ cell_type) + 
  theme_pretty(base_size = 6) +
  scale_shape_manual(name = NULL,values = c(1,2,6,7,10,11,15,17)) + 
  xlab("Proportion distribution in training data") + 
  ylab("Proportion distribution in validation data") + 
  scale_colour_manual(values = c("red4","darkorchid"),guide = F) + 
  theme(legend.key.size = unit(0,"cm")) +
  ggsave(sprintf("figures/%s/mile-vice-proportion-correlations-morph-disease-detection.pdf",output_str),
         width = 4,height = 2)

comparison_df %>%
  mutate(data_type = ifelse(grepl('\\.bc',model_name),'Morphology + B.C.',"Morphology"),
         model_name = decode_model_name(model_name)) %>%
  subset(data_type == "Morphology + B.C.") %>% 
  ggplot(aes(x = q50,y = q50_val, colour = cell_type)) + 
  geom_abline(slope = 1,size = 0.25,linetype = 2) +
  geom_point(aes(shape = fine_class),size = 0.5,alpha = 0.5) + 
  geom_linerange(aes(ymin = q05_val,ymax = q95_val),size = 0.25,
                 alpha = 0.5) +
  geom_linerange(aes(xmin = q05,xmax = q95),size = 0.25,
                 alpha = 0.5) +
  geom_text_repel(aes(label = label),colour = "black",size = 2,
                  segment.size = 0.25,min.segment.length = 0.1) + 
  facet_grid(cell_type ~ model_name) + 
  theme_pretty(base_size = 6) +
  scale_shape_manual(name = NULL,values = c(1,2,6,7,10,11,15,17)) + 
  xlab("Proportion distribution in training data") + 
  ylab("Proportion distribution in validation data") + 
  scale_colour_manual(values = c("red4","darkorchid"),guide = F) + 
  theme(legend.key.size = unit(0,"cm")) +
  ggsave(sprintf("figures/%s/mile-vice-proportion-correlations-morph-bc.pdf",output_str),
         width = 5,height = 2.5)

all_cell_proportions_val <- list()

for (model in unique(best_layers_subset_mo$model_name)) {
  path <- sprintf("../mile-vice/cell-proportions/adden_2/%s.csv",model)
  proportions_val <- read_proportions(path)
  proportions_val$model_name <- model
  all_cell_proportions_val[[model]] <- proportions_val
}

cell_proportions_val_df <- do.call(rbind,all_cell_proportions_val) %>%
  merge(all_conditions,by = "slide_id")

cell_prop_summary <- cell_proportions_mo_df %>%
  rowwise() %>%
  mutate(fine_class = conversion_list$`Disease detection`[as.character(fine_class)]) %>% 
  na.omit() %>% 
  group_by(model_name,fine_class,cell_type,virtual_cell_type,max_vc) %>%
  summarise(q05 = quantile(proportion,0.05),q50 = quantile(proportion,0.5),
            q95 = quantile(proportion,0.95))

cell_prop_summary_val <- cell_proportions_val_df %>%
  mutate(fine_class = ifelse(fine_class=="SF3B1-mutant","SF3B1 mutant",
                             as.character(fine_class))) %>%
  rowwise() %>%
  mutate(fine_class = conversion_list$`Disease detection`[fine_class]) %>% 
  na.omit() %>%
  group_by(model_name,fine_class,cell_type,virtual_cell_type,max_vc) %>%
  summarise(q05_val = quantile(proportion,0.05),q50_val = quantile(proportion,0.5),
            q95_val = quantile(proportion,0.95))

comparison_df <- merge(
  cell_prop_summary,cell_prop_summary_val,
  by = c("model_name","fine_class","virtual_cell_type","cell_type")) %>%
  group_by(cell_type,model_name) %>%
  mutate(label = ifelse(virtual_cell_type %in% virtual_cell_type[abs(q50-q50_val)>0.1],
                        virtual_cell_type,NA)) %>%
  mutate(fine_class = factor(fine_class,conversion_list_order))

Max <- max(comparison_df$q50,comparison_df$q50_val)

comparison_df %>%
  mutate(data_type = ifelse(grepl('\\.bc',model_name),'Morphology + B.C.',"Morphology"),
         model_name = decode_model_name(model_name)) %>%
  subset(data_type == "Morphology + B.C.") %>% 
  ggplot(aes(x = q50,y = q50_val, colour = cell_type)) + 
  geom_abline(slope = 1,size = 0.25,linetype = 2) +
  geom_point(aes(shape = fine_class),size = 1,alpha = 0.5) + 
  geom_linerange(aes(ymin = q05_val,ymax = q95_val),size = 0.25,
                 alpha = 0.5) +
  geom_linerange(aes(xmin = q05,xmax = q95),size = 0.25,
                 alpha = 0.5) +
  geom_text_repel(aes(label = label),colour = "black",size = 2,
                  segment.size = 0.25,min.segment.length = 0.1) + 
  facet_grid(. ~ cell_type) + 
  theme_pretty(base_size = 6) +
  scale_shape_manual(name = NULL,values = c(1,2,6,7,10,11,15,17)) + 
  xlab("Proportion distribution in training data") + 
  ylab("Proportion distribution in validation data") + 
  coord_cartesian(xlim = c(NA,Max),ylim = c(NA,Max)) +
  scale_colour_manual(values = c("red4","darkorchid"),guide = F) + 
  theme(legend.key.size = unit(0,"cm")) +
  ggsave(sprintf("figures/%s/mile-vice-proportion-correlations-mo-disease-detection.pdf",output_str),
         width = 4,height = 2)
