# setup -------------------------------------------------------------------

source("function-library.R")

args <- commandArgs(trailingOnly=TRUE)
output_str <- ifelse(length(args)>0,args[1],"subset")
dir.create(paste('figures',output_str,sep = "/"),showWarnings = F)

library(ggrepel)
library(umap)
library(cowplot)
library(glmnet)
library(MASS)
library(ggrepel)

select <- dplyr::select

multi_objective_matching <- c("anemia_binary","binary","mds_binary","disease_binary")

conversion_list <- list(
  `Disease detection` = c(
    Normal = "Normal",`SF3B1-wildtype` = "Disease",
    `SF3B1-mutant` = "Disease",`SRSF2-mutant` = "Disease",
    `RUNX1-mutant` = "Disease",`U2AF1-mutant` = "Disease",
    `Iron deficiency` = "Disease",
    `Megaloblastic` = "Disease"),
  `Disease classification` = c(
    `SF3B1-wildtype` = "MDS",`SF3B1-mutant` = "MDS",
    `SRSF2-mutant` = "MDS",`RUNX1-mutant` = "MDS",`U2AF1-mutant` = "MDS",
    `Iron deficiency` = "Anaemia",`Megaloblastic` = "Anaemia"),
  `SF3B1mut detection` = c(
    `SF3B1-mutant` = "SF3B1-mutant",`SRSF2-mutant` = "SF3B1-wildtype",
    `RUNX1-mutant` = "SF3B1-wildtype",`U2AF1-mutant` = "SF3B1-wildtype",
    `SF3B1-wildtype` = "SF3B1-wildtype"),
  `Anaemia classification` = c(
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"),
  `Multi-objective` = c(
    Normal = "Normal",`SF3B1-mutant` = "SF3B1-mutant",`SRSF2-mutant` = "SF3B1-wildtype",
    `RUNX1-mutant` = "SF3B1-wildtype",`U2AF1-mutant` = "SF3B1-wildtype",
    `SF3B1-wildtype` = "SF3B1-wildtype",
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

read_vcq_layer <- function(path,wbc_subset_path,rbc_subset_path) {
  wbc_feature_subset <- unlist(read_csv(wbc_subset_path,col_names = F))
  rbc_feature_subset <- unlist(read_csv(rbc_subset_path,col_names = F))
  df <- read_csv(path,col_names=F)
  colnames(df) <- c("model_name","fold","feature","cell_type","virtual_cell_type","coefficient")
  df$max_vc <- as.numeric(str_match(df$model_name,'[0-9]+'))
  
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

read_consensus <- function(path) {
  read_csv(path,col_names = F) %>%
    subset(grepl("ID_",X3)) %>%
    mutate(cell_type = ifelse(X3 == "ID_0","WBC","RBC"),
           cells = X4) %>%
    select(cell_type,cells) %>%
    group_by(cell_type) %>%
    summarise(virtual_cell_type = 1+as.numeric(unlist(str_split(cells,":")))) %>%
    return
}

coefficients_plot <- function(x) {
  x %>%
    ggplot(aes(x = virtual_cell_type_fctr,
               y = feature,
               fill = coefficient)) + 
    geom_tile() + 
    facet_wrap( ~ sprintf('%s (%s)',data_type,cell_type),scales = "free") + 
    theme_pretty(base_size = 6) + 
    theme(legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) +
    xlab("Virtual cell type") + 
    ylab("Feature") +
    scale_fill_gradient2(low = 'blue4',mid = "white",high = "red4",
                         name = "Coefficient") + 
    theme(legend.key.height = unit(0.1,"cm"),
          strip.text = element_text(margin = margin(t = 1,b = 2)),
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

relevant_classes <- function(x) {
  rc <- c(
    "Normal","SF3B1-mutant","Iron deficiency","U2AF1-mutant",
    "SRSF2-mutant","Megaloblastic","RUNX1-mutant")
  if (x == "binary") {
    return(rc)
  } else if (x == "disease_binary") {
    return(Filter(function(j) j != "Normal",rc))
  } else if (x == "mds_binary") {
    return(c("SF3B1-mutant","U2AF1-mutant","SRSF2-mutant","RUNX1-mutant"))
  } else if (x == "anemia_binary") {
    return(c("Megaloblastic","Iron deficiency"))
  }
}

decode_model_name <- function(model_names) {
  o <- ifelse(
    grepl("anemia_binary",model_names),"Anaemia classification",
    ifelse(grepl("mds_binary",model_names),"SF3B1mut detection",
           ifelse(grepl("disease_binary",model_names),
                  "Disease classification",
                  "Disease detection"))) %>%
    as.factor()
  levels(o) <- c("Disease detection","Disease classification",
                 "SF3B1mut detection","Anaemia classification")
  return(o)
}

# data loading and processing ---------------------------------------------

rbc_counts <-  rbind(
  read_csv(
    "datasets/rbc_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "MLLC"),
  read_csv(
    "datasets/rbc_adden_1_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double()))  %>%
    cbind(dataset = "AC1"),
  read_csv(
    "datasets/rbc_adden_2_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "AC2")
) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  group_by(slide_id,dataset) %>%
  summarise(n_rbc = sum(counts))

wbc_counts <-  rbind(
  read_csv(
    "datasets/wbc_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "MLLC"),
  read_csv(
    "datasets/wbc_adden_1_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double()))  %>%
    cbind(dataset = "AC1"),
  read_csv(
    "datasets/wbc_adden_2_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "AC2")
) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  group_by(slide_id,dataset) %>%
  summarise(n_wbc = sum(counts))

# select slides with > 50 RBC and > 50 WBC
representative_slides <- merge(rbc_counts,wbc_counts,by = c("slide_id","dataset")) %>% 
  subset(n_rbc > 50 & n_wbc > 50) %>% 
  select(slide_id) %>%
  unlist

so_layers <- paste(Filter(
  nchar,c("../mil-comori/best_models/best_layers",output_str)),collapse="_")
mo_layers <- paste(Filter(
  nchar,c("../mil-comori/best_models/best_layers_mo",output_str)),collapse="_")
vcq_layers <- paste(Filter(
  nchar,c("../mil-comori/best_models/best_vcq_layers",output_str)),collapse="_")
vcq_layers_mo <- paste(Filter(
  nchar,c("../mil-comori/best_models/best_vcq_layers_mo",output_str)),collapse="_")
wbc_subset_path <- NULL
rbc_subset_path <- NULL
if (grepl("subset",output_str)){
  wbc_subset_path <- "data_output/wbc_feature_subset"
  rbc_subset_path <- "data_output/rbc_feature_subset"
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
all_consensus <- list()
all_consensus_mo <- list()

for (model in unique(best_layers_subset$model_name)) {
  task_name <- str_split(model,"\\.")[[1]][2]
  path <- sprintf("../mil-comori/cell-proportions/mll/%s.csv",model)
  proportions <- read_proportions(path)
  proportions$model_name <- model
  all_cell_proportions[[model]] <- proportions
  path <- sprintf("../mil-comori/ev-scores-consensus/%s.%s.csv",model,task_name)
  consensus <- read_consensus(path)
  consensus$model_name <- model
  consensus$is_consensual <- T
  all_consensus[[model]] <- as.data.frame(consensus)
}

blmt <- best_layers_subset_mo %>%
  select(model_name,task_idx) %>%
  distinct

for (i in 1:nrow(blmt)) {
  model <- blmt$model_name[i]
  path <- sprintf("../mil-comori/cell-proportions/mll/%s.csv",model)
  proportions <- read_proportions(path) %>%
    distinct
  proportions$model_name <- model
  all_cell_proportions_mo[[model]] <- proportions
}

for (i in 1:nrow(blmt)) {
  model <- blmt$model_name[i]
  task_name <- multi_objective_matching[blmt$task_idx[i]+1]
  path <- sprintf("../mil-comori/ev-scores-consensus/mo_%s.%s.csv",model,task_name)
  consensus <- read_consensus(path)
  consensus$model_name <- model
  consensus$task_idx <- blmt$task_idx[i]
  consensus$is_consensual <- T
  all_consensus_mo[[paste(model,task_name)]] <- as.data.frame(consensus)
}

cell_proportions_df <- do.call(rbind,all_cell_proportions) %>%
  merge(all_conditions,by = "slide_id")
cell_proportions_mo_df <- do.call(rbind,all_cell_proportions_mo) %>%
  merge(all_conditions,by = "slide_id")
all_consensus_df <- do.call(rbind,all_consensus)
all_consensus_mo_df <- do.call(rbind,all_consensus_mo)

full_proportions_cell_type <- merge(
  cell_proportions_df,best_layers_subset,
  by=c("cell_type","virtual_cell_type","max_vc","model_name")) %>%
  merge(all_consensus_df,by = c("cell_type","virtual_cell_type","model_name"),all.x = T) %>%
  mutate(is_consensual = ifelse(is.na(is_consensual),F,is_consensual))
full_proportions_cell_type_mo <- merge(
  cell_proportions_mo_df,best_layers_subset_mo,
  by=c("cell_type","virtual_cell_type","max_vc","model_name")) %>%
  merge(all_consensus_mo_df,by = c("cell_type","virtual_cell_type","model_name","task_idx"),all.x = T) %>%
  mutate(is_consensual = ifelse(is.na(is_consensual),F,is_consensual))
full_proportions_cell_type_mo$task_idx <- multi_objective_matching[full_proportions_cell_type_mo$task_idx+1]

rbc_all_cells_summaries <- read_csv(
  "datasets/rbc_summaries.csv",
  col_names = c("slide_id",features_all,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>%
  as_tibble() %>% 
  mutate(fine_class = ifelse(coarse_class == "MDS",
                             ifelse(grepl("SF3B1",fine_class),
                                    "SF3B1-mutant",
                                    "SF3B1-wildtype"),
                             as.character(fine_class))) %>%
  mutate(fine_class = factor(
    as.character(fine_class),
    levels = fine_simple_levels)) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  subset(slide_id %in% representative_slides) 
rbc_means <- rbc_all_cells_summaries %>%
  subset(f == "mean") %>%
  mutate_if(is.numeric,scale) 
colnames(rbc_means)[2:ncol(rbc_means)] <- paste(
  colnames(rbc_means)[2:ncol(rbc_means)],"mean",sep = "_")
rbc_variances <- rbc_all_cells_summaries %>% 
  subset(f == "variance") %>%
  mutate_if(is.numeric,scale) 

wbc_all_cells_summaries <- read_csv(
  "datasets/wbc_summaries.csv",
  col_names = c("slide_id",features_all,features_nuclear,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all,features_nuclear)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>% 
  as_tibble() %>% 
  mutate(fine_class = ifelse(coarse_class == "MDS",
                             ifelse(grepl("SF3B1",fine_class),
                                    "SF3B1-mutant",
                                    "SF3B1-wildtype"),
                             as.character(fine_class))) %>%
  mutate(fine_class = factor(
    as.character(fine_class),
    levels = fine_simple_levels)) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  subset(slide_id %in% representative_slides)
wbc_means <- wbc_all_cells_summaries %>%
  subset(f == "mean") %>%
  mutate_if(is.numeric,scale) 
colnames(wbc_means)[2:ncol(wbc_means)] <- paste(
  colnames(wbc_means)[2:ncol(wbc_means)],"mean",sep = "_")
wbc_variances <- wbc_all_cells_summaries  %>%
  subset(f == "variance") %>%
  mutate_if(is.numeric,scale) 

vc_proportions_wide <- full_proportions_cell_type_mo %>%
  subset(data_type == "Morphology + B.C." & is_consensual == T) %>% 
  mutate(tmp = as.numeric(vct_conversion(virtual_cell_type,cell_type))) %>% 
  group_by(cell_type) %>%
  mutate(tmp = ifelse(is.na(tmp),max(tmp,na.rm=T) + as.numeric(as.factor(virtual_cell_type)),tmp)) %>%
  ungroup %>%
  mutate(virtual_cell_type = paste0(cell_type,tmp)) %>%
  select(-tmp) %>%
  select(-cell_type,-model_name,-`0`,-`1`,-class_difference,-fold) %>%
  distinct %>% 
  group_by(virtual_cell_type,task_idx) %>%
  #mutate(proportion = scale(proportion)) %>%
  spread(key = "virtual_cell_type",value = "proportion") 

rvc_proportions_wide <- vc_proportions_wide[,!grepl("WBC",colnames(vc_proportions_wide))] %>%
  mutate(fine_class = ifelse(
    coarse_class == "MDS",
    ifelse(fine_class == "SF3B1-mutant","SF3B1-mutant","SF3B1-wildtype"),
    as.character(fine_class)))
wvc_proportions_wide <- vc_proportions_wide[,!grepl("RBC",colnames(vc_proportions_wide))] %>%
  mutate(fine_class = ifelse(
    coarse_class == "MDS",
    ifelse(fine_class == "SF3B1-mutant","SF3B1-mutant","SF3B1-wildtype"),
    as.character(fine_class)))

glmnet_models <- readRDS("data_output/glmnet_models.rds")

task_conversion <- c(
  disease_detection = "binary",disease_classification = "disease_binary",
  sf3b1 = "mds_binary",anaemia_classification = "anemia_binary")

glmnet_coefficients <- names(glmnet_models$full) %>%
  lapply(
    function(n){
      x <- glmnet_models$full[[n]]
      best_model_idx<-which.max(sapply(x,function(y) y$metrics$Validation$AUC$auc))
      C <- as.matrix(coef(x[[best_model_idx]]$model))
      o <- data.frame(features = rownames(C),values = as.vector(C),task = n)
      return(o)
    }
  ) %>%
  do.call(what = rbind) %>%
  subset(grepl("WBC|RBC",features)) %>% 
  mutate(moment = ifelse(grepl("\\.vars\\.",features),"Variances","Means"),
         cell_type = ifelse(grepl("WBC",features),"WBC","RBC"),
         feature_name = str_match(features,"^[a-zA-Z0-9_]+\\.")) %>% 
  mutate(feature_name = gsub("\\.","",feature_name)) %>% 
  subset(values != 0) %>%
  subset(!is.na(feature_name)) %>%
  subset(moment == "Variances") %>%
  mutate(task = task_conversion[task])

rbc_merged <- merge(
  subset(rbc_all_cells_summaries,f=="variance"),
  subset(rbc_all_cells_summaries,f=="mean"),by = c("slide_id","fine_class","coarse_class"),
  suffixes = c(".variance",".mean")) %>%
  merge(rvc_proportions_wide,by = c("slide_id","fine_class"))
wbc_merged <- merge(
  subset(wbc_all_cells_summaries,f=="variance"),
  subset(wbc_all_cells_summaries,f=="mean"),by = c("slide_id","fine_class","coarse_class"),
  suffixes = c(".variance",".mean")) %>%
  merge(wvc_proportions_wide,by = c("slide_id","fine_class"))

# associations between vct prop and variance + examples -------------------


all_correlations <- list()
all_coefficient_associations <- list() 
for (tsk in unique(glmnet_coefficients$task)) {
  cl <- relevant_classes(tsk)
  s <- glmnet_coefficients %>% 
    subset(task == tsk)
  rbc_features <- unique(s$feature_name[s$cell_type == "RBC"])
  rbc_variance_subset <- rbc_variances[
    ,colnames(rbc_variances) %in% c("slide_id",rbc_features)]

  wbc_features <- unique(s$feature_name[s$cell_type == "WBC"])
  wbc_variance_subset <- wbc_variances[
    ,colnames(wbc_variances) %in% c("slide_id",wbc_features)]

  sub_rvc <- rvc_proportions_wide %>%
    subset(task_idx == tsk) %>% 
    merge(rbc_variance_subset,by = "slide_id") %>%
    select_if(function(x) any(!is.na(x))) %>%
    subset(fine_class %in% cl) %>%
    mutate_if(is.numeric,scale)
  sub_wvc <- wvc_proportions_wide %>%
    subset(task_idx == tsk) %>%
    merge(wbc_variance_subset,by = "slide_id") %>%
    select_if(function(x) any(!is.na(x))) %>%
    subset(fine_class %in% cl) %>%
    mutate_if(is.numeric,scale)
  
  rvc <- grep("RBC",colnames(sub_rvc),value=T)
  wvc <- grep("WBC",colnames(sub_wvc),value=T)
  
  rbc_correlations <- cor(sub_rvc[,rbc_features],sub_rvc[,rvc]) %>%
    as.data.frame()
  rbc_correlations$features <- rownames(rbc_correlations)
  rbc_correlations <- gather(rbc_correlations,"key","R",-features) %>%
    mutate(n = nrow(sub_rvc)) %>%
    mutate(t_val = R*sqrt(n-2)/sqrt(1-R^2)) %>%
    mutate(p.val = pt(0.05,n-2,abs(t_val))) %>%
    mutate(cell_type = "RBC") %>%
    mutate(task = tsk)

  wbc_correlations <- cor(sub_wvc[,wbc_features],sub_wvc[,wvc]) %>%
    as.data.frame()
  wbc_correlations$features <- rownames(wbc_correlations)
  wbc_correlations <- gather(wbc_correlations,"key","R",-features) %>%
    mutate(n = nrow(sub_rvc)) %>%
    mutate(t_val = R*sqrt(n-2)/sqrt(1-R^2)) %>%
    mutate(p.val = pt(0.05,n-2,abs(t_val))) %>%
    mutate(cell_type = "WBC") %>%
    mutate(task = tsk)
  
  all_correlations[[length(all_correlations)+1]] <- rbc_correlations
  all_correlations[[length(all_correlations)+1]] <- wbc_correlations
  
  rbc_coefficient_associations <- lapply(rbc_features,function(x) {
    mean_name <- paste(x,"mean",sep='_')
    fo <- as.formula(paste(x,paste(rvc,collapse="+"),sep = "~"))
    o <- as.data.frame(summary(glm(fo,data = sub_rvc))$coefficients)
    o$features <- x
    o$vct <- rownames(o)
    o$cell_type <- "RBC"
    rownames(o) <- NULL
    return(o)
  }) %>%
    do.call(what = rbind) %>%
    subset(vct != "(Intercept)") %>%
    mutate(task = tsk)
  wbc_coefficient_associations <- lapply(wbc_features,function(x) {
    mean_name <- paste(x,"mean",sep='_')
    fo <- as.formula(paste(x,paste(wvc,collapse="+"),sep = "~"))
    o <- as.data.frame(summary(glm(fo,data = sub_wvc))$coefficients)
    o$features <- x
    o$vct <- rownames(o)
    o$cell_type <- "WBC"
    rownames(o) <- NULL
    return(o) }) %>%
    do.call(what = rbind) %>% 
    subset(vct != "(Intercept)") %>%
    mutate(task = tsk)
  all_coefficient_associations[[
    length(all_coefficient_associations) + 1]] <- rbc_coefficient_associations
  all_coefficient_associations[[
    length(all_coefficient_associations) + 1]] <- wbc_coefficient_associations
}

do.call(rbind,all_correlations) %>%
  mutate(features = features_conversion[features]) %>% 
  mutate(features = factor(features,levels = rev(features_conversion))) %>% 
  mutate(task = decode_model_name(task)) %>% 
  mutate(sig.b = ifelse(p.val < (0.05/length(p.val)),"*",NA)) %>%
  mutate(vct_order = as.numeric(as.factor(cell_type))*1000+as.numeric(str_match(key,'[0-9]+'))) %>%
  ggplot(aes(x = reorder(key,vct_order),y = features,fill = R,label = sig.b)) + 
  geom_tile() + 
  geom_text(vjust = 0.75,colour = "white") +
  scale_fill_gradient2(low = "blue4",mid = "white",high = "red4",midpoint = 0) + 
  theme_pretty(base_size = 6) + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.key.width = unit(0.1,"cm"),axis.title = element_blank(),
        panel.border = element_rect(fill=NA)) + 
  facet_grid(task ~ .,scales = "free_y",space = "free_y")

do.call(rbind,all_coefficient_associations) %>%
  mutate(features = features_conversion[features]) %>% 
  mutate(features = factor(features,levels = rev(features_conversion))) %>% 
  mutate(task = decode_model_name(task)) %>% 
  mutate(p.val = `Pr(>|t|)`) %>% 
  mutate(sig.b = ifelse(p.adjust(p.val,"bonferroni") < 0.05,"*",NA)) %>%
  mutate(E = Estimate) %>% 
  mutate(vct = ifelse(grepl("mean",vct),"Mean",vct)) %>% 
  mutate(vct_order = as.numeric(as.factor(cell_type))*1000+as.numeric(str_match(vct,'[0-9]+'))) %>%
  ggplot(aes(x = reorder(vct,vct_order),y = features,fill = E,label = sig.b)) + 
  geom_tile() + 
  geom_text(vjust = 0.75,colour = "white",size = 4) +
  geom_text(vjust = 0.75,colour = "black",size = 2) +
  scale_fill_gradient2(low = "blue4",mid = "white",high = "red4",midpoint = 0,
                       name = "Coefficient") + 
  theme_pretty(base_size = 6) + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.key.width = unit(0.1,"cm"),axis.title = element_blank(),
        panel.border = element_rect(fill=NA)) + 
  facet_grid(task + cell_type ~ .,scales = "free_y",space = "free_y") + 
  ggsave(sprintf("figures/%s/mil-comori-vct-variance-association.pdf",output_str),
         width = 4,height = 5.5)

rbc_merged %>%
  subset(task_idx == "binary") %>% 
  mutate(coarse_class.x = ifelse(coarse_class.x == "Normal","Normal","Disease")) %>% 
  ggplot(aes(x = RBC1,y = RBC16,size = mass_displacement_green.variance,colour = coarse_class.x,
             fill = coarse_class.x)) +
  geom_smooth(method = "glm",formula = y ~ x,alpha = 0.2,size = 0.25) +
  geom_point(alpha = 0.2) +
  ggside::geom_xsidepoint(aes(size = mass_displacement_green.variance),
                          y = 0.5,alpha = 0.5) +
  ggside::geom_ysidepoint(aes(size = mass_displacement_green.variance),
                          x = 0.5,alpha = 0.5) +
  theme_pretty(base_size = 6) + 
  theme(legend.position = "bottom",legend.key.size = unit(0.1,"cm")) + 
  scale_size(name = NULL,range = c(0.2,5)) + 
  scale_colour_manual(values = fine_colours,breaks = c("Normal","Disease"),name = NULL) + 
  scale_fill_manual(values = fine_colours,breaks = c("Normal","Disease"),name = NULL) +
  xlab("RBC1 proportion") +
  ylab("RBC16 proportion") + 
  guides(colour = guide_legend(nrow = 2)) + 
  ggsave(sprintf("figures/%s/mil-comori-vct-variance-association-example-1.pdf",output_str),
         width = 2,height = 2)

wbc_merged %>%
  subset(task_idx == "binary") %>% 
  mutate(coarse_class.x = ifelse(coarse_class.x == "Normal","Normal","Disease")) %>% 
  ggplot(aes(x = WBC1,y = perimeter_separate_nuclear.variance,colour = coarse_class.x,
             fill = coarse_class.x)) +
  geom_smooth(method = "glm",formula = y ~ x,alpha = 0.2,size = 0.25) +
  geom_point(alpha = 0.5,size = 0.25) +
  theme_pretty(base_size = 6) + 
  theme(legend.position = "bottom",legend.key.size = unit(0.1,"cm")) + 
  scale_size(name = "Variance",range = c(1,7)) + 
  scale_colour_manual(values = fine_colours,breaks = c("Normal","Disease"),name = NULL) + 
  scale_fill_manual(values = fine_colours,breaks = c("Normal","Disease"),name = NULL) + 
  xlab("WBC1 proportion") + 
  ylab("var(WBC nuclear perimeter)") + 
  ggsave(sprintf("figures/%s/mil-comori-vct-variance-association-example-2.pdf",output_str),
         width = 1.7,height = 1.7)

wbc_merged %>%
  subset(task_idx == "disease_binary" & coarse_class.x != "Normal") %>% 
  subset(slide_id != "XVI_6") %>%
  ggplot(aes(x = WBC3,y = convexity.variance,colour = coarse_class.x,
             fill = coarse_class.x)) +
  geom_smooth(method = "glm",formula = y ~ x,alpha = 0.2,size = 0.25) +
  geom_point(alpha = 0.5,size = 0.25) +
  theme_pretty(base_size = 6) + 
  theme(legend.position = "bottom",legend.key.size = unit(0.1,"cm")) + 
  scale_size(name = "Variance",range = c(1,7),trans = 'log10') + 
  scale_colour_manual(values = fine_colours,
                      breaks = c("MDS","Anaemia"),name = NULL) + 
  scale_fill_manual(values = fine_colours,
                    breaks = c("MDS","Anaemia"),name = NULL) + 
  xlab("WBC3") + 
  ylab("var(WBC convexity)") + 
  ggsave(sprintf("figures/%s/mil-comori-vct-variance-association-example-3.pdf",output_str),
         width = 1.7,height = 1.7)

# feature distributions vs. vct -------------------------------------------

wbc_consensus <- all_consensus_mo_df %>%
  subset(cell_type == "WBC" & grepl("\\.bc",model_name)) %>% 
  select(virtual_cell_type,is_consensual,task_idx)
wbc_consensus$task <- str_split(rownames(wbc_consensus),pattern = " ") %>%
  sapply(function(x) str_match(x[2],'[a-z_]+'))
rbc_consensus <- all_consensus_mo_df %>%
  subset(cell_type == "RBC" & grepl("\\.bc",model_name)) %>% 
  select(virtual_cell_type,is_consensual,task_idx)
rbc_consensus$task <- str_split(rownames(rbc_consensus),pattern = " ") %>%
  sapply(function(x) str_match(x[2],'[a-z_]+'))

features_rbc <- features_all[
  unlist(read_csv("data_output/rbc_feature_subset",col_names = F))]
features_wbc <- c(features_all,features_nuclear)[
  unlist(read_csv("data_output/wbc_feature_subset",col_names = F))]

rbc_many_cells <- list.files("datasets/many-cells",pattern = "^rbc-cv_subset\\.multi.*bc",
                             full.names = T) %>%
  lapply(read_csv,col_names = c(
    "model_name","slide_id","cell_type","virtual_cell_type",features_rbc)) %>%
  do.call(what = rbind) %>%
  merge(all_conditions,by = "slide_id") %>%
  mutate(fine_class = ifelse((fine_class != "SF3B1-mutant") & (coarse_class == "MDS"),
                             "SF3B1-wildtype",as.character(fine_class))) %>%
  mutate(virtual_cell_type = virtual_cell_type + 1) %>%
  mutate(binary_class = ifelse(coarse_class == "Normal","Normal","Disease"))
wbc_many_cells <- list.files("datasets/many-cells",pattern = "^wbc-cv_subset\\.multi.*bc",
                             full.names = T) %>%
  lapply(read_csv,col_names = c(
    "model_name","slide_id","cell_type","virtual_cell_type",features_wbc)) %>%
  do.call(what = rbind) %>%
  merge(all_conditions,by = "slide_id") %>%
  mutate(fine_class = ifelse((fine_class != "SF3B1-mutant") & (coarse_class == "MDS"),
                            "SF3B1-wildtype",as.character(fine_class))) %>%
  mutate(virtual_cell_type = virtual_cell_type + 1) %>%
  mutate(binary_class = ifelse(coarse_class == "Normal","Normal","Disease"))

rbc_vct_centers <- rbc_many_cells %>%
  mutate(M = length(virtual_cell_type)) %>% 
  select(-cell_type) %>%
  group_by(virtual_cell_type,fine_class,coarse_class) %>%
  mutate(N = length(virtual_cell_type)) %>% 
  gather(key = "feature",value = "value",-slide_id,-model_name,
         -virtual_cell_type,-fine_class,-coarse_class,-binary_class,-N,-M) %>% 
  group_by(virtual_cell_type,feature) %>%
  mutate(center = mean(value),
         center_median = median(value)) %>% 
  ungroup %>%
  select(virtual_cell_type,feature,center,center_median,fine_class,coarse_class,binary_class,N,M) %>%
  distinct %>% 
  merge(rbc_consensus,by = "virtual_cell_type") %>%
  select(-task,-task_idx) %>%
  distinct

wbc_vct_centers <- wbc_many_cells %>%
  mutate(M = length(virtual_cell_type)) %>% 
  select(-cell_type) %>%
  group_by(virtual_cell_type,fine_class,coarse_class) %>%
  mutate(N = length(virtual_cell_type)) %>% 
  gather(key = "feature",value = "value",-slide_id,-model_name,
         -virtual_cell_type,-fine_class,-coarse_class,-binary_class,-N,-M) %>% 
  group_by(virtual_cell_type,feature) %>%
  mutate(center = median(value)) %>% 
  ungroup %>%
  select(virtual_cell_type,feature,center,fine_class,coarse_class,binary_class,N,M) %>%
  distinct %>% 
  merge(wbc_consensus,by = "virtual_cell_type") %>%
  select(-task,-task_idx) %>%
  distinct

# WBC binary classification plots - nuclear perimeter

center_data <- subset(wbc_vct_centers,feature == "perimeter_separate_nuclear") %>%
  group_by(binary_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(binary_class) %>%
  mutate(M = sum(N)) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"wbc")) %>%
  na.omit() %>% 
  group_by(virtual_cell_type) %>%
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 40
wbc_many_cells %>% 
  ggplot() +
  geom_density(aes(x = perimeter_separate_nuclear,colour = binary_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = binary_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lavender") + 
  xlab("WBC nuclear perimeter") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("Normal","Disease")) +
  scale_colour_manual(values = fine_colours,breaks = c("Normal","Disease"),guide=F) + 
  coord_cartesian(xlim = c(NA,450)) + 
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-wbc-nuclear-perimeter.pdf",output_str),
         width = 3.5,height = 2)

# RBC disease classification plots - cdf std

center_data <- subset(rbc_vct_centers,feature == "cdf_std") %>%
  subset(coarse_class %in% c("MDS","Anaemia")) %>% 
  group_by(coarse_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(coarse_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"rbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 1
rbc_many_cells %>% 
  subset(coarse_class %in% c("MDS","Anaemia")) %>% 
  ggplot() +
  geom_density(aes(x = cdf_std,colour = coarse_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = coarse_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lightpink") + 
  xlab("RBC std(CDF)") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("MDS","Anaemia")) +
  scale_colour_manual(values = fine_colours,breaks = c("MDS","Anaemia"),guide = F) +
  scale_alpha(guide = F) +
  #coord_cartesian(xlim = c(10,NA)) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-rbc-std-cdf.pdf",output_str),
         width = 3.5,height = 2)

# WBC disease classification plots - convexity

center_data <- subset(wbc_vct_centers,feature == "convexity") %>%
  subset(coarse_class %in% c("MDS","Anaemia")) %>% 
  group_by(coarse_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(coarse_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"wbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 0.003
wbc_many_cells %>% 
  subset(coarse_class %in% c("MDS","Anaemia")) %>% 
  ggplot() +
  geom_density(aes(x = convexity,colour = coarse_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = coarse_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lavender") + 
  xlab("WBC convexity") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("MDS","Anaemia")) +
  scale_colour_manual(values = fine_colours,breaks = c("MDS","Anaemia"),guide = F) + 
  scale_alpha(guide = F) +
  coord_cartesian(xlim = c(0.84,0.97)) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-wbc-convexity.pdf",output_str),
         width = 3.5,height = 2)

# RBC SF3B1 plots - mean cdf

center_data <- subset(rbc_vct_centers,feature == "cdf_mean") %>%
  subset(coarse_class %in% c("MDS")) %>% 
  group_by(fine_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(fine_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"rbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 80
rbc_many_cells %>% 
  subset(coarse_class %in% c("MDS")) %>% 
  ggplot() +
  geom_density(aes(x = cdf_mean,colour = fine_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = fine_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lightpink") + 
  xlab("RBC mean(CDF)") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype")) +
  scale_colour_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype"),guide=F) +
  scale_alpha(guide = F) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-rbc-cdf-mean.pdf",output_str),
         width = 3.5,height = 2)

# RBC SF3B1 plots - area

center_data <- subset(rbc_vct_centers,feature == "area") %>%
  subset(coarse_class %in% c("MDS")) %>% 
  group_by(fine_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(fine_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"rbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 120
rbc_many_cells %>% 
  subset(coarse_class %in% c("MDS")) %>% 
  ggplot() +
  geom_density(aes(x = area,colour = fine_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = fine_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lightpink") + 
  xlab("RBC area") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype")) +
  scale_colour_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype"),guide=F) +
  scale_alpha(guide = F) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-rbc-area.pdf",output_str),
         width = 3.5,height = 2)

# WBC SF3B1 plots - nuclear convexity

center_data <- subset(wbc_vct_centers,feature == "convexity_separate_nuclear") %>%
  subset(coarse_class %in% c("MDS")) %>% 
  group_by(fine_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(fine_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"wbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 0.035
wbc_many_cells %>% 
  subset(coarse_class %in% c("MDS")) %>% 
  ggplot() +
  geom_density(aes(x = convexity_separate_nuclear,colour = fine_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = fine_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lavender") + 
  xlab("WBC nuclear convexity") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype")) +
  scale_colour_manual(values = fine_colours,breaks = c("SF3B1-mutant","SF3B1-wildtype"),guide=F) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-wbc-nuclear-convexity.pdf",output_str),
         width = 3.5,height = 2)

# WBC anaemia plots - solidity

center_data <- subset(wbc_vct_centers,feature == "solidity") %>%
  subset(coarse_class %in% c("Anaemia")) %>% 
  group_by(fine_class,virtual_cell_type,center,M) %>%
  summarise(N = sum(N)) %>%
  group_by(fine_class) %>%
  mutate(M = sum(N)) %>%
  group_by(virtual_cell_type) %>%
  mutate(virtual_cell_type = vct_conversion(virtual_cell_type,"wbc")) %>%
  na.omit() %>% 
  mutate(plot_label = as.character(ifelse(N/M == max(N/M),virtual_cell_type,NA)))

d <- 0.005
wbc_many_cells %>% 
  subset(coarse_class %in% c("Anaemia")) %>% 
  ggplot() +
  geom_density(aes(x = solidity,colour = fine_class)) + 
  geom_segment(data = center_data,
               aes(x = center,y = 0,xend = center,yend = N/M/d),
               linetype = 3,size = 0.25) +
  geom_point(data = center_data,
             aes(x = center,y = N/M/d,colour = fine_class),
             stat = "identity",width = 10,
             position=position_dodge(width = 0),size = 1.5) + 
  geom_label(data = center_data,colour = "black",size = 2.1,alpha = 0.7,
             label.r = unit(0.1,"cm"),label.padding = unit(0.05,"cm"),
             aes(x = center,y = N/M/d,label = plot_label),vjust = -0.5,
             fill = "lightpink") + 
  xlab("WBC solidity") + 
  ylab("Density (line)") + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .*d,name = "VCT proportion (dots)"),
                     expand = c(0,0,0.01,0)) + 
  theme_pretty(base_size = 6) + 
  theme(legend.key.size = unit(0.2,"cm"),legend.title = element_blank(),
        legend.position = "bottom",legend.box.spacing = unit(0.05,"cm")) + 
  scale_fill_manual(values = c(`Iron deficiency` = "palevioletred1",Megaloblastic = "red4")) +
  scale_colour_manual(values = c(`Iron deficiency` = "palevioletred1",Megaloblastic = "red4"),
                      guide = F) + 
  scale_alpha(guide = F) +
  coord_cartesian(xlim = c(NA,1.3)) + 
  scale_x_continuous(breaks = c(1,1.1,1.2,1.3,1.4,1.5)) +
  ggsave(sprintf("figures/%s/mil-comori-density-vs-vct-wbc-solidty.pdf",output_str),
         width = 3.5,height = 2)

# feature wise variance values

remove_outliers <- function(x,m=1.5) {
  iqr <- IQR(x)
  iqr_range <- quantile(x,c(0.25,0.75))
  x <- (x > (iqr_range[1] - iqr*m)) & (x < (iqr_range[2] + iqr*m)) 
  return(x)
}

wbc_many_cells %>%
  filter(!(slide_id %in% poor_quality_slides)) %>% 
  ungroup %>% 
  group_by(binary_class) %>%
  filter(remove_outliers(perimeter_separate_nuclear,5)) %>%
  summarise(V = var(perimeter_separate_nuclear))

rbc_many_cells %>%
  filter(!(slide_id %in% poor_quality_slides)) %>% 
  ungroup %>% 
  group_by(coarse_class) %>%
  filter(remove_outliers(cdf_std,5)) %>%
  summarise(V = var(cdf_std))

wbc_many_cells %>%
  filter(!(slide_id %in% poor_quality_slides)) %>% 
  subset(coarse_class %in% c("MDS")) %>%
  ungroup %>% 
  group_by(fine_class) %>%
  filter(remove_outliers(convexity_separate_nuclear,5)) %>%
  summarise(V = var(convexity_separate_nuclear))

rbc_many_cells %>%
  filter(!(slide_id %in% poor_quality_slides)) %>% 
  subset(coarse_class %in% c("MDS")) %>%
  ungroup %>% 
  group_by(fine_class) %>%
  filter(remove_outliers(area,5)) %>%
  summarise(V = mean(area))
