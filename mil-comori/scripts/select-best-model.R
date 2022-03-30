# setup -------------------------------------------------------------------

library(tidyverse)
library(caret)
library(pROC)

get.coords.for.ggplot <- function(roc) {
  # from pROC source code
  df <- coords(roc, "all", transpose = FALSE)
  return(df[rev(seq(nrow(df))),])
}

select <- dplyr::select

set.seed(42)

read_prediction <- function(file_path) {
  read_csv(file_path,
           col_names = c("slide_id","labels","p1","p2")) %>%
    mutate(task = gsub(".csv","",gsub("labels/","",labels))) %>%
    return
}

task_conversion <- c(
  cv_binary = "Disease detection",
  cv_disease_binary = "Disease classification",
  cv_mds_binary = "SF3B1mut detection",
  cv_anemia_binary = "Anaemia classification",
  cv_multi_objective = "cv_multi_objective")
task_conversion_reverse <- names(task_conversion)
names(task_conversion_reverse) <- task_conversion

task_conversion_mo <- c("cv_anemia_binary","cv_binary","cv_mds_binary","cv_disease_binary")

args <- c("metrics/cv_subset_mll_adden_1_metrics_df.csv", "metrics/cv_subset_mll_adden_1_probs.csv",
          "best_models/best_models_subset","best_models/best_models_mo_subset",
          "cv_subset_mll_adden_1")
args <- commandArgs(trailingOnly=TRUE)
metrics_file <- args[1]
probs_file <- args[2]
out_binary <- args[3]
out_mo <- args[4]
cv_string <- args[5]

all_metrics <- read_csv(
    metrics_file,
    col_names = c("set","fold","metric","value","nvc","nc",
                  "task_number","task","dataset","feature_set")) 

all_metrics <- all_metrics %>%
  mutate(task_original = task) %>%
  mutate(multi_objective = ifelse(
    task == "cv_multi_objective","Multiple objective","Single objective")) %>% 
  mutate(multi_objective = factor(
    multi_objective,levels = c("Single objective","Multiple objective"))) %>%
  mutate(task = ifelse(task == "cv_multi_objective",task_conversion_mo[task_number+1],task)) %>%
  mutate(task = factor(task_conversion[task],levels = rev(task_conversion))) %>%
  mutate(dataset_pretty = c(cells = "Morphology",cells_bc = "Morphology + B.C.")[dataset])

probs_class <- read_csv(
    probs_file,
    col_names = c("X","fold","task_idx","prob","class",
                  "class_true","task","nvc","data_type","feature_set"))

probs_class_idx <- probs_class %>%
  select(task_idx,task,data_type,nvc,feature_set) %>%
  distinct

all_roc <- list()

for (i in 1:nrow(probs_class_idx)) {
  l_s <- length(all_roc)+1
  tmp <- probs_class_idx[i,] %>%
    unlist
  sub_probs_class <- probs_class %>%
    subset(task_idx==tmp[1] & task==tmp[2] & data_type==tmp[3] & nvc==tmp[4] & feature_set==tmp[5])
  full_roc <- roc(sub_probs_class$class_true,sub_probs_class$prob)
  roc_coords <- get.coords.for.ggplot(full_roc) 
  
  if (ncol(roc_coords) > 4) {
    roc_coords <- t(roc_coords)
  }
  roc_coords <- roc_coords %>%
    as.tibble %>%
    mutate(data_type = tmp[3])
  all_roc[[l_s]] <- roc_coords
  all_roc[[l_s]]$task_idx <- tmp[1]
  all_roc[[l_s]]$task <- tmp[2]
  all_roc[[l_s]]$nvc <- tmp[4]
  all_roc[[l_s]]$feature_set <- tmp[5]
  all_roc[[l_s]]$auc_value <- full_roc$auc
}

all_roc_df <- do.call(rbind,all_roc) %>%
  subset(is.finite(threshold)) %>%
  arrange(-threshold) %>%
  mutate(mo = task == "cv_multi_objective") %>% 
  mutate(task_idx = as.numeric(task_idx)) %>% 
  mutate(task = ifelse(mo,task_conversion_mo[task_idx+1],task)) %>%
  mutate(task = gsub("cv_","",task)) %>%
  mutate(task = factor(task,
                       levels = c("binary","disease_binary",
                                  "mds_binary","anemia_binary"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  mutate(data_type = factor(data_type,levels = c("bc","cells","cells_bc"),
                            labels = c("B.C.","Morphology","Morphology + B.C."))) %>%
  mutate(nvc = as.numeric(nvc))

# assess which hyperparameters are best -----------------------------------

best_models_so <- all_roc_df %>%
  select(nvc,task,multi_objective = mo,value = auc_value,dataset = data_type) %>%
  subset(multi_objective == F) %>%
  group_by(nvc,task,multi_objective,dataset) %>%
  summarise(value = mean(value)) %>%
  mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
  group_by(task,multi_objective,dataset) %>% 
  filter(value == max(value)) %>%
  filter(nvc == min(as.numeric(as.character(nvc)))) %>% 
  as_tibble() %>% 
  arrange(dataset,task)

best_models_mo <- all_roc_df %>%
  select(nvc,task,multi_objective = mo,value = auc_value,dataset = data_type) %>%
  subset(multi_objective == T) %>%
  group_by(nvc,task,multi_objective,dataset) %>%
  summarise(value = mean(value)) %>%
  mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
  group_by(nvc,multi_objective,dataset) %>% 
  mutate(best_average = mean(value)) %>%
  group_by(dataset) %>%
  filter(best_average == max(best_average)) %>%
  select(-best_average) %>%
  filter(nvc == min(as.numeric(as.character(nvc)))) %>% 
  as_tibble()

best_models_roc_curves <- list()
best_models_so_list <- list()
best_models_mo_list <- list()

for (i in 1:nrow(best_models_so)) {
  tmp <- best_models_so[i,]
  best_models_so_list[[i]] <- all_metrics %>%
    subset(set == "TEST" & metric == "AUC_0" &
             dataset_pretty == tmp$dataset & 
             multi_objective == "Single objective" &
             nvc == tmp$nvc & 
             as.character(task) == as.character(tmp$task)) %>%
    filter(fold == fold[which.max(value)[1]])
  best_models_roc_curves[[i]] <- all_roc_df %>%
    subset(data_type == tmp$dataset & 
             mo == F &
             nvc == tmp$nvc & 
             as.character(task) == as.character(tmp$task))
}

for (i in 1:nrow(distinct(select(best_models_mo,nvc,dataset)))) {
  tmp <- distinct(select(best_models_mo,nvc,dataset))[i,]
  best_models_mo_list[[i]] <- all_metrics %>%
    subset(set == "TEST" & metric == "AUC_0" &
             dataset_pretty == tmp$dataset & 
             multi_objective == "Multiple objective" &
             nvc == tmp$nvc) %>%
    group_by(fold) %>%
    mutate(M = mean(value)) %>%
    ungroup %>% 
    filter(fold == fold[which.max(M)[1]]) %>%
    select(fold,nvc,dataset)
  best_models_roc_curves[[length(best_models_roc_curves)+1]] <- all_roc_df %>%
    subset(data_type == tmp$dataset & 
             mo == T &
             nvc == tmp$nvc)
}

print(do.call(what=rbind,best_models_so_list))

do.call(what = rbind,best_models_so_list) %>%
  arrange() %>% 
  mutate(to = gsub("cv_","",task_original)) %>%
  transmute(a = sprintf("models/%s.%s.%s%s",
                        cv_string,to,nvc,
                        ifelse(dataset == "cells","",".bc")),
            b = fold,c = sprintf("labels/%s.csv",
                                 to)) %>%
  write.table(file = out_binary,
              row.names = F,col.names = F,quote = F,sep = ',')

do.call(what = rbind,best_models_mo_list) %>%
  distinct %>%
  transmute(a = sprintf("models/%s.%s.%s%s",
                        cv_string,"multi_objective",nvc,
                        ifelse(dataset == "cells","",".bc")),
            b = fold,
            c = sprintf(
              "labels/binary.csv",NA)) %>%
  write.table(file = out_mo,
              row.names = F,col.names = F,quote = F,sep = ',')
