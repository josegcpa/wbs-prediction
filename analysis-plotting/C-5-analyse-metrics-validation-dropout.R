# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(caret)
library(pROC)
library(ggrepel)

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
           coarse_class = factor(class_conversion[coarse_class],class_conversion))) %>%
  mutate(fine_class = ifelse(fine_class == "SRSF2-mutant","Non-SF3B1-mutant",as.character(fine_class)))

label_conversion <- list(
  binary = c(
    `Normal` = 0,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 1,
    `Iron deficiency` = 1,`Megaloblastic` = 1),
  disease_binary = c(
    `Normal` = NA,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 1,
    `Iron deficiency` = 0,`Megaloblastic` = 0),
  mds_binary = c(
    `Normal` = NA,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 0,
    `Iron deficiency` = NA,`Megaloblastic` = NA),
  anemia_binary = c(
    `Normal` = NA,
    `SF3B1-mutant` = NA,`Non-SF3B1-mutant` = NA,
    `Iron deficiency` = 0,`Megaloblastic` = 1)
)

all_metrics <- rbind(
  read_csv(
    "../mile-vice/metrics/cv_subset_dropout_metrics_df.csv",
    col_names = c("set","fold","metric","value","nvc","nc",
                  "task_number","task","dataset","feature_set"))) %>%
  mutate(task_original = task) %>%
  subset(task != "cv_multi_class") %>% 
  mutate(multi_objective = ifelse(
    task == "cv_multi_objective","Multiple objective","Single objective")) %>% 
  mutate(multi_objective = factor(
    multi_objective,levels = c("Single objective","Multiple objective"))) %>%
  mutate(task = ifelse(task == "cv_multi_objective",task_conversion_mo[task_number+1],task)) %>%
  mutate(task = factor(task_conversion[task],levels = rev(task_conversion))) %>%
  mutate(dataset_pretty = c(cells = "Morphology",cells_bc = "Morphology + B.C.")[dataset])

probs_class <- rbind(
  read_csv(
    "../mile-vice/metrics/cv_subset_dropout_probs.csv",
    col_names = c("X","fold","task_idx","prob","class","class_true","task","nvc","data_type","feature_set"))) %>%
  subset(feature_set == "subset_features")

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
  roc_coords <- get.coords.for.ggplot(full_roc) %>%
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

decode_model_name <- function(model_names) {
  ifelse(
    grepl("anemia_binary",model_names),"Anaemia classification",
    ifelse(grepl("mds_binary",model_names),"SF3B1mut detection",
           ifelse(grepl("disease_binary",model_names),
                  "Disease classification",
                  "Disease detection"))) %>%
    return
}

model_levels <- c("Disease detection","Disease classification","SF3B1mut detection","Anaemia classification")

metrics_cv <- best_models_so
all_auc_validation <- list()
all_auc_validation_mo <- list()
for (file_path in list.files("../mile-vice/ev-scores",pattern = "cv_subset_dropout\\.",full.names = T)) {
  if (!grepl("multi_obj",file_path)) {
    data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
    tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
    N <- tmp$value[tmp$metric == "N_0"]
    tmp <- tmp %>%
      subset(metric == "AUC_0")
    all_auc_validation[[file_path]] <- data.frame(
      task = decode_model_name(tmp$model_id),dataset = data_type,value = tmp$value,N=N)
  } else {
    data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
    tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
    N <- tmp$value[tmp$metric == "N_0"]
    tmp <- tmp %>%
      subset(metric == "AUC_0")
    all_auc_validation_mo[[file_path]] <- data.frame(
      task = decode_model_name(file_path),dataset = data_type,value = tmp$value,N=N)
  }
}

all_auc_validation_df <- mutate(do.call(rbind,all_auc_validation),set = "External validation") %>%
  mutate(dataset = factor(dataset,levels = c("Morphology","Morphology + B.C.")))
all_auc_validation_mo_df <- mutate(do.call(rbind,all_auc_validation_mo),set = "External validation") %>%
  mutate(dataset = factor(dataset,levels = c("Morphology","Morphology + B.C.")))

metrics_cv %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_df,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = ifelse(value + 1/sqrt(N) > 1,1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.2,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
  scale_color_discrete(guide = F) + 
  ggsave(filename = "figures/auc-bars-w-validation-mile-vice-dropout.pdf",height=1.7,width=3)

best_models_mo %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = ifelse(value + 1/sqrt(N) > 1,1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.2,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
  scale_color_discrete(guide = F) + 
  ggsave(filename = "figures/auc-bars-w-validation-mile-vice-dropout-mo.pdf",height=1.7,width=3)

metrics_cv <- best_models_so
all_auc_validation <- list()
all_auc_validation_mo <- list()
for (file_path in list.files("../mile-vice/ev-scores-consensus/",pattern = "cv_subset_dropout\\.",full.names = T)) {
  if (!grepl("multi_obj",file_path)) {
    data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
    tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
    N <- tmp$value[tmp$metric == "N_0"]
    tmp <- tmp %>%
      subset(metric == "AUC_0")
    all_auc_validation[[file_path]] <- data.frame(
      task = decode_model_name(tmp$model_id),dataset = data_type,value = tmp$value,N=N)
  } else {
    data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
    tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
    N <- tmp$value[tmp$metric == "N_0"]
    tmp <- tmp %>%
      subset(metric == "AUC_0")
    all_auc_validation_mo[[file_path]] <- data.frame(
      task = decode_model_name(file_path),dataset = data_type,value = tmp$value,N=N)
  }
}

all_auc_validation_df <- mutate(do.call(rbind,all_auc_validation),set = "External validation") %>%
  mutate(dataset = factor(dataset,levels = c("Morphology","Morphology + B.C."))) %>%
  mutate(N = as.numeric(as.character(N)),value = as.numeric(as.character(value)))
all_auc_validation_mo_df <- mutate(do.call(rbind,all_auc_validation_mo),set = "External validation") %>%
  mutate(dataset = factor(dataset,levels = c("Morphology","Morphology + B.C."))) %>%
  mutate(N = as.numeric(as.character(N)),value = as.numeric(as.character(value)))

metrics_cv %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_df,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = ifelse(value + 1/sqrt(N) > 1,1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.2,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
  scale_color_discrete(guide = F) + 
  ggsave(filename = "figures/auc-bars-w-validation-mile-vice-dropout-consensus.pdf",height=1.7,width=3)

best_models_mo %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = ifelse(value + 1/sqrt(N) > 1,1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_mo_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.2,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
  scale_color_discrete(guide = F) + 
  ggsave(filename = "figures/auc-bars-w-validation-mile-vice-dropout-mo-consensus.pdf",height=1.7,width=3)
