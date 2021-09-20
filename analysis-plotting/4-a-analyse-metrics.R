# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(caret)
library(pROC)

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

all_metrics <- read_csv(
  "../mile-vice/cv_metrics.csv",
  col_names = c("set","fold","metric","value","nvc","nc","task_number","task","dataset","feature_set")) %>%
  mutate(task_original = task) %>%
  mutate(multi_objective = ifelse(
    task == "cv_multi_objective","Multiple objective","Single objective")) %>% 
  mutate(multi_objective = factor(
    multi_objective,levels = c("Single objective","Multiple objective"))) %>%
  mutate(task = ifelse(task == "cv_multi_objective",task_conversion_mo[task_number+1],task)) %>%
  mutate(task = factor(task_conversion[task],levels = rev(task_conversion))) %>%
  mutate(dataset_pretty = c(cells = "Morphology",cells_bc = "Morphology + B.C.")[dataset]) %>%
  subset(feature_set == "subset_features")

probs_class <- read_csv(
  "../mile-vice/probs_class_train.csv",
  col_names = c("X","fold","task_idx","prob","class","class_true","task","nvc","data_type","feature_set")) %>%
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

# assess which hyperparameters are best -----------------------------------

all_roc_df %>%
  subset(data_type == "Morphology") %>%
  select(nvc,task,multi_objective = mo,value = auc_value,data_type) %>%
  mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>%
  ggplot(aes(x = task,y = value,
             fill = nvc,colour = multi_objective,
         group = reorder(paste(nvc,multi_objective),
                         as.numeric(nvc)+ 1000 * as.numeric(multi_objective)))) + 
  geom_bar(position = position_dodge(width = 0.95),
           stat = "identity",size = 0.25) +
  theme_pretty(base_size = 6) +
  coord_flip(ylim = c(0.85,1)) +
  scale_fill_manual(values = colorRampPalette(colors = c("gold","red4"))(length(unique(all_metrics$nvc))),
                    name = NULL) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2,"cm"),
        axis.title.y = element_blank()) + 
  ylab("AUC") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_colour_manual(values = c(NA,"black"),guide = F) +
  ggsave("figures/mile-vice-cv-performance.pdf",height = 1.8,width = 3)

all_roc_df %>%
  select(nvc,task,multi_objective = mo,value = auc_value,dataset = data_type) %>%
  mutate(multi_objective = factor(
    ifelse(multi_objective,"Multiple objective","Single objective"),
    levels = c("Single objective","Multiple objective"))) %>%
  mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
  group_by(task,multi_objective,dataset) %>% 
  mutate(MAX = ifelse(value == max(value),"X",NA)) %>% 
  ggplot(aes(x = task,y = nvc,fill = value)) + 
  geom_tile() +
  geom_tile(aes(colour = MAX),fill = NA,size = 0.5) +
  geom_text(aes(label = sprintf("%.3f",value)),colour = "black",size = 2) +
  theme_pretty(base_size = 6) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        axis.title.y = element_blank()) + 
  facet_grid(dataset ~ multi_objective,scales = "free_x",space = "free_x") +
  ylab("Number of virtual cells") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  scale_colour_manual(values = c("black",NA),guide = F) +
  scale_fill_gradient(low = "lightcyan1",high = "purple",name = "AUC") +
  ggsave("figures/mile-vice-cv-performance-tile.pdf",height = 1.7,width = 3.5)

all_roc_df %>%
  select(nvc,task,multi_objective = mo,value = auc_value,dataset = data_type) %>%
  mutate(multi_objective = factor(
    ifelse(multi_objective,"Multiple objective","Single objective"),
    levels = c("Single objective","Multiple objective"))) %>%
  mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
  group_by(task,multi_objective,dataset) %>% 
  mutate(MAX = ifelse(value == max(value),"X",NA)) %>%
  ggplot(aes(x = nvc,y = value,
             colour = task,group = paste(task,multi_objective),
             shape = multi_objective,linetype = multi_objective)) + 
  geom_point() +
  geom_line(alpha = 0.8) +
  theme_pretty(base_size = 6) +
  scale_colour_brewer(type = "qual",palette = 2,name = NULL) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0,"cm")) + 
  guides(shape = guide_legend(nrow = 4),
         colour = guide_legend(nrow = 4)) + 
  scale_shape_manual(name = NULL,values = c(16,4)) + 
  scale_linetype(name = NULL) +
  ylab("Median AUC") + 
  xlab("Number of virtual cells") +
  facet_wrap(~ dataset) +
  ggsave("figures/mile-vice-cv-performance-lines.pdf",height = 1.8,width = 2.5)

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
    subset(multi_objective == "Multiple objective") %>%
    subset(set == "TEST" & metric == "AUC_0" &
             dataset_pretty == tmp$dataset & 
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

do.call(what = rbind,best_models_so_list) %>%
  arrange() %>% 
  mutate(to = gsub("cv_","",task_original)) %>%
  transmute(a = sprintf("models/cv_subset.%s.%s%s",
                        to,nvc,
                        ifelse(dataset == "cells","",".bc")),
            b = fold,c = sprintf("labels/%s.csv",
                                 to)) %>%
  write.table(file = "../mile-vice/best_models",
              row.names = F,col.names = F,quote = F,sep = ',')

do.call(what = rbind,best_models_mo_list) %>%
  distinct %>%
  transmute(a = sprintf("models/cv_subset.%s.%s%s",
                        "multi_objective",nvc,
                        ifelse(dataset == "cells","",".bc")),
            b = fold,
            c = sprintf(
              "labels/binary.csv",NA)) %>%
  write.table(file = "../mile-vice/best_models_mo",
              row.names = F,col.names = F,quote = F,sep = ',')

# auroc analysis ----------------------------------------------------------

best_models_roc_curves_df <- do.call(rbind,best_models_roc_curves)

best_models_roc_curves_df %>%
  arrange(sensitivity,-specificity) %>% 
  subset(mo == F) %>%
  ggplot(aes(x = 1-specificity,y = sensitivity,colour = data_type)) + 
  geom_abline(slope = 1,intercept = 0,size = 0.25,linetype = 2,alpha = 0.5) +
  geom_line(size = 0.5,alpha = 0.7) +
  facet_wrap(~ task) +
  theme_pretty(base_size = 6) + 
  scale_colour_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm"),
        panel.spacing = unit(0.8,"lines")) +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) + 
  ggsave("figures/mile-vice-roc-curve.pdf",height = 2.3,width = 2.3)

best_models_roc_curves_df %>%
  arrange(sensitivity,-specificity) %>% 
  subset(mo == T) %>%
  ggplot(aes(x = 1-specificity,y = sensitivity,colour = data_type)) + 
  geom_abline(slope = 1,intercept = 0,size = 0.25,linetype = 2,alpha = 0.5) +
  geom_line(size = 0.5,alpha = 0.7) +
  facet_wrap(~ task) +
  theme_pretty(base_size = 6) + 
  scale_colour_manual(values = c("red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm"),
        panel.spacing = unit(0.8,"lines")) +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) + 
  ggsave("figures/mile-vice-roc-curve-mo.pdf",height = 2.3,width = 2.3)

# compare with glmnet -----------------------------------------------------

glmnet_scores <- read.csv("data_output/glmnet-auroc.csv") %>% 
  select(dataset = data_type,task,glmnet_auc = auc_value) %>%
  distinct %>%
  merge(rbind(best_models_so,best_models_mo),
        by = c("task","dataset")) %>%
  gather(key = "key",value = "value",value,glmnet_auc) %>%
  mutate(key = ifelse(key == "glmnet_auc","glmnet","MILe-ViCe")) %>%
  mutate(mo = multi_objective) %>%
  subset(!(mo == T & key == "glmnet")) %>%
  mutate(task = factor(
    task,
    levels = c("Disease detection","Disease classification",
               "SF3B1mut detection","Anaemia classification") %>% 
      rev))

ggplot(glmnet_scores,aes(x = task,y = value,fill = key)) + 
  geom_bar(stat = "identity",
           aes(group = reorder(paste(key,mo),mo*1000 + as.numeric(key))),
           position = position_dodge(width = 0.95)) + 
  geom_bar(stat = "identity",
           aes(group = reorder(paste(key,mo),mo*1000 + as.numeric(key)),
               colour = mo),
           fill = NA,
           position = position_dodge(width = 0.95)) +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%",x*100),expand = c(0.05,0,0,0)) + 
  ylab("AUC") +
  xlab("") +
  coord_flip(ylim = c(0.7,1)) +
  theme_pretty(base_size = 6) + 
  facet_wrap(~ dataset) +
  scale_colour_manual(guide = F,values = c(NA,"black")) +
  scale_fill_manual(name = NULL,values = c("lightblue","goldenrod")) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2,"cm"),
        panel.spacing = unit(1.5,"lines")) +
  ggsave("figures/mile-vice-vs-glmnet.pdf",height = 1.5,width = 4)
  
spread(glmnet_scores,key = "key",value = "value") %>%
  group_by(dataset,task) %>%
  mutate(glmnet = ifelse(is.na(glmnet),glmnet[!is.na(glmnet)][1],glmnet)) %>%
  ggplot(aes(x = glmnet,y = `MILe-ViCe`,colour = dataset,shape = mo)) + 
  geom_abline(slope = 1,linetype = 3,alpha = 0.5,size = 0.5) +
  geom_point(size = 1) + 
  scale_x_continuous(labels = function(x) sprintf("%.1f%%",x*100)) + 
  scale_y_continuous(labels = function(x) sprintf("%.1f%%",x*100)) + 
  ylab("MILe-ViCe AUC") +
  xlab("glmnet AUC") +
  coord_cartesian(xlim = c(0.8,1),ylim = c(0.8,1)) +
  theme_pretty(base_size = 6) + 
  scale_colour_manual(values = c("red4","orange"),name = NULL) + 
  scale_shape_manual(values = c(3,16),
                     labels = c("Single objective","Multiple objective"),name = NULL) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2,"cm"),
        panel.spacing = unit(1.5,"lines")) +
  guides(colour = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  ggsave("figures/mile-vice-vs-glmnet-scatter.pdf",height = 2.5,width = 2.5)
