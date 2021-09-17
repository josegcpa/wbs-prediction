# TODO: 

# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(WRS2)
library(umap)
library(ggpubr)
library(caret)
library(glmnet)
library(pROC)
library(RRF)
library(MLmetrics)

select <- dplyr::select

set.seed(42)

task_conversion <- c(
  cv_binary = "Disease detection",
  cv_disease_binary = "Disease classification",
  cv_mds_binary = "SF3B1mut detection",
  cv_anemia_binary = "Anaemia classification",
  cv_multi_objective = "cv_multi_objective"
)

task_conversion_mo <- c("cv_anemia_binary","cv_binary","cv_mds_binary","cv_disease_binary")

all_metrics <- read_csv(
  "../mile-vice/cv_metrics.csv",
  col_names = c("set","fold","metric","value","nvc","nc","task_number","task","dataset")) %>%
  mutate(multi_objective = ifelse(
    task == "cv_multi_objective","Multiple objective","Single objective")) %>% 
  mutate(multi_objective = factor(
    multi_objective,levels = c("Single objective","Multiple objective"))) %>%
  mutate(task = ifelse(task == "cv_multi_objective",task_conversion_mo[task_number+1],task)) %>%
  mutate(task = factor(task_conversion[task],levels = rev(task_conversion))) %>%
  mutate(dataset_pretty = c(cells = "Morphology",cells_bc = "Morphology + B.C.")[dataset])

# cv metrics --------------------------------------------------------------

all_metrics %>%
  subset(metric == "AUC_0" & dataset == "cells" & set == "TEST") %>%
  group_by(nvc,task,multi_objective) %>%
  summarise(value = mean(value)) %>%
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

all_metrics %>%
  subset(metric == "AUC_0" & set == "TEST") %>%
  group_by(nvc,task,multi_objective,dataset,dataset_pretty) %>%
  summarise(value = mean(value)) %>%
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
  facet_grid(dataset_pretty ~ multi_objective,scales = "free_x",space = "free_x") +
  ylab("Number of virtual cells") + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  scale_colour_manual(values = c("black",NA),guide = F) +
  scale_fill_gradient(low = "lightcyan1",high = "purple",name = "AUC") +
  ggsave("figures/mile-vice-cv-performance-tile.pdf",height = 1.7,width = 3.5)

all_metrics %>%
  subset(metric == "AUC_0" & dataset == "cells" & set == "TEST") %>%
  group_by(nvc,task,multi_objective) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = nvc,y = value,
             colour = task,group = paste(task,multi_objective),
             shape = multi_objective,linetype = multi_objective)) + 
  geom_point(position = position_dodge(width = 0.7)) +
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
  ggsave("figures/mile-vice-cv-performance-lines.pdf",height = 1.8,width = 2.5)

best_models <- rbind(
  all_metrics %>%
    subset(metric == "AUC_0" & set == "TEST"& multi_objective == "Single objective") %>%
    group_by(nvc,task,multi_objective,dataset,dataset_pretty) %>%
    summarise(value = mean(value)) %>%
    mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
    group_by(task,multi_objective,dataset) %>% 
    filter(value == max(value)) %>%
    as_tibble(),
  all_metrics %>%
    subset(metric == "AUC_0" & set == "TEST"& multi_objective == "Multiple objective") %>%
    group_by(nvc,task,multi_objective,dataset,dataset_pretty) %>%
    summarise(value = mean(value)) %>%
    mutate(nvc = factor(nvc,levels = sort(unique(nvc)))) %>% 
    group_by(nvc,multi_objective,dataset) %>% 
    mutate(best_average = mean(value)) %>%
    group_by(dataset) %>%
    filter(best_average == max(best_average)) %>%
    select(-best_average) %>%
    as_tibble()
)

# compare with glmnet -----------------------------------------------------

glmnet_scores <- read.csv("data_output/glmnet-auroc.csv") %>% 
  select(dataset = data_type,task,glmnet_auc = auc_value) %>%
  distinct %>%
  merge(select(best_models,-dataset),
        by.y = c("task","dataset_pretty"),by.x = c("task","dataset")) %>%
  gather(key = "key",value = "value",value,glmnet_auc) %>%
  mutate(key = ifelse(key == "glmnet_auc","glmnet","MILe-ViCe")) %>%
  mutate(mo = multi_objective == "Multiple objective") %>%
  subset(!(mo == T & key == "glmnet")) %>%
  mutate(task = factor(task,
                       levels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification") %>% 
                         rev))

ggplot(glmnet_scores,
       aes(x = task,y = value,
           fill = key,colour = mo)) + 
  geom_bar(stat = "identity",
           aes(group = reorder(paste(key,mo),mo*1000 + as.numeric(key))),
           position = position_dodge(width = 0.9)) + 
  scale_y_continuous(labels = function(x) sprintf("%.1f%%",x*100)) + 
  ylab("AUC") +
  xlab("") +
  coord_flip(ylim = c(0.8,1)) +
  theme_pretty(base_size = 6) + 
  facet_wrap(~ dataset) +
  scale_colour_manual(guide = F,values = c(NA,"black")) +
  scale_fill_manual(name = NULL,values = c("lightblue","goldenrod"))
