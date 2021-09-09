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
  cv_mds_binary = "SF3B1mut detection",
  cv_disease_binary = "Disease classification",
  cv_binary = "Disease detection",
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
  mutate(task = task_conversion[task])


# cv metrics --------------------------------------------------------------

all_metrics %>%
  subset(metric == "AUC_0" & dataset == "cells" & set == "TEST") %>%
  group_by(nvc,task,multi_objective) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(x = task,y = value,
             colour = as.factor(nvc),group = paste(multi_objective,nvc),
             shape = multi_objective)) + 
  geom_point(position = position_dodge(width = 0.7),size = 0.5) +
  theme_pretty(base_size = 6) +
  coord_flip() +
  scale_colour_brewer(type = "div",name = NULL) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0,"cm"),
        axis.title.y = element_blank()) + 
  ylab("AUC") + 
  scale_shape_manual(name = NULL,values = c(16,4)) + 
  guides(shape = guide_legend(nrow = 2)) + 
  ggsave("figures/mile-vice-cv-performance.pdf",height = 1.8,width = 3)

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

