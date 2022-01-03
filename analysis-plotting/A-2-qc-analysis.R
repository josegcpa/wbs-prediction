# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(WRS2)
library(ggpubr)
library(progress)
library(glmnet)

quality_dataset_labels <- read_csv(
  "datasets/blur-quality-net.csv",col_names = c("tile_id","class","blurriness",paste("HIST",seq(0,255))))

quality_dataset_labels_train <- quality_dataset_labels[1:8050,]
quality_dataset_labels_test <- quality_dataset_labels[8051:nrow(quality_dataset_labels),]
X <- as.matrix(quality_dataset_labels_train[,3:ncol(quality_dataset_labels_train)])
X_norm_params <- list(
  mean = colMeans(X),std = apply(X,2,sd)
)

simple_features_model <- cv.glmnet(
  (X - X_norm_params$mean)/(X_norm_params$std),
  quality_dataset_labels_train$class,family = "binomial",trace.it = T)

X_test <- as.matrix(quality_dataset_labels_test[,3:ncol(quality_dataset_labels_test)])
assess.glmnet(simple_features_model,
              newx = (X_test - X_norm_params$mean)/(X_norm_params$std),
              newy = quality_dataset_labels_test$class)

quality_dataset_labels_train %>%
  mutate(class = c("Poor","Good")[class+1]) %>%
  ggplot(aes(y = class,x = blurriness)) + 
  geom_jitter(size = 0.25,alpha = 0.5) + 
  geom_boxplot(alpha = 0.9,outlier.alpha = 0,size = 0.25) + 
  scale_x_continuous(trans = 'log10') +
  xlab("Tile sharpness") + 
  ylab("Tile quality") + 
  theme_pretty(base_size = 6) + 
  ggsave("figures/sharpness-distribution.pdf",height = 1,width = 3)
 
qc_summaries <- read_csv(
  "datasets/qc_summary.csv",
  col_names = c("slide_id","dataset","good_quality_tiles","total_tiles"))

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

blood_parameters <- read_csv(
  "data/blood_count_data.csv",
  col_types = cols_only(
    id = col_character(),wbc_ul = col_double(),
    hb_g_dl = col_double(),plt_ul = col_double())) %>%
  select(slide_id = id,wbc_ul,hb_g_dl,plt_ul)

fine_to_not_so_fine <- c(
  `Normal` = "Normal",`SF3B1-mutant` = "SF3B1-mutant",
  `RUNX1-mutant` = "Non-SF3B1-mutant",
  `SRSF2-mutant` = "Non-SF3B1-mutant",
  `U2AF1-mutant` = "Non-SF3B1-mutant",
  `Iron deficiency` = "Iron deficiency",
  `Megaloblastic` = "Megaloblastic")

qc_conditions <- merge(qc_summaries,all_conditions,by = "slide_id") %>%
  mutate(p = good_quality_tiles/total_tiles) %>%
  mutate(fine_class = fine_to_not_so_fine[fine_class])

qc_conditions %>%  
  mutate(dataset = factor(dataset,
                          levels = c("ADDEN2","MLL","ADDEN1"),
                          labels = c("AC2","MLLC","AC1"))) %>% 
  ggplot(aes(x = coarse_class,y = p,colour = fine_class,
             shape = dataset,group = dataset)) + 
  geom_point(position = position_jitterdodge(
    jitter.width = 0.4,seed = 42,jitter.height = 0,dodge.width = 0.6),
    size = 0.25,alpha = 0.3) + 
  stat_summary(geom = "errorbar",
               fun.data = function(x) return(data.frame(ymin = quantile(x,0.25),ymax = quantile(x,0.75))),
               size = 0.4,colour = "black",position = position_dodge(0.6),width = 0.2,
               alpha = 0.9) +
  stat_summary(geom = "point",
               fun.data = function(x) return(c(y = median(x))),
               size = 1,
               colour = "black",position = position_dodge(0.6),
               alpha = 0.9) +
  scale_y_continuous(labels = function(x) sprintf("%s%%",x*100)) + 
  theme_pretty(base_size = 6) + 
  ylab("Percentage of good quality tiles per slide") + 
  theme(axis.title.y = element_blank(),
        legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank(),
        legend.margin = margin()) + 
  scale_colour_manual(values = fine_colours) + 
  scale_shape(breaks = c("AC1","MLLC","AC2")) +
  coord_flip() +
  ggsave("figures/qc-slide.pdf",height = 1,width = 3.5) 

glm(p ~ dataset + fine_class,data = qc_conditions) %>% 
  summary

glm(p ~ dataset - 1,data = qc_conditions[qc_conditions$coarse_class == "MDS",]) %>% 
  summary
