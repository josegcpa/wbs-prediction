# TODO: independent validation cohort

# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(caret)
library(pROC)

read_prediction <- function(file_path) {
  read_csv(file_path,
           col_names = c("slide_id","labels","p1","p2")) %>%
    mutate(task = gsub(".csv","",gsub("labels/","",labels))) %>%
    return
}

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

all_roc_validation <- list()
par(mfrow=c(2,2))
for (file_path in list.files("../mile-vice/predictions/",full.names = T,pattern = "adden_2")) {
  tmp <- merge(read_prediction(file_path),all_conditions,all=F,by="slide_id") %>%
    mutate(labels_binary = label_conversion[[task[1]]][as.character(fine_class)]) %>%
    subset(!is.na(labels_binary))
  print(file_path)
  print(plot(roc(tmp$labels_binary,tmp$p2)))
}
