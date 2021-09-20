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

all_roc_validation <- list()
for (file_path in list.files("../mile-vice/predictions/",full.names = T,pattern = "adden_2")) {
  tmp <- merge(read_prediction(file_path),all_conditions,all=F,by="slide_id") %>%
    mutate(labels_binary = label_conversion[[task[1]]][as.character(fine_class)]) %>%
    subset(!is.na(labels_binary))
  plot(roc(tmp$labels_binary,tmp$p2),main = file_path)
}