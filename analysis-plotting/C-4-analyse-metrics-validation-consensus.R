# setup -------------------------------------------------------------------

source("function-library.R")

args <- commandArgs(trailingOnly=TRUE)
output_str <- ifelse(length(args)>0,args[1],"subset")
dir.create(paste('figures',output_str,sep = "/"),showWarnings = F)

library(tidyverse)
library(cowplot)
library(ggpubr)
library(caret)

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

pattern <- ifelse(
  "full" == output_str,
  '^cv\\.',sprintf('^cv_%s\\.',gsub("full_","",output_str)))
all_auc_validation <- list()
for (file_path in list.files("../mile-vice/ev-scores-consensus/",pattern = pattern,full.names = T)) {
  data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
  tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
  N <- as.numeric(tmp$value[tmp$metric == "N_0"])
  tmp <- tmp %>%
    subset(metric == "AUC_0")
  all_auc_validation[[file_path]] <- tibble(
    task = decode_model_name(tmp$model_id),dataset = data_type,value = tmp$value,N=N,obj="Single")
}

for (file_path in list.files("../mile-vice/ev-scores-consensus/",pattern = "^mo_cv_subset\\.",full.names = T)) {
  data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
  tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
  N <- tmp$value[tmp$metric == "N_0"]
  tmp <- tmp %>%
    subset(metric == "AUC_0")
  x <- str_split(file_path,'\\.')[[1]]
  x <- x[length(x)-1]
  all_auc_validation[[file_path]] <- tibble(
    task = decode_model_name(x),dataset=data_type,value = tmp$value,N=N,obj="Multiple")
}

all_auc_validation_df <- all_auc_validation %>%
  do.call(what = rbind) %>%
  mutate(set = "External validation") %>%
  mutate(N = as.numeric(N),value = as.numeric(as.character(value)))
all_auc_validation_df_so <- all_auc_validation_df %>%
  subset(obj == 'Single')
all_auc_validation_df_mo <- all_auc_validation_df %>%
  subset(obj == 'Multiple')

metrics_cv <- read_csv(sprintf("data_output/best_models_so_%s.csv",output_str)) %>%
  select(-`...1`,-nvc,-multi_objective)
metrics_cv %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_df_so,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = pmin(1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_df_so,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_df_so,position = position_dodge(width = 0.9),
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
  ggsave(filename = sprintf("figures/%s/auc-bars-w-validation-mile-vice-consensus.pdf",output_str),
         height=1.7,width=3)

metrics_cv <- read_csv(sprintf("data_output/best_models_so_%s.csv",output_str)) %>%
  select(-`...1`,-nvc,-multi_objective)
metrics_cv %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_df_mo,position = position_dodge(width = 0.9),
                 aes(ymin = value - 1/sqrt(N),
                     ymax = pmin(1,value + 1/sqrt(N))),
                 alpha = 1,size = 0.5,shape = 5) + 
  geom_point(data = all_auc_validation_df_mo,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = all_auc_validation_df_mo,position = position_dodge(width = 0.9),
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
  ggsave(filename = sprintf("figures/%s/auc-bars-w-validation-mile-vice-mo-consensus.pdf",output_str),
         height=1.7,width=3)
