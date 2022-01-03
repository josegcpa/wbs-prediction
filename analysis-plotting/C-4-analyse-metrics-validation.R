# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(caret)
library(pROC)

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

metrics_cv <- read_csv("data_output/best_models_so.csv") %>%
  select(-X1,-nvc,-multi_objective)
all_auc_validation <- list()
par(mfrow=c(2,2))
for (file_path in list.files("../mile-vice/ev-scores/",pattern = "cv_subset\\.",full.names = T)) {
  data_type <- ifelse(grepl('bc',file_path),"Morphology + B.C.","Morphology")
  tmp <- read_csv(file_path,col_names = c("model_id","fold","metric","value"))
  N <- tmp$value[tmp$metric == "N_0"]
  tmp <- tmp %>%
    subset(metric == "AUC_0")
  all_auc_validation[[file_path]] <- data.frame(
    task = decode_model_name(tmp$model_id),dataset = data_type,value = tmp$value,N=N)
}

all_auc_validation_df <- mutate(do.call(rbind,all_auc_validation),set = "External validation")

ggplot() + 
  geom_bar(data = metrics_cv,aes(x = value,y = task),stat = "identity") + 
  geom_point()

metrics_cv %>% 
  ggplot(aes(x = factor(task,rev(model_levels)),y = value, fill = dataset)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) +
  geom_linerange(data = all_auc_validation_df,position = position_dodge(width = 0.9),
                 aes(ymin = auc_lower_ci(value,N/2,N/2),
                     ymax = auc_upper_ci(value,N/2,N/2)),
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
  ggsave(filename = "figures/auc-bars-w-validation-mile-vice.pdf",height=1.7,width=3)
