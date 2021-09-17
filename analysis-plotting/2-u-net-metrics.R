# setup -------------------------------------------------------------------


source("function-library.R")

library(tidyverse)
library(ggplot2)
library(ggsci)

margin <- ggplot2::margin

dataset_conversion <- list(
  train = c(
    ADDENBROOKES = "AC1",
    MUNICH = "MLLC"
  ),
  test = c(
    ADDENBROOKES = "AC1",
    MUNICH = "MLLC"
  )
)

# plotting ----------------------------------------------------------------

data <- read_csv("data/test_results_unet_full_.csv") %>%
  subset(SAE == FALSE & TRANSFORMED == F) %>%
  select(-SAE,-TRANSFORMED)

data_ac2 <- rbind(
  read_csv("data/test_results_unet_adden_2.csv",
           col_names = c("Set","Metric","Metric_desc","value")) %>%
    mutate(model = "Trained on AC1"),
  read_csv("data/test_results_unet_adden_2_ft.csv",
           col_names = c("Set","Metric","Metric_desc","value")) %>%
    mutate(model = "Trained on AC1,\nfinetuned on AC2")
)

data %>%
  mutate(TRAIN_DATASET = dataset_conversion$train[TRAIN_DATASET],
         TEST_DATASET = dataset_conversion$test[TEST_DATASET]) %>%
  mutate(FINETUNE = ifelse(FINETUNE,"Trained on AC1,\nfinetuned on MLLC","Trained on AC1")) %>% 
  mutate(TUMBLE = ifelse(TUMBLE,"Yes","No")) %>% 
  mutate(TUMBLE = factor(TUMBLE,levels = c("Yes","No"))) %>% 
  subset(TRAIN_DATASET == "AC1") %>% 
  mutate(DEPTH = sprintf("Depth mult. = %.2f",DEPTH)) %>% 
  mutate(DEPTH = factor(DEPTH,levels = sprintf("Depth mult. = %.2f",c(0.1,0.25,0.5,1.0)))) %>%
  ggplot(aes(x = TEST_DATASET,y = MeanIOU,fill = TUMBLE)) + 
  geom_bar(stat = "identity",position = "dodge",colour = "black",size = 0.25) + 
  facet_wrap(~ FINETUNE + DEPTH,ncol = 4) + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by = 0.1)) + 
  coord_cartesian(ylim = c(0.6,1.0)) +
  theme_pretty(base_size = 6) + 
  xlab("Testing set") +
  ylab("Jaccard index") + 
  scale_fill_aaas(name = "Test-time data\naugmentation") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.1,"cm"),
        strip.text = element_text(margin = margin())) +
  ggsave("figures/u-net-metrics-full.pdf",height = 4.2,width = 3.8) 

data %>%
  mutate(TRAIN_DATASET = dataset_conversion$train[TRAIN_DATASET],
         TEST_DATASET = dataset_conversion$test[TEST_DATASET]) %>%
  mutate(FINETUNE = ifelse(FINETUNE,"Trained on AC1,\nfinetuned on MLLC","Trained on AC1")) %>% 
  mutate(TUMBLE = ifelse(TUMBLE,"Yes","No")) %>% 
  mutate(TUMBLE = factor(TUMBLE,levels = c("Yes","No"))) %>% 
  subset(TRAIN_DATASET == "AC1") %>% 
  mutate(DEPTH = sprintf("Depth mult. = %.2f",DEPTH)) %>% 
  mutate(DEPTH = factor(DEPTH,levels = sprintf("Depth mult. = %.2f",c(0.1,0.25,0.5,1.0)))) %>%
  ggplot(aes(x = `AUC`,y = MeanIOU,colour = TUMBLE,group = paste(TUMBLE,DEPTH),shape = TEST_DATASET)) + 
  geom_line(size = 0.25,aes(linetype = DEPTH)) +
  geom_point(stat = "identity",position = "dodge",size = 1.0) + 
  facet_wrap(~ FINETUNE,scales = "free") + 
  scale_y_continuous(trans = 'log10') + 
  scale_x_continuous(trans = 'log10') + 
  #coord_cartesian(ylim = c(0.8,1.0)) +
  theme_pretty(base_size = 6) + 
  xlab("AUC") +
  ylab("Jaccard index") + 
  scale_colour_aaas(name = "Test-time data\naugmentation") + 
  scale_shape(name = "Testing\ndataset") + 
  scale_linetype_manual(name = NULL,values = c(3,4,2,1)) +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),
        strip.text = element_text(margin = margin())) +
  guides(linetype = guide_legend(nrow = 4,keywidth = 1),
         shape = guide_legend(nrow = 2,keywidth = 0.5),
         colour = guide_legend(nrow = 2,keywidth = 0.5)) +
  ggsave("figures/u-net-metrics-scatter-full.pdf",height = 2,width = 3.2) 

data %>%
  mutate(TRAIN_DATASET = dataset_conversion$train[TRAIN_DATASET],
         TEST_DATASET = dataset_conversion$test[TEST_DATASET]) %>%
  mutate(FINETUNE = ifelse(FINETUNE,"Trained on AC1,\nfinetuned on MLLC","Trained on AC1")) %>% 
  mutate(TUMBLE = ifelse(TUMBLE,"Yes","No")) %>% 
  mutate(TUMBLE = factor(TUMBLE,levels = c("Yes","No"))) %>% 
  subset(DEPTH == 0.25 & TRAIN_DATASET == "AC1") %>% 
  ggplot(aes(x = TEST_DATASET,y = MeanIOU,fill = TUMBLE)) + 
  geom_bar(stat = "identity",position = "dodge",colour = "black",size = 0.25) + 
  facet_wrap(~ FINETUNE) + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by = 0.1)) + 
  coord_cartesian(ylim = c(0.7,1.0)) +
  theme_pretty(base_size = 6) + 
  xlab("Testing set") +
  ylab("Jaccard index") + 
  scale_fill_aaas(name = "Test-time data\naugmentation") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.1,"cm"),
        strip.text = element_text(margin = margin(t = 2))) +
  ggsave("figures/u-net-metrics.pdf",height = 2.1,width = 1.6) 

data_ac2 %>% 
  subset(Metric == "IOU" & Metric_desc == "global") %>%
  ggplot(aes(x = model,y = value)) + 
  geom_bar(stat = "identity",position = "dodge",fill = "black",colour = NA,size = 0.25) + 
  scale_y_continuous(expand = c(0,0),labels = function(x) sprintf("%.1f%%",x*100),
                     breaks = c(0.845,0.85)) + 
  theme_pretty(base_size = 6) + 
  ylab("Jaccard index") + 
  scale_fill_aaas(name = "Test-time data\naugmentation") + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.1,"cm"),
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(t = 2))) +
  coord_flip(ylim = c(0.845,0.85)) +
  ggsave("figures/u-net-metrics-ac2.pdf",height = 0.7,width = 1.25) 

