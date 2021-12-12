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

df <- read_csv("data/test_results_unet_tf2.csv",
                 col_names = c("dataset","epoch","depth","TTA","refine","blank","metric","blank2","value")) %>%
  select(-blank,-blank2) %>%
  mutate(dataset = c(adden_1 = "Adden1",adden_2 = "Adden2","mll" = "MLL")[dataset]) %>%
  mutate(dataset = factor(dataset,levels = c("Adden1","MLL","Adden2"))) %>% 
  mutate(depth = factor(depth,levels = c(0.25,0.5,1.0),labels = c("Depth mult. = 0.25","Depth mult. = 0.5","Depth mult. = 1.0"))) %>%
  mutate(epoch = factor(epoch,levels = c("50","100"),labels = c("Early stopping\n(at 50 epochs)","Trained\nfor 100 epochs"))) %>%
  mutate(TTA = ifelse(TTA == 1, "TTA","No TTA"))

MAX <- max(df$value[df$metric == "IOU"])
MIN <- 0.5
B <- c(0,0.2,0.5,0.6,0.8,0.9)

# basic 
df %>%
  subset(TTA == "No TTA" & refine == F & metric == "IOU" & 
           epoch == "Trained\nfor 100 epochs") %>% 
  ggplot(aes(x = dataset,y = value)) +
  geom_bar(stat = "identity",position = "dodge") + 
  facet_wrap(~ depth) + 
  theme_pretty(base_size = 6) +
  theme(legend.key.width = unit(0.1,"cm"),
        legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),legend.margin = margin()) +
  scale_y_continuous(expand = c(0,0,0.04,0),breaks = B) +
  coord_cartesian(ylim = c(MIN,MAX)) +
  xlab("Cohort") + 
  ylab("Mean IoU") + 
  ggsave("figures/u-net-metrics-basic.pdf",height=1.3,width=3)

# early vs late

df %>%
  subset(TTA == "No TTA" & refine == F & metric == "IOU") %>% 
  ggplot(aes(x = dataset,y = value,fill = epoch)) +
  geom_bar(stat = "identity",position = "dodge") + 
  facet_wrap(~ depth) + 
  theme_pretty(base_size = 6) +
  theme(legend.key.width = unit(0.1,"cm"),
        legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),legend.margin = margin()) +
  scale_y_continuous(expand = c(0,0,0.04,0),breaks = B) +
  scale_fill_aaas(name = NULL) + 
  coord_cartesian(ylim = c(MIN,MAX)) +
  xlab("Cohort") + 
  ylab("Mean IoU") + 
  ggsave("figures/u-net-metrics-early-late.pdf",height=1.5,width=3)

# TTA vs. no TTA

df %>%
  subset(refine == F & metric == "IOU" & epoch == "Trained\nfor 100 epochs") %>% 
  ggplot(aes(x = dataset,y = value,fill = TTA)) +
  geom_bar(stat = "identity",position = "dodge") + 
  facet_wrap(~ depth) + 
  theme_pretty(base_size = 6) +
  theme(legend.key.width = unit(0.1,"cm"),
        legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),legend.margin = margin()) +
  scale_y_continuous(expand = c(0,0,0.04,0),breaks = B) +
  scale_fill_jama(name = NULL) + 
  coord_cartesian(ylim = c(MIN,MAX)) +
  xlab("Cohort") + 
  ylab("Mean IoU") + 
  ggsave("figures/u-net-metrics-tta.pdf",height=1.5,width=3)

# Refine vs. no refine

df %>%
  subset(TTA == "TTA" & metric == "IOU" & epoch == "Trained\nfor 100 epochs") %>% 
  mutate(refine = factor(refine,c(0,1),labels = c("No prediction\npost-processing","Prediction\npost-processing"))) %>% 
  ggplot(aes(x = dataset,y = value,fill = refine)) +
  geom_bar(stat = "identity",position = "dodge") + 
  facet_wrap(~ depth) + 
  theme_pretty(base_size = 6) +
  theme(legend.key.width = unit(0.1,"cm"),
        legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),legend.margin = margin()) +
  scale_y_continuous(expand = c(0,0,0.04,0),breaks = B) +
  scale_fill_simpsons(name = NULL) + 
  coord_cartesian(ylim = c(MIN,MAX)) +
  xlab("Cohort") + 
  ylab("Mean IoU") + 
  ggsave("figures/u-net-metrics-refine.pdf",height=1.5,width=3)
