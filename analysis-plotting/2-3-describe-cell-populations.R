# setup -------------------------------------------------------------------


source("function-library.R")

library(tidyverse)
library(cowplot)
library(WRS2)
library(umap)
library(ggpubr)

poor_quality_slides <- c(
  "SRSF2_10", # the whole slide is practically blurred
  "II_26","VII_11","VIII_10",
  "VIII_21","X_9","X_11","XV_11",
  "XV_19","XVII_1","XVII_3","XVII_9"
)

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




# cell counts -------------------------------------------------------------

qc_count <- read_csv("../mile-vice/data_output/qc_count.csv",
                     col_names = c("slide_id","good_quality_tiles","total_tiles"),
                     col_types = c(col_character(),col_integer(),col_integer()))

rbc_counts <-  rbind(
  read_csv(
    "../mile-vice/data_output/rbc_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "MLLC"),
  read_csv(
    "../mile-vice/data_output/rbc_counts_adden_1.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double()))  %>%
    cbind(dataset = "AC1"),
  read_csv(
    "../mile-vice/data_output/rbc_counts_adden_2.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "AC2")
) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  group_by(slide_id,dataset) %>%
  summarise(counts = sum(counts))

wbc_counts <-  rbind(
  read_csv(
    "../mile-vice/data_output/wbc_counts.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "MLLC"),
  read_csv(
    "../mile-vice/data_output/wbc_counts_adden_1.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double()))  %>%
    cbind(dataset = "AC1"),
  read_csv(
    "../mile-vice/data_output/wbc_counts_adden_2.csv",
    col_names = c("slide_id","counts"),
    col_types = c(col_character(),col_double())) %>%
    cbind(dataset = "AC2")
) %>%
  subset(!(slide_id %in% poor_quality_slides)) %>%
  group_by(slide_id,dataset) %>%
  summarise(counts = sum(counts))

range(wbc_counts$counts)
range(rbc_counts$counts)

rbind(cbind(rbc_counts,cell_type = "RBC"),
      cbind(wbc_counts,cell_type = "WBC")) %>% 
  merge(all_conditions,by = "slide_id") %>% 
  ggplot(aes(x = coarse_class,y = counts,colour = fine_class,
             shape = dataset)) + 
  geom_point(position = position_jitter(width = 0.3,seed = 42,height = 0),size = 0.25,alpha = 0.8) + 
  stat_summary(geom = "linerange",fun.data = median_hilow,size = 0.5,
               colour = "goldenrod",position = position_dodge(0.3)) +
  stat_summary(geom = "point",fun.data = function(x) return(c(y = median(x))),size = 1,
               colour = "goldenrod",position = position_dodge(0.3)) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(10,100,1000,10000,100000),
                     labels = scientific) + 
  theme_pretty(base_size = 6) + 
  ylab("Number of cells") + 
  theme(axis.title.x = element_blank(),
        legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank()) + 
  scale_colour_manual(values = fine_colours) + 
  scale_shape() +
  facet_wrap(~ cell_type) +
  ggsave("figures/no-of-cells.pdf",height = 1.5,width = 3.5) 

merge(select(rbc_counts,slide_id,dataset,rbc = counts),
      select(wbc_counts,slide_id,wbc = counts),
      by = "slide_id") %>% 
  merge(all_conditions,by = "slide_id") %>% 
  subset(!(slide_id %in% poor_quality_slides)) %>%
  ggplot(aes(x = rbc,y = wbc,shape = dataset)) + 
  geom_point(position = position_jitter(width = 0.3,seed = 42,height = 0),size = 0.25,alpha = 0.8) + 
  scale_y_continuous(trans = 'log10',
                     breaks = c(10,100,1000,10000,100000),
                     labels = scientific) + 
  scale_x_continuous(trans = 'log10',
                     breaks = c(10,100,1000,10000,100000),
                     labels = scientific) + 
  theme_pretty(base_size = 6) + 
  ylab("Number of cells") + 
  theme(axis.title.x = element_blank(),
        legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank()) + 
  ggsave("figures/no-of-cells-scatter.pdf",height = 1.5,width = 2.5) 




# density to feature ------------------------------------------------------

tmp <- merge(wbc_counts,qc_count,by = "slide_id") %>%
  merge(select(blood_parameters,slide_id,wbc_ul),by = "slide_id") %>%
  mutate(cellular_density = counts/good_quality_tiles) %>%
  na.omit()

rob_cor_est <- pbcor(tmp$wbc_ul,tmp$cellular_density,ci = T,seed = 42)

cat(sprintf("robust R2 = %.4f (CI = [%.4f,%.4f],p-value = %f)",
           rob_cor_est$cor^2,rob_cor_est$cor_ci[1]^2,rob_cor_est$cor_ci[2]^2,
           rob_cor_est$p.value))

tmp %>%
  ggplot(aes(x = cellular_density,y = wbc_ul)) +
  geom_point(size = 0.25) + 
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') + 
  geom_smooth(method = "lm",formula = y ~ x,colour = "goldenrod",size = 0.25) + 
  theme_pretty(base_size = 6) +
  xlab("Detected WBC (per tile)") + 
  ylab("WBC (g/dL)") +
  ggsave("figures/wbc-density-concentration.pdf",height = 1.5,width = 1.7) 



# feature distribution per condition (wbc) --------------------------------

wbc_all_cells_summaries <- read_csv(
  "../mile-vice/data_output/wbc_summaries.csv",
  col_names = c("slide_id",features_all,features_nuclear,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all,features_nuclear)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>%
  gather(key = "key",value = "value",-slide_id,-coarse_class,-fine_class,-f) %>%
  spread(key = "f",value = "value") %>%
  as_tibble() %>% 
  mutate(fine_class = ifelse(coarse_class == "MDS",
                             ifelse(grepl("SF3B1",fine_class),
                                    "SF3B1-mutant",
                                    "Non-SF3B1-mutant"),
                             as.character(fine_class))) %>%
  mutate(fine_class = factor(
    as.character(fine_class),
    levels = fine_simple_levels)) %>%
  subset(!(slide_id %in% poor_quality_slides))
  
wbc_all_cells_summaries$features <- factor(
  features_conversion[wbc_all_cells_summaries$key],
  levels = features_conversion)

wbc_all_cells_summaries %>%
  subset(features %in% features_conversion[1:30]) %>%
  group_by(features) %>%
  mutate(mean = scale(mean)) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised mean") +
  ggsave("figures/wbc-mean-feature-distribution-1.pdf",height = 6,width = 6) 

wbc_all_cells_summaries %>%
  subset(features %in% features_conversion[31:length(features_conversion)]) %>%
  group_by(features) %>%
  mutate(mean = scale(mean)) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised mean") +
  ggsave("figures/wbc-mean-feature-distribution-2.pdf",height = 6,width = 6) 

wbc_all_cells_summaries %>%
  subset(features %in% features_conversion[1:30]) %>%
  group_by(features) %>%
  mutate(variance = scale(variance)) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised variance") +
  ggsave("figures/wbc-variance-feature-distribution-1.pdf",height = 6,width = 6) 

wbc_all_cells_summaries %>%
  subset(features %in% features_conversion[31:length(features_conversion)]) %>%
  group_by(features) %>%
  mutate(variance = scale(variance)) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised variance") +
  ggsave("figures/wbc-variance-feature-distribution-2.pdf",height = 6,width = 6) 

kt_values_mean <- list()
kt_values_var <- list()
for (k in unique(wbc_all_cells_summaries$key)) {
  tmp_data <- subset(wbc_all_cells_summaries,key == k)
  kt <- kruskal.test(mean ~ fine_class,data = tmp_data)
  kt_values_mean[[k]] <- data.frame(
    feature = k,
    statistic = kt$statistic,
    p.value = kt$p.value)
  kt <- kruskal.test(variance ~ fine_class,data = tmp_data)
  kt_values_var[[k]] <- data.frame(
    feature = k,
    statistic = kt$statistic,
    p.value = kt$p.value)
  
}

kt_values_df <- rbind(
  kt_values_mean %>%
    do.call(what = rbind) %>%
    cbind(moment = "mean"),
  kt_values_var %>%
    do.call(what = rbind) %>%
    cbind(moment = "var")) %>%
  subset(p.adjust(p.value,method = "BH") < 0.05) %>%
  arrange(-statistic)

kt_values_df_mean <- subset(kt_values_df,moment == "mean")
kt_values_df_var <- subset(kt_values_df,moment == "var")

dt_values_mean <- list()
dt_values_var <- list()
for (feature in kt_values_df$feature) {
  tmp_data <- subset(wbc_all_cells_summaries,key == feature) %>%
    arrange(fine_class)
  x <- dunn.test::dunn.test(tmp_data$mean,tmp_data$fine_class,
                            method = "none",kw = F)
  dt_values_mean[[feature]] <- data.frame(
    chi2 = x$Z,
    p.value = x$P,
    comparisons = x$comparisons,
    feature = feature
  )
  x <- dunn.test::dunn.test(tmp_data$variance,tmp_data$fine_class,
                            method = "none",kw = F)
  dt_values_var[[feature]] <- data.frame(
    chi2 = x$Z,
    p.value = x$P,
    comparisons = x$comparisons,
    feature = feature
  )
}

dt_values_mean_df <- do.call(rbind,dt_values_mean)
dt_values_var_df <- do.call(rbind,dt_values_var)

rbind(dt_values_mean_df,
      dt_values_var_df) %>% 
  mutate(adj.p.value = p.adjust(p.value,method = "BH")) %>%
  group_by(feature,comparisons) %>% 
  filter(adj.p.value == min(adj.p.value)) %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(adj.p.value < 0.05)) %>%
  arrange(N_features) %>%
  mutate(Proportion = N_features / 60)

dt_values_mean_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05)) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features,label = N_features)) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),
             size = 2.5) +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "No. of stat.\nsign. features") + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.3,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/wbc-feature-significant-dunn-heatmap-mean.pdf",height = 2,width = 2.5) 

dt_values_var_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05)) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features,label = N_features)) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),
             size = 2.5) +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "No. of stat.\nsign. features") + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.3,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/wbc-feature-significant-dunn-heatmap-var.pdf",height = 2,width = 2.5) 

wbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_mean$feature[1:4]) %>%
  mutate(key = factor(key,levels = kt_values_df_mean$feature[1:3])) %>% 
  group_by(features) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = fine_class)) + 
  geom_boxplot(size = 0.25,position = position_dodge(width = 0.85),
               outlier.size = 0.25) +
  facet_wrap(~ reorder(features,key),scales = "free",ncol = 4) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Mean") +
  rotate_x_text(angle = 40) +
  ggsave("figures/wbc-feature-distribution-subset.pdf",height = 1.7,width = 4) 

wbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_var$feature[1:4]) %>%
  mutate(key = factor(key,levels = kt_values_df_var$feature[1:3])) %>% 
  group_by(features) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = fine_class)) + 
  geom_boxplot(size = 0.25,position = position_dodge(width = 0.85),
               outlier.size = 0.25) +
  facet_wrap(~ reorder(features,key),scales = "free",ncol = 4) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Variance") +
  rotate_x_text(angle = 40) +
  ggsave("figures/wbc-feature-distribution-subset-var.pdf",height = 1.7,width = 4) 


# feature distribution per condition (rbc) --------------------------------

rbc_all_cells_summaries <- read_csv(
  "../mile-vice/data_output/rbc_summaries.csv",
  col_names = c("slide_id",features_all,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>%
  gather(key = "key",value = "value",-slide_id,-coarse_class,-fine_class,-f) %>%
  spread(key = "f",value = "value") %>%
  as_tibble() %>% 
  mutate(fine_class = ifelse(coarse_class == "MDS",
                             ifelse(grepl("SF3B1",fine_class),
                                    "SF3B1-mutant",
                                    "Non-SF3B1-mutant"),
                             as.character(fine_class))) %>%
  mutate(fine_class = factor(
    as.character(fine_class),
    levels = fine_simple_levels)) %>%
  subset(!(slide_id %in% poor_quality_slides))

rbc_all_cells_summaries$features <- factor(
  features_conversion[rbc_all_cells_summaries$key],
  levels = features_conversion)

rbc_all_cells_summaries %>%
  subset(features %in% features_conversion[1:30]) %>%
  group_by(features) %>%
  mutate(mean = scale(mean)) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised mean") +
  ggsave("figures/rbc-mean-feature-distribution-1.pdf",height = 6,width = 6) 

rbc_all_cells_summaries %>%
  subset(features %in% features_conversion[31:length(features_conversion)]) %>%
  group_by(features) %>%
  mutate(mean = scale(mean)) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised mean") +
  ggsave("figures/rbc-mean-feature-distribution-2.pdf",height = 6,width = 6) 

rbc_all_cells_summaries %>%
  subset(features %in% features_conversion[1:30]) %>%
  group_by(features) %>%
  mutate(variance = scale(variance)) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised variance") +
  ggsave("figures/rbc-variance-feature-distribution-1.pdf",height = 6,width = 6) 

rbc_all_cells_summaries %>%
  subset(features %in% features_conversion[31:length(features_conversion)]) %>%
  group_by(features) %>%
  mutate(variance = scale(variance)) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = coarse_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.25) +
  facet_wrap(~ features,scales = "free",ncol = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Standardised variance") +
  ggsave("figures/rbc-variance-feature-distribution-2.pdf",height = 6,width = 6) 

kt_values_mean <- list()
kt_values_var <- list()
for (k in unique(rbc_all_cells_summaries$key)) {
  tmp_data <- subset(rbc_all_cells_summaries,key == k)
  kt <- kruskal.test(mean ~ fine_class,data = tmp_data)
  kt_values_mean[[k]] <- data.frame(
    feature = k,
    statistic = kt$statistic,
    p.value = kt$p.value)
  kt <- kruskal.test(variance ~ fine_class,data = tmp_data)
  kt_values_var[[k]] <- data.frame(
    feature = k,
    statistic = kt$statistic,
    p.value = kt$p.value)
  
}

kt_values_df <- rbind(
  kt_values_mean %>%
    do.call(what = rbind) %>%
    cbind(moment = "mean"),
  kt_values_var %>%
    do.call(what = rbind) %>%
    cbind(moment = "var")) %>%
  subset(p.adjust(p.value,method = "BH") < 0.05) %>%
  arrange(-statistic)

kt_values_df_mean <- subset(kt_values_df,moment == "mean")
kt_values_df_var <- subset(kt_values_df,moment == "var")

dt_values_mean <- list()
dt_values_var <- list()
for (feature in kt_values_df$feature) {
  tmp_data <- subset(rbc_all_cells_summaries,key == feature) %>%
    arrange(fine_class)
  x <- dunn.test::dunn.test(tmp_data$mean,tmp_data$fine_class,
                            method = "none",kw = F)
  dt_values_mean[[feature]] <- data.frame(
    chi2 = x$Z,
    p.value = x$P,
    comparisons = x$comparisons,
    feature = feature
  )
  x <- dunn.test::dunn.test(tmp_data$variance,tmp_data$fine_class,
                            method = "none",kw = F)
  dt_values_var[[feature]] <- data.frame(
    chi2 = x$Z,
    p.value = x$P,
    comparisons = x$comparisons,
    feature = feature
  )
}

dt_values_mean_df <- do.call(rbind,dt_values_mean)
dt_values_var_df <- do.call(rbind,dt_values_var)

rbind(dt_values_mean_df,
      dt_values_var_df) %>% 
  mutate(adj.p.value = p.adjust(p.value,method = "BH")) %>%
  group_by(feature,comparisons) %>% 
  filter(adj.p.value == min(adj.p.value)) %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(adj.p.value < 0.05)) %>%
  arrange(N_features) %>%
  mutate(Proportion = N_features / 60)

dt_values_mean_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05)) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features,label = N_features)) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),
             size = 2.5) +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "No. of stat.\nsign. features") + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.3,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/rbc-feature-significant-dunn-heatmap-mean.pdf",height = 2,width = 2.5) 

dt_values_var_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05)) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features,label = N_features)) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),
             size = 2.5) +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "No. of stat.\nsign. features") + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.3,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/rbc-feature-significant-dunn-heatmap-var.pdf",height = 2,width = 2.5) 

rbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_mean$feature[1:4]) %>%
  mutate(key = factor(key,levels = kt_values_df_mean$feature[1:3])) %>% 
  group_by(features) %>%
  ggplot(aes(x = coarse_class,y = mean,
             fill = fine_class)) + 
  geom_boxplot(size = 0.25,position = position_dodge(width = 0.85),
               outlier.size = 0.25) +
  facet_wrap(~ reorder(features,key),scales = "free",ncol = 4) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Mean") +
  rotate_x_text(angle = 40) +
  ggsave("figures/rbc-feature-distribution-subset.pdf",height = 1.7,width = 4) 

rbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_var$feature[1:4]) %>%
  mutate(key = factor(key,levels = kt_values_df_var$feature[1:3])) %>% 
  group_by(features) %>%
  ggplot(aes(x = coarse_class,y = variance,
             fill = fine_class)) + 
  geom_boxplot(size = 0.25,position = position_dodge(width = 0.85),
               outlier.size = 0.25) +
  facet_wrap(~ reorder(features,key),scales = "free",ncol = 4) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL) +
  scale_y_continuous(labels = function(x) format(x,digits = 4)) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        axis.title.x = element_blank()) + 
  ylab("Variance") +
  rotate_x_text(angle = 40) +
  ggsave("figures/rbc-feature-distribution-subset-var.pdf",height = 1.7,width = 4) 









# u-map (not used for now) ------------------------------------------------

# wbc_all_cells <- read_csv(
#   "../mile-vice/data_output/wbc_all_cells.csv",
#   col_names = c("slide_id",features_all,features_nuclear),
#   col_types = c(list(col_character()),replicate(length(c(features_all,features_nuclear)),col_double()))) %>%
#   merge(all_conditions,by = "slide_id")
# 
# rbc_all_cells <- read_csv(
#   "../mile-vice/data_output/rbc_all_cells.csv",
#   col_names = c("slide_id",features_all),
#   col_types = c(list(col_character()),replicate(length(features_all),col_double()))) %>%
#   merge(all_conditions,by = "slide_id")
# 
# set.seed(42)
# wbc_subsample <- sample.int(nrow(wbc_all_cells),size = 10000,replace = F)
# wbc_tsne <- umap(wbc_all_cells[wbc_subsample,-c(1,62,63)],
#                  verbose=T,
#                  random_state=42)
# 
# cbind(data.frame(wbc_tsne$layout),
#       coarse_class = wbc_all_cells$coarse_class[wbc_subsample],
#       fine_class = wbc_all_cells$fine_class[wbc_subsample]) %>%
#   ggplot(aes(x = X1,y = X2)) + 
#   geom_density2d(colour = "black") +
#   geom_point(alpha = 0.1,size = 0.25,aes(colour = fine_class)) + 
#   scale_colour_manual(values = fine_colours,name = NULL) +
#   theme_pretty(base_size = 6) +
#   xlab("UMAP 1") +
#   ylab("UMAP 2") +
#   facet_wrap(~ coarse_class) +
#   theme(legend.position = "bottom",
#         legend.key.size = unit(0,"cm")) + 
#   ggsave("figures/tsne-wbc.pdf",height = 2,width = 5) 
# 
# rbc_subsample <- sample.int(nrow(wbc_all_cells),size = 10000,replace = F)
# rbc_tsne <- umap(rbc_all_cells[rbc_subsample,-c(1,44,45)],
#                  verbose=T,
#                  random_state=42,
#                  spread = 0.2)
# 
# cbind(data.frame(rbc_tsne$layout),
#       coarse_class = rbc_all_cells$coarse_class[rbc_subsample],
#       fine_class = rbc_all_cells$fine_class[rbc_subsample]) %>%
#   ggplot(aes(x = X1,y = X2)) + 
#   geom_density2d(colour = "black") +
#   geom_point(alpha = 0.1,size = 0.25,aes(colour = fine_class)) + 
#   scale_colour_manual(values = fine_colours,name = NULL) +
#   theme_pretty(base_size = 6) +
#   xlab("UMAP 1") +
#   ylab("UMAP 2") +
#   facet_wrap(~ coarse_class) +
#   theme(legend.position = "bottom",
#         legend.key.size = unit(0,"cm")) + 
#   ggsave("figures/tsne-rbc.pdf",height = 2,width = 5) 
# 
# 
