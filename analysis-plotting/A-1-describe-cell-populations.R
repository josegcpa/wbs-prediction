# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(WRS2)
library(umap)
library(ggpubr)

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
  `Normal` = "Normal",`SF3B1-mutant` = "SF3B1-mutant MDS",
  `RUNX1-mutant` = "non-SF3B1-mutant MDS",
  `SRSF2-mutant` = "non-SF3B1-mutant MDS",
  `U2AF1-mutant` = "non-SF3B1-mutant MDS",
  `Iron deficiency` = "Iron deficiency anaemia",
  `Megaloblastic` = "Megaloblastic anaemia")

# cell counts -------------------------------------------------------------

qc_count <- read_csv("datasets/qc_summary.csv",
                     col_names = c("slide_id","dataset",
                                   "good_quality_tiles","total_tiles"))

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
  summarise(counts = sum(counts)) %>%
  subset(!(slide_id %in% poor_quality_slides))

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
  summarise(counts = sum(counts)) %>%
  subset(!(slide_id %in% poor_quality_slides))

rbind(cbind(rbc_counts,cell_type = "RBC"),
      cbind(wbc_counts,cell_type = "WBC")) %>% 
  merge(all_conditions,by = "slide_id") %>% 
  mutate(dataset = factor(dataset,levels = rev(c("AC1","MLLC","AC2")))) %>% 
  ggplot(aes(x = coarse_class,y = counts,colour = fine_class,
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
  scale_y_continuous(trans = 'log10',
                     breaks = c(10,100,1000,10000,100000),
                     labels = scientific) + 
  theme_pretty(base_size = 6) + 
  ylab("Number of cells") + 
  theme(axis.title.y = element_blank(),
        legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank()) + 
  scale_colour_manual(values = fine_colours) + 
  scale_shape(breaks = c("AC1","MLLC","AC2")) +
  facet_wrap(~ cell_type) +
  coord_flip() +
  ggsave("figures/no-of-cells.pdf",height = 1,width = 3.5) 

rbind(cbind(rbc_counts,cell_type = "RBC"),
      cbind(wbc_counts,cell_type = "WBC")) %>% 
  merge(all_conditions,by = "slide_id") %>% 
  mutate(dataset = factor(dataset,levels = rev(c("AC1","MLLC","AC2")))) %>% 
  mutate(fine_class = fine_to_not_so_fine[fine_class]) %>%
  mutate(fine_class = factor(fine_class,levels = rev(unique(fine_to_not_so_fine)))) %>%
  ggplot(aes(x = fine_class,y = counts,
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
  scale_y_continuous(trans = 'log10',
                     breaks = c(10,100,1000,10000,100000),
                     labels = scientific) + 
  theme_pretty(base_size = 6) + 
  ylab("Number of cells") + 
  theme(axis.title.y = element_blank(),
        legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank()) + 
  scale_shape(breaks = c("AC1","MLLC","AC2")) +
  facet_wrap(~ cell_type) +
  coord_flip() +
  ggsave("figures/no-of-cells-fine.pdf",height = 1.5,width = 3.5) 

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
  xlab("Number of RBC") + 
  ylab("Number of WBC") + 
  theme(legend.key.height = unit(0,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_blank()) + 
  ggsave("figures/no-of-cells-scatter.pdf",height = 1.5,width = 2.5) 

quantile(wbc_counts$counts,c(0,0.25,0.5,0.75,1))
quantile(rbc_counts$counts,c(0,0.25,0.5,0.75,1))

# density to feature ------------------------------------------------------

tmp <- merge(wbc_counts,qc_count,by = "slide_id") %>%
  merge(select(blood_parameters,slide_id,wbc_ul),by = "slide_id") %>%
  mutate(cellular_density = counts/good_quality_tiles) %>%
  na.omit()

tmp_rbc <- merge(rbc_counts,qc_count,by = "slide_id") %>%
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
  xlab("WBC density (WBC/good quality tile)") + 
  ylab("WBC counts (WBC/dL)") +
  ggsave("figures/wbc-density-concentration.pdf",height = 1.5,width = 1.9) 

tmp %>%
  ggplot(aes(x = good_quality_tiles,y = counts)) +
  geom_point(size = 0.25) + 
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') + 
  geom_smooth(method = "lm",formula = y ~ x,colour = "goldenrod",size = 0.25) + 
  theme_pretty(base_size = 6) +
  xlab("Detected WBC") + 
  ylab("Good quality slides") +
  ggsave("figures/wbc-good-quality.pdf",height = 1.5,width = 1.9) 

tmp_rbc %>%
  ggplot(aes(x = good_quality_tiles,y = counts)) +
  geom_point(size = 0.25) + 
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') + 
  geom_smooth(method = "lm",formula = y ~ x,colour = "goldenrod",size = 0.25) + 
  theme_pretty(base_size = 6) +
  xlab("Detected RBC") + 
  ylab("Good quality slides") +
  ggsave("figures/rbc-good-quality.pdf",height = 1.5,width = 1.9) 

cor_wbc <- pbcor(tmp$good_quality_tiles,tmp$counts,ci = T,seed = 42)
cor_rbc <- pbcor(tmp_rbc$good_quality_tiles,tmp_rbc$counts,ci = T,seed = 42)

cat(sprintf("robust R2 (Detected WBC vs. good quality tiles) = %.4f (CI = [%.4f,%.4f],p-value = %f)",
            cor_wbc$cor^2,cor_wbc$cor_ci[1]^2,cor_wbc$cor_ci[2]^2,
            cor_wbc$p.value))

cat(sprintf("robust R2 (Detected RBC vs. good quality tiles) = %.4f (CI = [%.4f,%.4f],p-value = %f)",
            cor_rbc$cor^2,cor_rbc$cor_ci[1]^2,cor_rbc$cor_ci[2]^2,
            cor_rbc$p.value))

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

dt_values_mean_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05),
            Total = length(unique(wbc_all_cells_summaries$features))) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features/Total,
             label = sprintf("%s/%s",N_features,Total))) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),size = 2,
             fill = "white") +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "% of stat.\nsign. features",
                      labels = function(x) sprintf("%.0f%%",x * 100),
                      limits = c(0,1),
                      breaks = c(0.25,0.5,0.75,1)) + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.45,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/wbc-feature-significant-dunn-heatmap-mean.pdf",height = 2,width = 2.5) 

dt_values_var_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05),
            Total = length(unique(wbc_all_cells_summaries$features))) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features/Total,
             label = sprintf("%s/%s",N_features,Total))) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),size = 2,
             fill = "white") +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "% of stat.\nsign. features",
                      labels = function(x) sprintf("%.0f%%",x * 100),
                      limits = c(0,1),
                      breaks = c(0.25,0.5,0.75,1)) + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.45,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/wbc-feature-significant-dunn-heatmap-var.pdf",height = 2,width = 2.5) 

wbc_mean_plot <- wbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_mean$feature[1:2]) %>%
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

wbc_var_plot <- wbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_var$feature[1:2]) %>%
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

dt_wbc <- rbind(dt_values_mean_df,dt_values_var_df) %>%
  as_tibble()

dt_wbc %>% 
  mutate(adj.p.value = p.adjust(p.value,method = "BH")) %>%
  group_by(feature,comparisons) %>% 
  filter(adj.p.value == min(adj.p.value)) %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(adj.p.value < 0.05)) %>%
  arrange(N_features) %>%
  mutate(Proportion = N_features / length(unique(dt_values_mean_df$feature)))

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

dt_values_mean_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05),
            Total = length(unique(rbc_all_cells_summaries$features))) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features/Total,
             label = sprintf("%s/%s",N_features,Total))) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),size = 2,
             fill = "white") +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "% of stat.\nsign. features",
                      labels = function(x) sprintf("%.0f%%",x * 100),
                      limits = c(0,1),
                      breaks = c(0.25,0.5,0.75,1)) + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.45,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/rbc-feature-significant-dunn-heatmap-mean.pdf",height = 2,width = 2.5) 

dt_values_var_df %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(p.adjust(p.value,method = "BH") < 0.05),
            Total = length(unique(rbc_all_cells_summaries$features))) %>%
  mutate(A = gsub(" -","",str_match(comparisons,'[A-za-z0-9 -]+ -')),
         B = gsub("- ","",str_match(comparisons,'- [A-za-z0-9 -]+'))) %>%
  ggplot(aes(x = A,y = B,fill = N_features/Total,
             label = sprintf("%s/%s",N_features,Total))) + 
  geom_tile() + 
  geom_label(label.padding = unit(0.05,"cm"),
             label.r = unit(0,"cm"),
             label.size = unit(0,"cm"),size = 2,
             fill = "white") +
  scale_fill_gradient(low = "goldenrod",high = "red4",
                      name = "% of stat.\nsign. features",
                      labels = function(x) sprintf("%.0f%%",x * 100),
                      limits = c(0,1),
                      breaks = c(0.25,0.5,0.75,1)) + 
  theme_pretty(base_size = 6) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.45,"cm")) +
  rotate_x_text(angle = 20) +
  ggsave("figures/rbc-feature-significant-dunn-heatmap-var.pdf",height = 2,width = 2.5) 

rbc_mean_plot <- rbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_mean$feature[1:2]) %>%
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

rbc_var_plot <- rbc_all_cells_summaries %>%
  subset(key %in% kt_values_df_var$feature[1:2]) %>%
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

 dt_rbc <- rbind(dt_values_mean_df,dt_values_var_df)
 
 dt_rbc %>% 
  mutate(adj.p.value = p.adjust(p.value,method = "BH")) %>%
  group_by(feature,comparisons) %>% 
  filter(adj.p.value == min(adj.p.value)) %>%
  group_by(comparisons) %>%
  summarise(N_features = sum(adj.p.value < 0.05)) %>%
  arrange(N_features) %>%
  mutate(Proportion = N_features / length(unique(dt_values_mean_df$feature)))

# big image (wbc + rbc) ---------------------------------------------------

rbind(
  mutate(dt_wbc,cell_type="WBC",Total = length(unique(wbc_all_cells_summaries$features))),
  mutate(dt_rbc,cell_type="RBC",Total = length(unique(rbc_all_cells_summaries$features)))) %>%
   mutate(adj.p.value = p.adjust(p.value,method = "BH")) %>%
   group_by(feature,comparisons,cell_type,Total) %>% 
   filter(adj.p.value == min(adj.p.value)) %>%
   group_by(comparisons,cell_type,Total) %>%
   summarise(N_features = sum(adj.p.value < 0.05)) %>%
   group_by(comparisons) %>%
   mutate(N_features_both = sum(N_features)) %>%
   arrange(N_features_both,comparisons,cell_type) %>%
   mutate(Proportion = N_features/Total)
 
LL <- get_legend(wbc_mean_plot + guides(fill = guide_legend(nrow = 3)))

plot_grid(
  wbc_mean_plot + theme(legend.position = "none"),
  rbc_mean_plot + theme(legend.position = "none"),
  LL,
  ncol = 1,
  rel_heights = c(1,1,0.3)
) + 
  ggsave("figures/dunn-bonferroni-wbc-rbc-mean.pdf",height = 3,width = 2.5) 

wilcox.test(
  formula = mean ~ ifelse(coarse_class=="Normal","Normal","Not"),
  data = subset(wbc_all_cells_summaries,features == "Ellipse variance")
)

wilcox.test(
  formula = mean ~ ifelse(coarse_class=="Normal","Normal","Not"),
  data = subset(wbc_all_cells_summaries,features == "CDF 1st noiseless mom.")
)

wilcox.test(
  formula = mean ~ fine_class,
  data = subset(wbc_all_cells_summaries,features == "Ellipse variance" & coarse_class == "Anaemia")
)

wilcox.test(
  formula = mean ~ fine_class,
  data = subset(wbc_all_cells_summaries,features == "CDF 1st noiseless mom." & coarse_class == "Anaemia")
)

wilcox.test(
  formula = mean ~ fine_class,
  data = subset(rbc_all_cells_summaries,features == "Perimeter" & coarse_class == "MDS")
)

wilcox.test(
  formula = mean ~ fine_class,
  data = subset(rbc_all_cells_summaries,features == "Min(CDF)" & coarse_class == "MDS")
)

plot_grid(
  wbc_var_plot + theme(legend.position = "none"),
  rbc_var_plot + theme(legend.position = "none"),
  LL,
  ncol = 1,
  rel_heights = c(1,1,0.2)
) + 
  ggsave("figures/dunn-bonferroni-wbc-rbc-var.pdf",height = 3,width = 2.5) 
