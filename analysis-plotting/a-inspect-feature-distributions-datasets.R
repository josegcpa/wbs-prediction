source("function-library.R")

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

RBC_data <- rbind(
  read_csv("../mile-vice/data_output/rbc_adden_1_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "AC1"),
  read_csv("../mile-vice/data_output/rbc_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "MLL"),
  read_csv("../mile-vice/data_output/rbc_adden_2_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "AC2")) %>%
  gather(key = "key",value = "value",-slide_id,-moment,-dataset) %>%
  merge(all_conditions,by = "slide_id")

WBC_data <- rbind(
  read_csv("../mile-vice/data_output/wbc_adden_1_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "AC1"),
  read_csv("../mile-vice/data_output/wbc_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "MLL"),
  read_csv("../mile-vice/data_output/wbc_adden_2_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "AC2")) %>%
  gather(key = "key",value = "value",-slide_id,-moment,-dataset) %>%
  merge(all_conditions,by = "slide_id")

RBC_data %>%
  subset(dataset %in% c("AC2","MLL") & fine_class == "Normal") %>% 
  subset(moment %in% c("mean","variance")) %>%
  mutate(key = factor(key,levels = rev(c(features_all,features_nuclear)))) %>%
  group_by(key,moment) %>% 
  mutate(value = scale(value)) %>% 
  group_by(key,moment,dataset) %>%
  summarise(m = median(value),q05 = quantile(value,0.05),q95 = quantile(value,0.95)) %>%
  group_by(key,moment) %>%
  mutate(d = abs(diff(m))) %>% 
  ggplot(aes(x = reorder(key,d),y = m,ymin = q05,ymax = q95,colour  = dataset)) + 
  geom_point(position = position_dodge(width = 0.9),shape = 3) +
  geom_linerange(position = position_dodge(width = 0.9)) +
  theme_pretty(base_size = 6) + 
  coord_flip() +
  facet_wrap(~ moment) +
  xlab("Fetures") + 
  ylab("Standardized mean") + 
  ggsave("figures/feture_density.pdf",height = 5,width = 5)

wide_WBC_data <- WBC_data %>%
  subset(moment == "mean" & dataset %in% c("MLL","AC2")) %>%
  spread(key = "key",value = "value")
WBC_PCA <- cbind(wide_WBC_data[,c(1:5)],prcomp(wide_WBC_data[,-c(1:5)])$x) 
ggplot(WBC_PCA,aes(x=PC2,y=PC3,colour=dataset)) + 
  geom_point(alpha=0.5) + facet_wrap(~fine_class) + 
  theme_pretty(base_size = 6)

auc_list <- list()
ratio_list <- list()
for (k in unique(WBC_data$key)) {
  X <- WBC_data %>%
    subset(key == k) %>% 
    subset(dataset != "AC2")
  
  n1 <- sum(subset(X,moment == "mean")$dataset == "MLL")
  n2 <- sum(subset(X,moment == "mean")$dataset == "AC1")
  a <- wilcox.test(value ~ dataset,data = subset(X,moment == "mean"),
                   paired=F)
  b <- wilcox.test(sqrt(value) ~ dataset,data = subset(X,moment == "variance"),
                   paired=F)
  
  a_auc <- a$statistic/(n1*n2)
  b_auc <- b$statistic/(n1*n2)
  
  auc_list[[length(auc_list)+1]] <- data.frame(
    feature = k,
    mean_auc = a_auc,
    var_auc = b_auc,
    dataset_pair = "MLL vs. AC1")
  ratio_list[[length(ratio_list)+1]] <-   X %>%
    subset(moment == "mean" | moment == "variance") %>% 
    group_by(dataset,key,moment) %>% 
    summarise(mean = mean(value),.groups = "drop")
}
for (k in unique(WBC_data$key)) {
  X <- WBC_data %>%
    subset(key == k) %>% 
    subset(dataset != "AC1")
  
  n1 <- sum(subset(X,moment == "mean")$dataset == "MLL")
  n2 <- sum(subset(X,moment == "mean")$dataset == "AC2")
  a <- wilcox.test(value ~ dataset,data = subset(X,moment == "mean"),
                   paired=F)
  b <- wilcox.test(sqrt(value) ~ dataset,data = subset(X,moment == "variance"),
                   paired=F)
  
  a_auc <- a$statistic/(n1*n2)
  b_auc <- b$statistic/(n1*n2)
  
  auc_list[[length(auc_list)+1]] <- data.frame(
    feature = k,
    mean_auc = a_auc,
    var_auc = b_auc,
    dataset_pair = "MLL vs. AC2")
  ratio_list[[k]] <-   X %>%
    subset(moment == "mean" | moment == "variance") %>% 
    group_by(dataset,key,moment) %>% 
    summarise(mean = mean(value),.groups = "drop")
}

do.call(rbind,auc_list) %>%
  mutate(I = grepl("mass_disp|invariant",feature)) %>% 
  subset(dataset_pair == "MLL vs. AC2") %>%
  ggplot(aes(x = mean_auc,y = var_auc)) + 
  geom_point() + 
  geom_hline(yintercept = 0.5,size = 0.25,linetype = 2) +
  geom_vline(xintercept = 0.5,size = 0.25,linetype = 2) + 
  geom_abline(slope = 1,linetype = 3,alpha = 0.5, size = 0.5) +
  theme_pretty(base_size = 6) + 
  scale_x_continuous(expand = c(0.01,0)) + 
  scale_y_continuous(expand = c(0.01,0)) + 
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) +
  ylab("AUC for variances") +
  xlab("AUC for means")

do.call(rbind,ratio_list) %>%
  spread(key = "dataset",value = "mean") %>%
  subset(dataset_pair == "MLL vs. AC2") %>%
  ggplot(aes(x = MLL,y = AC2)) + 
  geom_point() + 
  geom_abline(slope = 1,linetype = 3,alpha = 0.5, size = 0.5) +
  theme_pretty(base_size = 6) + 
  facet_wrap(~ moment, scales = "free") +
  ylab("AC2") +
  xlab("MLL") + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')

do.call(rbind,ratio_list) %>%
  spread(key = "dataset",value = "mean") %>%
  mutate(ratio = MLL/AC2) %>% 
  mutate(ratio = ifelse(ratio > 10,10,ratio)) %>%
  mutate(key = factor(key,levels = c(features_all,features_nuclear))) %>% 
  mutate(moment = c(mean = "Means",variance = "Variances")[moment]) %>% 
  ggplot(aes(x = key,y = ratio)) + 
  geom_hline(yintercept = 1,size = 0.5,linetype = 3, alpha = 0.5) +
  geom_point(shape = 5) + 
  theme_pretty(base_size = 6) + 
  facet_wrap(~ moment, scales = "free") +
  ylab("Ratio of averages") +
  xlab("") + 
  coord_flip() +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.1,1,10),
                     labels = c(0.1,"1.0",">10"))
