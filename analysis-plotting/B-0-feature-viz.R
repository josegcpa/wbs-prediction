source("function-library.R")

library(umap)
library(cowplot)

set.seed(42)

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

RBC_data_all <- rbind(
  read_csv("datasets/rbc_adden_1_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "AC1"),
  read_csv("datasets/rbc_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "MLL"),
  read_csv("datasets/rbc_adden_2_summaries.csv",
           col_names = c("slide_id",features_all,"moment")) %>%
    mutate(dataset = "AC2")) %>%
  merge(all_conditions,by = "slide_id")
WBC_data_all <- rbind(
  read_csv("datasets/wbc_adden_1_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "AC1"),
  read_csv("datasets/wbc_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "MLL"),
  read_csv("datasets/wbc_adden_2_summaries.csv",
           col_names = c("slide_id",features_all,features_nuclear,"moment")) %>%
    mutate(dataset = "AC2")) %>%
  merge(all_conditions,by = "slide_id") 

RBC_cells <- rbind(
  read_csv("datasets/rbc_adden_1_all_cells.csv",
           col_names = c("slide_id",features_all)) %>%
    mutate(dataset = "AC1"),
  read_csv("datasets/rbc_all_cells.csv",
           col_names = c("slide_id",features_all)) %>%
    mutate(dataset = "MLL"),
  read_csv("datasets/rbc_adden_2_all_cells.csv",
           col_names = c("slide_id",features_all)) %>%
    mutate(dataset = "AC2")) %>%
  merge(all_conditions,by = "slide_id")
WBC_cells <- rbind(
  read_csv("datasets/wbc_adden_1_all_cells.csv",
           col_names = c("slide_id",features_all,features_nuclear)) %>%
    mutate(dataset = "AC1"),
  read_csv("datasets/wbc_all_cells.csv",
           col_names = c("slide_id",features_all,features_nuclear)) %>%
    mutate(dataset = "MLL"),
  read_csv("datasets/wbc_adden_2_all_cells.csv",
           col_names = c("slide_id",features_all,features_nuclear)) %>%
    mutate(dataset = "AC2")) %>%
  merge(all_conditions,by = "slide_id")

umap_cells <- function(RBC_data,WBC_data,m=0.5,Merge=T) {
  WBC_umap <- WBC_data[,unlist(summarise_all(WBC_data,is.numeric))] %>%
    scale %>%
    umap(min_dist = m)
  WBC_umap_df <- WBC_umap$layout %>%
    as.tibble() %>% 
    mutate(dataset = WBC_data$dataset,
           coarse_class = WBC_data$coarse_class,
           fine_class = WBC_data$fine_class)
  
  RBC_umap <- RBC_data[,unlist(summarise_all(RBC_data,is.numeric))] %>%
    scale %>%
    umap(min_dist = m)
  RBC_umap_df <- RBC_umap$layout %>%
    as.tibble() %>% 
    mutate(dataset = RBC_data$dataset,
           coarse_class = RBC_data$coarse_class,
           fine_class = RBC_data$fine_class)
  if (Merge == T) {
    WBC_RBC_data <- merge(RBC_data,WBC_data,by = c("slide_id","coarse_class","fine_class","dataset"))
    WBC_RBC_umap <- WBC_RBC_data[,unlist(summarise_all(WBC_RBC_data,is.numeric))] %>%
      scale %>% 
      umap(min_dist = m)
    WBC_RBC_umap_df <- WBC_RBC_umap$layout %>%
      as.tibble() %>% 
      mutate(dataset = WBC_RBC_data$dataset,
             coarse_class = WBC_RBC_data$coarse_class,
             fine_class = WBC_RBC_data$fine_class)
  } else {
    WBC_RBC_umap_df <- NULL
  }
  return(list(wbc = WBC_umap_df,rbc = RBC_umap_df,wbc_rbc = WBC_RBC_umap_df))
}

RBC_data <- RBC_data_all %>% 
  subset(moment == "mean")
WBC_data <- WBC_data_all %>%
  subset(moment == "mean")
mean_umap_output <- umap_cells(RBC_data,WBC_data,m=0.9)

RBC_data <- RBC_data_all %>% 
  subset(moment == "variance")
WBC_data <- WBC_data_all %>%
  subset(moment == "variance")
variance_umap_output <- umap_cells(RBC_data,WBC_data,m=0.9)

plot_grid(
  ggplot(mean_umap_output$wbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: WBC mean"),
  ggplot(mean_umap_output$rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: RBC mean"),
  ggplot(mean_umap_output$wbc_rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: WBC mean"),
  ncol = 1) + 
  ggsave(filename = "figures/umap-datasets-mean.pdf",height=4.5,width=4.5)

plot_grid(
  ggplot(variance_umap_output$wbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: WBC variance"),
  ggplot(variance_umap_output$rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: RBC variance"),
  ggplot(variance_umap_output$wbc_rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.5,size = 1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) + 
    ggtitle("UMAP input: RBC and WBC variance"),
  ncol = 1) + 
  ggsave(filename = "figures/umap-datasets-variance.pdf",height=4.5,width=4.5)

RBC_subset <- RBC_cells[sample(nrow(RBC_cells),12500),]
WBC_subset <- WBC_cells[sample(nrow(WBC_cells),12500),]

cells_umap <- umap_cells(
  RBC_subset,WBC_subset,m = 0.99,Merge = F)

wbc_lims <- list(x = quantile(cells_umap$wbc$V1,c(0.01,0.99)),
                 y = quantile(cells_umap$wbc$V2,c(0.01,0.99)))
rbc_lims <- list(x = quantile(cells_umap$rbc$V1,c(0.01,0.99)),
                 y = quantile(cells_umap$rbc$V2,c(0.01,0.99)))

plot_grid(
  ggplot(cells_umap$wbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.4,size = 0.1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) +
    coord_cartesian(xlim = wbc_lims$x,ylim = wbc_lims$y) + 
    scale_x_continuous(expand = c(0.2,0)) + 
    scale_y_continuous(expand = c(0.2,0)) + 
    ggtitle("UMAP input: WBC features"),
  ggplot(cells_umap$rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_point(alpha = 0.4,size = 0.1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) +
    theme(legend.key.size = unit(0,"cm")) +
    coord_cartesian(xlim = rbc_lims$x,ylim = rbc_lims$y) + 
    scale_x_continuous(expand = c(0.2,0)) + 
    scale_y_continuous(expand = c(0.2,0)) + 
    ggtitle("UMAP input: RBC features"),
  ncol = 1) + 
  ggsave(filename = "figures/umap-datasets-cells.pdf",height = 6,width = 6)

plot_grid(
  ggplot(cells_umap$wbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_density2d(size = 0.25) + 
    geom_point(alpha = 0.4,size = 0.1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) +
    coord_cartesian(xlim = wbc_lims$x,ylim = wbc_lims$y) + 
    scale_x_continuous(expand = c(0.3,0)) + 
    scale_y_continuous(expand = c(0.3,0)) + 
    ggtitle("UMAP input: WBC features"),
  ggplot(cells_umap$rbc,aes(x = V1,y = V2,colour = dataset)) + 
    geom_density2d(size = 0.25) + 
    geom_point(alpha = 0.4,size = 0.1) + 
    facet_wrap(~ coarse_class) +
    theme_pretty(base_size = 6) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_jama(name = NULL) + 
    theme(legend.key.size = unit(0,"cm")) +
    theme(legend.key.size = unit(0,"cm")) +
    coord_cartesian(xlim = rbc_lims$x,ylim = rbc_lims$y) + 
    scale_x_continuous(expand = c(0.3,0)) + 
    scale_y_continuous(expand = c(0.3,0)) + 
    ggtitle("UMAP input: RBC features"),
  ncol = 1) + 
  ggsave(filename = "figures/umap-datasets-cells-density.pdf",height = 6,width = 6)