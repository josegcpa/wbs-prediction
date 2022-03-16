# setup -------------------------------------------------------------------

source("function-library.R")

library(ggrepel)
library(umap)
library(cowplot)
library(MASS)

select <- dplyr::select

multi_objective_matching <- c("anemia_binary","binary","mds_binary","disease_binary")

N_points <- 50000
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

conversion_list <- list(
  `Disease detection` = c(
    Normal = "Normal",`SF3B1-wildtype` = "Disease",
    `SF3B1-mutant` = "Disease",`Iron deficiency` = "Disease",
    `Megaloblastic` = "Disease"),
  `Disease classification` = c(
    `SF3B1-wildtype` = "MDS",
    `SF3B1-mutant` = "MDS",`Iron deficiency` = "Anaemia",
    `Megaloblastic` = "Anaemia"),
  `SF3B1mut detection` = c(
    `SF3B1-wildtype` = "SF3B1-wildtype",`SF3B1-mutant` = "SF3B1-mutant"),
  `Anaemia classification` = c(
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"),
  `Multi-objective` = c(
    Normal = "Normal",`SF3B1-wildtype` = "SF3B1-wildtype",`SF3B1-mutant` = "SF3B1-mutant",
    `Iron deficiency` = "Iron deficiency",
    `Megaloblastic` = "Megaloblastic"))

fine_class_order <- c(
  "Normal","Disease","MDS","Anaemia","SF3B1-mutant","SF3B1-wildtype",
  "Iron deficiency","Megaloblastic")

read_final_layer <- function(path) {
  df <- read_csv(path,col_names=F)
  colnames(df) <- c("model_name","fold","class_id","task_idx","virtual_cell_type","coefficient")
  df$max_vc <- as.numeric(str_match(df$model_name,'[0-9]+'))
  
  df <- df %>%
    na.omit %>%
    mutate(model_name = gsub("models/","",model_name),
           virtual_cell_type = virtual_cell_type + 1) %>% 
    group_by(model_name) %>%
    filter(virtual_cell_type <= (max_vc*2)) %>%
    mutate(cell_type = ifelse(virtual_cell_type <= max_vc,"WBC","RBC")) %>% 
    mutate(
      virtual_cell_type = ifelse(
        virtual_cell_type > max_vc,virtual_cell_type - max_vc,virtual_cell_type)) %>% 
    mutate(data_type = ifelse(grepl("bc$",model_name),"Morphology + B.C.","Morphology")) %>%
    ungroup %>% 
    spread(key = "class_id",value = "coefficient") %>% 
    mutate(class_difference = as.numeric(`1`) - as.numeric(`0`))
  
  return(df)} 

read_vcq_layer <- function(path,wbc_subset_path,rbc_subset_path) {
  wbc_feature_subset <- unlist(read_csv(wbc_subset_path,col_names = F))
  rbc_feature_subset <- unlist(read_csv(rbc_subset_path,col_names = F))
  df <- read_csv(path,col_names=F)
  colnames(df) <- c("model_name","fold","feature","cell_type","virtual_cell_type","coefficient")
  df$max_vc <- as.numeric(str_match(df$model_name,'[0-9]+'))
  
  df <- df %>%
    na.omit %>%
    mutate(model_name = gsub("models/","",model_name),
           virtual_cell_type = 1+virtual_cell_type) %>% 
    group_by(model_name) %>%
    mutate(cell_type = ifelse(cell_type == 'dataset_0',"WBC","RBC")) %>% 
    mutate(feature = ifelse(
      cell_type == 'WBC',
      c(features_all,features_nuclear)[wbc_feature_subset][feature+1],
      features_all[rbc_feature_subset][feature+1]
    )) %>%
    mutate(
      virtual_cell_type = ifelse(
        virtual_cell_type > max_vc,virtual_cell_type - max_vc,virtual_cell_type)) %>% 
    mutate(data_type = ifelse(grepl("bc$",model_name),"Morphology + B.C.","Morphology")) %>%
    ungroup 
  
  return(df)}

read_proportions <- function(path) {
  df <- read_csv(path,col_names = F)
  colnames(df)[c(1,2)] <- c("slide_id","cell_type") 
  
  df <- df %>%
    mutate(cell_type = ifelse(cell_type == "cell_type_0","WBC","RBC")) %>%
    gather(key = "virtual_cell_type",value = "proportion",-slide_id,-cell_type) %>%
    mutate(virtual_cell_type = gsub("X","",virtual_cell_type)) %>%
    mutate(virtual_cell_type = as.numeric(virtual_cell_type)) %>%
    mutate(virtual_cell_type = virtual_cell_type - min(virtual_cell_type)+1) %>%
    mutate(max_vc = max(virtual_cell_type))
  
  return(df)}

coefficients_plot <- function(x) {
  x %>%
    ggplot(aes(x = virtual_cell_type_fctr,
               y = feature,
               fill = coefficient)) + 
    geom_tile() + 
    facet_wrap( ~ sprintf('%s (%s)',data_type,cell_type),scales = "free") + 
    theme_pretty(base_size = 6) + 
    theme(legend.position = "bottom") +
    xlab("Virtual cell type") + 
    ylab("Feature") +
    scale_fill_gradient2(low = 'blue4',mid = "white",high = "red4",
                         name = "Coefficient") + 
    theme(legend.key.height = unit(0.1,"cm"),
          strip.text = element_text(margin = margin(t = 1,b = 2)),
          axis.text = element_text(colour = "black")) +
    scale_x_discrete(expand = c(0,0),labels = function(x) str_match(x,'[0-9]+$')) +
    scale_y_discrete(expand = c(0,0)) %>%
    return
}

subset_vcq <- function(best_vcq_subset,S) {
  best_vcq_subset %>%
    subset(virtual_cell_type_fctr %in% S) %>%
    mutate(virtual_cell_type_fctr = factor(
      virtual_cell_type_fctr,levels=S,ordered = T)) %>%
    return
}

# general data-processing -------------------------------------------------

feat_conv <- rev(gsub('\n',' ',features_conversion))

best_vcq_subset <- rbind(
  read_vcq_layer(
    path = "../mile-vice/best_models/best_vcq_layers_subset",
    wbc_subset_path = "data_output/wbc_feature_subset",
    rbc_subset_path = "data_output/rbc_feature_subset") %>%
    mutate(virtual_cell_type_fctr = paste(
      decode_model_name(model_name),data_type,cell_type,virtual_cell_type)) %>%
    mutate(feature = factor(feat_conv[feature],levels = feat_conv)),
  read_vcq_layer(
    path = "../mile-vice/best_models/best_vcq_layers_mo_subset",
    wbc_subset_path = "data_output/wbc_feature_subset",
    rbc_subset_path = "data_output/rbc_feature_subset") %>%
    mutate(virtual_cell_type_fctr = paste(data_type,cell_type,virtual_cell_type)) %>%
    mutate(feature = factor(feat_conv[feature],levels = feat_conv)))

# rbc data-processing -----------------------------------------------------

rbc_cells <- list.files("datasets/many-cells/",pattern = "^rbc-cv_subset\\..*",
                        full.names = T) %>%
  lapply(read_csv,col_names = F) %>%
  do.call(what = rbind)
rbc_cells_subset <- rbc_cells[sample(nrow(rbc_cells),size = N_points,replace = F),]
colnames(rbc_cells_subset)[1:4] <- c("model_name","slide_id",
                                     "cell_type","virtual_cell_type")
colnames(rbc_cells_subset)[5:ncol(rbc_cells_subset)] <- features_all[
  unlist(read.csv("data_output/rbc_feature_subset",header=F))]
rbc_cells_subset$unique_idx <- 1:nrow(rbc_cells_subset)
rbc_cells_subset_long <- rbc_cells_subset %>%
  mutate(cell_type = toupper(cell_type)) %>%
  gather(key = "feature",value = "value",-model_name,-slide_id,
         -cell_type,-virtual_cell_type,-unique_idx) %>%
  mutate(feature = factor(feat_conv[feature],levels = feat_conv)) %>%
  group_by(feature) %>%
  merge(all_conditions,by = "slide_id") %>%
  mutate(fine_class = ifelse(
    coarse_class == "MDS",
    ifelse(grepl('SF3B1',fine_class),"SF3B1-mutant","SF3B1-wildtype"),
    as.character(fine_class))) %>% 
  subset(!grepl("axis",feature)) %>% 
  subset(!grepl("GLCM",feature))

rbc_feature_matrix <- rbc_cells_subset_long %>%
  select(model_name,feature,unique_idx,virtual_cell_type,slide_id,fine_class,
         coarse_class,value) %>%
  group_by(feature) %>%
  filter(value <= quantile(value,0.975) & value >= quantile(value,0.025)) %>%
  mutate(value = (value-mean(value))/sd(value)) %>%
  mutate(value = (value - min(value))/(max(value)-min(value))) %>%
  spread(key = "feature",value = "value") %>%
  subset(grepl("\\.bc",model_name)) %>%
  na.omit()

M <- rbc_feature_matrix[,7:ncol(rbc_feature_matrix)]
rbc_feature_umap <- umap(M,min_dist = 0.9,verbose=T)
rbc_feature_umap_data <- cbind(rbc_feature_matrix[,1:6],rbc_feature_umap$layout)

rbc_feature_matrix_ <- cbind(rbc_feature_matrix[,1:6],M)
rbc_vc_means <- rbc_feature_matrix_ %>%
  select(-unique_idx,-slide_id,-fine_class,-coarse_class) %>% 
  group_by(model_name,virtual_cell_type) %>% 
  summarise_all(.funs = mean)

rbc_vc_centers <- predict(rbc_feature_umap,rbc_vc_means[,3:ncol(rbc_vc_means)]) %>%
  as.data.frame
colnames(rbc_vc_centers) <- c("UMAP1","UMAP2")

rbc_vc_centers_data <- cbind(
  as.data.frame(select(rbc_vc_means,model_name,virtual_cell_type)),rbc_vc_centers)

# wbc data-processing -----------------------------------------------------

wbc_cells <- list.files("datasets/many-cells/",pattern = "^wbc-cv_subset\\..*",
                        full.names = T) %>%
  lapply(read_csv,col_names = F) %>%
  do.call(what = rbind)
wbc_cells_subset <- wbc_cells[sample(nrow(wbc_cells),size = N_points,replace = F),]
colnames(wbc_cells_subset)[1:4] <- c("model_name","slide_id","cell_type",
                                     "virtual_cell_type")
colnames(wbc_cells_subset)[5:ncol(wbc_cells_subset)] <- c(features_all,features_nuclear)[
  unlist(read.csv("data_output/scripts/wbc_feature_subset",header=F))]
wbc_cells_subset$unique_idx <- 1:nrow(wbc_cells_subset)
wbc_cells_subset_long <- wbc_cells_subset %>%
  mutate(cell_type = toupper(cell_type)) %>%
  gather(key = "feature",value = "value",-model_name,-slide_id,-cell_type,
         -virtual_cell_type,-unique_idx) %>%
  mutate(feature = factor(feat_conv[feature],levels = feat_conv)) %>%
  group_by(feature) %>%
  merge(all_conditions,by = "slide_id") %>%
  mutate(fine_class = ifelse(
    coarse_class == "MDS",
    ifelse(grepl('SF3B1',fine_class),"SF3B1-mutant","SF3B1-wildtype"),
    as.character(fine_class)))

wbc_feature_matrix <- wbc_cells_subset_long %>%
  select(model_name,feature,unique_idx,virtual_cell_type,slide_id,fine_class,
         coarse_class,value) %>%
  group_by(feature) %>%
  filter(value <= quantile(value,0.99) & value >= quantile(value,0.01)) %>%
  mutate(value = (value-mean(value))/sd(value)) %>%
  mutate(value = (value - min(value))/(max(value)-min(value))) %>%
  spread(key = "feature",value = "value") %>%
  subset(grepl("\\.bc",model_name)) %>%
  na.omit()

M <- wbc_feature_matrix[,7:ncol(wbc_feature_matrix)]
wbc_feature_umap <- umap(M,min_dist = 0.1,verbose=T,n_neighbors=200)
wbc_feature_umap_data <- cbind(wbc_feature_matrix[,1:6],wbc_feature_umap$layout)

wbc_feature_matrix_ <- cbind(wbc_feature_matrix[,1:6],M)
wbc_vc_means <- wbc_feature_matrix_ %>%
  select(-unique_idx,-slide_id,-fine_class,-coarse_class) %>% 
  group_by(model_name,virtual_cell_type) %>% 
  summarise_all(.funs = mean)

wbc_vc_centers <- predict(wbc_feature_umap,wbc_vc_means[,3:ncol(wbc_vc_means)]) %>%
  as.data.frame
colnames(wbc_vc_centers) <- c("UMAP1","UMAP2")

wbc_vc_centers_data <- cbind(
  as.data.frame(select(wbc_vc_means,model_name,virtual_cell_type)),wbc_vc_centers)

# plotting simple u-map ---------------------------------------------------

bw <- bandwidth.nrd(c(rbc_feature_umap_data$`1`,rbc_feature_umap_data$`2`))
bw_rbc <- bw
L <- c(range(rbc_feature_umap_data$`1`),range(rbc_feature_umap_data$`2`))
L_rbc <- L
long_kde2d <- function(x,y,h,n=25,lims = c(range(x), range(y))) {
  X <- kde2d(x,y,h,n,lims)
  data.frame(x = X$x,X$z) %>%
    gather(key = "y",value = "z",-x) %>%
    mutate(y = X$y[as.numeric(gsub("X","",y))]) %>%
    return
}

plot_grid(
  ggplot(rbc_feature_umap_data,aes(x = `1`,y = `2`)) + 
    geom_point(size = 0.25,colour = "grey",alpha = 0.5) + 
    geom_point(
      mapping = aes(UMAP1,UMAP2,colour = as.factor(virtual_cell_type)),
      data = rbc_vc_centers_data,
      size = 1) +
    theme_pretty(base_size = 6) + 
    theme(legend.key.size = unit(0,"cm"),
          legend.position = "bottom",
          legend.box.spacing = unit(0,"cm")) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_discrete(guide = F) +
    facet_wrap(~ decode_model_name(model_name),nrow = 1) +
    scale_x_continuous(limits = L[1:2]) + 
    scale_y_continuous(breaks = seq(-10,10,by = 2),limits = L[3:4]),
  NULL,
  rbc_feature_umap_data %>%
    rowwise() %>%
    mutate(fine_class = conversion_list[[decode_model_name(model_name)]][fine_class]) %>%
    group_by(model_name,fine_class) %>%
    summarise(long_kde2d(`1`,`2`,bw,100,L)) %>%
    subset(z >= 0.0005) %>%
    group_by(model_name,x,y) %>%
    summarise(P = max(z)/sum(sort(z,decreasing = T)[1:length(z)]),
              cl = fine_class[which.max(z)]) %>%
    mutate(cl = factor(cl,levels = fine_class_order)) %>% 
    na.omit %>% 
    ggplot(aes(x = x,y = y,fill = cl,alpha = P)) +
    geom_tile() + 
    scale_fill_manual(values = fine_colours,name = "Predominant class") + 
    facet_wrap(~ decode_model_name(model_name),nrow = 1) + 
    theme_pretty(base_size = 6) + 
    theme(axis.text = element_blank(),legend.position = "bottom",
          legend.key.size = unit(0.2,"cm")) + 
    scale_alpha(name = "density ratio",range = c(1e-8,1)) +
    xlab("UMAP1") +
    ylab("UMAP2") + 
    scale_x_continuous(limits = L[1:2]) + 
    scale_y_continuous(limits = L[3:4]) + 
    guides(alpha = guide_legend(nrow = 3)),
  align = "hv",axis = "tblr",ncol = 1,rel_heights = c(1,-0.14,1)) + 
  ggsave("figures/u-map-rbc.pdf",width=6,height=3)

bw <- bandwidth.nrd(c(wbc_feature_umap_data$`1`,wbc_feature_umap_data$`2`))
L <- c(range(wbc_feature_umap_data$`1`),range(wbc_feature_umap_data$`2`))
long_kde2d <- function(x,y,h,n=25,lims = c(range(x), range(y))) {
  X <- kde2d(x,y,h,n,lims)
  data.frame(x = X$x,X$z) %>%
    gather(key = "y",value = "z",-x) %>%
    mutate(y = X$y[as.numeric(gsub("X","",y))]) %>%
    return
}

plot_grid(
  ggplot(wbc_feature_umap_data,aes(x = `1`,y = `2`)) + 
    geom_point(size = 0.25,colour = "grey",alpha = 0.5) + 
    geom_point(
      mapping = aes(UMAP1,UMAP2,colour = as.factor(virtual_cell_type)),
      data = wbc_vc_centers_data,alpha = 0.8,size = 1.5) +
    theme_pretty(base_size = 6) + 
    theme(legend.key.size = unit(0,"cm"),
          legend.position = "bottom",
          legend.box.spacing = unit(0,"cm")) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_discrete(guide = F) +
    facet_wrap(~ decode_model_name(model_name),nrow = 1) +
    scale_x_continuous(limits = L[1:2]) + 
    scale_y_continuous(breaks = seq(-9,9,by = 3),limits = L[3:4]),
  ggplot(rbc_feature_umap_data,aes(x = `1`,y = `2`)) + 
    geom_point(size = 0.25,colour = "grey",alpha = 0.5) + 
    geom_point(
      mapping = aes(UMAP1,UMAP2,colour = as.factor(virtual_cell_type)),
      data = rbc_vc_centers_data,alpha = 0.8,size = 1.5) +
    theme_pretty(base_size = 6) + 
    theme(legend.key.size = unit(0,"cm"),
          legend.position = "bottom",
          legend.box.spacing = unit(0,"cm")) + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    scale_colour_discrete(guide = F) +
    facet_wrap(~ decode_model_name(model_name),nrow = 1) +
    scale_x_continuous(limits = L_rbc[1:2]) + 
    scale_y_continuous(breaks = seq(-9,9,by = 3),limits = L_rbc[3:4]),
  align = "hv",axis = "tblr",ncol = 1,rel_heights = c(1,1)) + 
  ggsave("figures/u-map-with-vc.pdf",width=5,height=2.5)

plot_grid(
  wbc_feature_umap_data %>%
    rowwise() %>%
    mutate(fine_class = conversion_list[[decode_model_name(model_name)]][fine_class]) %>%
    group_by(model_name,fine_class) %>%
    summarise(long_kde2d(`1`,`2`,bw,100,L)) %>%
    subset(z >= 0.01) %>%
    group_by(model_name,x,y) %>%
    summarise(P = max(z)/sum(sort(z,decreasing = T)[1:length(z)]),
              cl = fine_class[which.max(z)]) %>%
    mutate(cl = factor(cl,levels = fine_class_order)) %>% 
    ggplot(aes(x = x,y = y,fill = cl,
               alpha = P)) +
    geom_tile() + 
    scale_fill_manual(values = fine_colours,guide=F) + 
    facet_wrap(~ decode_model_name(model_name),nrow = 1) + 
    theme_pretty(base_size = 6) + 
    theme(legend.position = "bottom",legend.key.size = unit(0.2,"cm"),
          legend.box = "vertical",legend.margin = margin()) + 
    scale_alpha(name = "density ratio",trans = 'log10',range = c(1e-8,1)) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_x_continuous(limits = L[1:2]) + 
    scale_y_continuous(limits = L[3:4]),
  rbc_feature_umap_data %>%
    rowwise() %>%
    mutate(fine_class = conversion_list[[decode_model_name(model_name)]][fine_class]) %>%
    group_by(model_name,fine_class) %>%
    summarise(long_kde2d(`1`,`2`,bw_rbc,100,L_rbc)) %>%
    subset(z >= 0.0005) %>%
    group_by(model_name,x,y) %>%
    summarise(P = max(z)/sum(sort(z,decreasing = T)[1:length(z)]),
              cl = fine_class[which.max(z)]) %>%
    mutate(cl = factor(cl,levels = fine_class_order)) %>% 
    na.omit %>% 
    ggplot(aes(x = x,y = y,fill = cl,
               alpha = P)) +
    geom_tile() + 
    scale_fill_manual(values = fine_colours,name = "Predominant class") + 
    facet_wrap(~ decode_model_name(model_name),nrow = 1) + 
    theme_pretty(base_size = 6) + 
    theme(legend.position = "bottom",legend.key.size = unit(0.2,"cm"),
          legend.box = "vertical",legend.margin = margin()) + 
    scale_alpha(name = "density ratio",trans = 'log10',range = c(1e-8,1)) +
    xlab("UMAP1") +
    ylab("UMAP2"),
  align = "hv",axis = "tblr",ncol = 1,rel_heights = c(1,1)) + 
  ggsave("figures/u-map-density-ratio.pdf",width=5,height=3.5)

