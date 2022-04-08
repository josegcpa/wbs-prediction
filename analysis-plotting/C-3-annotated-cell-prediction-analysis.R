# setup -------------------------------------------------------------------

source("function-library.R")

library(cowplot)

set.seed(42)

threshold <- 0.1

wbc_type_conversion <- c(
  gn  ="Neutrophil",gb  ="Basophil",ge = "Eosinophil",bc = "Band cell/Pelger-Huët",
  ly = "Lymphocyte",mo = "Monocyte",bl = "Blast",
  poor = "Incomplete WBC",multi = "Multiple WBC",null = "No WBC",
  lowres = "Poor resolution",
  rbc = "Red blood cell",imma = "Nucleated RBC",reti = "Reticulocyte",
  plt = "Platelet clump",dead = "Fragmented/unrec.")

rbc_type_conversion <- c(
  no = "Normal",imma = "Nucleated RBC",sp = "Spherocyte",tc = "Target cell",
  ir = "Irregularly contracted", da = "Dacrocyte (tear-shaped)",
  ke = "Keratocyte",ec = "Echinocyte", el = "Eliptocyte",ac = "Acanthocyte",
  poor = "Incomplete RBC",multi = "Multiple RBC",null = "No RBC",
  pl = "Platelet/dye",npl = "Nearby platelet",una = "Fragmented/unrec."
)

model_levels <- c("Disease detection","Disease classification","SF3B1mut detection","Anaemia classification")

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

filter_label <- function(x,seed=42) {
  set.seed(42)
  if (length(x) == 1) {
    return(x)
  } else {
    tbl <- sort(table(x))  
    tbl <- tbl[tbl == max(tbl)]
    return(sample(names(tbl),1))
  }
  }

decode_model_name <- function(model_names) {
  ifelse(
    grepl("anemia_binary",model_names),"Anaemia classification",
    ifelse(grepl("mds_binary",model_names),"SF3B1mut detection",
           ifelse(grepl("disease_binary",model_names),
                  "Disease classification",
                  "Disease detection"))) %>%
    return
}

read_annotations <- function(path) {
  model_id <- str_split(gsub('\\.csv','',path),'/')[[1]]
  model_id <- gsub('[wr]bc_cv_subset\\.','',model_id[length(model_id)])
  nvc <- str_match(model_id,'[0-9]+')
  data_type <- ifelse(grepl("bc",model_id),"Morphology + B.C.","Morphology")
  model_id <- decode_model_name(str_split(model_id,'\\.')[[1]][1])
  read_csv(
    path,
    col_names = c("center_idx","slide_id","x","y","virtual_cell_type","p","label","user")) %>%
    mutate(model_id = model_id,nvc = nvc,data_type = data_type) %>%
    return
}

process_data <- function(x) {
  x %>% group_by(center_idx) %>%
    filter(label == filter_label(label)) %>%
    group_by(virtual_cell_type) %>%
    mutate(total = length(unique(center_idx))) %>%
    group_by(label) %>%
    mutate(total_annotated = length(unique(center_idx))) %>%
    group_by(model_id,data_type) %>%
    mutate(big_total = length(unique(center_idx))) %>% 
    group_by(label,virtual_cell_type,total,model_id,data_type,nvc,big_total,total_annotated) %>%
    summarise(N = length(unique(center_idx))) %>%
    mutate(model_id = factor(model_id,levels = model_levels)) %>%
    return
}

generate_plots <- function(tmp_data) {
  entropy_data <- tmp_data %>%
    mutate(P = N/total_annotated) %>%
    group_by(label,data_type,model_id) %>%
    summarise(Entropy = -sum(P * log(P)))
  LL <- c(0,max(entropy_data$Entropy))
  entropy_plot <- entropy_data %>%
    ggplot(aes(x = Entropy,y = label)) +
    geom_bar(stat = "identity",colour = "black",size = 0.25,fill = "grey90") + 
    facet_wrap(~model_id,ncol=1,scales = "free_x") + 
    theme_pretty(base_size = 6) + 
    theme(axis.text.y = element_blank(),axis.title.y = element_blank()) + 
    scale_x_continuous(expand = c(0,0),limits = LL,breaks = c(0,1,2))

  count_plot <- tmp_data %>%
    ungroup %>% 
    select(label,total = total_annotated,model_id) %>%
    distinct %>%
    ggplot(aes(x = total,y = label)) +
    geom_bar(stat = "identity",colour = "black",size = 0.25,fill = "grey90") + 
    facet_wrap(~model_id,ncol=1,scales = "free_x") + 
    theme_pretty(base_size = 6) + 
    theme(axis.text.y = element_blank(),axis.title.y = element_blank()) + 
    scale_x_continuous(expand = c(0,0,0.05,0.05),trans = 'log10') + 
    xlab("Count") 
    
  heatmap_plot <- tmp_data %>%
    group_by(virtual_cell_type,model_id) %>%
    mutate(expected_proportion = total_annotated/big_total,
           observed_proportion = N/total) %>%
    mutate(Enrichment = observed_proportion /expected_proportion) %>% 
    rowwise() %>%
    mutate(p.val = fisher.test(matrix(
      c(N,total-N,total_annotated,big_total-total_annotated),nrow = 2),
      alternative = "greater")$p.val) %>%
    mutate(sig = ifelse(p.val < 0.05,"*",NA)) %>% 
    mutate(vcf = paste(model_id,virtual_cell_type)) %>%
    group_by(vcf) %>%
    filter(any(sig == "*")) %>%
    mutate(S = -as.numeric(label)[which.max(Enrichment)]) %>%
    ungroup %>%
    mutate(vcf = reorder(vcf,S)) %>%
    group_by(vcf) %>%
    ggplot(aes(y = label,x = vcf,
               fill = Enrichment)) + 
    geom_tile() + 
    geom_text(aes(label = sig),vjust = 0.75) +
    scale_fill_gradient2(name = "Enrichment",midpoint = 0, trans = 'log10',
                         low = "blue4",mid = "grey95",high = "red4") + 
    facet_wrap(~ model_id,scales="free_x",ncol = 1) + 
    theme_pretty(base_size = 6) + 
    xlab("Virtual cell type") + 
    ylab("Expert annotation") + 
    theme(legend.position = "bottom",
          legend.key.height = unit(0.2,"cm")) + 
    guides(fill = guide_legend(nrow = 1)) + 
    scale_size(limits = c(0,1)) + 
    scale_x_discrete(labels = function(x) str_match(x,'[CM0-9]+$'))
  
  return(list(heatmap_plot = heatmap_plot,entropy_plot = entropy_plot,
              count_plot = count_plot))
}


# inter-individual agreement ----------------------------------------------

wbc_labels <- lapply(
  list.files(path = "../mil-comori/labelled-cells/output/",
             pattern = "^wbc.+csv",full.names = T),
  read_annotations
) %>%
  do.call(what = rbind) %>%
  mutate(cell_type = "wbc",virtual_cell_type = virtual_cell_type + 1) %>%
  subset(!(user %in% c(1,2))) %>%
  mutate(label = wbc_type_conversion[label]) %>% 
  mutate(label = factor(label,levels = rev(wbc_type_conversion)))

wbc_labels_mo <- lapply(
  list.files(path = "../mil-comori/labelled-cells/output/",
             pattern = "^mo_wbc.+csv",full.names = T),
  read_annotations
) %>%
  do.call(what = rbind) %>%
  mutate(cell_type = "wbc",virtual_cell_type = virtual_cell_type + 1) %>%
  subset(!(user %in% c(1,2))) %>%
  mutate(label = wbc_type_conversion[label]) %>% 
  mutate(label = factor(label,levels = rev(wbc_type_conversion)))

rbc_labels <- lapply(
  list.files(path = "../mil-comori/labelled-cells/output/",
             pattern = "^rbc.+csv",full.names = T),
  read_annotations
) %>%
  do.call(what = rbind) %>%
  mutate(cell_type = "rbc",virtual_cell_type = virtual_cell_type + 1) %>%
  subset(!(user %in% c(1,2))) %>%
  mutate(label = rbc_type_conversion[label]) %>% 
  mutate(label = factor(label,levels = rev(rbc_type_conversion)))

rbc_labels_mo <- lapply(
  list.files(path = "../mil-comori/labelled-cells/output/",
             pattern = "^mo_rbc.+csv",full.names = T),
  read_annotations
) %>%
  do.call(what = rbind) %>%
  mutate(cell_type = "rbc",virtual_cell_type = virtual_cell_type + 1) %>%
  subset(!(user %in% c(1,2))) %>%
  mutate(label = rbc_type_conversion[label]) %>% 
  mutate(label = factor(label,levels = rev(rbc_type_conversion)))

wbc_labels %>% 
  select(center_idx,user,label) %>% 
  distinct %>%
  group_by(center_idx,user) %>% 
  filter(label == label[length(label)]) %>% 
  group_by(center_idx) %>% 
  filter(length(user)>1) %>%
  summarise(all_same = length(unique(label)) == 1,label = filter_label(label)) %>% 
  group_by(label) %>% 
  summarise(N = sum(all_same),Total = length(unique(center_idx))) %>% 
  mutate(P = N/Total) %>% 
  arrange(P) %>% 
  rowwise() %>%
  mutate(ymin = prop.test(N,Total)$conf.int[1],
         ymax = prop.test(N,Total)$conf.int[2]) %>%
  ggplot(aes(x = P,y = reorder(label,P),
             xmin = ymin,xmax = ymax)) + 
  geom_bar(stat = "identity") + 
  geom_linerange() +
  geom_text(aes(x = 1,label = sprintf("%s/%s",N,Total)),size = 2,hjust = -0.04) + 
  xlab("Proportion of labels with no disagreements") +
  ylab("Expert annotations") + 
  theme_pretty(base_size = 6) +
  theme(axis.text = element_text(colour = "black"),panel.grid.major.y = element_line(size = 0.25,colour = "grey70")) +
  scale_x_continuous(expand = c(0,0,0.2,0),breaks = c(0,0.25,0.5,0.75,1)) +
  ggsave("figures/mil-comori-annotated-cells-wbc-agreement.pdf",height=1.7,width=3)

wbc_labels %>% 
  select(center_idx,user,label) %>% 
  distinct %>%
  group_by(center_idx,user) %>% 
  filter(label == label[length(label)]) %>%
  group_by(center_idx) %>% 
  filter(length(unique(user))>1) %>%
  distinct %>% 
  mutate(major_label = filter_label(label)) %>%
  group_by(major_label) %>%
  mutate(total = length(center_idx)) %>%
  group_by(major_label,label,total) %>%
  summarise(N = length(center_idx)) %>%
  mutate(major_label = factor(major_label,levels = wbc_type_conversion)) %>% 
  ggplot(aes(x = major_label,y = label,fill = N/total)) +
  geom_tile() + 
  theme_pretty(base_size = 6) +
  ggpubr::rotate_x_text(angle = 45) + 
  xlab("Majority label") +
  ylab("Expert annotation") +
  scale_fill_material(name = "Proportion") +
  theme(legend.key.width = unit(0.2,"cm"),axis.text = element_text(colour = "black")) +
  ggsave("figures/mil-comori-annotated-cells-wbc-agreement-heatmap.pdf",height=2.5,width=3.5)

rbc_labels %>% 
  select(center_idx,user,label) %>% 
  distinct %>%
  group_by(center_idx) %>% 
  filter(length(user)>1) %>%
  summarise(all_same = length(unique(label)) == 1,label = filter_label(label)) %>% 
  group_by(label) %>% 
  summarise(N = sum(all_same),Total = length(unique(center_idx))) %>% 
  mutate(P = N/Total) %>% 
  arrange(P) %>% 
  ggplot(aes(x = P,y = reorder(label,P))) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = sprintf("%s/%s",N,Total)),size = 2,hjust = -0.04) + 
  xlab("Proportion of labels with no disagreements") +
  ylab("Expert annotations") + 
  theme_pretty(base_size = 6) +
  theme(axis.text = element_text(colour = "black")) +
  scale_x_continuous(expand = c(0,0,0.2,0),breaks = c(0,0.25,0.5,0.75,1)) +
  ggsave("figures/mil-comori-annotated-cells-rbc-agreement.pdf",height=1.7,width=3)

# proportions per class ---------------------------------------------------
X <- wbc_labels %>%
  merge(all_conditions,by = "slide_id") %>%
  select(fine_class = coarse_class,slide_id,center_idx,label) %>%
  distinct %>%
  group_by(label) %>%
  mutate(total_label = length(unique(center_idx))) %>%
  group_by(fine_class) %>%
  mutate(total_class = length(unique(slide_id))) %>%
  group_by(fine_class,label,total_label,total_class) %>%
  summarise(N = length(unique(center_idx)))

all_poisson_tests <- list()
for (l in unique(X$label)) {
  Y <- X %>%
    subset(label == l)
  if (length(Y$fine_class) > 1) {
    all_c <- t(combn(Y$fine_class,2))
    for (i in 1:nrow(all_c)) {
      CC <- sort(all_c[i,])
      Z <- subset(Y,fine_class %in% CC)
      pt <- poisson.test(Z$N,Z$total_class)
      pt <- data.frame(ratio_estimate = pt$estimate,p.val = pt$p.value,
                       A = CC[1],B = CC[2],label = l)
      all_poisson_tests[[length(all_poisson_tests)+1]] <- pt
    }
  }
}

do.call(rbind,all_poisson_tests) %>% 
  mutate(colour_code = ifelse(
    p.val < 0.05,
    ifelse(ratio_estimate > 1,"Yes (larger in\ncondition 1)","Yes (larger in\ncondition 2)"),
    "No")) %>%
  subset(grepl("Pelg|Blast",label)) %>%
  ggplot(aes(x = A,y = B,fill = colour_code)) +
  geom_tile() + 
  facet_wrap(~ label) + 
  theme_pretty(base_size = 6) + 
  scale_fill_lancet(name = "p-value for\nratio test < 0.05") + 
  xlab("Condition 1") + 
  ylab("Condition 2") + 
  theme(legend.key.size = unit(0.4,"cm")) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  ggsave("figures/annotated-wbc-rate-test.pdf",height = 1,width = 3.5)

wbc_labels %>%
  merge(all_conditions,by = "slide_id") %>%
  select(fine_class = coarse_class,slide_id,center_idx,label) %>%
  distinct %>%
  group_by(label) %>%
  mutate(total_label = length(unique(center_idx))) %>%
  group_by(fine_class) %>%
  mutate(total_class = length(unique(slide_id))) %>%
  group_by(fine_class,label,total_label,total_class) %>%
  summarise(N = length(unique(center_idx))) %>%
  ggplot(aes(y = label,x = N/total_class,fill = fine_class)) + 
  geom_bar(stat = "identity",position = "dodge") + 
  scale_fill_manual(values = fine_colours,name = NULL) + 
  theme_pretty(base_size = 6)  + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.y = element_blank(),legend.key.size = unit(0.2,"cm")) +
  xlab("Instances per condition") +
  ggsave("figures/annotated-wbc-rate.pdf",height = 1.8,width = 3.5) 
  
# wbc (morphology) --------------------------------------------------------

tmp_data <- wbc_labels %>%
  subset(data_type == "Morphology") %>% 
  process_data()

plots <- generate_plots(tmp_data)

plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.25),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-wbc-morphology.pdf",width=4,height=6.5)

# wbc (morphology + bc) ---------------------------------------------------

tmp_data <- wbc_labels %>%
  subset(data_type == "Morphology + B.C.") %>% 
  process_data()

plots <- generate_plots(tmp_data)

plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.25),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-wbc-morphology-bc.pdf",width=4,height=6.5)

tmp_data <- wbc_labels_mo %>%
  subset(data_type == "Morphology + B.C.") %>%
  process_data()

plots <- tmp_data %>%
  mutate(vct_l = vct_conversion(virtual_cell_type,"wbc")) %>%
  mutate(virtual_cell_type = ifelse(
    is.na(vct_l),as.character(virtual_cell_type),paste0("CM",vct_l))) %>%
  select(-vct_l) %>%
  generate_plots()
plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.17),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-mo-wbc-morphology-bc.pdf",
         height=2,width = 6)

# rbc (morphology) --------------------------------------------------------

tmp_data <- rbc_labels %>%
  subset(data_type == "Morphology") %>% 
  process_data()

plots <- generate_plots(tmp_data)

plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.25),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-rbc-morphology.pdf",width=4,height=6.5)


# rbc (morphology + bc) ---------------------------------------------------

tmp_data <- rbc_labels %>%
  subset(data_type == "Morphology + B.C.") %>% 
  process_data()

plots <- generate_plots(tmp_data)

plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.25),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-rbc-morphology-bc.pdf",width=4,height=6.5)

tmp_data <- rbc_labels_mo %>%
  subset(data_type == "Morphology + B.C.") %>%
  process_data()

plots <- tmp_data %>%
  mutate(vct_l = vct_conversion(virtual_cell_type,"rbc")) %>%
  mutate(virtual_cell_type = ifelse(
    is.na(vct_l),as.character(virtual_cell_type),paste0("CM",vct_l))) %>%
  select(-vct_l) %>%
  generate_plots()
plot_grid(plots$heatmap_plot,plots$count_plot,
          align = "h",axis = "lrtb",rel_widths = c(1,0.17),nrow = 1) + 
  ggsave("figures/mil-comori-annotated-cells-mo-rbc-morphology-bc.pdf",
         height=2,width = 6)

