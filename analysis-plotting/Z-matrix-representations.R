library(tidyverse)
library(gtools)

n_classes <- 2
n_features_rbc <- 10
n_features_wbc <- 14
nvc <- 10
n_slides <- 7
feature_centers_rbc <- lapply(
  1:n_classes,
  function(x) replicate(nvc,rnorm(n_features_rbc)))
feature_centers_wbc <- lapply(
  1:n_classes,
  function(x) replicate(nvc,rnorm(n_features_wbc)))
rvc_prob <- rdirichlet(n_classes,rep(1,nvc)) + rdirichlet(n_classes,rep(1,nvc))
wvc_prob <- rdirichlet(n_classes,rep(1,nvc)) + rdirichlet(n_classes,rep(1,nvc))

all_n_rbc <- round(runif(n_slides,10,20))
all_n_wbc <- round(runif(n_slides,8,all_n_rbc))

all_feature_maps <- list()
all_slide_probs <- list()
all_slide_class <- c()
for (i in 1:n_slides) {
  slide_class <- sample.int(n_classes,size = 1)
  all_slide_class <- c(all_slide_class,slide_class)
  n_rbc <- all_n_rbc[i]
  n_wbc <- all_n_wbc[i]
  rvc_slide_probs <- rdirichlet(n_rbc,rvc_prob[slide_class,])
  wvc_slide_probs <- rdirichlet(n_wbc,wvc_prob[slide_class,])
  rvc <- apply(rvc_slide_probs,1,which.max) 
  wvc <- apply(wvc_slide_probs,1,which.max)
  all_slide_probs[[length(all_slide_probs) + 1]] <- rvc_slide_probs %>%
    as.tibble() %>%
    mutate(slide_id = i,slide_class = slide_class,
           cell_type = "rbc",vc = rvc,cell_id = 1:n_rbc)
  all_slide_probs[[length(all_slide_probs) + 1]] <- wvc_slide_probs %>%
    as.tibble() %>%
    mutate(slide_id = i,slide_class = slide_class,
           cell_type = "wbc",vc = wvc,cell_id = 1:n_wbc)
  all_feature_maps[[length(all_feature_maps) + 1]] <- sapply(rvc,function(x) {
    rnorm(length(feature_centers_rbc[[slide_class]][,x]),
          feature_centers_rbc[[slide_class]][,x])}) %>%
    t %>% 
    as_tibble() %>%
    mutate(cell_id = 1:n_rbc) %>% gather(key = "feature",value = "value",-cell_id) %>% 
    mutate(slide_id = i,slide_class = slide_class,cell_type = "rbc")
  all_feature_maps[[length(all_feature_maps) + 1]] <- sapply(wvc,function(x) {
    rnorm(length(feature_centers_wbc[[slide_class]][,x]),
          feature_centers_wbc[[slide_class]][,x])}) %>%
    t %>% 
    as_tibble() %>%
    mutate(cell_id = 1:n_wbc) %>% gather(key = "feature",value = "value",-cell_id) %>% 
    mutate(slide_id = i,slide_class = slide_class,cell_type = "wbc")
}

all_feature_maps <- do.call(rbind,all_feature_maps)
all_slide_probs <- do.call(rbind,all_slide_probs) %>%
  gather(key = "virtual_cell_type",value = "value",-slide_id,
         -slide_class,-cell_type,-vc,-cell_id)

all_feature_maps %>%
  ggplot(aes(x = as.factor(cell_id),y = feature,alpha = value,fill = cell_type)) + 
  geom_tile() + 
  theme_minimal() +
  theme(axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = "none",
        panel.grid = element_blank(),strip.text = element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c(rbc = "red4",wbc = "purple4")) +
  scale_alpha(range = c(0.2,1)) +
  facet_grid(cell_type ~ slide_id,scales = "free",space = "free") + 
  ggsave("figures/representation-feature-maps.pdf",height=2,width = 5)

all_slide_probs %>%
  ggplot(aes(x = as.factor(cell_id),y = virtual_cell_type,
             alpha = value,fill = cell_type)) + 
  geom_tile() + 
  theme_minimal() +
  theme(axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = "none",
        panel.grid = element_blank(),strip.text = element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c(rbc = "red4",wbc = "purple4")) +
  scale_alpha(range = c(0.2,1)) +
  facet_grid(cell_type ~ slide_id,scales = "free",space = "free") + 
  ggsave("figures/representation-vc-prob.pdf",
         height=2*((2*nvc)/(n_features_rbc + n_features_wbc)),width = 5)

all_slide_probs %>%
  group_by(slide_id,cell_type,virtual_cell_type) %>%
  summarise(proportion = mean(value)) %>% 
  ggplot(aes(x = 1,y = virtual_cell_type,alpha = proportion,fill = cell_type)) +
  geom_tile() + 
  theme_minimal() +
  theme(axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),legend.position = "none",
        panel.grid = element_blank(),strip.text = element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c(rbc = "red4",wbc = "purple4")) +
  scale_alpha(range = c(0.2,1)) +
  facet_grid(cell_type ~ slide_id,scales = "free",space = "free") +
  ggsave("figures/representation-vc-prop.pdf",
         height=2*((2*nvc)/(n_features_rbc + n_features_wbc)),width = 2)

all_slide_class
