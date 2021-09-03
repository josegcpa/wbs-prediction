# setup -------------------------------------------------------------------

source("function-library.R")

library(tidyverse)
library(cowplot)
library(WRS2)
library(umap)
library(ggpubr)
library(caret)
library(glmnet)
library(pROC)
library(RRF)

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

blood_parameters <- read_csv(
  "data/blood_count_data.csv",
  col_types = cols_only(
    id = col_character(),wbc_ul = col_double(),
    hb_g_dl = col_double(),plt_ul = col_double())) %>%
  select(slide_id = id,wbc_ul,hb_g_dl,plt_ul)

wbc_all_cells_summaries <- read_csv(
  "../mile-vice/data_output/wbc_summaries.csv",
  col_names = c("slide_id",features_all,features_nuclear,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all,features_nuclear)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>% 
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

rbc_all_cells_summaries <- read_csv(
  "../mile-vice/data_output/rbc_summaries.csv",
  col_names = c("slide_id",features_all,"f"),
  col_types = c(list(col_character()),
                replicate(length(c(features_all)),col_double()),
                list(col_character()))) %>%
  merge(all_conditions,by = "slide_id") %>%
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

glmnet_models <- list()
rf_models <- list()

IDs <- c("disease_detection","disease_classification","sf3b1","anaemia_classification")

K <- 5

label_conversion <- list(
  disease_detection = c(
    `Normal` = 0,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 1,
    `Iron deficiency` = 1,`Megaloblastic` = 1),
  disease_classification = c(
    `Normal` = NA,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 1,
    `Iron deficiency` = 0,`Megaloblastic` = 0),
  sf3b1 = c(
    `Normal` = NA,
    `SF3B1-mutant` = 1,`Non-SF3B1-mutant` = 0,
    `Iron deficiency` = NA,`Megaloblastic` = NA),
  anaemia_classification = c(
    `Normal` = NA,
    `SF3B1-mutant` = NA,`Non-SF3B1-mutant` = NA,
    `Iron deficiency` = 0,`Megaloblastic` = 1)
)

reverse_label_conversion <- list(
  disease_detection = c("Normal","Disease"),
  disease_classification = c("Anaemia","MDS"),
  sf3b1 = c("nonSF3B1","SF3B1"),
  anaemia_classification = c("Iron deficiency","Megaloblastic")
)

make_folds <- function(f,k=5) {
  f_idx <- seq(1,length(f))
  output_folds <- lapply(1:k,function(x) list(train = c(),test = c()))
  f_unique <- unique(f)
  for (f_element in f_unique) {
    sub_f_idx <- f_idx[f == f_element]
    fold_assignment <- sample.int(k,length(sub_f_idx),replace = T)
    for (fold in sort(unique(fold_assignment))) {
      test_idx <- sub_f_idx[fold_assignment == fold]
      train_idx <- sub_f_idx[!(sub_f_idx %in% test_idx)]
      output_folds[[fold]]$train <- c(output_folds[[fold]]$train,train_idx)
      output_folds[[fold]]$test <- c(output_folds[[fold]]$test,test_idx)
    }
  }
  return(output_folds)
}

make_class_weights <- function(f) {
  weights_f <- 1 - table(f) / length(f)
  for (n in names(weights_f)) {
    f[f == n] <- weights_f[n]  
  }
  return(f)
}

make_class_weights_rf <- function(f) {
  weights_f <- 1 - table(f) / length(f)
  return(weights_f)
}


# data splitting ----------------------------------------------------------

wbc_feat_cols <- 2:61
rbc_feat_cols <- 2:43

wbc_means <- wbc_all_cells_summaries %>%
  subset(f == "mean") %>%
  select(-f)
colnames(wbc_means)[wbc_feat_cols] <- paste(
  colnames(wbc_means)[wbc_feat_cols],"means","WBC",sep = ".")
wbc_vars <- wbc_all_cells_summaries %>%
  subset(f == "variance") %>%
  select(-f)
colnames(wbc_vars)[wbc_feat_cols] <- paste(
  colnames(wbc_vars)[wbc_feat_cols],"vars","WBC",sep = ".")

rbc_means <- rbc_all_cells_summaries %>%
  subset(f == "mean") %>%
  select(-f)
colnames(rbc_means)[rbc_feat_cols] <- paste(
  colnames(rbc_means)[rbc_feat_cols],"means","RBC",sep = ".")
rbc_vars <- rbc_all_cells_summaries %>%
  subset(f == "variance") %>%
  select(-f)
colnames(rbc_vars)[rbc_feat_cols] <- paste(
  colnames(rbc_vars)[rbc_feat_cols],"vars","RBC",sep = ".")


full_dataset_morphology <- merge(
  merge(
    wbc_means,
    wbc_vars,
    by = c("slide_id","coarse_class","fine_class")),
  merge(
    rbc_means,
    rbc_vars,
    by = c("slide_id","coarse_class","fine_class")),
  by = c("slide_id","coarse_class","fine_class"))

full_dataset_morphology <- list(
  data = select(full_dataset_morphology,-slide_id,-coarse_class,-fine_class),
  coarse_labels = full_dataset_morphology$coarse_class,
  fine_labels = full_dataset_morphology$fine_class,
  slide_id = full_dataset_morphology$slide_id
)

stratified_folds <- make_folds(full_dataset_morphology$fine_labels,k = K)



# model training (glmnet; morphology) -------------------------------------


for (ID in IDs) {
  glmnet_models[[ID]] <- list()
  
  for (i in 1:K) {
    print(sprintf("Training fold = %s (glmnet; %s)",i,ID))
    # define folds and training/validation sets
    glmnet_models[[ID]][[i]] <- list()
    fold <- stratified_folds[[i]]
    training_X <- full_dataset_morphology$data[fold$train,] %>%
      as.matrix
    training_y <- label_conversion[[ID]][
      full_dataset_morphology$fine_labels[fold$train]]
    NA_train <- is.na(training_y)
    training_X <- training_X[!NA_train,]
    training_y <- training_y[!NA_train]
    
    testing_X <- full_dataset_morphology$data[fold$test,] %>%
      as.matrix
    testing_y <- label_conversion[[ID]][
      full_dataset_morphology$fine_labels[fold$test]]
    NA_test <- is.na(testing_y)
    testing_X <- testing_X[!NA_test,]
    testing_y <- testing_y[!NA_test]
    
    # calculate class weights
    class_weights <- make_class_weights(training_y)
    
    # normalise data
    #scale_params <- list(means = colMeans(training_X),
    #                     stds = apply(training_X,2,sd))
    #training_X <- (training_X - scale_params$means)
    #testing_X <- (testing_X - scale_params$means)
    
    # training and inference for training and validation sets
    glmnet_models[[ID]][[i]]$model <- cv.glmnet(
      training_X,training_y,family = "binomial",measure = "deviance",
      weights = class_weights)
    
    train_proba_pred <- predict(glmnet_models[[ID]][[i]]$model,
                                newx = training_X,type = "response")
    train_class_pred <- predict(glmnet_models[[ID]][[i]]$model,
                                newx = training_X,type = "class")
    
    test_proba_pred <- predict(glmnet_models[[ID]][[i]]$model,
                               newx = testing_X,type = "response")
    test_class_pred <- predict(glmnet_models[[ID]][[i]]$model,
                               newx = testing_X,type = "class")
    
    glmnet_models[[ID]][[i]]$metrics <- list(
      Training = list(
        ConfusionMatrix = confusionMatrix(table(
          reverse_label_conversion[[ID]][as.numeric(train_class_pred)+1],
          reverse_label_conversion[[ID]][training_y+1])),
        AUC = roc(reverse_label_conversion[[ID]][training_y+1],
                  train_proba_pred,quiet = T,
                  levels = reverse_label_conversion[[ID]])),
      Validation = list(
        ConfusionMatrix = confusionMatrix(table(
          reverse_label_conversion[[ID]][as.numeric(test_class_pred)+1],
          reverse_label_conversion[[ID]][testing_y+1])),
        AUC = roc(reverse_label_conversion[[ID]][testing_y+1],
                  test_proba_pred,quiet = T,
                  levels = reverse_label_conversion[[ID]]))
    )
  }
}


# model training (rrf; morphology) ----------------------------------------



for (ID in IDs) {
  rf_models[[ID]] <- list()
  
  for (i in 1:K) {
    print(sprintf("Training fold = %s (rrf; %s)",i,ID))
    # define folds and training/validation sets
    rf_models[[ID]][[i]] <- list()
    fold <- stratified_folds[[i]]
    training_X <- full_dataset_morphology$data[fold$train,] %>%
      as.matrix
    training_y <- label_conversion[[ID]][
      full_dataset_morphology$fine_labels[fold$train]]
    NA_train <- is.na(training_y)
    training_X <- training_X[!NA_train,]
    training_y <- training_y[!NA_train]
    
    testing_X <- full_dataset_morphology$data[fold$test,] %>%
      as.matrix
    testing_y <- label_conversion[[ID]][
      full_dataset_morphology$fine_labels[fold$test]]
    NA_test <- is.na(testing_y)
    testing_X <- testing_X[!NA_test,]
    testing_y <- testing_y[!NA_test]
    
    # calculate class weights
    class_weights <- make_class_weights_rf(training_y)
    
    # normalise data
    # scale_params <- list(means = colMeans(training_X),
    #                      stds = apply(training_X,2,sd))
    # training_X <- (training_X - scale_params$means)/scale_params$stds
    # testing_X <- (testing_X - scale_params$means)/scale_params$stds
    
    # training and inference for training and validation sets
    rf_models[[ID]][[i]]$model <- RRF(
      training_X,as.factor(training_y),importance = T,
      classwt = class_weights,ntree = 500)
    
    train_proba_pred <- predict(rf_models[[ID]][[i]]$model,
                                newx = training_X,type = "prob")[,1]
    train_class_pred <- predict(rf_models[[ID]][[i]]$model,
                                newx = training_X,type = "response") %>%
      as.character
    test_proba_pred <- predict(rf_models[[ID]][[i]]$model,
                               newdata = testing_X,type = "prob")[,1]
    test_class_pred <- predict(rf_models[[ID]][[i]]$model,
                               newdata = testing_X,type = "response") %>%
      as.character
    
    rf_models[[ID]][[i]]$metrics <- list(
      Training = list(
        ConfusionMatrix = confusionMatrix(table(
          reverse_label_conversion[[ID]][as.numeric(train_class_pred)+1],
          reverse_label_conversion[[ID]][training_y+1])),
        AUC = roc(reverse_label_conversion[[ID]][training_y+1],
                  train_proba_pred,quiet = T,
                  levels = reverse_label_conversion[[ID]])),
      Validation = list(
        ConfusionMatrix = confusionMatrix(table(
          reverse_label_conversion[[ID]][as.numeric(test_class_pred)+1],
          reverse_label_conversion[[ID]][testing_y+1])),
        AUC = roc(reverse_label_conversion[[ID]][testing_y+1],
                  test_proba_pred,quiet = T,
                  levels = reverse_label_conversion[[ID]]))
    )
  }
}

# comparing performance ---------------------------------------------------

all_metrics <- list()

for (ID in IDs) {
  tmp <- list()
  for (s in c("Training","Validation")) {
    tmp[[s]] <- tibble(
      AUC = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$AUC$auc) %>% 
        as.numeric,
      Accuracy = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$overall["Accuracy"]) %>% 
        as.numeric,
      Kappa = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$overall["Kappa"]) %>% 
        as.numeric,
      Sensitivity = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Sensitivity"]) %>% 
        as.numeric,
      Specificity = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Specificity"]) %>% 
        as.numeric,
      Precision = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Precision"]) %>% 
        as.numeric,
      Recall = glmnet_models[[ID]] %>% 
        lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Recall"]) %>% 
        as.numeric,
      Set = s,
      Fold = c(1:K))
  }
  all_metrics[[ID]] <- do.call(rbind,tmp)
  all_metrics[[ID]]$task <- ID
}

all_metrics_df <- do.call(rbind,all_metrics) 

all_metrics_df_long <- all_metrics_df %>%
  gather(key = "key",value = "value",-Set,-task,-Fold)

all_metrics_df_long %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  group_by(Set,task,key) %>%
  summarise(m = mean(value),se = sd(value),.groups = "drop") %>% 
  ggplot(aes(x = task,y = m,ymin = m-se,ymax = m+se,shape = Set,fill = Set,group = Set)) +
  geom_bar(position = position_dodge(width = 0.9),stat = "identity") + 
  geom_linerange(position = position_dodge(width = 0.9)) +
  facet_wrap(~ key,scales = "free_y",ncol = 4) + 
  theme_pretty(base_size = 6) +
  rotate_x_text(angle = 40) + 
  theme(axis.title = element_blank()) +
  scale_shape(name = NULL) +
  scale_fill_manual(values = c("green4","red4"),name = NULL) + 
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(ylim = c(0.4,1))
