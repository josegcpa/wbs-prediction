# TODO: look at feature specific importances, get some examples, see how it goes

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
library(MLmetrics)
library(RANN)

select <- dplyr::select
set.seed(42)

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

rearrange_data <- function(dataset,id="WBC") {
  means <- dataset %>%
    subset(f == "mean") %>%
    select(-f) %>%
    mutate(key = paste(key,"means",id,sep='.')) %>%
    spread(key = key,value = value)
  vars <- dataset %>%
    subset(f == "variance") %>%
    select(-f) %>%
    mutate(key = paste(key,"vars",id,sep='.')) %>%
    spread(key = key,value = value)
  output <- merge(means,vars,by = c(
    "slide_id","coarse_class","fine_class"))
  return(output)
}

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

demographic_data <- read_csv(
  "data/demographic_data.csv",
  col_types = cols_only(
    id = col_character(),sex = col_character(),age = col_double())) %>%
  select(slide_id = id,sex,age)

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

glmnet_models <- list(
  morphology = list(),full = list(),bcdem = list())
rf_models <- list(
  morphology = list(),full = list(),bcdem = list())

# data splitting ----------------------------------------------------------

wbc_feat_cols <- c(1:19,27:50,58:64) #1:ncol(wbc_all_cells_summaries)
rbc_feat_cols <- 1:ncol(rbc_all_cells_summaries)

wbc_all_cells_summaries <- wbc_all_cells_summaries[,wbc_feat_cols] %>%
  gather(key = "key",value = "value",-slide_id,-f,-fine_class,-coarse_class)

rbc_all_cells_summaries <- rbc_all_cells_summaries[,rbc_feat_cols] %>%
  gather(key = "key",value = "value",-slide_id,-f,-fine_class,-coarse_class)

wbc_dataset_morphology <- rearrange_data(
  wbc_all_cells_summaries,id = "WBC")
rbc_dataset_morphology <- rearrange_data(
  rbc_all_cells_summaries,id = "RBC")

full_dataset_morphology <- merge(
  wbc_dataset_morphology,
  rbc_dataset_morphology,
  by = c("slide_id","coarse_class","fine_class"))

full_dataset <- merge(
  full_dataset_morphology,blood_parameters,by = "slide_id") %>%
  merge(demographic_data,by = "slide_id") %>%
  mutate(sex = ifelse(sex == "male",0,1)) %>% 
  mutate(sex.BCDEM = sex,
         age.BCDEM = age,
         wbc_ul.BCDEM = wbc_ul,
         hb_g_dl.BCDEM = hb_g_dl,
         plt_ul.BCDEM = plt_ul) %>%
  select(-sex,-age,-wbc_ul,-hb_g_dl,-plt_ul)

full_dataset_collection <- list(
  data = list(
    morphology = select(full_dataset,-slide_id,-coarse_class,-fine_class,
                        -wbc_ul.BCDEM,-hb_g_dl.BCDEM,-plt_ul.BCDEM,
                        -sex.BCDEM,-age.BCDEM),
    full = select(full_dataset,-slide_id,-coarse_class,-fine_class,-sex.BCDEM,-age.BCDEM),
    bcdem = select(full_dataset,wbc_ul.BCDEM,hb_g_dl.BCDEM,plt_ul.BCDEM)
  ),
  coarse_labels = full_dataset$coarse_class,
  fine_labels = full_dataset$fine_class,
  slide_id = full_dataset$slide_id) 

if (file.exists("data_output/folds_preproc.RDS")) {
  tmp <- readRDS("data_output/folds_preproc.RDS")
  stratified_folds <- tmp$stratified_folds
  scales <- tmp$scales
} else {
  stratified_folds <- make_folds(full_dataset_collection$fine_labels,k = K)
  
  scales <- list(
    full = lapply(
      stratified_folds,
      function(x) {
        pp <- preProcess(full_dataset_collection$data$full[x$train,],
                         method = c("center","scale","corr","nzv","medianImpute"))
        return(pp)}),
    morphology = lapply(
      stratified_folds,
      function(x) {
        pp <- preProcess(full_dataset_collection$data$morphology[x$train,],
                         method = c("center","scale","corr","nzv","medianImpute"))
        return(pp)}),
    bcdem =   lapply(
      stratified_folds,
      function(x) {
        pp <- preProcess(full_dataset_collection$data$bcdem[x$train,],
                         method = c("center","scale","corr","nzv","medianImpute"))
        return(pp)})
  )
  saveRDS(list(stratified_folds = stratified_folds,scales = scales),
          file = "data_output/folds_preproc.RDS")
}

# model training (glmnet) -------------------------------------------------

if (file.exists("data_output/glmnet_models.rds")) {
  glmnet_models <- readRDS("data_output/glmnet_models.rds")
} else {
  for (data_type in c("bcdem","morphology","full")) {
    for (ID in IDs) {
      glmnet_models[[data_type]][[ID]] <- list()
      
      for (i in 1:K) {
        print(sprintf("Training fold = %s (glmnet; %s with %s)",i,ID,data_type))
        
        # define folds and training/validation sets
        glmnet_models[[data_type]][[ID]][[i]] <- list()
        fold <- stratified_folds[[i]]
        training_X <- full_dataset_collection$data[[data_type]][fold$train,] %>%
          as.matrix
        training_y <- label_conversion[[ID]][
          full_dataset_collection$fine_labels[fold$train]]
        NA_train <- is.na(training_y)
        training_X <- training_X[!NA_train,]
        training_y <- training_y[!NA_train]
        
        testing_X <- full_dataset_collection$data[[data_type]][fold$test,] %>%
          as.matrix
        testing_y <- label_conversion[[ID]][
          full_dataset_collection$fine_labels[fold$test]]
        NA_test <- is.na(testing_y)
        testing_X <- testing_X[!NA_test,]
        testing_y <- testing_y[!NA_test]
        
        # calculate class weights
        class_weights <- make_class_weights(training_y)
        
        # scaling
        training_X <- predict(scales[[data_type]][[i]],training_X)
        testing_X <- predict(scales[[data_type]][[i]],testing_X)
        
        # training and inference for training and validation sets
        alphas <- c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)
        models <- list()
        for (alpha in alphas) {
          m <- cv.glmnet(
            training_X,training_y,family = "binomial",type.measure = "auc",
            nfolds = 5,weights = class_weights,alpha = alpha)
          l <- m$lambda
          l.min <- m$lambda.min
          models[[as.character(alpha)]]$model <- m
          models[[as.character(alpha)]]$metric <- m$cvm[l == l.min]
        }
        best_model_idx <- lapply(models,function(x) x$metric) %>%
          unlist
        if (m$name == "Binomial Deviance") {
          best_model_idx <- which.min(best_model_idx)
        } else {
          best_model_idx <- which.max(best_model_idx)
        }
        trained_model <- models[[best_model_idx]]$model
        
        train_proba_pred <- predict(trained_model,
                                    newx = training_X,type = "response")
        train_class_pred <- predict(trained_model,
                                    newx = training_X,type = "class")
        
        test_proba_pred <- predict(trained_model,
                                   newx = testing_X,type = "response")
        test_class_pred <- predict(trained_model,
                                   newx = testing_X,type = "class")
        
        glmnet_models[[data_type]][[ID]][[i]]$model <- trained_model
        glmnet_models[[data_type]][[ID]][[i]]$metrics <- list(
          Training = list(
            ConfusionMatrix = confusionMatrix(table(
              reverse_label_conversion[[ID]][as.numeric(train_class_pred)+1],
              reverse_label_conversion[[ID]][training_y+1])),
            AUC = roc(reverse_label_conversion[[ID]][training_y+1],
                      train_proba_pred[,1],quiet = T,
                      levels = reverse_label_conversion[[ID]])),
          Validation = list(
            ConfusionMatrix = confusionMatrix(table(
              reverse_label_conversion[[ID]][as.numeric(test_class_pred)+1],
              reverse_label_conversion[[ID]][testing_y+1])),
            AUC = roc(reverse_label_conversion[[ID]][testing_y+1],
                      test_proba_pred[,1],quiet = T,
                      levels = reverse_label_conversion[[ID]])))
      }
    }
  }
  saveRDS(glmnet_models,"data_output/glmnet_models.rds")
}

# model training (rrf) ----------------------------------------------------

if (file.exists("data_output/rf_models.rds")) {
  rf_models <- readRDS("data_output/rf_models.rds")
} else {
  for (data_type in c("bcdem","morphology","full")) {
    for (ID in IDs) {
      rf_models[[data_type]][[ID]] <- list()
      
      for (i in 1:K) {
        print(sprintf("Training fold = %s (rrf; %s with %s)",i,ID,data_type))
        
        # define folds and training/validation sets
        rf_models[[data_type]][[ID]][[i]] <- list()
        fold <- stratified_folds[[i]]
        training_X <- full_dataset_collection$data[[data_type]][fold$train,] %>%
          as.matrix
        training_y <- label_conversion[[ID]][
          full_dataset_collection$fine_labels[fold$train]]
        NA_train <- is.na(training_y)
        training_X <- training_X[!NA_train,]
        training_y <- training_y[!NA_train]
        
        testing_X <- full_dataset_collection$data[[data_type]][fold$test,] %>%
          as.matrix
        testing_y <- label_conversion[[ID]][
          full_dataset_collection$fine_labels[fold$test]]
        NA_test <- is.na(testing_y)
        testing_X <- testing_X[!NA_test,]
        testing_y <- gsub(" ","_",reverse_label_conversion[[ID]][testing_y[!NA_test]+1])
        training_y <- gsub(" ","_",reverse_label_conversion[[ID]][training_y+1])
        
        # calculate class weights
        class_weights <- make_class_weights_rf(training_y)
        
        # scaling
        training_X <- predict(scales[[data_type]][[i]],training_X)
        testing_X <- predict(scales[[data_type]][[i]],testing_X)
        
        # training and inference for training and validation sets
        # training and inference for training and validation sets
        # using cv here to estimate the ideal parameters:
        # mtry - number of candidate variables per split
        # coefReg - regularization coefficients (strength of regul.)
        train_control <- trainControl(method="boot", number=1,
                                      search = "grid",
                                      summaryFunction = twoClassSummary,
                                      classProbs = T)
        trained_model <- train(
          x = training_X, y = training_y,
          classwt = class_weights,method = "RRFglobal",
          metric = "ROC",trControl = train_control)
        
        train_proba_pred <- predict(trained_model,
                                    newdata = training_X,type = "prob")[,1]
        train_class_pred <- predict(trained_model,
                                    newdata = training_X,type = "raw") %>%
          as.character
        
        test_proba_pred <- predict(trained_model,
                                   newdata = testing_X,type = "prob")[,1]
        test_class_pred <- predict(trained_model,
                                   newdata = testing_X,type = "raw") %>%
          as.character
        
        rf_models[[data_type]][[ID]][[i]]$model <- trained_model
        rf_models[[data_type]][[ID]][[i]]$metrics <- list(
          Training = list(
            ConfusionMatrix = confusionMatrix(table(train_class_pred,training_y)),
            AUC = roc(training_y,train_proba_pred,quiet = T)),
          Validation = list(
            ConfusionMatrix = confusionMatrix(table(test_class_pred,testing_y)),
            AUC = roc(testing_y,test_proba_pred,quiet = T)))
      }
    }
  }
  saveRDS(rf_models,"data_output/rf_models.rds")
}

# comparing performance ---------------------------------------------------

all_metrics <- list()

for (data_type in c("bcdem","morphology","full")) {
  for (ID in IDs) {
    tmp <- list()
    for (s in c("Training","Validation")) {
      tmp[[s]] <- rbind(
        tibble(
          AUC = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$AUC$auc) %>% 
            as.numeric,
          Accuracy = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$overall["Accuracy"]) %>% 
            as.numeric,
          Sensitivity = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Sensitivity"]) %>% 
            as.numeric,
          Specificity = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Specificity"]) %>% 
            as.numeric,
          Precision = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Precision"]) %>% 
            as.numeric,
          Recall = glmnet_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Recall"]) %>% 
            as.numeric,
          Set = s,Fold = c(1:K),Model = "glmnet",data_type = data_type),
        tibble(
          AUC = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$AUC$auc) %>% 
            as.numeric,
          Accuracy = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$overall["Accuracy"]) %>% 
            as.numeric,
          Sensitivity = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Sensitivity"]) %>% 
            as.numeric,
          Specificity = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Specificity"]) %>% 
            as.numeric,
          Precision = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Precision"]) %>% 
            as.numeric,
          Recall = rf_models[[data_type]][[ID]] %>% 
            lapply(function(x) x$metrics[[s]]$ConfusionMatrix$byClass["Recall"]) %>% 
            as.numeric,
          Set = s,Fold = c(1:K),Model = "RRF",data_type = data_type)
      )
    }
    all_metrics[[paste(data_type,ID,sep="_")]] <- do.call(rbind,tmp)
    all_metrics[[paste(data_type,ID,sep="_")]]$task <- ID
  }
}

all_metrics_df <- do.call(rbind,all_metrics) 

all_metrics_df_long <- all_metrics_df %>%
  gather(key = "key",value = "value",-Set,-task,-Fold,-Model,-data_type)

all_metrics_df_long %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification") %>%
                         rev,
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification") %>%
                         rev)) %>%
  group_by(Set,task,key,Model,data_type) %>%
  summarise(m = mean(value),se = sd(value),.groups = "drop") %>% 
  subset(Set == "Validation") %>% 
  subset(data_type == "full") %>% 
  mutate(data_type = factor(data_type,levels = c("bcdem","morphology","full"),
                            labels = c("B.C.","Morphology","Morphology + B.C."))) %>%
  ggplot(aes(x = task,y = m,ymin = m-se,ymax = m+se,
             colour = Model,group = paste(Model,Set,data_type))) +
  geom_point(position = position_dodge(width = 0.5),shape = 3) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  facet_wrap(~ key,scales = "free_y",ncol = 2) + 
  theme_pretty(base_size = 6) +
  theme(axis.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = c("green4","red4"),name = NULL) + 
  scale_y_continuous(expand = c(0,0,0,0.01)) + 
  coord_flip(ylim = c(0,1)) + 
  ggsave("figures/rf-glmnet-all-metrics.pdf",height = 4,width = 4)

all_metrics_df_long %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification") %>%
                         rev,
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification") %>%
                         rev)) %>%
  group_by(Set,task,key,Model,data_type) %>%
  summarise(m = mean(value),se = sd(value),.groups = "drop") %>% 
  subset(Set == "Validation") %>% 
  subset(Model == "glmnet") %>% 
  mutate(data_type = factor(data_type,levels = c("bcdem","morphology","full"),
                            labels = c("B.C.","Morphology","Morphology + B.C."))) %>% 
  ggplot(aes(x = task,y = m,ymin = m-se,ymax = m+se,
             colour = data_type,group = paste(Model,Set,data_type))) +
  geom_point(position = position_dodge(width = 0.7),size = 0.5,shape = 3) + 
  geom_linerange(position = position_dodge(width = 0.7),size = 0.25) +
  facet_wrap(~ key,scales = "free_y",ncol = 2) + 
  theme_pretty(base_size = 6) +
  theme(axis.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values = c("blue4","red4","orange"),name = NULL) + 
  scale_y_continuous(expand = c(0,0,0,0.01)) + 
  coord_flip(ylim = c(0.2,1)) + 
  ggsave("figures/glmnet-all-metrics.pdf",height = 4,width = 4)

all_roc <- list()
roc_list <- list()
  
for (data_type in c("bcdem","morphology","full")) {
  for (ID in IDs) {
    all_probs_preds <- glmnet_models[[data_type]][[ID]] %>%
      lapply(function(x) {
        x$metrics$Validation$AUC[c("original.predictor","original.response")]}) %>%
      lapply(function(x) return(do.call(cbind,x))) %>%
      do.call(what = rbind) %>%
      as.data.frame
    all_probs_preds$original.predictor <- all_probs_preds$original.predictor %>%
      as.character %>% as.numeric()
    full_roc <- roc(all_probs_preds$original.response,
                    all_probs_preds$original.predictor,quiet = T)
    tmp <- get.coords.for.ggplot(full_roc) %>%
      as.tibble %>%
      mutate(data_type = data_type)
    all_roc[[paste(data_type,ID,sep="_")]] <- tmp
    all_roc[[paste(data_type,ID,sep="_")]]$task <- ID
    all_roc[[paste(data_type,ID,sep="_")]]$auc_value <- full_roc$auc
    roc_list[[paste(data_type,ID)]] <- full_roc
  }
}

all_roc_df <- do.call(rbind,all_roc) %>%
  subset(is.finite(threshold)) %>%
  arrange(-threshold) %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  mutate(data_type = factor(data_type,levels = c("bcdem","morphology","full"),
                            labels = c("B.C.","Morphology","Morphology + B.C.")))

ggplot(all_roc_df,aes(colour = data_type)) + 
  geom_abline(slope = 1,intercept = 0,size = 0.25,linetype = 2,alpha = 0.5) +
  geom_line(aes(x = 1-specificity,y = sensitivity),alpha = 0.7,size = 0.5) + 
  facet_wrap(~ task) + 
  theme_pretty(base_size = 6) + 
  scale_colour_manual(values = c("blue4","red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm"),
        panel.spacing = unit(0.8,"lines")) +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  ggsave("figures/roc-curves.pdf",height=2.3,width=2.3)

all_roc_df %>%
  select(task,data_type,auc_value) %>%
  distinct %>%
  ggplot(aes(x = reorder(task,desc(task)),y = auc_value, fill = data_type)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_text(aes(label = sprintf("%.1f%%",auc_value*100)),
            position = position_dodge(width = 0.9),
            hjust = -0.1,size = 2) +
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("blue4","red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.5,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
    ggsave("figures/auc-bars.pdf",height=1.7,width=3)

write.csv(x = subset(all_metrics_df_long,Model == "glmnet"),
          file = "data_output/glmnet-metrics.csv")
write.csv(x = all_roc_df,file = "data_output/glmnet-auroc.csv")

# feature importance ------------------------------------------------------

feature_colours <- c(RBC = "red4",`RBC means` = "red4",`RBC variances` = "red1",
                     WBC = "purple4",`WBC means` = "purple4",`WBC variances` = "purple1",
                     `B.C.` = "green4")
feature_shapes <- c(Mean = 16,Variance = 17,X = 15)
  
cat("\nconsistency within glmnet models\n")
# - compare how often a specific feature is selected
# - compare average ranks of different features

for (data_type in c("bcdem","morphology","full")) {
  for (ID in IDs) {
    all_coefs <- c(1:K) %>% 
      lapply(function(x) {
        C <- as.matrix(coef(glmnet_models[[data_type]][[ID]][[x]]$model))
        o <- data.frame(features = rownames(C),coef = unlist(C[,1]),fold = x)
        rownames(o) <- NULL
        return(o)}) %>%
      do.call(what = rbind)
    all_coef_ranks <- all_coefs %>% 
      group_by(fold) %>%
      mutate(feature_rank = rank(coef,ties.method = "last")) %>%
      group_by(features) %>%
      summarise(average_rank = mean(feature_rank),.groups = "drop") %>%
      arrange(average_rank)
    all_imp_coefs <- all_coefs %>%
      group_by(features) %>%
      subset(coef > 0) %>% 
      summarise(N = length(coef > 0),.groups = "drop") 
    cat(sprintf("\nTask: %s (%s)",ID,data_type))
    print(table(all_imp_coefs$N))
    cat("\nTop 10 features with highest average rank:\n")
    cat(paste(apply(head(all_coef_ranks,10),1,
                    function(x) sprintf(
                      "%s = %.2f (s.d. = %.3f)",
                      x[1],as.numeric(x[2]),as.numeric(x[3]))),
              collapse='\n'))
    cat("\n------------------------------------\n")
  }
}

cat("\nconsistency within rrf models\n")
# - compare average ranks

for (data_type in c("bcdem","morphology","full")) {
  for (ID in IDs) {
    all_coefs <- c(1:K) %>% 
      lapply(function(x) {
        C <- as.matrix(varImp(rf_models[[data_type]][[ID]][[x]]$model)$importance)
        o <- data.frame(features = rownames(C),coef = unlist(C[,1]),fold = x)
        rownames(o) <- NULL
        return(o)}) %>%
      do.call(what = rbind)
    all_coef_ranks <- all_coefs %>% 
      group_by(fold) %>%
      mutate(feature_rank = rank(coef,ties.method = "last")) %>%
      group_by(features) %>%
      summarise(average_rank = mean(feature_rank),
                sd_rank = sd(feature_rank),.groups = "drop") %>%
      arrange(average_rank)
    cat(sprintf("\nTask: %s (%s)",ID,data_type))
    cat("\nTop 10 features with highest average rank:\n")
    cat(paste(apply(head(all_coef_ranks,10),1,
                    function(x) sprintf(
                      "%s = %.2f (s.d. = %.3f)",
                      x[1],as.numeric(x[2]),as.numeric(x[3]))),
              collapse='\n'))
    cat("\n------------------------------------\n")
  }
}

cat("\nconsistency between glmnet and rrf models\n")
# measure kendall's tau between both coefficient estimates

feature_associations <- list()

for (data_type in c("bcdem","morphology","full")) {
  cat("\n")
  for (ID in IDs) {
    all_coefs_glmnet <- c(1:K) %>% 
      lapply(function(x) {
        C <- as.matrix(coef(glmnet_models[[data_type]][[ID]][[x]]$model))
        o <- data.frame(features = rownames(C),coef = unlist(C[,1]),fold = x)
        rownames(o) <- NULL
        return(o)}) %>%
      do.call(what = rbind) %>%
      mutate(model = "glmnet")
    all_coefs_rrf <- c(1:K) %>% 
      lapply(function(x) {
        C <- as.matrix(caret::varImp(rf_models[[data_type]][[ID]][[x]]$model)$importance)
        o <- data.frame(features = rownames(C),coef = unlist(C[,1]),fold = x)
        rownames(o) <- NULL
        return(o)}) %>%
      do.call(what = rbind) %>%
      mutate(model = "RRF")
    non_zero_feature_association_df <- rbind(all_coefs_glmnet,all_coefs_rrf) %>%
      group_by(model,features) %>%
      summarise(f = mean(abs(coef),na.rm=T),.groups = "drop") %>%
      spread(key = "model",value = "f") %>%
      na.omit()
    
    C <- cor.test(non_zero_feature_association_df$glmnet,
                  non_zero_feature_association_df$RRF,
                  method = "kendall",exact = F)
    
    cat(sprintf("Task: %s with %s\n",ID,data_type))
    cat(sprintf("Kendall's tau between RRF and glmnet coefficients = %s (%s)\n",
                format(C$estimate,digits = 2),format(C$p.value,digits = 2)))
    
    feature_associations[[paste(ID,data_type)]] <- non_zero_feature_association_df %>%
      mutate(task = ID,data_type = data_type)
  }
}

feature_associations_df <- feature_associations %>%
  do.call(what = rbind) %>%
  subset(glmnet > 0 & RRF > 0) %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  mutate(WBC = grepl("WBC",features),
         RBC = grepl("RBC",features),
         means = grepl("means",features),
         variances = grepl("vars",features)) %>% 
  mutate(cell_type = ifelse(WBC,"WBC",ifelse(RBC,"RBC","B.C.")),
         moment = ifelse(means,"Mean",ifelse(variances,"Variance","X"))) %>%
  mutate(features_raw = gsub('\\.','',str_match(features,'[a-z_A-Z0-9]+\\.')))

feature_associations_df %>% 
  subset(data_type == "bcdem") %>% 
  ggplot(aes(x = glmnet,y = RRF,colour = cell_type,shape = moment)) +
  geom_point() + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  theme_pretty(base_size = 6) +
  xlab("glmnet coefficient mangnitude") +
  ylab("RRF feature importance (%)") +
  facet_wrap(~ task) + 
  scale_shape(name = NULL) + 
  scale_colour_manual(values = feature_colours,name = NULL) +
  theme(legend.key.size = unit(0,"cm")) + 
  ggsave("figures/compare-feature-imp-bcdem.pdf",height=2.5,width=3)

feature_associations_df %>% 
  subset(data_type == "morphology") %>% 
  ggplot(aes(x = glmnet,y = RRF,colour = cell_type,shape = moment)) +
  geom_point() + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  theme_pretty(base_size = 6) +
  xlab("glmnet coefficient mangnitude") +
  ylab("RRF feature importance (%)") +
  facet_wrap(~ task) + 
  scale_shape(name = NULL) + 
  scale_colour_manual(values = feature_colours,name = NULL) +
  theme(legend.key.size = unit(0,"cm")) + 
  ggsave("figures/compare-feature-imp-morphology.pdf",height=2.5,width=3)

feature_associations_df %>% 
    subset(data_type == "full") %>% 
    ggplot(aes(x = glmnet,y = RRF,colour = cell_type,shape = moment)) +
    geom_point() + 
    scale_y_continuous(trans = 'log10') + 
    theme_pretty(base_size = 6) +
    xlab("glmnet coefficient mangnitude (/1000)") +
    ylab("RRF feature importance (%)") +
    facet_wrap(~ task) + 
    scale_shape(name = NULL) + 
    scale_colour_manual(values = feature_colours,name = NULL) +
    theme(legend.key.size = unit(0,"cm")) + 
    scale_x_continuous(breaks = c(0.0001,0.001,0.01,0.1,1),labels = c(0.1,1,10,100,1000),trans = 'log10') +
    ggsave("figures/compare-feature-imp-full.pdf",height=2.5,width=3)

feature_no <- seq(1,length(features_conversion))
names(feature_no) <- names(features_conversion)

file_connection <- file("data_output/wbc_feature_subset")
feature_associations_df %>% 
  subset(glmnet > 0.001) %>%
  subset(data_type == "full" & cell_type == "WBC") %>% 
  select(features_raw) %>%
  distinct %>%
  mutate(features_raw = feature_no[features_raw]) %>%
  unlist %>%
  sort %>%
  paste(collapse = ',') %>%
  writeLines(file_connection)
close(file_connection)

file_connection <- file("data_output/rbc_feature_subset")
feature_associations_df %>% 
  subset(glmnet > 0.001) %>%
  subset(data_type == "full" & cell_type == "RBC") %>% 
  select(features_raw) %>%
  distinct %>%
  mutate(features_raw = feature_no[features_raw]) %>%
  unlist %>%
  sort %>%
  paste(collapse = ',') %>%
  writeLines(file_connection)
close(file_connection)

# feature group importance

d <- "full"
var_group_importance <- list()

for (ID in IDs) {
  idx <- lapply(glmnet_models[[d]][[ID]],
                           function(x) x$metrics$Validation$AUC$auc) %>%
    which.max
  D <- as.matrix(predict(scales[[d]][[idx]],full_dataset_collection$data[[d]]))
  coefficients <- as.matrix(coef(glmnet_models[[d]][[ID]][[idx]]$model))
  coefficients <- coefficients[2:nrow(coefficients),]
  fx <- D * coefficients
  var_group_importance[[ID]] <- data.frame(
    WBCMeans = rowSums(fx[,grepl("means\\.WBC",colnames(fx))]),
    WBCVars = rowSums(fx[,grepl("vars\\.WBC",colnames(fx))]),
    RBCMeans = rowSums(fx[,grepl("means\\.RBC",colnames(fx))]),
    RBCVars = rowSums(fx[,grepl("vars\\.RBC",colnames(fx))]),
    BCDem = rowSums(fx[,grepl("BCDEM",colnames(fx))])) %>%
    cov %>%
    rowSums %>%
    t %>%
    as_tibble() %>%
    mutate(task = ID)
}

var_group_importance %>% 
  do.call(what = rbind) %>% 
  gather(key = "key",value = "value",-task) %>%
  group_by(task) %>%
  mutate(S = sum(value)) %>%
  mutate(Proportions = value / sum(value)) %>%
  mutate(S = sum(value)^(1/6)) %>%
  mutate(key = factor(key,
                      levels = c("WBCMeans","WBCVars","RBCMeans","RBCVars","BCDem"),
                      labels = c("WBC means","WBC variances","RBC means","RBC variances",
                                 "B.C."))) %>% 
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  ggplot(aes(x = S/2,y = Proportions,fill = key,width = S)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_polar(theta = "y") + 
  facet_wrap(~ task) + 
  theme_pretty(base_size = 6) +
  scale_fill_manual(values = feature_colours,name = NULL) +
  theme(panel.grid = element_blank(),axis.title = element_blank(),
        axis.line = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(),legend.position = "bottom",
        legend.key.size = unit(0.2,"cm"),
        strip.text = element_text(margin = ggplot2::margin()),
        legend.box.margin = ggplot2::margin()) + 
  guides(fill = guide_legend(nrow = 2)) + 
  ggsave("figures/glmnet-feature-group-importance.pdf",height = 3,width = 2.5)

# out-of-sample validation ------------------------------------------------

wbc_cells_summaries_oos <- read_csv(
  "../mile-vice/data_output/wbc_adden_2_summaries.csv",
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

rbc_cells_summaries_oos <- read_csv(
  "../mile-vice/data_output/rbc_adden_2_summaries.csv",
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

wbc_cells_summaries_oos <- wbc_cells_summaries_oos[,wbc_feat_cols] %>%
  gather(key = "key",value = "value",-slide_id,-f,-fine_class,-coarse_class)

rbc_cells_summaries_oos <- rbc_cells_summaries_oos[,rbc_feat_cols] %>%
  gather(key = "key",value = "value",-slide_id,-f,-fine_class,-coarse_class)

wbc_dataset_morphology_oos <- rearrange_data(
  wbc_cells_summaries_oos,"WBC")
rbc_dataset_morphology_oos <- rearrange_data(
  rbc_cells_summaries_oos,"RBC")

full_dataset_morphology_oos <- merge(
  wbc_dataset_morphology_oos,
  rbc_dataset_morphology_oos,
  by = c("slide_id","coarse_class","fine_class"))

full_dataset_morphology_oos <- list(
  data = select(full_dataset_morphology_oos,-slide_id,-coarse_class,-fine_class),
  coarse_labels = full_dataset_morphology_oos$coarse_class,
  fine_labels = full_dataset_morphology_oos$fine_class,
  slide_id = full_dataset_morphology_oos$slide_id
)

data_type <- "morphology"

validation_aucs <- list()
all_roc_val <- list()

par(mfrow = c(2,4),mar = c(2,0,2,0))
for (ID in IDs) {
  best_model_idx <- lapply(
    glmnet_models[[data_type]][[ID]],
    function(x) c(
      as.numeric(x$metrics$Validation$AUC$auc))) %>%
    which.max
  best_model <- glmnet_models[[data_type]][[ID]][[best_model_idx]]$model
  ground_truth <- label_conversion[[ID]][full_dataset_morphology_oos$fine_labels]
  d <- as.matrix(predict(scales[[data_type]][[best_model_idx]],full_dataset_morphology_oos$data))
  
  class_prediction <- predict(best_model,type = "class",newx = d)[,1]
  prob_prediction <- predict(best_model,type = "response",newx = d)[,1]
  dat <- data.frame(
    pred = reverse_label_conversion[[ID]][as.numeric(class_prediction)+1],
    obs = reverse_label_conversion[[ID]][ground_truth+1],
    prob = prob_prediction
  ) %>%
    na.omit
  dat$pred <- factor(dat$pred,levels = reverse_label_conversion[[ID]])
  dat$obs <- factor(dat$obs,levels = reverse_label_conversion[[ID]])
  roc_curve <- roc(dat$obs,dat$prob)
  plot(dat$obs,dat$prob)
  plot(roc_curve)
  validation_aucs[[ID]] <- data.frame(
    auc_value = roc_curve$auc,
    task = c(disease_detection = "Disease detection",
             disease_classification = "Disease classification",
             sf3b1 = "SF3B1mut detection",
             anaemia_classification = "Anaemia classification")[ID],
    data_type = "Morphology"
  )
  tmp <- get.coords.for.ggplot(roc_curve) %>%
    as.tibble %>%
    mutate(data_type = data_type)
  all_roc_val[[paste(data_type,ID,sep="_")]] <- tmp
  all_roc_val[[paste(data_type,ID,sep="_")]]$task <- ID
  all_roc_val[[paste(data_type,ID,sep="_")]]$auc_value <- full_roc$auc
}

all_roc_val_df <- do.call(rbind,all_roc_val) %>%
  subset(is.finite(threshold)) %>%
  arrange(-threshold) %>%
  mutate(task = factor(task,
                       levels = c("disease_detection","disease_classification",
                                  "sf3b1","anaemia_classification"),
                       labels = c("Disease detection","Disease classification",
                                  "SF3B1mut detection","Anaemia classification"))) %>%
  mutate(data_type = factor(data_type,levels = c("bcdem","morphology","full"),
                            labels = c("B.C.","Morphology","Morphology + B.C.")))

validation_aucs_df <- do.call(what = rbind,validation_aucs)

rbind(
  mutate(all_roc_val_df,type = "Independent validation set"),
  mutate(all_roc_df,type = "Cross-validated")
  ) %>%
  ggplot(aes(colour = data_type)) + 
  geom_abline(slope = 1,intercept = 0,size = 0.25,linetype = 3,alpha = 0.5) +
  geom_line(aes(x = 1-specificity,y = sensitivity,linetype = type),alpha = 0.7,size = 0.25) + 
  facet_wrap(~ task) + 
  theme_pretty(base_size = 6) + 
  scale_colour_manual(values = c("blue4","red4","orange"),name = NULL,guide = F) + 
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm"),
        panel.spacing = unit(0.8,"lines")) +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1)) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_linetype_manual(name = NULL,values = c(1,2)) +
  guides(linetype = guide_legend(nrow = 2)) +
  ggsave("figures/roc-curves-val.pdf",height=2.5,width=2.5)

all_roc_df %>%
  select(task,data_type,auc_value) %>%
  distinct %>%
  ggplot(aes(x = reorder(task,desc(task)),y = auc_value, fill = data_type)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = validation_aucs_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 18,colour = "white") + 
  geom_point(data = validation_aucs_df,position = position_dodge(width = 0.9),
             alpha = 1,size = 2,shape = 5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = c("blue4","red4","orange"),name = NULL) + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.1,"cm"),
        legend.key.width = unit(0.2,"cm")) +
  scale_x_discrete() +
  ylab("Area under the curve") +
  coord_flip(ylim = c(0.5,1)) + 
  scale_y_continuous(expand = c(0,0,0.02,0.02)) + 
  scale_color_discrete(guide = F) + 
  ggsave(filename = "figures/auc-bars-w-validation.pdf",height=1.7,width=3)
