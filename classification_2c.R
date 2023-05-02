
#title: Zircon classification-2
#author: Keita Itano
#date: 2023/02/20

# library
library(tidyverse)
library(caret) # for classification model
library(compositions) # for clr
library(rsample) # data splitting
library(rattle)  # for plotting decision trees
library(progress) # for progress bar
library(ggrepel) # for data visualization

# function
IQR_screen <- function(x){
  tmp_index <-c()
  for (i in 1:(ncol(x)-1)) {
    element <- x %>% select(-index) %>% colnames() %>% .[i]
    tmp <-  x %>% select(element) %>% quantile(., na.rm=TRUE)
    IQR <- tmp[[4]] - tmp[[2]]
    upper <- tmp[[4]] + IQR * 1.5
    lower <- tmp[[2]] - IQR * 1.5
    index <- which(select(x, element) > upper | select(x, element) < lower)
    tmp_index <-c(tmp_index, x$index[index]) 
  }
  tmp_index <- sort(tmp_index) %>% unique(.)
  return(tmp_index)
}

est_performance <- function(x){
  precision <-c()
  recall <-c()
  accuracy <-c()
  for ( i in 1:ncol(x)) {
    precision[i] <- x[i,i]/sum(x[i,], na.rm=T)
    recall[i] <- x[i,i]/sum(x[,i], na.rm=T)
    accuracy[i] <- x[i,i]
  }
  precision <- precision  %>% replace(is.na(.), 0) %>% mean(.)
  recall<-mean(recall, na.rm = T)
  F1 <- 2*precision*recall/(precision+recall)
  accuracy <- sum(accuracy)/sum(x)
  performance <- c(precision=precision, recall=recall, F1=F1, accuracy=accuracy)
  return(performance)
}

est_performance_class<- function(x){
  class <-colnames(x)
  precision <-c()
  recall <-c()
  accuracy <-c()
  for ( i in 1:ncol(x)) {
    precision[i] <- x[i,i]/sum(x[i,], na.rm=T)
    recall[i] <- x[i,i]/sum(x[,i])
    accuracy[i] <- x[i,i]
  }
  precision <- precision  %>% replace(is.na(.), 0)
  names(precision) <- class
  recall<- recall %>% replace(is.na(.), 0)
  names(precision) <- class
  F1 <- 2*precision*recall/(precision+recall)
  F1 <- F1 %>% replace(is.na(.), 0)
  performance_class <- c(precision=precision, recall=recall, F1=F1)
}

# data import
df <- read_csv("zircon_trace_imputation.csv",
               col_types="fffffcfffffcfnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn")
df <- data.frame(index=c(1:nrow(df)),df)
df<- df %>% filter(Label2 != "NA")

# Data Preprocess-1 
## feature selection
df_new <-
  df %>%
  select(-Ref:-Author, -SiO2:-Methods) %>%
  select(Label2:Label3, Y:Nb, La:Lu, Hf, Th, U, index)
element_order <- c("Y","Nb", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er","Tm", "Yb","Lu", "Hf", "Th","U")

# Data Preprocess-2
## remove noise from mineral inclusion: La > 10 ppm
df_screened <- df_new %>% filter(., La < 10)

# Data Preprocess-3
## remove outliers
Index_removed <-c()
## igneous (except for SIAM type granite)
tmp <- subset(df_screened, is.na(df_screened$Label3))
### Basic
Index_removed <-
  tmp %>%
  filter(Label2 == "Basic") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### Intermediate
Index_removed <-
  tmp %>%
  filter(Label2 == "Intermediate") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### Acidic
Index_removed <-
  tmp %>%
  filter(Label2 == "Acidic") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### Kimberlite
Index_removed <-
  tmp %>%
  filter(Label2 == "Kimberlite") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### Carbonatite
Index_removed <-
  tmp %>%
  filter(Label2 == "Carbonatite") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### Ne-syenite
Index_removed <-
  tmp %>%
  filter(Label2 == "Ne-syenite") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)

### igneous (SIAM type granite)
tmp <- df_screened %>% filter(Label3 != "NA")
### S-type
Index_removed <-
  tmp %>%
  filter(Label3 == "S") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### I-type
Index_removed <-
  tmp %>%
  filter(Label3 == "I") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### A-type
Index_removed <-
  tmp %>%
  filter(Label3 == "A") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)
### M-type
Index_removed <-
  tmp %>%
  filter(Label3 == "M") %>% 
  select(-Label2, -Label3) %>% 
  IQR_screen(.) %>% 
  append(Index_removed, .)

df_screened2 <-  df_screened %>%  subset(., !index %in% Index_removed)

# Data Preprocess-4
## remove missing values
index_missing <-
  df_screened2 %>% 
  select(-Label2:-Label3) %>% 
  na.omit() %>% 
  select(index) %>%
  unlist()

df_screened3 <- subset(df_screened2, index %in% index_missing)

# Data Preprocess-5
## log-ratio conversion
df_clr <- 
  df_screened3 %>% 
  select(-Label2:-Label3, -index) %>%
  clr(.) %>% 
  data.frame(index=df_screened3$index, Label2=df_screened3$Label2, Label3=df_screened3$Label3, .)

# Data Preprocess-6
## PCA
pca <- df_clr %>% select(-index:-Label3) %>% prcomp(., scale=T)
df_pca <- df_clr %>% select(index:Label3) %>% data.frame(., pca$x)

# Export the data used in the analysis
## sample size
tmp<- rbind(table(df$Label2), table(df_clr$Label2)) 
rownames(tmp) <- c("Raw data", "Screened")
write.csv(tmp, "Output/classification_2c/data.csv")
## correlation plot
png("Output/classification_2c/correlation.png")
df_clr %>% select(-index:-Label3) %>% cor() %>% corrplot::corrplot(method="shade", shade.col=NA, tl.col="black", col = colorRampPalette(c("#f39800","white", "#4393C3"))(10)) 
dev.off()
## PCA contribution
sink("Output/classification_2c/pca_contribution.txt")
summary(pca)
closeAllConnections()
## PCA loading
pca$rotation %>% data.frame(.) %>% mutate(element=rownames(.)) %>% 
  ggplot(aes(x=PC1, y=PC2))+ geom_point()+
  geom_text_repel(aes(label=element))+
  geom_hline(yintercept=0, size = 0.2, col="blue")+
  geom_vline(xintercept=0, size = 0.2, col="blue")+
  coord_fixed(ratio=1) +
  theme_bw()
ggsave("Output/classification_2c/PC1_PC2.png")
pca$rotation %>% data.frame(.) %>% mutate(element=rownames(.)) %>% 
  ggplot(aes(x=PC1, y=PC3))+ geom_point()+
  geom_text_repel(aes(label=element))+
  geom_hline(yintercept=0, size = 0.2, col="blue")+
  geom_vline(xintercept=0, size = 0.2, col="blue")+
  coord_fixed(ratio=0.5) +
  theme_bw()
ggsave("Output/classification_2c/PC1_PC3.png")
pca$rotation %>% data.frame(.) %>% mutate(element=rownames(.)) %>% 
  ggplot(aes(x=PC1, y=PC4))+ geom_point()+
  geom_text_repel(aes(label=element))+
  geom_hline(yintercept=0, size = 0.2, col="blue")+
  geom_vline(xintercept=0, size = 0.2, col="blue")+
  coord_fixed(ratio=1) +
  theme_bw()
ggsave("Output/classification_2c/PC1_PC4.png")

# split data into training and test sets
k_fold <- 10
## stratified k-fold cross-validation
df_cv<-c()
df_cv <- vfold_cv(df_clr, strata = "Label2", v = k_fold) 
df_cv_class <-  df_clr$Label2 %>% table() %>%　t()
set.seed(50)
for (i in 1: k_fold){
  df_cv$splits[[i]] %>% assessment() %>% select(Label2) %>% table() %>% t() -> tmp
  df_cv_class <- rbind(df_cv_class, tmp)
}
rownames(df_cv_class) <- c("Pre-splitting", 1:k_fold)
## stratified k-fold cross-validation (PCA data)
df_pca_cv <- vfold_cv(df_pca, strata = "Label2", v = k_fold)
df_pca_cv_class <-  df_pca$Label2 %>% table() %>%　t()
set.seed(32)
for (i in 1: k_fold){
  df_pca_cv$splits[[i]] %>% assessment() %>% select(Label2) %>% table() %>% t() -> tmp
  df_pca_cv_class <- rbind(df_pca_cv_class, tmp)
}
rownames(df_pca_cv_class) <- c("Pre-splitting", 1:k_fold)

print("splitted data")
print(df_pca_cv_class)

# CART
pb <- progress_bar$new(total = k_fold)
hp_cart <- c()
p_train_cart <- c()
p_test_cart <- c()
vi_cart <- data.frame(element=factor(element_order, levels=element_order))
for (i in 1:k_fold) {
  pb$tick()
  Sys.sleep(1 / 100)
  ### dataset
  df_cv$splits[[i]] %>% analysis() %>% select(-index, -Label3) -> df_train
  df_cv$splits[[i]] %>% assessment() %>% select(-index, -Label3) -> df_test
  ### model training
  set.seed(32)
  model_cart <- train(data = df_train, Label2 ~ ., method = "rpart", trControl = trainControl("cv", number = 10))
  ### results of CART
  # fancyRpartPlot(model_cart$finalModel)
  hp_tmp <- model_cart$bestTune  %>% unlist()
  p_train_cart_tmp <- confusionMatrix(model_cart) %>% .$table
  p_test_cart_tmp <- model_cart %>% predict(df_test) %>% confusionMatrix(., df_test$Label2) %>% .$table
  vi_cart_tmp <- varImp(model_cart) %>% .$importance 
  hp_cart[i] <- hp_tmp
  p_train_cart[[i]] <-  p_train_cart_tmp
  p_test_cart[[i]] <-  p_test_cart_tmp
  vi_cart <- vi_cart_tmp %>% mutate(element=factor(rownames(.), levels=element_order)) %>% arrange(., element) %>% left_join(vi_cart, ., by=c("element"="element"))
}
vi_cart <- vi_cart %>% select(-element) %>% t()
colnames(vi_cart) <- element_order
rownames(vi_cart) <- c(1:k_fold)
result_cart<- tibble(hp=hp_cart, cm_train=p_train_cart, cm_test=p_test_cart, vi=vi_cart)

## Summary of analytical results
### hyperperameter
cart_hyperparameter <- c(average=mean(result_cart$hp), sd=sd(result_cart$hp))
### variable importance
#variable_importance
cart_variableimportance <- result_cart$vi %>%  apply(., 2, mean) %>% data.frame(CART=.)
### estimation for training set
performance_cart_train<-
  result_cart$cm_train　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_cart_test<-
  result_cart$cm_test　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_cart <- data.frame(Training=performance_cart_train,Test=performance_cart_test)

performance_cart_train_class<-
  result_cart$cm_train　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_cart_train_class) <- c("precision", "recall", "F1")
colnames(performance_cart_train_class) <- c(colnames(p_train_cart[[1]]))
performance_cart_test_class<-
  result_cart$cm_test　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_cart_test_class) <- c("precision", "recall", "F1")
colnames(performance_cart_test_class) <- c(colnames(p_train_cart[[1]]))

print("DONE: CART")

# Random Forest
pb <- progress_bar$new(total = k_fold)
hp_rf <- c()
p_train_rf <- c()
p_test_rf <- c()
vi_rf <- c()
for (i in 1:k_fold) {
  pb$tick()
  df_cv$splits[[i]] %>% analysis() %>% select(-index, -Label3) -> df_train
  df_cv$splits[[i]] %>% assessment() %>% select(-index, -Label3) -> df_test
  ### model training
  # tune mtry minimizing OOB error.
  # mtry:the number of variables to randomly sample as candidates at each split.
  # OOB: Out-Of-Bag
  fitControl <-trainControl(method = "repeatedcv",number = 10,repeats=3, selectionFunction = "oneSE") # 10-fold, 3 times
  grid <- expand.grid(mtry=1:8)
  set.seed(32)
  model_rf <- train(Label2~.,
                    data=df_train,
                    mothod="rf",
                    trControl=fitControl,
                    tuneGrid=grid)
  hp_rf_tmp <- model_rf$bestTune  %>% unlist()
  p_train_rf_tmp <- confusionMatrix(model_rf) %>% .$table
  p_test_rf_tmp <- model_rf %>% predict(df_test) %>% confusionMatrix(., df_test$Label2) %>% .$table
  vi_rf_tmp <- varImp(model_rf) %>% .$importance 
  hp_rf[i] <- hp_rf_tmp
  p_train_rf[[i]] <-  p_train_rf_tmp
  p_test_rf[[i]] <-  p_test_rf_tmp
  vi_rf[[i]] <- vi_rf_tmp
}
vi_rf_n <- vi_rf %>% unlist() %>% matrix(.,nrow=length(element_order),ncol=k_fold) %>%  t(.) %>% data.frame(.)
colnames(vi_rf_n) <- colnames(df_test) %>% .[-1]
result_rf<- tibble(hp=hp_rf, cm_train=p_train_rf, cm_test=p_test_rf, vi=vi_rf_n)

## Summary of analytical results
### hyperperameter
rf_hyperparameter <- c(average=mean(result_rf$hp), sd=sd(result_rf$hp))
### variable importance
#variable_importance <- 1
rf_variableimportance <- result_rf$vi %>% apply(., 2, mean) %>% data.frame(RF=.)
### estimation for training set
performance_rf_train<-
  result_rf$cm_train　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_rf_test<-
  result_rf$cm_test　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_rf <- data.frame(Training=performance_rf_train,Test=performance_rf_test)

performance_rf_train_class<-
  result_rf$cm_train　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_rf_train_class) <- c("precision", "recall", "F1")
colnames(performance_rf_train_class) <- c(colnames(p_train_cart[[1]]))
performance_rf_test_class<-
  result_rf$cm_test　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_rf_test_class) <- c("precision", "recall", "F1")
colnames(performance_rf_test_class) <- c(colnames(p_train_cart[[1]]))
print("DONE: RF")

# SVM
pb <- progress_bar$new(total = k_fold)
hp_svm <- data.frame(matrix(ncol=2, nrow=0)) %>% rename(sigma =X1, C=X2)
p_train_svm <- c()
p_test_svm <- c()
vi_svm <- c()
for (i in 1:k_fold) {
  pb$tick()
  df_pca_cv$splits[[i]] %>% analysis() %>% select(-index, -Label3) -> df_train
  df_pca_cv$splits[[i]] %>% assessment() %>% select(-index, -Label3) -> df_test
  ### model training
  set.seed(32)
  model_svm <- train(
    Label2 ~.,
    data = df_train,
    method = "svmRadial", #svmLinear #"svmPoly"
    trControl = trainControl("cv", number = 10),
    tuneLength = 10 #tuneGrid = expand.grid(C = seq(0, 3, length = 20))
  )
  p_train_svm_tmp <- confusionMatrix(model_svm) %>% .$table
  p_test_svm_tmp <- model_svm %>% predict(df_test) %>% confusionMatrix(., df_test$Label2) %>% .$table
  vi_svm_tmp <- varImp(model_svm) %>% .$importance
  hp_svm[i,] <-model_svm$bestTune
  p_train_svm[[i]] <-  p_train_svm_tmp
  p_test_svm[[i]] <-  p_test_svm_tmp
  vi_svm[[i]] <- vi_svm_tmp
}
result_svm <- tibble(hp=hp_svm, cm_train=p_train_svm, cm_test=p_test_svm, vi=vi_svm)
print("DONE: SVM")

## Summary of analytical results
### hyperperameter
svm_hyperparameter <- c(average=apply(result_svm$hp, 2, mean), sd=apply(result_svm$hp, 2, sd))
### variable importance
#variable_importance <- 1
svm_variableimportance <-result_svm$vi %>% 
  sapply(.,function(x){apply(x, 1, mean)}) %>% 
  apply(., 1, mean)
### estimation for training set
performance_svm_train<-
  result_svm$cm_train　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_svm_test<-
  result_svm$cm_test　%>% 
  sapply(., est_performance) %>% 
  apply(.,1,mean)
performance_svm<- data.frame(Training=performance_svm_train,Test=performance_svm_test)

performance_svm_train_class<-
  result_svm$cm_train　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_svm_train_class) <- c("precision", "recall", "F1")
colnames(performance_svm_train_class) <- c(colnames(p_train_cart[[1]]))
performance_svm_test_class<-
  result_svm$cm_test　%>% 
  sapply(., est_performance_class) %>% 
  apply(., 1, mean, na.rm=T) %>%
  matrix(., nrow=ncol(p_train_cart[[1]])) %>% t() %>% 
  data.frame()
rownames(performance_svm_test_class) <- c("precision", "recall", "F1")
colnames(performance_svm_test_class) <- c(colnames(p_train_cart[[1]]))

# Summary
predict_peformances <- data.frame(performance_cart, performance_rf,performance_svm) %>%
  rename(CART_training=Training, CART_test=Test,RF_training=Training.1, RF_test=Test.1,SVM_training=Training.2, SVM_test=Test.2)　%>% 
  format(digits = 3)

# export results
## precision, recall, F1 score, accuracy
write.csv(predict_peformances, "Output/classification_2c/predict_all.csv")

##confusion matrix (training data)
file_path <- "Output/classification_2c/ConfusionMatrix_train.txt"
sink(file_path)
"Confusion matrix: CART"
sink(file_path,  append = TRUE)
p_train_cart
sink(file_path,  append = TRUE)
"Confusion matrix: RF"
sink(file_path,  append = TRUE)
p_train_rf
sink(file_path,  append = TRUE)
"Confusion matrix: SVM"
sink(file_path,  append = TRUE)
p_train_svm
sink()
closeAllConnections() #Sink does not release file --> release the file

##confusion matrix (test data)
file_path <- "Output/classification_2c/ConfusionMatrix_test.txt"
sink(file_path)
"Confusion matrix: CART"
sink(file_path,  append = TRUE)
p_test_cart
sink(file_path,  append = TRUE)
"Confusion matrix: RF"
sink(file_path,  append = TRUE)
p_test_rf
sink(file_path,  append = TRUE)
"Confusion matrix: SVM"
sink(file_path,  append = TRUE)
p_test_svm
closeAllConnections()

## Prediction performance indices (train)
file_path <- "Output/classification_2c/prediction_class_train.txt"
sink(file_path)
"CART" 
sink(file_path,  append = TRUE)
performance_cart_train_class
sink(file_path,  append = TRUE)
"RF"
sink(file_path,  append = TRUE)
performance_rf_train_class
sink(file_path,  append = TRUE)
"SVM"
sink(file_path,  append = TRUE)
performance_svm_train_class
closeAllConnections()

## Prediction performance indices (test)
file_path <- "Output/classification_2c/prediction_class_test.txt"
sink(file_path)
"CART" 
sink(file_path,  append = TRUE)
performance_cart_test_class
sink(file_path,  append = TRUE)
"RF"
sink(file_path,  append = TRUE)
performance_rf_test_class
sink(file_path,  append = TRUE)
"SVM"
sink(file_path,  append = TRUE)
performance_svm_test_class
closeAllConnections()

## hyper parameter
file_path <- "Output/classification_2c/hyper perameter.txt"
sink(file_path)
"CART: average_cp, sd_cp" 
sink(file_path,  append = TRUE)
cart_hyperparameter
sink(file_path,  append = TRUE)
"RF"
sink(file_path,  append = TRUE)
rf_hyperparameter
sink(file_path,  append = TRUE)
"SVM"
sink(file_path,  append = TRUE)
svm_hyperparameter
closeAllConnections()

## variable importance
file_path <- "Output/classification_2c/variable importance.txt"
sink(file_path )
"CART" 
sink(file_path,  append = TRUE)
cart_variableimportance
sink(file_path,  append = TRUE)
"RF"
sink(file_path,  append = TRUE)
rf_variableimportance

closeAllConnections()

