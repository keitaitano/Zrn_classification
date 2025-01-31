#title: Zircon classification: Igneous rocks
#author: Keita Itano
#date: 2025/1/29

#### library ####
library(tidyverse)
library(caret) # for classification model
library(compositions) # for clr
library(rsample) # data splitting
library(rattle)  # for plotting decision trees
library(progress) # for progress bar
library(ggrepel) # for data visualization
library(ROSE) # for Undersampling
library(randomForest)

#### function ####
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




#### data import ####
df <- read_csv("zircon_trace_imputation.csv",
               col_types="fffffcfffffcfnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn")
df <- data.frame(index=c(1:nrow(df)),df)
df<- df %>% filter(Label2 != "NA")

#### Pre-processing ####
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
  data.frame(Label2=df_screened3$Label2, .) #index=df_screened3$index,


# Export the data used in the analysis
## sample size
tmp<- rbind(table(df$Label2), table(df_clr$Label2)) 
dir.create("Output")
dir.create("Output/ModelEvaluation")
rownames(tmp) <- c("Raw data", "Screened")
write.csv(tmp, "Output/ModelEvaluation/used_data.csv")

closeAllConnections()

#### Model Training ####

custom_rf <- function(data, target, majority_classes, reduction_ratios, ntree = 100, replace = FALSE) {
  classes <- unique(data[[target]])
  if (!all(majority_classes %in% classes)) {
    stop("The specified majority class is not present in the data")
  }
  # Split into majority and minority classes.
  minority_classes <- setdiff(classes, majority_classes)
  
  # Majority-class data reduction
  majority_data_list <- lapply(majority_classes, function(cls) data[data[[target]] == cls, ])
  new_majority_sizes <- mapply(function(df, ratio) round(nrow(df) * ratio), 
                               majority_data_list, reduction_ratios, SIMPLIFY = TRUE)
  minority_data <- data[data[[target]] %in% minority_classes, ]
  
  trees <- list()
  
  for (i in 1:ntree) {
    # Undersampling: Only reduce the majority classes
    sampled_majorities <- mapply(function(df, size) df[sample(nrow(df), size, replace = FALSE), ],
                                 majority_data_list, new_majority_sizes, SIMPLIFY = FALSE)
    
    # Use all data from minority classes without changes
    sampled_data <- do.call(rbind, c(sampled_majorities, list(minority_data)))
    
    # Train a random forest tree
    model <- randomForest(as.formula(paste(target, "~ .")), data = sampled_data, ntree = 1)
    trees[[i]] <- model
  }
  
  return(trees)
}


set.seed(123)
custom_models <- custom_rf(
  data = df_clr, 
  target = "Label2", 
  majority_classes = c("Acidic", "Intermediate"), 
  reduction_ratios = c(0.1, 0.5), 
  ntree = 100, 
  replace = FALSE
)


# split data
set.seed(123)
train_index <- createDataPartition(df_clr$Label2, p = 0.8, list = FALSE)
train_data <- df_clr[train_index, ]
test_data <- df_clr[-train_index, ]

# Prediction (bagging)
predict_custom_rf <- function(models, test_data) {
  predictions <- sapply(models, function(model) predict(model, test_data, type = "class"))
  majority_vote <- apply(predictions, 1, function(x) names(sort(table(x), decreasing = TRUE))[1])
  return(majority_vote)
}

# Evaluation
predictions <- predict_custom_rf(custom_models, test_data)
conf_matrix <- confusionMatrix(factor(predictions, levels = levels(test_data$Label2)), test_data$Label2)

conf_matrix_text <- capture.output(print(conf_matrix))
writeLines(conf_matrix_text, "Output/ModelEvaluation/Evaluation_TrainedModel.txt")

# Save the trained model
save(rf_model, file = "rf_model.RData")



#### Predict Unknown data ####
unknown_data <- read.csv("test_data.csv")
unknown_data_clr <- unknown_data %>% select(-sampleID) %>% clr(.) %>% data.frame(ID=unknown_data$sampleID, .)
predictions_unknowndata <- predict_custom_rf(custom_models, unknown_data_clr)
data.frame(prediction=predictions_unknowndata, unknown_data) %>% write.csv(., "Output/prediction result.csv")
 
