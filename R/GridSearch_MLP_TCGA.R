#required libraries
set.seed(1234)
library(AMoNet)
library(keras)
library(tensorflow)
library(ggplot2)
library(survival)


# Generate AMoNet model

## Initiate
histo<-"lusc|coad|luad|read|brca|stad|ucec|hnsc|blca|cesc|esca|chol"
GENEman <- c("AKT1", "APC", "BRAF", "BRCA2", "CDKN2A", "CTNNB1", "ERBB2", "ERBB3", "FBXW7", "FGFR3", "HRAS", "IDH1", "KRAS", "MAP2K4", "MAP3K1", "NF1", "NFE2L2", "NRAS", "PIK3CA", "PPP2R1A", "PTEN", "SMAD4", "SMARCB1", "STK11", "TP53", "MTORC1", "MTOR")
net<-AMoNet(GENESman = GENEman)

## Build
net<-build.AMoNet(net, RestrictedBuilding = T, MeanWinit = 0.001, SdWinit = 0.1, Phenotypes=GenesSelecHall, MECA=AMoNet::names_MECA, RestrictedBase = T, Interval = 10, Optimizer = "Adam", KeepPhenotypes=T, no_cores = 3, NameProj = "AMoNet_TCGA")


## Load TCGA data
Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)
net<-LoadCleanTCGA(net, Species = Species, Param = "", RestrictUnique = T, organ = histo)
net<-RandomiStates(net)

## Split data into Train and Validation set
net <- split(net, 0.8)

## Input data format for Keras MLP
input_tensor <- net$Data$MUTa
input_tensor[is.na(input_tensor)]<-0
Target<-t(net$Data$y)
Train<-net$TrainSplit$Train
Val<-net$TrainSplit$Val

input_train = as.matrix(input_tensor[Train,])
Target_train = as.matrix(Target[Train,])
input_val = as.matrix(input_tensor[Val,])
Target_val = as.matrix(Target[Val,])

## Model building
build_model <- function(dropout = 0.4, learning_rate = 0.003, l2_reg = 0.01) {
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = 682, activation = 'relu', input_shape = c(208)) %>%  
    layer_dropout(rate = dropout) %>% 
    layer_dense(units = 3862, activation = 'relu') %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 490, activation = 'relu') %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 10, activation = 'sigmoid', kernel_regularizer = regularizer_l2(l = l2_reg))
  
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_rmsprop(learning_rate = learning_rate), 
    metrics = "mse"
  )
  
  return(model)
}

## Setting up the parameters for grid search
param_grid <- expand.grid(
  dropout_rate = c(0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9),
  learning_rate = c(0.001, 0.01, 0.1),
  l2_reg = c(0.0001, 0.001, 0.01)
)

## Initializing the data.frame to store the performances
performance_data <- data.frame(dropout_rate = numeric(),
                               learning_rate = numeric(),
                               l2_reg = numeric(),
                               c_index_train = numeric(),
                               mse_train = numeric(),
                               c_index_val = numeric(),
                               mse_val = numeric())

best_c_index_val <- 0

## Iteration over each set of parameters in the grid
for(i in 1:nrow(param_grid)) {
  current_params <- param_grid[i,]
  
  # Retrieve the current values of the parameters
  dropout_rate <- current_params$dropout_rate
  learning_rate <- current_params$learning_rate
  l2_reg <- current_params$l2_reg
  
  # Build and train the model
  model <- build_model(dropout = dropout_rate, learning_rate = learning_rate, l2_reg = l2_reg)
  history <- model %>% fit(input_train, Target_train, epochs = 5, validation_data = list(input_val, Target_val))
  
  # Evaluate on the training set
  pred_train <- model %>% predict(input_train)
  pred_train_mean = rowMeans(pred_train)
  c_index_train <- concordance(Surv(net$Data$SurvData[1, Train], net$Data$SurvData[2, Train])~pred_train_mean)
  c_index_train <- round(c_index_train$concordance, 2)
  mse_train <- mean((Target_train - pred_train)^2)
  
  # Evaluate on the validation set
  pred_val <- model %>% predict(input_val)
  pred_val_mean = rowMeans(pred_val)
  c_index_val <- concordance(Surv(net$Data$SurvData[1, Val], net$Data$SurvData[2, Val])~pred_val_mean)
  c_index_val <- round(c_index_val$concordance, 2)
  mse_val <- mean((Target_val - pred_val)^2)
  
  # If the current model is better than the previous best model, save the weights
  if(c_index_val > best_c_index_val) {
    best_c_index_val <- c_index_val
    model %>% save_model_weights_hdf5("best_model_weights.h5")
  }
  
  # Store the metrics and parameters for this iteration
  performance_data <- rbind(performance_data,
                            data.frame(dropout_rate = dropout_rate,
                                       learning_rate = learning_rate,
                                       l2_reg = l2_reg,
                                       c_index_train = c_index_train,
                                       mse_train = mse_train,
                                       c_index_val = c_index_val,
                                       mse_val = mse_val))
}


## Sort the performance dataframe by validation c-index and select the first row
best_model_cindex_val <- performance_data[order(-performance_data$c_index_val),][1,]
best_model_cindex_val

## Retrieve the best parameters
best_dropout_rate <- best_model_cindex_val$dropout_rate
best_learning_rate <- best_model_cindex_val$learning_rate
best_l2_reg <- best_model_cindex_val$l2_reg

## Retrain the model with the best parameters on the complete training data set
best_model <- build_model(dropout = best_dropout_rate, learning_rate = best_learning_rate, l2_reg = best_l2_reg)
best_model %>% load_model_weights_hdf5("best_model_weights.h5")
history <- best_model %>% fit(input_train, Target_train, epochs = 10, validation_data = list(input_val, Target_val))

## Extract the metrics from the history object
train_loss <- history$metrics$loss
val_loss <- history$metrics$val_loss

## Create a dataframe for plotting
df <- data.frame(
  epoch = seq_len(length(train_loss)),
  train_loss = train_loss,
  val_loss = val_loss
)

## Plot the learning curve
ggplot(df, aes(x = epoch)) +
  geom_line(aes(y = train_loss, color = 'Train Loss')) +
  geom_line(aes(y = val_loss, color = 'Validation Loss')) +
  ggtitle("Courbe d'apprentissage") +
  xlab("Époque") +
  ylab("Loss") +
  scale_color_manual("", 
                     breaks = c("Train Loss", "Validation Loss"),
                     values = c("blue", "red"))


# Set the plotting areas
par(mfrow=c(3,2)) # 3 rows, 2 columns of sub-plots

# For each hyperparameter, create a subplot
for(hyperparam in names(current_params)){
  # Sub-plot for c_index_train
  plot(performance_data[[hyperparam]], performance_data$c_index_train, main=paste("c_index_train according to", hyperparam), xlab=hyperparam, ylab="c_index_train", col="blue", pch=20)
  
  # Sub-plot for c_index_val
  plot(performance_data[[hyperparam]], performance_data$c_index_val, main=paste("c_index_val according to", hyperparam), xlab=hyperparam, ylab="c_index_val", col="red", pch=20)
  
  # Sub-plot for mse_train
  plot(performance_data[[hyperparam]], performance_data$mse_train, main=paste("mse_train according to", hyperparam), xlab=hyperparam, ylab="mse_train", col="blue", pch=20)
  
  # Sub-plot for mse_val
  plot(performance_data[[hyperparam]], performance_data$mse_val, main=paste("mse_val according to", hyperparam), xlab=hyperparam, ylab="mse_val", col="red", pch=20)
}

## Reset the plotting areas
par(mfrow=c(1,1))

## Prediction of the best model on the training set
pred_train_mat <- best_model %>% predict((input_train))

## Training performance (C-index and MSE)
pred_train = rowMeans(pred_train_mat)
cindex <- concordance(Surv(net$Data$SurvData[1,Train],net$Data$SurvData[2,Train])~pred_train)
cindex_train <- round(cindex$concordance,2)
mse <- mean((Target_train - pred_train_mat)^2)
cat("c-index = ", cindex_train, "\n")
cat("mse = ", mse)

## Plot for the training set
plot(pred_train, net$Data$SurvData[1,Train], pch=ifelse(net$Data$SurvData[2,Train]==0,2,1), main="Training cohort", xlab = "Predicted Survival Probability",ylab = "True survival (months)", xlim=c(0,1))
legend("topleft", legend = c(paste("C-index =", cindex_train )),cex=0.5)


## Prediction of the best model on the validation set
pred_val_mat <- best_model %>% predict((input_val))

## Validation performance (C-index et MSE)
pred_val = rowMeans(pred_val_mat)
cindex <- concordance(Surv(net$Data$SurvData[1,Val],net$Data$SurvData[2,Val])~pred_val)
cindex_val <- round(cindex$concordance,2)
mse <- mean((Target_val - pred_val_mat)^2)
cat("c-index = ", cindex_val, "\n")
cat("mse = ", mse)

## Plot for the validation set
plot(pred_val, net$Data$SurvData[1,Val], pch=ifelse(net$Data$SurvData[2,Val]==0,2,1), main="Validation cohort", xlab = "Predicted Survival Probability",ylab = "True survival(months)", xlim=c(0,1))
legend("topleft", legend = c(paste("C-index =", cindex_val)), cex=0.5)



