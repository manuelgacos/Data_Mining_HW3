require(MASS)
require(glmnet)

# Function to obtain the correlation between two variables (i,j) from the data
correlation <- function(data, i, j, indices){
  temp <- data[indices, ] # Allows boot() to select the bootstrap sample
  return( cor(temp[,i],temp[,j]) )
}

# Shuffles 'column' from 'data' 
shuffle <-function(data = NULL,column){
  temp <- sample(1:nrow(data), replace = FALSE)
  output <- data[ ,column]
  output <- output[temp]
  return( output )
}

# Computes the correlation between 'fix' feature and 'change' from 'data'. 
#   'change' feature is shuffled
shuffled.cor <- function(data,fix,change){
  temp1 <- data[,fix]
  temp2 <- shuffle(data = data, column = change)
  return( cor(temp1, temp2) )
}

# Changes the 'i' column from the list 'columns' from 'data' ordering it
#   using 'index'
shuffle.column <- function(i, data, columns, index){
  output <- data[ , columns[i] ]
  output <- output[index]
  return( output )
}

# Shuffles the features indicated by the vector 'columns' in 'data'
shuffle.features <- function(data, columns){
  index <- sample(1:nrow(data), replace = FALSE)
  for (i in 1:length(columns)) {
    temp <- shuffle.column(i, data = data, columns = columns, index = index)
    data[ , columns[i] ] <- temp
  }
  return(data)
}

# Computes the test error for a linear model while shuffling the columns 
#   indicated by the vector 'features' in the training dataset 'data'. 'test'
#   is the test dataset
shuffled.mse <- function(data,features,test){
  train <- shuffle.features(data = data, columns = features)
  fit <- lm(y~., data = train)
  fitted.test <- predict.lm(fit, newdata = subset(test, select = -y) )
  mse <- mean( (test$y - fitted.test )^2 )
  return( mse )
}

# Generates a covariance matrix with a given number of 'features' and a given
#   sample correlation 'sample.corr'
# NOTE: All features will have the same correlation with the others
# NOTE: All features have variance 1
covariance.matrix <- function(features, sample.corr){
  temp <- matrix(0, nrow = features, ncol = features)
  i.upr <- which(upper.tri(temp, diag = FALSE), arr.ind=TRUE)
  temp[i.upr] <- sample.corr
  l.upr <- which(lower.tri(temp, diag = FALSE), arr.ind=TRUE)
  temp[l.upr] <- sample.corr
  output <- temp + diag(x=1, nrow = features, ncol = features)
  return(output)
}

# Creates a set of 'n' observations for 'features' different normal random 
#   varaibles which have a given sample correlation 'sample.corr'. If
#   'empirical' is set to TRUE, the sample correlation will be exactly 
#   'sample.corr'
# NOTE: All features will have mean and variance 1
simulate.normal <- function(n, features, sample.corr, empirical = TRUE){
  myCov <- covariance.matrix(features = features, sample.corr =  sample.corr)
  output <- mvrnorm(n = n, mu = rep(1, features), Sigma = myCov,
                  empirical = empirical)
  return(output)
}

# Used after regsubset(). Indicates which varaibles are present in the best model
#   with size 'size'
best.model <- function(size, summ.models = NULL){
  temp <- summ.models$which
  output <- temp[size,]
  return(as.numeric(output))
}

# Only works if length(x)/k is an integer
# Returns k folds with the indexes of the observations
k.folds <- function(x, k = 10){
  if( is.null(dim(x)) ){
    temp <- sample(1:length(x), size = length(x), replace = FALSE)
    split(temp, ceiling(seq_along(temp)/( length(x)/k )))
  } else {
    temp <- sample(1:dim(x)[1], size = dim(x)[1], replace = FALSE)
    split(temp, ceiling(seq_along(temp)/( dim(x)[1]/k )))
  }
}

# For a vector of length 21 returns the string 'formula' for a linear model
get.formula <- function(vector){
  intercept <- vector[1]
  ind <-which(vector[2:21] != 0)
  temp.variables <- paste0('x', ind)
  temp.variables <- paste(temp.variables, collapse = '+')
  if (intercept == 0){
    temp.variables <- paste0(temp.variables,'-1')
  }
  return(temp.variables)
}

# For a given fold 'i' from the set of folds 'folds', gets the test/train MSE
#   of linear regression
fit.fold.lm <- function(i, folds = NULL, features = NULL, data = NULL){
  temp <- folds[-i]
  ind.train <- do.call(c,temp)
  ind.test <- folds[[i]]
  train <- data[ ind.train, ]
  test <- data[ ind.test, ]
  temp.variables <- get.formula(vector = features)
  fit <- lm( as.formula(paste0('y','~',temp.variables)), 
             data = train)
  mse.train <- mean( (fit$residuals)^2 )
  fitted.test <- predict.lm(fit, newdata = subset(test, select = -y) )
  mse.test <- mean( (test$y - fitted.test)^2 )
  return(c(mse.train, mse.test))
}

# Runs k-fold cross validation (OLS) once and reports train/test errors
single.cv.lm <- function(features, data, k = 10){
  folds <- k.folds(data, k = k)
  info <- lapply(1:k, fit.fold.lm, folds = folds, features = features, 
                 data = data)
  temp <- do.call(rbind, info)
  colnames(temp) <- c('mse.train','mse.test')
  output <- colMeans(temp)
  return(output)
}

# Runs k-fold cross validation (OLS) as many times as indicated by 'repeats' and 
#   reports train/test errors
cv.lm <- function(features, data = NULL, k = 10, repeats = 1){
  info <- lapply(1:repeats, function(z){
    single.cv.lm(features = features, data = data, k = k)
  })
  temp <- do.call(rbind,info)
  output <- colMeans(temp)
  return(output)
}

##############################################################################
# NOTE: This functions were originally used in HW2

# For a glmnet model, 'model', training dataset ('x.train', 'y.train') and test
#   dataset ('x.test', 'y.test') obtains the test/train error
single.lambda.mse <- function(lambda, x.train = NULL, y.train = NULL, x.test = NULL, 
                              y.test = NULL, model = NULL){
  pred.train <- predict(model, s = lambda, newx = x.train)
  mse.train <- mean( (pred.train - y.train)^2 )
  pred.test <- predict(model, s = lambda, newx = x.test)
  mse.test <- mean( (pred.test - y.test)^2 )
  c(mse.train,mse.test, lambda)
}

# For a training dataset ('x.train', 'y.train') and test dataset 
#   ('x.test', 'y.test') obtains the test/train error while running Lasso
#   for a set of tunning parameters 'lambda'
get.mse.lasso <- function(x.train, y.train, x.test, y.test, alpha = 1, 
                          lambda = NULL){
  fit <- glmnet(x = x.train, y = y.train, alpha = alpha, lambda = lambda)
  info <- lapply(lambda, single.lambda.mse, x.train = x.train, y.train = y.train,
                 x.test = x.test, y.test = y.test, model = fit)
  output <- do.call(rbind,info)
  colnames(output) <- c('train','test','lambda')
  output
}

# For a given fold 'i' from the set of folds 'folds', gets the test/train MSE
#   of Lasso regression
fit.fold.lasso <- function(i, folds = NULL, x = NULL, y = NULL, alpha = 1,
                           lambda = NULL){
  train <- folds[-i]
  ind.train <- do.call(c,train)
  ind.test <- folds[[i]]
  x.train <- x[ind.train, ]
  y.train <- y[ind.train]
  x.test <- x[ind.test, ]
  y.test <- y[ind.test]
  output <- get.mse.lasso(x.train = x.train, y.train = y.train, x.test = x.test,
                          y.test = y.test, alpha = alpha, lambda = lambda)
  output
}

# Reports the test/train CV error for the values of lambda
cv.error.lasso <- function(matrix){
  df <- as.data.frame(matrix)
  df.train <- aggregate(train~lambda, data = df , FUN = mean)
  df.test <- aggregate(test~lambda, data = df , FUN = mean)
  merge(df.train,df.test)
}

# Runs k-fold cross validation once and reports train/test errors for the set
#   of tunning parameters 'lambda'
single.cv.lasso <- function(x,y,alpha = 1,lambda = NULL, k = 5){
  folds <- k.folds(x, k = k)
  info <- lapply(1:k, fit.fold.lasso, folds = folds, x = x, y = y, alpha = alpha,
                 lambda = lambda)
  temp <- do.call(rbind,info)
  output <- cv.error.lasso(temp)
  output
}

# Runs k-fold cross validation as many times as indicated by 'repeats' and reports
#   train/test errors for the set of tunning parameters 'lambda'
cv.lasso <- function(x, y, alpha = 1, lambda = NULL, k = 5, repeats = 1){
  info <- lapply(1:repeats, function(z){single.cv.lasso(x = x, y = y ,
                                                        alpha = alpha,
                                                        lambda = lambda, k = k)})
  temp <- do.call(rbind,info)
  output <- cv.error.lasso(temp)
  output
}

# Gets the best tunning parameter for a model ran with cv.lasso() using MSE train
best.lambda.train <- function(model){
  model$lambda[which.min(model$train)]
}

# Gets the best tunning parameter for a model ran with cv.lasso() using MSE test
best.lambda.test <- function(model){
  model$lambda[which.min(model$test)]
}

###############################################################################

# boots.sample.index() Returns the indexes for a bootrstrap sample of 'data'
boots.sample.index <- function(data){
  if( is.null(dim(data)) ){
    temp <- sample(1:length(data), size = length(data), replace = TRUE)
    temp
  } else{
    temp <- sample(1:dim(data)[1], size = dim(data)[1], replace = TRUE)
    temp
  }
}

# boots.sample() Returns a bootstrap sample of 'data'  
boots.sample <- function(data){
  index <- boots.sample.index(data = data)
  as.data.frame(data[index,])
}

# Fits lasso using 'data' for a given vector of values for lambda 'lambda'.
#   Returns a 0-1 vector indicating which features were assigned non-zero
#   coefficient estimates
lasso.non.zero <- function(data = NULL, alpha = 1, lambda = NULL){
  x <- as.matrix(subset(data, select = -y))
  y <- data$y
  fit <- glmnet(x = x, y = y, alpha = alpha, lambda = lambda)
  output <- as.numeric(fit$beta != 0)
  return(output)
}

# Fits group lasso using 'data' for a given vector of values for lambda 'lambda'.
#   Returns a 0-1 vector indicating which features were assigned non-zero
#   coefficient estimates
group.lasso.non.zero <- function(data, group = NULL, lambda = NULL){
  x <- as.matrix(subset(data, select = -y))
  y <- data$y
  fit <- gglasso(x = x, y = y, group = group, loss = 'ls', lambda = lambda)
  output <- as.numeric(fit$beta != 0)
  return(output)
}

# Fot a bootstrap sample 'data' fits a linear model and returns which of the 
#   two variables was selected first. Reports also the correlation between 
#   X1 and Y and the correlation between X2 and Y.
bot.proportion <- function(data){
  fit <- lm(y~., data = data)
  temp <- summary(fit)
  temp1 <- temp$coefficients[,4]
  temp1 <- temp1[2:3]
  temp1 <- as.numeric(temp1 == min(temp1))
  cor1 <- cor(data$x1,data$y)
  cor2 <- cor(data$x2,data$y)
  output <- c(temp1,cor1,cor2)
  return(output)
}

# Fot a dataset 'data' fits creates 500 boostrap samples and applies 
#   the bot.proportion() function to them. Reports the proortion of times the 
#   features X1 and X2 were selected first along with the maximum of the two.
bot.proportions <- function(data){
  bot <- lapply(1:500, function(z){
    boots.sample(data = data)
  })
  results <- lapply(bot, bot.proportion)
  resultados <- do.call(rbind,results)[,1:2]
  resumen <- colSums(resultados)/500
  output <- c(resumen, max(resumen))
  return(output)
}
