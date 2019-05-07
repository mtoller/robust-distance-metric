switch(Sys.info()[['sysname']],
       Windows= {pathSymbol <- '\\'},
       Linux  = {pathSymbol <- '/'},
       Darwin = {pathSymbol <- '/'})
computeDistanceResults <- function(pathToUCR="../data",measure=ensembleMetric){
  source('metrics.R')
  
  delimiter <- nchar(pathToUCR)+nchar('/UCRArchive_2018/')+1
  data_names <- c()
  
  trainXs <- list()
  trainYs <- list()
  testXs <- list()
  testYs <- list()
  
  #load data
  
  for (path in list.dirs(pathToUCR)){
    if (path == pathToUCR) next
    if (path == paste0(pathToUCR,pathSymbol,'UCRArchive_2018')) next
    dataset_name <- substr(path,start = delimiter,stop = nchar(path))
    print(paste0('Loading ',dataset_name,'...'))
    pathTrain <- paste0(path,pathSymbol,dataset_name,'_TRAIN.tsv')
    pathTest <- paste0(path,pathSymbol,dataset_name,'_TEST.tsv')

    train <- read.table(pathTrain)
    test <- read.table(pathTest)
    
    if (anyNA(train) | anyNA(test)){
      #print('Dataset contains NAs, skipping...')
      next
    }
    
    if (nrow(test) > 1000 | nrow(train) > 1000){
      #print('Dataset too large, skipping...')
      next
    }
    data_names <- c(data_names,dataset_name)
    
    trainX <- train[,2:ncol(train)]
    trainY <- train[,1]
    
    testX <- test[,2:ncol(test)]
    testY <- test[,1]
    
    trainXs <- append(trainXs,list(trainX))
    trainYs <- append(trainYs,list(trainY))
    testXs <- append(testXs,list(testX))
    testYs <- append(testYs,list(testY))
  }
  
  #classification accuracy error
  accuracy <- c()
  for (j in 1:length(trainXs)){
    result <- Dknn(trainXs[[j]], trainYs[[j]], testXs[[j]], 1, dis = measure)
    accuracy <- c(accuracy, 1-length(which(result == testYs[[j]]))/length(testYs[[j]]))
  }
  #contamination tolerance
  cont_tolerance <- c()
  for (j in 1:length(testXs)){
    print(data_names[j])
    cont_tolerance <- c(cont_tolerance,contaminationTolerance(testXs[[j]],testYs[[j]],measure = measure))
  }
  
  #imprecision invariance
  invar <- c()
  for (j in 1:length(testXs)){
    print(data_names[j])
    invar <- c(invar,imprecisionInvariance(testXs[[j]],testYs[[j]],measure = measure))
  }
  
  #accuracy_results <- as.list(ensemble_accuracy)
  #names(accuracy_results) <- data_names
  
  return(list(accuracy,cont_tolerance,invar))
  
}

Dknn <- function(trainX, trainY, test, k, dis=ensembleMetric,...){
  i <- 0
  n <- nrow(test)
  apply(test, 1,function(testSeries){
    i <<- i + 1
    if (n > 1000){
      if (i %% 100 == 0){
        cat(paste0("Progress: ",i, " out of ", n, "\n"))
      }
    }
    distanceProfile <- apply(trainX,1,function(trainSeries)dis(trainSeries,testSeries,...))
    return(trainY[order(distanceProfile,decreasing = F)[k]])
  })
}

contaminationTolerance <- function(testX,testY,iter.max=3,measure=ensembleMetric){
  if (iter.max > 0)
  n <- iter.max
  else
    n <- length(testY)
  performance <- rep(F,n)
  classes <- unique(testY)
  n <- length(classes)
  for (i in 1:n){
    x <- as.numeric(testX[which(testY==classes[i])[1],])
    other_classes <- which(testY != classes[i][1])
    results <- rep(F,length(other_classes))
    for(y in 1:length(other_classes)){
      results[y] <- testBreakDown(x,as.numeric(testX[other_classes[y],]),measure=measure,FUN=contaminate)
    }
    performance[i] <- length(which(results))/length(other_classes)
  }
  
  return(mean(performance))
}
imprecisionInvariance <- function(testX,testY,iter.max=3,measure=ensembleMetric){
  if (iter.max > 0)
    n <- iter.max
  else
    n <- length(testY)
  performance <- rep(F,n)
  classes <- unique(testY)
  n <- length(classes)
  for (i in 1:n){
    x <- as.numeric(testX[which(testY==classes[i])[1],])
    other_classes <- which(testY != classes[i][1])
    results <- rep(F,length(other_classes))
    for(y in 1:length(other_classes)){
      results[y] <- testBreakDown(x,as.numeric(testX[other_classes[y],]),measure=measure,FUN=addEpsilon)
    }
    performance[i] <- length(which(results))/length(other_classes)
  }
  
  return(mean(performance))
}
testBreakDown <- function(c1,c2,n=1,measure,FUN=contaminate){
  results <- rep(F,n)
  for (i in 1:n){
    results[i] <- (scaleMetric(measure(c1,c2)) - scaleMetric(measure(c1,FUN(c1)))) > 0
  }
  return(all(results))
}
contaminate <- function(x,k=0.05){
  
  n <- length(x)
  s <- max(c(floor(n*k)-1,1))
  o <- sample(1:n)[1:s]
  x[o] <- 1e+100*rnorm(1)
  return(x)
}
addEpsilon <- function(x,rand=rnorm,...){
  n <- length(x)
  s <- sd(x)
  return(x+rand(n,sd=s/n))
}