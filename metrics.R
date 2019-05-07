require(dplyr)
distance.mean <- function(a,b){
  abs(mean(a)-mean(b)) %>% return
}

testIfMetric <- function(f,iter.max=1000,rand.gen=rnorm,n=1000,epsilon=1e-10,verbose=F,ultra=F){
  #set.seed(0815)
  pass <- T
  for (i in 1:iter.max){
    if (verbose) print(i)
    x <- rand.gen(n)
    y <- rand.gen(n)
    z <- rand.gen(n)
    
    if (f(x,x) > epsilon | f(y,y) > epsilon | f(z,z) > epsilon |
        (f(x,y)<=epsilon & all(x!=y)) | (f(y,z)<=epsilon & all(y!=z)) | 
        (f(x,z)<=epsilon & all(x!=z))){
      print('Identity of indiscernables does not hold (f(x,x)!=0)')
      pass <- F
    }
    
    
    if (f(x,y)!=f(y,x)){
      print('Symmetry does not hold (f(x,y)!=f(y,x))')
      pass <- F
    }
    
    if (f(x,z) > (f(x,y) + f(y,z)+epsilon)){
      print('Triangle inequality does not hold (f(x,z) > f(x,y) + f(y,z))')
      #print(c(x,y,z))
      #print(c(f(x,z),f(x,y),f(y,z)))
      pass <- F
    }
    if (ultra){
      if (f(x,z) > max(c(f(x,y), f(y,z)))){
        print('Strong triangle inequality does not hold (f(x,z) > max(f(x,y), f(y,z)))')
        pass <- F
      }
    }
    if (!pass) return(F)
  }
  return(T)
}

mse <- function(x,y)
{
  return(sd(x-y))
}
smother <- function(x,rand=rnorm,...){
  s <- sd(x)
  x %>% length() %>% rand(sd=s/length(x),...) %>% {.+x} %>% return
}
discreteMetric <- function(x,y){
  if (all(x != y)) return(1)
  else return(0)
}

divAbs <- function(x){
  smalls <- x %>% {which(abs(.) < 1)}
  x[smalls] <- x[smalls]^-1
  return(x)
}

logDistance <- function(x,y){
  return(sum(log((1+abs(x-y)))))
  #{x-y} %>% abs %>% {.+1} %>% log %>% sum %>% {./1} %>% return
}
smotherOnce <- function(x,k=0.1){
  
  n <- length(x)
  s <- floor(n*k)-1
  o <- sample(1:n)[1:s]
  x[o] <- 1e+100*rnorm(1)
  return(x)
}
breakDownTest <- function(c1,c2,measure,n=100,smothering=smotherOnce,...){
  results <- rep(F,n)
  for (i in 1:n){
    results[i] <- (scaleMetric(measure(c1,c2,...)) - scaleMetric(measure(c1,smothering(c1),...))) > 0
  }
  return(all(results))
}
breakDownTestOverDataset <- function(pathData, measure,n=1,iter.max=0,smothering=smotherOnce,...){
  data <- read.table(pathData)
  print(dim(data))
  if (anyNA(data)){
    print("NAs in data")
    return()
  }
  X <- data[,2:ncol(data)]
  Y <- as.numeric(data[,1])
  if (iter.max > 0)
      n <- iter.max
  else
    n <- length(Y)
  performance <- rep(F,n)
  for (i in 1:n){
    #print(n-i+1)
    other_classes <- which(Y != Y[i])
    c1 <- as.numeric(X[i,])
    results <- rep(F,length(other_classes))
    for(y in 1:length(other_classes)){
      results[y] <- breakDownTest(c1,as.numeric(X[other_classes[y],]),measure,n = n,smothering = smothering,...)
    }
    performance[i] <- all(results)
  }
  
  return(length(which(performance))/n)
}
medianCentric <- function(x,y){
  return(euclid(median(x),median(y)))
}

ensembleMetric <- function(x,y,l=0.1){
  m <- rep(0,6)
  if (length(x) != length(y)) stop('x and y must have the same length')
  n <- length(x)
  w <- trunc(n*l)
  
  if (w %% 2 == 0) w <- w + 1
  
  omits <- floor(w/2)
  
  medX <- runmed(x,k = w,endrule = 'constant')[(omits+1):(n-w+1+omits)]
  medY <- runmed(y,k = w,endrule = 'constant')[(omits+1):(n-w+1+omits)]
  
  m[1] <- euclid(x,y)
  m[2] <- fastEDR(x,y)
  m[3] <- logDistance(x,y)
  m[4] <- euclid(medX,medY)
  m[5] <- fastEDR(medX,medY)
  m[6] <- logDistance(medX,medY)
  
  m <- unlist(lapply(m,scaleMetric))
  
  return(sqrt(sum(m^2)))
}

scaleMetric <- function(d){
  return(1-(1/(1+d)))
}
fastEDR <- function(x,y,g=0){
  if (g!=0) return(length(which(abs(x-y) > g)))
  return(length(which(x!=y)))
}
lpnorm <- function(x,y,p){
  return(sum(abs(x-y)^p) ^(1/p))
}