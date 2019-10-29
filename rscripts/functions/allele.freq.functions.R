getFreq <- function(data, loc){
  results <- matrix(,nrow(data), 1)
  colnames(results) <- c("allele.freq")
  for(i in 1:nrow(data)){
    if(loc=="auto"){
      auto <- data[i, 3]
      allele.freq <- mean(auto)
    }
  results[i, 1] <- allele.freq
  }
  return(as.matrix(results))
}
