extractwindows1 <- function(values,startindices){
  startindices = startindices[!is.na(startindices)]
  startindices = startindices[order(startindices)]
  startindices = startindices[startindices>=0]
  print(startindices)
  result = list();
  if(length(startindices)==0){
      return(NA)
  }
  if(startindices[1]>length(values)){
      return(NA)
  }
  for(i in 1:(length(startindices))){
    print(i)
    if(startindices[i]>length(values)){
      names(result) = paste('bi:',1:length(result),sep='')
      return(result)
    }
    if(i == length(startindices)){
      result[[i]] = values[startindices[i]:length(values)]
      names(result) = paste('bi:',1:length(result),sep='')
      return(result)
    }
    if(startindices[i+1]>=length(values)+1){
      result[[i]] = values[startindices[i]:length(values)]
      names(result) = paste('bi:',1:length(result),sep='')
      return(result)
    }
    result[[i]] = values[startindices[i]:(startindices[i+1]-1)]
  }
}