## z score function
# John M, Fall 2023
# provide vector as x. outputs zscore as y

z_score = function(x){
  
  z = (x-mean(x)) / sd(x)
  return(z)  
  
}