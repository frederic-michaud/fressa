
mod <- function(i,j){
  value <- i%%j
  if(value>0) return(value) else return(j)
}
