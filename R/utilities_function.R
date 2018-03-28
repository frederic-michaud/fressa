
mod <- function(value,modulo){
  new.value <- value%%modulo
  if(new.value>0) return(new.value) else return(modulo)
}
