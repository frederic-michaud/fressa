#' give the modulus of a number
#'
#' This function return the modulus of a number. The difference with
#' the R built-in function is that the value of n\%n is n instead of 0.
#' @return The remainder of the division of value by modulo
#' @param value A number from which you want to compute the modulo
#' @param modulo The value of the modulus.
#' @examples
#' mod(15,12)
#' mod(12,12)
mod <- function(value,modulo){
  new.value <- value%%modulo
  if(new.value>0) return(new.value) else return(modulo)
}
