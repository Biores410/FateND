
#' @name  get_genepair
#' @title get all genepairs by pairing the g1 and g2
#' @description  function to get gene pair
#' @param A refer to g1
#' @param B refer to g2
#' @return  all gene pairs
#' @export
get_genepair <- function(A,B){
  
  overall_genepair=merge(A,B)
  overall_genepair=as.matrix(overall_genepair)
  overall_genepair=unique(overall_genepair)
  overall_genepair <- overall_genepair[!apply(overall_genepair, 1, function(row) all(row == row[1])), ]
  
  U <- c(A,B)
  results=list(overall_genepair, U)
  return(results)
}