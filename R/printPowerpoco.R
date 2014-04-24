print.Powerpoco <- function(x, ...){
  if(x[[11]]==1 & x[[12]]=="Dunnett"){
    cat("Any pair power for simultaneous test \nfor non-inferiority \n  \n")
    cat("Power:", x[[1]], "\n")
    cat("per group sample size:", x[[2]],  " \n")
    cat("relative non-inferiority margin:" , x[[9]])
  }else{
    cat("Any pair power for simultaneous test \nfor treatment-by-subgroup interactions \n  \n")
    cat("Power:", x[[1]], "\n")
    cat("treatment-by-subgroup sample size:", x[[2]]," \n")
    cat("relative inconsistency margin:" , x[[9]])
  }
}