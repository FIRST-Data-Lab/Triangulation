#' @importFrom pracma inpolygon
del_tri <- function(Pt,V,Tr,order=NULL){
  n <- dim(Tr)[1]
  i <- 1
  # if (is.null(order)){
  #   order = polyord(Pt)
  # }
  while (i <= n) {
    C <- (V[Tr[i,1],]+V[Tr[i,2],]+V[Tr[i,3],])/3
    if ((!inpolygon(C[1],C[2],Pt[,1],Pt[,2],boundary = TRUE ) & order > 0) |
        (inpolygon(C[1],C[2],Pt[,1],Pt[,2],boundary = TRUE) & order < 0)) {
      if ((i-1) > 0) {
        if ((i+1) <= n) {
             Tr <- rbind(Tr[1:(i-1),],Tr[(i+1):n,])
        }
        else {
          Tr <- Tr[1:(i-1),]
        }
      }
      else {
        Tr <- Tr[(i+1):n,]
      }
      n <- n - 1
    }
    else {
      i <- i+1
    }
  }
  return(list(NewTr =Tr))
}
