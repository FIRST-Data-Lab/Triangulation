#' Triangulation Plot
#'
#' This function plots a 2D triangulation.
#' @importFrom graphics plot
#' @importFrom pracma numel isempty meshgrid
#' @param V an \code{N} by two matrix that lists vertices with the \code{i}th
#' row storing in Cartesian coordinates for the
#' \code{i}th vertex. \code{N} is the number of vertices.
#' \cr
#' @param Tr a \code{K} by three matrix that each row represents one triangle.
#' All the elements are the integers that stand for the indices of vertices in \code{V}.
#' \cr
#' @param col A specification for the plotting color, defaulting to 1.
#' \cr
#' @param lwd The line width, a positive number, defaulting to 1. The interpretation is device-specific,
#' and some devices do not implement line widths less than one.
#' @return A triangulation plot of a 2D region.
#'
#' @details
#'
#' @examples
#' # rectangular domain
#' bb=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
#' VT=TriMesh(bb,2)
#' TriPlot(VT$V,VT$Tr)
#'
#' # irregular domains
#' data("horseshoe")
#' VT=TriMesh(horseshoe,n=8)
#' TriPlot(VT$V,VT$Tr)
#'
#' data('shape')
#' VT=TriMesh(shape,15)
#' TriPlot(VT$V,VT$Tr)
#'
#' data('weird')
#' VT=TriMesh(weird,25)
#' TriPlot(VT$V,VT$Tr)
#'
#' # region with holes
#' data("BMP")
#' VT=TriMesh(BMP$bound,25,list(as.matrix(BMP$H1),as.matrix(BMP$H2)))
#' TriPlot(VT$V,VT$Tr)
#'
#' data("mymontreal")
#' VT=TriMesh(mymontreal$bound,25,list(mymontreal$H1,mymontreal$H2))
#' TriPlot(VT$V,VT$Tr)
#' @export

TriPlot<-function(V, Tr, col=1, lwd=1){
  if(ncol(V)!=2){
    stop("Matrix of vectice should have 2 columns.")
  }
  Vmin=apply(V,2,min); Vmax=apply(V,2,max); Vplot=rbind(Vmin,Vmax);
  tri.plot=plot(x=Vplot,col="white",xlab="",ylab="",axes=F)
  apply(Tr,1,tplot,V0=V,col=col,lwd=lwd)
  invisible()
}

tplot<-function(V0, Tr0, col=1, lwd=1){
  V1=V0[Tr0[1],]
  V2=V0[Tr0[2],]
  V3=V0[Tr0[3],]
  lines(x=rbind(V1,V2),col=col,lwd=lwd)
  lines(x=rbind(V2,V3),col=col,lwd=lwd)
  lines(x=rbind(V3,V1),col=col,lwd=lwd)
}
