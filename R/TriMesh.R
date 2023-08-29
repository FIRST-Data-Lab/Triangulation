#' Mesh generator and Delaunay triangulator for 2D/3D Domains
#'
#' This function triangulates the polygonal domain by using Delaunay Triangulation.
#' @importFrom graphics plot
#' @importFrom tripack tri.mesh triangles
#' @importFrom pracma numel isempty meshgrid
#' @importFrom proxy dist
#' @param bdy A two by \code{N} matrix which indicates the outer boundry points of a 2D region.
#' @param n An integer parameter controlling the fineness of the triangulation
#' and subsequent triangulation. As n increases the fineness increases. Usually, \code{n = 8} seems to be a
#' good choice.
#' @param holes A list of vertices that are the inner boundary points,
#' default set to '\code{NULL}' if there is no holes.
#'
#' @return
#' \item{V}{an \code{N} by two matrix that lists vertices with the \code{i}th
#' row storing in Cartesian coordinates for the
#' \code{i}th vertex. \code{N} is the number of vertices.}
#' \item{Tr}{a \code{K} by three matrix that each row represents one triangle.
#' All the elements are the integers that stand for the indices of vertices in \code{V}.}
#' @details In the function, we firstly get grid points inside and on the boundary of
#' the polygon with extreme points \code{bdy} and interior holes defined by \code{holes}. Then delaunay triangulation
#' is used to generate triangulations by using the grid points.
#' And lastly we delete triangles within the holes or outside the boundary of the region.
#'
#' @examples
#' # rectangular domain
#' bb = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' VT = TriMesh(bb, 3)
#' typeof(VT$V)
#' typeof(VT$Tr)
#'
#' # irregular domains
#' data("horseshoe")
#' TriMesh(hs, n = 9)
#'
#' pt = rbind(c(-0.475, -0.5), c(0.020, -0.5), c(0.435, -0.5), c(0.890, -0.5),
#' c(1.345, -0.5), c(1.800, -0.5), c(2.255, -0.5), c(2.710, -0.5), c(3.165, -0.5),
#' c(-0.475, 0.5), c(0.020, 0.5), c(0.435, 0.5), c(0.890, 0.5), c(1.345, 0.5),
#' c(1.800, 0.5), c(2.255, 0.5), c(2.710, 0.5), c(3.165, 0.5))
#' VT = TriMesh(hs, n = 4, pt = pt)
#'
#' data("shape")
#' TriMesh(shape, 15)
#'
#' data("weird")
#' TriMesh(weird, 30)
#'
#' data("USbb")
#' TriMesh(USbb, 15)
#'
#' # domains with holes
#' data("BMP")
#' TriMesh(BMP$bound, 25, holes = list(as.matrix(BMP$H1), as.matrix(BMP$H2)))
#'
#' data("mymontreal")
#' TriMesh(mymontreal$bound, 25, holes = list(mymontreal$H1, mymontreal$H2))
#' @export

TriMesh <- function(bdy, n, pt = NULL, holes = NULL) {
  X <- gridpoly(bdy, n, holes)$X
  Y <- gridpoly(bdy, n, holes)$Y
  X[abs(X) < 1e-12] <- 0
  Y[abs(Y) < 1e-12] <- 0
  tmp <- cbind(X, Y)
  tmp <- unique(tmp)
  if(!is.null(pt)){
    sx = abs(tmp[nrow(tmp), 1] - tmp[nrow(tmp)-1, 1])
    sy = abs(tmp[nrow(tmp), 2] - tmp[nrow(tmp)-1, 2])
    dat1 = data.frame(X = c(pt[, 1]), Y = c(pt[, 2]))
    dat2 = data.frame(X = c(tmp[, 1]), Y = c(tmp[, 2]))
    dis = dist(dat1, dat2, method = "euclidean")
    # ind = which(dis < sqrt(sx^2 + sy^2), arr.ind = TRUE)
    ind = which(dis < max(sx, sy), arr.ind = TRUE)
    ind.del = unique(ind[, 2])
    # points(tmp[ind.del, ], pch = 20, col = 2)
    X <- c(pt[, 1], tmp[-ind.del, 1])
    Y <- c(pt[, 2], tmp[-ind.del, 2])
  }else{
    X <- tmp[, 1]
    Y <- tmp[, 2]
  }

  tt <- tri.mesh(X, Y)
  Tr <- triangles(tt)[, 1:3]
  Tr <- as.matrix(Tr)
  V <- cbind(X, Y)
  Tr <- del_tri(bdy, V, Tr, 1)$NewTr
  Tr <- as.matrix(Tr)
  if (!isempty(holes)) {
    for (i in 1:numel(holes)) {
      Tr <- del_tri(holes[[i]], V, Tr, -1)$NewTr
    }
  }

  area <- c()
  # Tr <- as.data.frame(Tr)
  for (i in 1:nrow(Tr)) {
    area <- c(area, triarea(V[Tr[i, 1], ], V[Tr[i, 2], ], V[Tr[i, 3], ]))
  }
  dT <- which(area < 1e-12)
  if (length(dT) > 0) {
    Tr <- Tr[-dT, ]
  }
  TriPlot(V, Tr)
  return(list(V = V, Tr = Tr))
}
