#' @importFrom pracma numel inpolygon rem meshgrid
gridpoly <- function(Pt,n,H) {
  Pt <- unique(Pt)
  a <- min(Pt[,1])
  b <- max(Pt[,1])
  c <- min(Pt[,2])
  d <- max(Pt[,2])
  m <- dim(Pt)[1]
  delta_x <- (b-a)/n
  delta_y <- (d-c)/n
  #delta <- norm([delta_x,delta_y])
  x <- seq(a,b,by=delta_x)
  y <- seq(c,d,by=delta_y)
  delta <- min(delta_x,delta_y)
  x <- seq(a,b,by=delta)
  y <- seq(c,d,by=delta)
  X <- meshgrid(x,y)$X
  Y <- meshgrid(x,y)$Y
  X <- matrix(X,ncol=1)
  Y <- matrix(Y,ncol=1)

  ### I do this because of the machine error of R
  Pt2 <- which(abs(Pt) < 2.220446e-16,arr.ind = TRUE)
  Pt3 <- Pt
  Pt3[Pt2] <- 0
  Pt0 <- Pt3

  # We now delete all points which are not inside the polygon
  k <- 1
  #kk <- c()
  while (k <= length(X)) {
    #print(k)
    #kk <- c(kk,k)
    inHole <- c()
    if (numel(H) > 0) {
      for (i in 1:numel(H)) {
        inHole <- cbind(inHole,inpolygon(X[k],Y[k],H[[i]][,1],H[[i]][,2],boundary = TRUE)) ## NEED CHANGE
      }
    }
    inHoleTest <- ifelse(sum(inHole)>=1,1,0)
    if (!inpolygon(X[k],Y[k],Pt0[,1],Pt0[,2],boundary = TRUE)|inHoleTest) {
      if ((k-1) > 0) {
        if ((k+1) <= length(X)) {
          X = c(X[1:(k-1)],X[(k+1):length(X)])
          Y = c(Y[1:(k-1)],Y[(k+1):length(Y)])
        }
        else {
          X = c(X[1:(k-1)])
          Y = c(Y[1:(k-1)])
        }
      }
      else {
        X = c(X[(k+1):length(X)])
        Y = c(Y[(k+1):length(Y)])
      }
    }
    else {
      k = k + 1
    }
  }
  m <- length(X)
  p <- dim(Pt0)[1]
  # q <- dim(H)[1]
  # q <- 0
  v <- cbind(Pt[c(seq(2,p),1),1]-Pt[seq(1,p),1],Pt[c(seq(2,p),1),2]-Pt[seq(1,p),2])
  L <- sqrt(v[,1]^2 + v[,2]^2)
  v <- cbind(v[,1]/L,v[,2]/L)
  q <- c() # q vector for holes

  # For each hole
  if (numel(H) > 0) {
    w <- list()
    LL <- list()
    for (i in 1:numel(H)) {
      H[[i]] <- unique(H[[i]])
      q  <-  c(q,dim(H[[i]])[1])
      if (q[i]>0) {
        w[[i]] <- cbind(H[[i]][c(seq(2,q[i]),1),1]-H[[i]][seq(1,q[i]),1],
                     H[[i]][c(seq(2,q[i]),1),2]-H[[i]][seq(1,q[i]),2])
        LL[[i]] <- sqrt(w[[i]][,1]^2 + w[[i]][,2]^2)
        w[[i]] <- cbind(w[[i]][,1]/LL[[i]],w[[i]][,2]/LL[[i]])
      }
    }
  }
  newX <- c()
  newY <- c()
  for (k in 1:m) {
    IP <- (X[k]-Pt[,1])*v[,1] + (Y[k]-Pt[,2])*v[,2]
    Q <- cbind(Pt[,1]+IP*(v[,1]),Pt[,2]+IP*(v[,2]))
    Indx1 <- which(IP < 0)
    Indx2 <- which(IP > L)
    Q <- as.matrix(Q)
    Pt <- as.matrix(Pt)
    Q[Indx1,] <- Pt[Indx1,]
    Q[Indx2,] <- Pt[rem(Indx2,p)+1,]
    D <- sqrt((X[k]-Q[,1])^2 + (Y[k]-Q[,2])^2)

    # For each hole
    dH=c()
    if (numel(H) > 0) {
      for (i in 1:numel(H)) {
        if (q[i] > 0) {
          H[[i]] <- unique(H[[i]])
          IP2 <- (X[k]-H[[i]][,1])*w[[i]][,1] + (Y[k]-H[[i]][,2])*w[[i]][,2]
          Q2 <- cbind(H[[i]][,1]+IP2*(w[[i]][,1]),H[[i]][,2]+IP2*(w[[i]][,2]))
          Indx1 <- which(IP2 < 0)
          Indx2 <- which(IP2 > LL[[i]])
          Q2[Indx1,] <- H[[i]][Indx1,]
          Q2[Indx2,] <- H[[i]][rem(Indx2,q[i])+1,]
          D2 <- sqrt((X[k]-Q2[,1])^2 + (Y[k]-Q2[,2])^2)
        }
        else {
          D2=Inf
        }
        dH <- c(dH,D2)
      }
    }
    if (min(c(D,dH)) > delta*1/3) { #delta/2 is the thresh hold
      newX <- c(newX,X[k])
      newY <- c(newY,Y[k])
    }
  }
  X <- newX
  Y <- newY

  # Next, we add on bdr points
  # delta is the largest distance allowed between two point on the boundary
  newPts <- c()
  k <- 1
  numofpt=0 #length of newPts
  while (k <= p) {
    j=rem(k,p)+1
    newPts=rbind(newPts,Pt[k,]) #adds the current pt;
    numofpt=numofpt+1
    D=norm(as.matrix(Pt[j,]-Pt[k,]),type = "2") #distance between two adjacent pts;
    numpts=floor(D/(delta)-0.01)     #number of points should be added;
    if (numpts > 0) {
      for (i in 1:numpts){
        newPt <- (i*Pt[j,]+(numpts+1-i)*Pt[k,])/(numpts+1);
        newPts <- rbind(newPts,newPt)
        numofpt <- numofpt+1
      }
    }
    k=k+1
  }
  X <- c(X,newPts[,1])
  Y <- c(Y,newPts[,2])

  # For each hole
  if (numel(H) > 0) {
    for (nH in 1:numel(H)) {
      if (q[nH] > 0) {
        newPts <- c()
        k <- 1
        numofpt <- 0 # length of newPts
        while (k <= q[nH]) {
          j <- rem(k,q[nH])+1
          H[[nH]] <- unique(H[[nH]])
          newPts <- rbind(newPts,H[[nH]][k,]) #adds the current pt;
          numofpt <- numofpt+1
          D <- norm(as.matrix(H[[nH]][j,]-H[[nH]][k,]),type = "2") #distance between two adjacent pts;
          numpts <- floor(D/(delta)-0.01)     #number of points should be added;
          if (numpts>0){
            for (i in 1:numpts){
              newPt <- (i*H[[nH]][j,]+(numpts+1-i)*H[[nH]][k,])/(numpts+1);
              newPts <- rbind(newPts,newPt)
              numofpt <- numofpt+1
            }
          }
          k <- k+1
        }
        X <- c(X,newPts[,1])
        Y <- c(Y,newPts[,2])
      }
    }
  }
  names(X)<-NULL
  names(Y)<-NULL
  return(list(X=X,Y=Y))
}
