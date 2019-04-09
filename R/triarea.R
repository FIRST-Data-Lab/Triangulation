triarea <- function(V1,V2,V3) {
	x <- V1[1]
	y <- V1[2]
	a <- V2[1]
	b <- V2[2]
	c <- V3[1]
	d <- V3[2]
	A <- ((a-x)*(d-y)-(c-x)*(b-y))/2
	return(A)
}
