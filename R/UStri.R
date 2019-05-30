install_github("funstatpackages/Triangulation")
library(Triangulation)
bb=read.csv("boundusa.csv",head=FALSE)
VT=TriMesh(bb,6) # rough triangulation
VT=TriMesh(bb,9) # fine triangulation