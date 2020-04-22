

gSigma <- function(b,v,loca){
  n<-nrow(loca)
  R<-exp(-b*(as.matrix(dist(loca))))
  mat<-v*R
  mat
}