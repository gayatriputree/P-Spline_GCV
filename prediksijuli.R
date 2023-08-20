x=matrix(c(15114.95),
         1,
         1)
x
library(splines)
x=c(15114.95)
prediksi<-function(m,k,betatopi,lambda){
  #knot<-as.matrix(k)
  #k<-length(k)
  #betatopi<-as.matrix(betatopi)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  betatopi = as.matrix(betatopi)
  yhat = B %*% betatopi
  prediksi = cbind(x,yhat)
  prediksi_J = as.matrix(prediksi)
  cat("Nilai prediksi adalah \n \n")
  print(prediksi_J)
}
betatopi=c(13707.95, 15307.64, 14001.73, 15396.86, 15389.91)
k=c(13844.25, 14342.79)
lambda=0.01
m=3

prediksi(m,k,betatopi,lambda)
prediksi
