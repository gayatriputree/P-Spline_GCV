install.packages("splines")
library(splines)

MPL<-function(x,eps=1e-009)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1]) %*% t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)] %*% diag(1/diago) %*% t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

#memanggil data variabel prediktor dan respon
data=read.csv("datakurs.csv", header=TRUE, sep = ";")
data
x= data$x[1:42]
x
y=data$y[1:42]
y

#memanggil gcv dengan knot 1
gcv1<-function(x,y,m,k,lamda)
{
  n=length(y) #panjang y
  gcv=10^10 #untuk mendefinisikan data
  t1=seq(min(x),max(x),length.out=100) #maksimum dari x dan min x
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,1,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4) #untuk orde
  {
    k1=c1[,1] #
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=k1[j1], degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      #for (k in 1:pord) D <-diff(D)
      BtB = t(B) %*% B #buat mencari parameter
      deter=det(BtB)
      z = MPL(BtB)
      beta= MPL(BtB + lamda * t(D) %*% D) %*% (t(B) %*% y)
      #yhat <- B %*% beta
      L= B %*% MPL(BtB + lamda * t(D) %*% D) %*% t(B)
      yhat = B %*% beta
      MSE = t(y-yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-L))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        det=deter
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}
gcv1(x,y,m,k,0.01)
#gcv1(x,y,m,k,0.02)
#gcv1(x,y,m,k,0.03)
#gcv1(x,y,m,k,0.04)
#gcv1(x,y,m,k,0.05)
#gcv1(x,y,m,k,0.06)
#gcv1(x,y,m,k,0.07)
#gcv1(x,y,m,k,0.08)
#gcv1(x,y,m,k,0.09)
#gcv1(x,y,m,k,0.1)

#memanggil gcv dengan knot 2
gcv2<-function(x,y,m,k,lamda)
{
  n=length(y)
  gcv=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,2,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB)
      z = MPL(BtB)
      beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
      L = B %*% MPL(BtB + lamda*t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y-yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-L))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}
gcv2(x,y,m,k,0.01)
#gcv2(x,y,m,k,0.02)
#gcv2(x,y,m,k,0.03)
#gcv2(x,y,m,k,0.04)
#gcv2(x,y,m,k,0.05)
#gcv2(x,y,m,k,0.06)
#gcv2(x,y,m,k,0.07)
#gcv2(x,y,m,k,0.08)
#gcv2(x,y,m,k,0.09)
#gcv2(x,y,m,k,0.1)

#memanggil gcv dengan knot 3
gcv3<-function(x,y,m,k,lamda)
{
  n=length(y)
  gcv=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,3,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1) #membuat matriks B
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB)
      z = MPL(BtB)
      beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
      L = B %*% MPL(BtB + lamda*t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y-yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-L))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}
gcv3(x,y,m,k,0.01)
#gcv3(x,y,m,k,0.02)
#gcv3(x,y,m,k,0.03)
#gcv3(x,y,m,k,0.04)
#gcv3(x,y,m,k,0.05)
#gcv3(x,y,m,k,0.06)
#gcv3(x,y,m,k,0.07)
#gcv3(x,y,m,k,0.08)
#gcv3(x,y,m,k,0.09)
#gcv3(x,y,m,k,0.1)

#memanggil gcv dengan knot 4
gcv4<-function(x,y,m,k,lamda)
{
  n=length(y)
  gcv=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,4,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    k4=c1[,4]
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1],k4[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB)
      z = MPL(BtB)
      beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
      L = B %*% MPL(BtB + lamda*t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y-yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-L))
      #    AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(4+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        knot4=k4[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "knot 4 = ", knot4, "GCV = ", GCV, "determinan = ", deter,"\n")
    gcv=10^10
  }
}
gcv4(x,y,m,k,0.01)

#estimasi parameter
pspline<-function(x,y,m,k,lamda){
  n<-length(y)
  print(n)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  bs1 = bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
  cat("Nilai parameter adalah", "\n", beta, "\n")
  L = B %*% MPL(BtB + lamda*t(D) %*% D) %*% t(B)
  yhat = (B %*% beta)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  I <- matrix(0, ncol=n, nrow = n)
  for (j in 1:n) 
    I[j,j] = 1
  l = sum(diag(I-L))
  GCV=(MSE/(l/n)^2)
  cat("NIlai GCV adalah ","\n",GCV,"\n")
}
pspline(x,y,3,k=c(13844.25, 14342.79),0.01)

#pengujian parameter
Uji_Parameter<-function(x,y,m,k,lamda){
  n<-length(y)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  orde = (m-1)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
  Beta = as.matrix(beta)
  cat("Nilai parameter adalah", "\n", "\n")
  print(Beta)
  L = B %*% MPL(BtB + lamda*t(D) %*% D) %*% t(B)
  yhat = (B %*% beta)
  ybar <- mean(y)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  MSE <- as.numeric(MSE)
  SSR <- sum((yhat-ybar)^2)
  MSR <- SSR/((orde+k1)-1)
  SSE <- sum((y-yhat)^2)
  MSER <- SSE/(n-(orde+k1))
  JKT = sum((y-ybar)^2)
  Fhit_1 = MSR/MSER
  cat("Nilai F hitung pengujian serentak adalah","\n", Fhit_1,"\n")
  cat("Kesimpulan hasil uji serentak","\n")
  cat("-----------------------------------","\n")
  cat("Analysis Of Variance (ANOVA)","\n")
  cat("===================================================","\n")
  cat("Sumber df   SS    MS      Fhitung","\n")
  cat("Regresi",((orde+k1)-1), "",SSR,"",MSR,"",Fhit_1,"\n")
  cat("Error", (n-(orde+k1)), "", SSE, "", MSER, "", "\n")
  cat("Total ", (n-1), "", JKT, "", "\n")
  print(MSE)
  cat("===================================================", "\n")
}
Uji_Parameter(x,y,3,k=c(13844.25, 14342.79),0.01)

#Melihat f tabel
qf(p = 0.05,df1 =1,df2 =102,lower.tail = FALSE)

#Menghitung nilai prediksi
prediksi<-function(x,y,m,k,lamda){
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL(BtB + lamda*t(D) %*% D) %*% (t(B) %*% y)
  Beta = as.matrix(beta)
  yhat = (B %*% beta)
  prediksi = cbind(y,yhat)
  prediksi_B = as.matrix(prediksi)
  cat("Nilai prediksi adalah \n \n")
  print(prediksi_B)
}
prediksi(x,y,3,k=c(13844.25, 14342.79),0.01)
yhat = prediksi(x,y,3,k=c(13844.25, 14342.79),0.01)
MAPE = mean(abs(y-yhat)/y)*100
MAPE

#uji individu #F tidak terpenuhi
SEbeta <- matrix(0, nrow = (orde+k), ncol = (orde+k))
Ai <- solve(wtw)
SEbeta <- sqrt(diag(MSE * Ai))
SEbeta <- as.matrix(SEbeta)
Thit <- bet/SEbeta
Thit
wtw