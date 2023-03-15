library(gamlss.dist)
library(nleqslv)
###### parameters of the Simplex process ###########
mu0<-0.2
S0<-0.5 # sigma parameter
sigma0<-S0
########################
# calculate mean and variance of the simplex distribution
funX<-function(q1){q1*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX<-integrate(funX,0,1)$value
funX2<-function(q1){q1^2*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX2<-integrate(funX2,0,1)$value
VARX<-MUX2-MUX^2
cat(' mean of X:',MUX," variance of X:",VARX," stdev of X:",sqrt(VARX),'\n')
##############
# define a function for the IC ARL in terms of L
# here is denotes as Ly
# default value for lambda is 0.20
# change for any other value, e.g. 0.05 or 0.10
ARLinfun<-function(Ly,lambda=0.20){
m<-101 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
# control limits
LCL<-MUX-Ly*sqrt(VARX*lambda/(2-lambda))
UCL<-MUX+Ly*sqrt(VARX*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
Q0<-matrix(0,ncol=M,nrow=M)
###############################
# fill the IC Transition Probability Matrix 
  for(i in 1:M){
    for(j in 1:M){
upperL<-LCL+(delta*(j-(1-lambda)*(i-0.5)))/lambda
lowerL<-LCL+(delta*(j-1-(1-lambda)*(i-0.5)))/lambda
####
if(upperL<0||lowerL<0){Q0[i,j]<-0}else{
if(upperL>0&upperL<1&lowerL<0){Q0[i,j]<-pSIMPLEX(upperL,mu=mu0,sigma=S0)}else{
if(upperL>0&upperL<1&lowerL<1&lowerL>0&upperL>lowerL){Q0[i,j]<-pSIMPLEX(upperL,mu=mu0,sigma=S0)-pSIMPLEX(lowerL,mu=mu0,sigma=S0)}
}
}
####
}
}
avec<-rep(0,M);avec[m+1]<-1 # initial probabilities vector
l1<-rep(1,M) # vector of 1s
ID<-diag(M) # Identity Matrix
ARLin<-avec%*%solve(ID-Q0)%*%l1 # matrix formula for ARL
# solve the following equation and determine Ly
as.vector(ARLin)-370.4
}
#### apply nleqslv function for non-linear equation solving
Obj1<-nleqslv(c(2.45),ARLinfun,method="Newton",global="qline",control=list(btol=0.00001,maxit=100))
cat(" L:",Obj1$x,'\n')
# clear everything at the end
rm(list=ls(all=TRUE))

