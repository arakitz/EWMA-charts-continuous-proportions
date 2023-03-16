library(gamlss.dist)
library(nleqslv)
############### parameters of the beta process
mu0<-0.2
phi0<-31
alpha0<-mu0*phi0 # in-control value of parameter shape1
beta0<-(1-mu0)*phi0 # in-control value of parameter shape2
###################### Design state, Markov chain method
m<-500 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
# for each value lambda, we use the nleqslv function in R to find the unique value for L
# we define a function, named myfun2, for this purpose
# the function is for the case ARL0=370.4. It needs a slight modification to work for any other ARL0 value
for(lambda in seq(0.20,0.20,by=0.05)){
  myfun2<-function(y){
    #### below we determine the entries of the Transition Probability Matrix Q  
    mu0X<-alpha0/(alpha0+beta0)
    sigma0X2<-alpha0*beta0/(((alpha0+beta0)^2)*(alpha0+beta0+1))
    LCLy<-mu0X-y*sqrt(sigma0X2*lambda/(2-lambda))
    UCLy<-mu0X+y*sqrt(sigma0X2*lambda/(2-lambda))
    deltay<-(UCLy-LCLy)/(M)
    Q0y<-matrix(0,ncol=M,nrow=M)
    for(i in 1:M){
      for(j in 1:M){
        lowerly<-LCLy+deltay*(j-1-(1-lambda)*(i-0.5))/lambda
        upperly<-LCLy+deltay*(j-(1-lambda)*(i-0.5))/lambda
        Q0y[i,j]<-pbeta(upperly,shape1 = alpha0,shape2 = beta0)-pbeta(lowerly,shape1 = alpha0,shape2 = beta0)
      }
    }
    avecy<-rep(0,M);avecy[m+1]<-1 # initial probabilities vector
    l1y<-rep(1,M) # vector of ones
    IDy<-diag(M) # Identity Matrix
    ARLyin<-avecy%*%solve(IDy-Q0y)%*%l1y # ARL
    ARLyin-370.4 # we define ARL-ARL0 and with nleqslv we solve the equation ARL-ARL0=0
  }
  obj2<-nleqslv(c(2.81),myfun2,method="Newton",global="qline",control=list(btol=0.00001,maxit=100)) # solution of the equation ARL-ARL0=0
  L<-obj2$x # the L value
######## IN Control
# verification, for the L value, the ARL0 should be 370.4
  mu0X<-alpha0/(alpha0+beta0) # IC mean of the beta distribution
  sigma0X2<-alpha0*beta0/(((alpha0+beta0)^2)*(alpha0+beta0+1)) # IC variance of the beta distribution
  LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
  UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
  delta<-(UCL-LCL)/(M)
  Q0<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerl<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperl<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q0[i,j]<-pbeta(upperl,shape1 = alpha0,shape2 = beta0)-pbeta(lowerl,shape1 = alpha0,shape2 = beta0)
    }
  }
  avec<-rep(0,M);avec[m+1]<-1
  l1<-rep(1,M)
  ID<-diag(M)
  ARLin<-avec%*%solve(ID-Q0)%*%l1
  print(c(lambda,L,as.vector(ARLin))) # here the output is lambda, L and the ARL0 value
# if everything is OK, the ARL0 should be 370.4  
}
############## OOCm shifts in mu0
for(mu1 in seq(0.12,0.28,by=0.02)){
phi1<-phi0
alpha1<-mu1*phi1 # in-control value of parameter shape1
beta1<-(1-mu1)*phi1 # in-control value of parameter shape2
#### we evaluate the OOC performance of the beta-EWMA chart for various shifts in process parameters
LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
mu1X<-alpha1/(alpha1+beta1) # OOC mean of the beta distribution
sigma1X2<-alpha1*beta1/(((alpha1+beta1)^2)*(alpha1+beta1+1)) # OOC variance of the beta distribution
##################
Q1u<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerlu<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperlu<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q1u[i,j]<-pbeta(upperlu,shape1 = alpha1,shape2 = beta1)-pbeta(lowerlu,shape1 = alpha1,shape2 = beta1)
    }
  }
  ARLout<-avec%*%solve(ID-Q1u)%*%l1
  Mout<-solve(ID-Q1u)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1u%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2)
  m1out<-ARLout
  m2out<-E2RLout
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1u%*%Q1u%*%l1
  m3out<-ERL3out+3*E2RLout-2*m1out
  MM1<-m1out
  MM2<-m2out-m1out^2
  MM3<-m3out-3*m2out*m1out+2*(m1out^3)
  A1<-4*(MM2^3)/(MM3^2)
  B1<-MM3/(2*MM2)
  threspar<-MM1-(2*MM2^2)/(MM3)
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # the MRL value
  RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # the 0.95-percentile value of the run-length distribution
#####
  cat("L:",L," a0:",alpha0," b0:",beta0," a1:",alpha1," b1:",beta1," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",mu0X," varX0:",sigma0X2," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
####### OOC, shifts in the dispersion parameters
for(d1 in seq(0.6,1.4,by=0.1)){
  mu1<-0.2
  phi1<-d1*phi0
  alpha1<-mu1*phi1 # in-control value of parameter shape1
  beta1<-(1-mu1)*phi1 # in-control value of parameter shape2
  #### we evaluate the OOC performance of the beta-EWMA chart for various shifts in process parameters
  LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
  UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
  delta<-(UCL-LCL)/(M)
  mu1X<-alpha1/(alpha1+beta1) # OOC mean of the beta distribution
  sigma1X2<-alpha1*beta1/(((alpha1+beta1)^2)*(alpha1+beta1+1)) # OOC variance of the beta distribution
#################
  Q1u<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerlu<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperlu<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q1u[i,j]<-pbeta(upperlu,shape1 = alpha1,shape2 = beta1)-pbeta(lowerlu,shape1 = alpha1,shape2 = beta1)
    }
  }
  ARLout<-avec%*%solve(ID-Q1u)%*%l1
  Mout<-solve(ID-Q1u)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1u%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2)
  m1out<-ARLout
  m2out<-E2RLout
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1u%*%Q1u%*%l1
  m3out<-ERL3out+3*E2RLout-2*m1out
  MM1<-m1out
  MM2<-m2out-m1out^2
  MM3<-m3out-3*m2out*m1out+2*(m1out^3)
  A1<-4*(MM2^3)/(MM3^2)
  B1<-MM3/(2*MM2)
  threspar<-MM1-(2*MM2^2)/(MM3)
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # the MRL value
  RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # the 0.95-percentile value of the run-length distribution
  #####
  cat("L:",L," a0:",alpha0," b0:",beta0," a1:",alpha1," b1:",beta1," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",mu0X," varX0:",sigma0X2," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
###########################################3
rm(list=ls(all=TRUE))
