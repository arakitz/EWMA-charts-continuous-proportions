library(gamlss.dist)
library(nleqslv)
############### parameters of the uGA process
mu0<-0.2
tau0<-20
alpha0<-tau0 # in-control value of shape parameter
beta0<-(mu0^(1/alpha0))/(1-(mu0^(1/alpha0))) # in-control value of the rate parameter
##################
lambda<-0.20 # lambda value
# define the ARLin as a function of L (here Ly)
ARLiny<-function(Ly){
###################### Design state, Markov chain method
  m<-200 # states for the Markov chain
  M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
  mu0X<-mu0 # IC mean of the beta distribution
  sigma0X2<-mu0*(-mu0+1/(2-mu0^(1/tau0))^tau0) # IC variance of the uGA distribution
  LCLy<-max(0,mu0X-Ly*sqrt(sigma0X2*lambda/(2-lambda))) # lower control limit
  UCLy<-mu0X+Ly*sqrt(sigma0X2*lambda/(2-lambda)) # upper control limit
  deltay<-(UCLy-LCLy)/(M)
  # fill the transition probability matrix
  Q0<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerl<-max(0,LCLy+deltay*(j-1-(1-lambda)*(i-0.5))/lambda)
      upperl<-max(0,LCLy+deltay*(j-(1-lambda)*(i-0.5))/lambda)
      Q0[i,j]<-(1-pgamma(-log(upperl),shape = alpha0,rate = beta0))-(1-pgamma(-log(lowerl),shape = alpha0,rate = beta0))
    }
  }
  avec<-rep(0,M);avec[m+1]<-1 # initial probabilities vector
  l1<-rep(1,M) # vector of 1s
  ID<-diag(M) # identity matrix
  ARLin<-avec%*%solve(ID-Q0)%*%l1 # matrix formula for the ARL 
  # define the equation that must be solved to determine L
  as.vector(ARLin)-370.4
}
# use the nleqslv function for non-linear equation solving
Obj1<-nleqslv(c(2.7),ARLiny,method="Newton",global="qline",control=list(btol=0.00001,maxit=100))
L<-Obj1$x # L value
cat(" L:",L,'\n')
################### OOC performance, shifts in mu0

#alpha0<-tau0 # in-control value of parameter shape
#beta0<-(mu0^(1/alpha0))/(1-(mu0^(1/alpha0))) # in-control value of parameter rate
for(mu1 in seq(0.12,0.28,by=0.02)){ # calculate the ARLout for various mu1 values
tau1<-tau0
alpha1<-tau1 # out-of-control value of shape parameter
beta1<-(mu1^(1/alpha1))/(1-(mu1^(1/alpha1))) # out-of-control value of rate parameter
###################### Markov chain method
m<-200 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
mu0X<-mu0 # IC mean of the beta distribution
sigma0X2<-mu0*(-mu0+1/(2-mu0^(1/tau0))^tau0) # IC variance of the uGA distribution
mu1X<-mu1 # OOC mean of the uGA distribution
sigma1X2<-mu1*(-mu1+1/(2-mu1^(1/tau1))^tau1) # OOC variance of the uGA distribution
# control limits
LCL<-max(0,mu0X-L*sqrt(sigma0X2*lambda/(2-lambda)))
UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
# Transition Probabilities Matrix
Q1<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerl<-max(0,LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda)
      upperl<-max(0,LCL+delta*(j-(1-lambda)*(i-0.5))/lambda)
      Q1[i,j]<-(1-pgamma(-log(upperl),shape = alpha1,rate = beta1))-(1-pgamma(-log(lowerl),shape = alpha1,rate = beta1))
    }
  }
  avec<-rep(0,M);avec[m+1]<-1 # initial probabilities vector
  l1<-rep(1,M) # vector of 1s
  ID<-diag(M) # identity matrix
  ARLout<-avec%*%solve(ID-Q1)%*%l1 # matrix formula for ARL1
  # calculation of the standard deviation run length, the median run length (MRL) and the 95-th percentile point
  Mout<-solve(ID-Q1) 
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2) #the SDRL
  m1out<-ARLout
  m2out<-E2RLout
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1%*%Q1%*%l1
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
  cat("L:",L," a0:",alpha0," b0:",beta0," a1:",alpha1," b1:",beta1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",mu0X," varX0:",sigma0X2," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
######################## shifts in tau0
# the code is the same as for the case of shifts in mu0
# the OOC ARL values are calculated for shifts in the dispersion parameter tau0

#alpha0<-tau0 # in-control value of parameter shape
#beta0<-(mu0^(1/alpha0))/(1-(mu0^(1/alpha0))) # in-control value of parameter rate
# shifts in tau0
for(d1 in seq(0.6,1.4,by=0.1)){
mu1<-mu0
tau1<-d1*tau0
alpha1<-tau1 # out-of-control value of shape parameter
beta1<-(mu1^(1/alpha1))/(1-(mu1^(1/alpha1))) # out-of-control value of rate parameter
###################### Markov chain method
m<-200 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
mu0X<-mu0 # IC mean of the uGA distribution
  sigma0X2<-mu0*(-mu0+1/(2-mu0^(1/tau0))^tau0) # IC variance of the uGA distribution
mu1X<-mu1 # OOC mean of the uGA distribution
  sigma1X2<-mu1*(-mu1+1/(2-mu1^(1/tau1))^tau1) # OOC variance of the uGA distribution
  LCL<-max(0,mu0X-L*sqrt(sigma0X2*lambda/(2-lambda)))
  UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
  delta<-(UCL-LCL)/(M)
  Q1<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerl<-max(0,LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda)
      upperl<-max(0,LCL+delta*(j-(1-lambda)*(i-0.5))/lambda)
      Q1[i,j]<-(1-pgamma(-log(upperl),shape = alpha1,rate = beta1))-(1-pgamma(-log(lowerl),shape = alpha1,rate = beta1))
    }
  }
  avec<-rep(0,M);avec[m+1]<-1
  l1<-rep(1,M)
  ID<-diag(M)
  ARLout<-avec%*%solve(ID-Q1)%*%l1
  Mout<-solve(ID-Q1)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2) # the SDRL
  m1out<-ARLout
  m2out<-E2RLout
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1%*%Q1%*%l1
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
  cat("L:",L," a0:",alpha0," b0:",beta0," a1:",alpha1," b1:",beta1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",mu0X," varX0:",sigma0X2," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
#####################
rm(list=ls(all=TRUE))

