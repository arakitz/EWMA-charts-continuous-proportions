library(gamlss.dist)
library(nleqslv)
#library(simplexreg)
#library(bazar)
#library(pracma)
######## in-control parameters of the Simplex process ######
mu0<-0.2
S0<-0.50
############### parameters of the two-sided EWMA chart
lambda<-0.20 # smoothing parameter
L<-2.872233 # L
# evaluate the OOC ARL for shifts in mu0
for(mu1 in seq(0.12,0.28,by=0.02)){
S1<-S0
sigma0<-S0
sigma1<-S1
######### calculate the mean the variance of the simplex distribution, IC
funX<-function(q1){q1*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX<-integrate(funX,0,1)$value
funX2<-function(q1){q1^2*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX2<-integrate(funX2,0,1)$value
VARX<-MUX2-MUX^2
############## calculate the mean the variance of the simplex distribution, OOC
funX1<-function(q1){q1*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX1<-integrate(funX1,0,1)$value
funX21<-function(q1){q1^2*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX21<-integrate(funX21,0,1)$value
VARX1<-MUX21-MUX1^2
########### ARL calculation via a Markov chain method
m<-101 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
LCL<-MUX-L*sqrt(VARX*lambda/(2-lambda))
UCL<-MUX+L*sqrt(VARX*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
Q1<-matrix(0,ncol=M,nrow=M)
###############################
# fill the transition probabilities matrix
  for(i in 1:M){
    for(j in 1:M){
upperL<-LCL+(delta*(j-(1-lambda)*(i-0.5)))/lambda
lowerL<-LCL+(delta*(j-1-(1-lambda)*(i-0.5)))/lambda
if(upperL<0||lowerL<0){Q1[i,j]<-0}else{
if(upperL>0&upperL<1&lowerL<0){Q1[i,j]<-pSIMPLEX(upperL,mu=mu1,sigma=S1)}else{
if(upperL>0&upperL<1&lowerL<1&lowerL>0&upperL>lowerL){Q1[i,j]<-pSIMPLEX(upperL,mu=mu1,sigma=S1)-pSIMPLEX(lowerL,mu=mu1,sigma=S1)}
}
}
}
}
avec<-rep(0,M);avec[m+1]<-1 # initial probabilities vector
l1<-rep(1,M) # vector of 1s
ID<-diag(M) # identity matrix
ARLout<-avec%*%solve(ID-Q1)%*%l1 # matrix formula for ARL
Mout<-solve(ID-Q1)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2) # standard deviation of RL (SDRL)
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
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # the median run length (MRL) value
  RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # the 0.95-percentile value of the run-length distribution
#####
  cat("L:",L," mu0:",mu0," sigma0:",S0," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",MUX," varX0:",VARX," mu1:",MUX1," varX1:",VARX1,"\n")
}
############### calculate the OOC ARL for shifts in dispersion parameter
for(d1 in seq(0.6,1.4,by=0.1)){
mu1<-mu0
S1<-d1*S0
sigma0<-S0
sigma1<-S1
######### calculation IC mean and variance, simplex distribution
funX<-function(q1){q1*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX<-integrate(funX,0,1)$value
funX2<-function(q1){q1^2*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX2<-integrate(funX2,0,1)$value
VARX<-MUX2-MUX^2
############## calculation OOC mean and variance, simplex distribution
funX1<-function(q1){q1*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX1<-integrate(funX1,0,1)$value
funX21<-function(q1){q1^2*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX21<-integrate(funX21,0,1)$value
VARX1<-MUX21-MUX1^2
###########
# below we apply the Markov chain method for the calculation of OOC ARL
# the code is the same as previously (for shifts in mu0) but now the shifts 
# are in dispersion parameter
m<-101 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
LCL<-MUX-L*sqrt(VARX*lambda/(2-lambda))
UCL<-MUX+L*sqrt(VARX*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
Q1<-matrix(0,ncol=M,nrow=M)
###############################
  for(i in 1:M){
    for(j in 1:M){
upperL<-LCL+(delta*(j-(1-lambda)*(i-0.5)))/lambda
lowerL<-LCL+(delta*(j-1-(1-lambda)*(i-0.5)))/lambda
if(upperL<0||lowerL<0){Q1[i,j]<-0}else{
if(upperL>0&upperL<1&lowerL<0){Q1[i,j]<-pSIMPLEX(upperL,mu=mu1,sigma=S1)}else{
if(upperL>0&upperL<1&lowerL<1&lowerL>0&upperL>lowerL){Q1[i,j]<-pSIMPLEX(upperL,mu=mu1,sigma=S1)-pSIMPLEX(lowerL,mu=mu1,sigma=S1)}
}
}
}
}
avec<-rep(0,M);avec[m+1]<-1
l1<-rep(1,M)
ID<-diag(M)
ARLout<-avec%*%solve(ID-Q1)%*%l1
Mout<-solve(ID-Q1)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2) # SDRL
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
  cat("L:",L," mu0:",mu0," sigma0:",S0," d1:",d1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu0:",MUX," varX0:",VARX," mu1:",MUX1," varX1:",VARX1,"\n")
}
#############################
rm(list=ls(all=TRUE))
