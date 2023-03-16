library(gamlss.dist)
library(nleqslv)
############### parameters
phi0<-31
mu0<-0.2
mu0X<-mu0 # IC mean
sigma0X2<-mu0*(1-mu0)/(1+phi0) # IC variance
########
a0<-mu0*phi0 # shape1 parameter for beta distribution
b0<-(1-mu0)*phi0 # shape2 parameter for beta distribution
LCL<-qbeta(0.00135,shape1=a0,shape2=b0) # lower probability limit 
UCL<-qbeta(1-0.00135,shape1=a0,shape2=b0) # upper probability limit
# both limits are based on equal tail probability limits for ARL0=370.4
alpha<-1-pbeta(UCL,shape1=a0,shape2=b0)+pbeta(LCL,shape1=a0,shape2=b0)
alpha # False Alarm Rate (for verification)
##################### OOC shifts in mu0
for(mu1 in seq(0.12,0.28,by=0.02)){
    phi1<-phi0
    a1<-mu1*phi1 # out-of-control shape1 parameter
    b1<-(1-mu1)*phi1 # out-of-control shape2 parameter
    ##########################
beta1<-1-pbeta(UCL,shape1=a1,shape2=b1)+pbeta(LCL,shape1=a1,shape2=b1) # OOC probability
ARLout<-1/beta1 # OOC ARL
SDRLout<-sqrt(1-beta1)/beta1 # OOC SDRL
MRL<-ceiling(log(1-0.5)/log(1-beta1)) # OOC median run length
RL95<-ceiling(log(1-0.95)/log(1-beta1)) # OOC 95th-percentile point
    cat(" mu:",mu1," phi:",phi1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
        ,"\n")    
  }
##################### OOC shifts in phi0
# the code is almost the same as for the case of shifts in mu0
for(d1 in seq(0.6,1.4,by=0.1)){
mu1<-mu0
phi1<-d1*phi0
a1<-mu1*phi1
b1<-(1-mu1)*phi1
##########################
beta1<-1-pbeta(UCL,shape1=a1,shape2=b1)+pbeta(LCL,shape1=a1,shape2=b1)
ARLout<-1/beta1
SDRLout<-sqrt(1-beta1)/beta1
MRL<-ceiling(log(1-0.5)/log(1-beta1))
RL95<-ceiling(log(1-0.95)/log(1-beta1))
cat(" mu:",mu1," phi:",phi1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
        ,"\n")    
  }
##################
# remove everything
rm(list=ls(all=TRUE))

