library(gamlss.dist)
library(nleqslv)
##### IC values of Unit-Gamma process
mu0<-0.2
tau0<-20
alpha0<-tau0 # in-control value of shape parameter 
beta0<-(mu0^(1/alpha0))/(1-(mu0^(1/alpha0))) # in-control value of rate parameter
# equal tail probability limits
# the FAR is 0.0027
# the IC ARL for the two-sided Shewhart chart is 370.4
UCL<-exp(-qgamma(0.00135,shape = alpha0,rate = beta0))
LCL<-exp(-qgamma(1-0.00135,shape = alpha0,rate = beta0))
### calculate OOC ARL for shifts in mu0
for(mu1 in seq(0.12,0.28,by=0.02)){
tau1<-tau0 
alpha1<-tau1 # out-of-control value of parameter shape
beta1<-(mu1^(1/alpha1))/(1-(mu1^(1/alpha1))) # out-of-control value of parameter rate
# below is the OOC probability
pOOC<-1-((1-pgamma(-log(UCL),shape = alpha1,rate = beta1)))+((1-pgamma(-log(LCL),shape = alpha1,rate = beta1)))
ARLout<-1/pOOC # out-of-control ARL
SDRLout<-sqrt(1-pOOC)/pOOC # out-of-control SDRL
MRL<-ceiling(log(1-0.5)/log(1-pOOC)) # out-of-control MRL
RL95<-ceiling(log(1-0.95)/log(1-pOOC)) # out-of-control 95th percentile point
cat(" mu:",mu1," tau:",tau1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5),"\n")   
}
### calculate OOC ARL for shifts in dispersion parameter
# the code is the same as for the case of shifts in mu0
for(d1 in seq(0.6,1.4,by=0.1)){
mu1<-mu0
tau1<-d1*tau0
alpha1<-tau1 # out-of-control value of parameter shape
beta1<-(mu1^(1/alpha1))/(1-(mu1^(1/alpha1))) # out-of-control-control value of parameter rate
pOOC<-1-((1-pgamma(-log(UCL),shape = alpha1,rate = beta1)))+((1-pgamma(-log(LCL),shape = alpha1,rate = beta1)))
ARLout<-1/pOOC
SDRLout<-sqrt(1-pOOC)/pOOC
MRL<-ceiling(log(1-0.5)/log(1-pOOC))
RL95<-ceiling(log(1-0.95)/log(1-pOOC))
cat(" mu:",mu1," tau:",tau1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5),"\n")
}
## remove everything
rm(list=ls(all=TRUE))