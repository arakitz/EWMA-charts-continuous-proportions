library(gamlss.dist)
library(nleqslv)
#library(simplexreg)
#library(bazar)
#library(pracma)
##### IC parameters of the Simplex process
mu0<-0.2
S0<-0.5
############### OOC ARL values
# two-sided Shewhart chart
# shifts in mu0
for(mu1 in seq(0.12,0.28,by=0.02)){
S1<-S0
sigma0<-S0
sigma1<-S1
######### calculate the IC process parameters, Simplex case
funX<-function(q1){q1*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX<-integrate(funX,0,1)$value # mean
funX2<-function(q1){q1^2*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX2<-integrate(funX2,0,1)$value
VARX<-MUX2-MUX^2 # variance
############## calculate the OOC process parameters, Simplex case
funX1<-function(q1){q1*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX1<-integrate(funX1,0,1)$value # mean
funX21<-function(q1){q1^2*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX21<-integrate(funX21,0,1)$value
VARX1<-MUX21-MUX1^2 # variance
###########
# equal tail probability limits
# the IC ARL equals 370.4
LCL<-qSIMPLEX(0.00135,mu=mu0,sigma=S0)
UCL<-qSIMPLEX(1-0.00135,mu=mu0,sigma=S0)
alpha<-1-pSIMPLEX(UCL,mu=mu0,sigma=S0)+pSIMPLEX(LCL,mu=mu0,sigma=S0) # false alarm rate
#################
beta1<-1-pSIMPLEX(UCL,mu=mu1,sigma=S1)+pSIMPLEX(LCL,mu=mu1,sigma=S1) # OOC probability
ARLout<-1/beta1 # OOC ARL
SDRLout<-sqrt(1-beta1)/beta1 # OOC SDRL
MRL<-ceiling(log(1-0.5)/log(1-beta1)) # OOC MRL
RL95<-ceiling(log(1-0.95)/log(1-beta1)) # OOC 95-th percentile point
cat(" mu:",mu1," sigma:",S1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5),"\n")    
}
############### 
# below we calculate the OOC ARL values
# shifts in the dispersion parameter
# the code is exactly the same as for the case of shifts in mu0
for(d1 in seq(0.6,1.4,by=0.1)){
mu1<-mu0
S1<-d1*S0
sigma0<-S0
sigma1<-S1
######### IC process parameters
funX<-function(q1){q1*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX<-integrate(funX,0,1)$value
funX2<-function(q1){q1^2*dSIMPLEX(q1,mu=mu0,sigma=S0)}
MUX2<-integrate(funX2,0,1)$value
VARX<-MUX2-MUX^2
############## OOC process parameters
funX1<-function(q1){q1*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX1<-integrate(funX1,0,1)$value
funX21<-function(q1){q1^2*dSIMPLEX(q1,mu=mu1,sigma=S1)}
MUX21<-integrate(funX21,0,1)$value
VARX1<-MUX21-MUX1^2
########### control limits
LCL<-qSIMPLEX(0.00135,mu=mu0,sigma=S0)
UCL<-qSIMPLEX(1-0.00135,mu=mu0,sigma=S0)
alpha<-1-pSIMPLEX(UCL,mu=mu0,sigma=S0)+pSIMPLEX(LCL,mu=mu0,sigma=S0) # FAR
#################
beta1<-1-pSIMPLEX(UCL,mu=mu1,sigma=S1)+pSIMPLEX(LCL,mu=mu1,sigma=S1)
ARLout<-1/beta1
SDRLout<-sqrt(1-beta1)/beta1
MRL<-ceiling(log(1-0.5)/log(1-beta1))
RL95<-ceiling(log(1-0.95)/log(1-beta1))
cat(" mu:",mu0," d1:",d1," sigma:",S1," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95,
        " LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5),"\n")    
}
# remove everything
rm(list=ls(all=TRUE))