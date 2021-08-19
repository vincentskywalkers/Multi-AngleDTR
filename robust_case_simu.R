#setwd("/Users/mac/Desktop/newnewfunc/Latest")
#setwd("/Users/mac/Desktop/simunewkm")
library(personalized)
library(survival)
#source("all_func13.R")
#source("all_funcnew13.R")

Run = 150;
N = 500;
M = 5;
k=3
dimco = 25;
givent = 1.4;
mt = c(0,0.5,1,1.5,2)
censrate=NULL
truesurv=NULL
estsurv=NULL
estsurv_q=NULL
estsurv_sub=NULL
MCR=matrix(0,Run,1)
MCR_q=matrix(0,Run,1)
MCR_sub=matrix(0,Run,1)
V=XI.gen(k)


beta_set1 = c(-1,0.5,-1,rep(0,dimco-3))
beta_set2 = c(0.5,-1,-1,rep(0,dimco-3))
beta1 = c(beta_set1,beta_set2)
beta2 = c(-0.5*beta_set1,1*beta_set1,-0.3,-0.5*beta_set2,1*beta_set2,-0.3)
beta3 = c(0.25*beta_set1,-0.25*beta_set1,1.5*beta_set1,-0.05,-0.1,0.25*beta_set2,-0.25*beta_set2,1.5*beta_set2,-0.05,-0.1)
beta4 = c(0.25*beta_set1,-0.25*beta_set1,-0.25*beta_set1,1*beta_set1,0.05,0.05,-0.15,0.25*beta_set2,-0.25*beta_set2,-0.25*beta_set2,1*beta_set2,0.05,0.05,-0.15)
beta5 = c(0.05*beta_set1,-0.1*beta_set1,-0.25*beta_set1,-0.5*beta_set1,1*beta_set1,0.05,-0.05,0.05,-0.1,0.05*beta_set2,-0.1*beta_set2,-0.25*beta_set2,-0.5*beta_set2,1*beta_set2,0.05,-0.05,0.05,-0.1)
#####

for (rr in 1:Run) {
  #######################---DATA---############################
  set.seed(((rr+121)*555)+123)
  N=500
  #C# = exp(-1*X[,1,1]+0.5*X[,1,8]+rnorm(N,0,1)-1.5)
  #C# = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9)
  #C# = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.5)
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])
  ####################
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  }
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  
  #########################----proposed method---###############################
  
  lambda = 0.0001
  a=1; b=1;
  grid=length(S)
  l0_qt =matrix(0,grid,1);
  for (ss in 1:grid) {
    l0_qt[ss]=1-sum((Y==S[ss])*delta)/sum((Y>=S[ss]))
  }
  f_norm = -sum(log(l0_qt))/lambda
  C_x = max(max(X[,1:Nstep,]),max(A[,1:Nstep]))
  C_H = C_x*sqrt(dimco*Nstep+Nstep-1)
  c = -(C_H*sqrt(f_norm))
  
  betahat = optim(beta,Qfunc,lambda=lambda,X=X,V=V,A=A,a=a,b=b,c=c,
                  phat=phat,TT=S,Y=Y,delta=delta,dimco=dimco,mt=mt,indexstep=indexstep,method="BFGS",
                  control = list(maxit = 5, trace = TRUE, REPORT = 500))$par
  
  set.seed(((rr+121)*555)+321)
  N=2000
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  # Tf # = exp(0.4*apply((A==d),1,sum)-2*abs(X[,1,1])+X[,2,2])# # #
  # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2]) #
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  ################################
  
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  }
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  KMtc=KMcurve_wtc(grid,S,Y,delta,A,d,phat,mt,indexstep)
  if (sum(KMtc == "NaN") > 60){
    truesurv=c(truesurv,NA)
  } else {
    truesurv=c(truesurv,min(KMtc,na.rm = T))
  }
  if(rr==1){
    plot(c(0,S),KMtc,main="KM",type = "l",xlab = "t",ylim=c(0,1),ylab = "Survival Proportion",col="red")
  }else{
    points(c(0,S),KMtc,type = "l",col="red")
  }
  
  fhat1=decifunc(phat,X,A,betahat,dimco)
  fhat2=decifunc2(phat,X,A,betahat,dimco)
  
  dhat = matrix(0,N,Nstep)
  for (tt in 1:Nstep){
    inner.matrix=matrix(0,N,3)
    for (ii in 1:k){
      inner.matrix[,ii] = apply(V[,ii]*cbind(fhat1[,tt],fhat2[,tt]),1,sum)
    }
    dhat[,tt]=apply(inner.matrix,1,pred)
  }
  
  
  KMtce=KMcurve_wtc(grid,S,Y,delta,A,dhat,phat,mt,indexstep)
  if (sum(KMtce == "NaN") > 60){
    estsurv=c(estsurv,NA)
  } else {
    estsurv=c(estsurv,min(KMtce,na.rm = T))
  }
  points(c(0,S),KMtce,type = "l",col="blue")
  MCR[rr] = mean((d[,1:Nstep]!=dhat))
  
  ########################---subgroup identification--######################
  
  dhat_sub = matrix(0,N,Nstep)
  
  
  for (tt in 1:Nstep) {
    
    #
    set.seed(((rr+121)*555)+123)
    N=500
    X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
    
    A=f=d=matrix(0,N,M)
    
    Pr=matrix(1/3,N,M)
    
    for (i in 1:N) {
      A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    
    f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
    f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
    
    d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
    d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
    f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
    f_21 = (f_21)^3
    f_22 = (f_22)^3
    
    d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
    d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    
    for (i in 1:N) {
      A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
    f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
    f_31 = (f_31)^3
    f_32 = (f_32)^3
    
    d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
    d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
    f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
    f_41 = (f_41)^3
    f_42 = (f_42)^3
    
    d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
    d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
    f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
    f_51 = (f_51)^3
    f_52 = (f_52)^3
    
    d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
    d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
    # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
    
    Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
    
    Y = apply(cbind(Tf,C),1,min)
    delta = (Tf<=C)
    censrate = c(censrate,mean(1-delta))
    ################################
    indexstep = c()
    for (j in 1:N){
      realsurv = min(Y[j],givent)
      if(realsurv>=mt[1]&realsurv<mt[2]){
        indexstep[j] = 1
      }
      if(realsurv>=mt[2]&realsurv<mt[3]){
        indexstep[j] = 2
      }
      if(realsurv>=mt[3]&realsurv<mt[4]){
        indexstep[j] = 3
      }
      if(realsurv>=mt[4]&realsurv<mt[5]){
        indexstep[j] = 4
      }
      if(realsurv>=mt[5]){
        indexstep[j] = 5
      }
    }
    
    #############################---get true survival curve---##########################
    
    S=timegrid(Y,delta,givent)
    grid=length(S)
    if(givent>=mt[1]&givent<mt[2]){
      phat=matrix(1/3,N,1)
      beta = beta1
    }
    if(givent>=mt[2]&givent<mt[3]){
      phat=matrix(1/3,N,2)
      beta = c(beta1,beta2)
    }
    if(givent>=mt[3]&givent<mt[4]){
      phat=matrix(1/3,N,3)
      beta = rep(0,9*dimco+4+3*dimco+2)
      #beta = c(beta1,beta2,beta3)
    }
    if(givent>=mt[4]&givent<mt[5]){
      phat=matrix(1/3,N,4)
      beta =c(beta1,beta2,beta3)
    }
    if(givent>=mt[5]){
      phat=matrix(1/3,N,5)
      beta = rep(0,25*dimco+16+5*dimco+4)
    }
    phat=as.matrix(phat)
    Nstep=length(phat[1,])
    
    A1  = A[,tt]
    A1[A1==3] = 0
    A1[A1==2] = 0
    
    A2   = A[,tt]
    A2[A2==1] = 0
    A2[A2==3] = 0
    A2[A2==2] = 1
    
    
    A3   = A[,tt]
    A3[A3==1] = 0
    A3[A3==2] = 0
    A3[A3==3] = 1
    
    n = N
    nstage = Nstep
    y =  matrix(0, nrow=n, ncol=nstage)
    u =  matrix(0, nrow=n, ncol=nstage)
    delta_mat = matrix(0, nrow=n, ncol=nstage)
    
    ##--------------------calculate interval survival outcome---------------------
    y[,1] = exp(0.4*((A[,1]==d[,1]))-3.5*abs(X[,1,3])+X[,2,2])
    
    y[,2] =  exp(0.4*apply((A[,1:2]==d[,1:2]),1,sum)
                 -3.5*abs(X[,1,3])+X[,2,2])-  y[,1]
    
    for (i in 3:nstage){
      y[,i] =  exp(0.4*apply((A[,1:i]==d[,1:i]),1,sum)
                   -3.5*abs(X[,1,3])+X[,2,2])- exp(0.4*apply((A[,1:(i-1)]==d[,1:(i-1)]),1,sum)
                                                   -3.5*abs(X[,1,3])+X[,2,2])
    }
    
    
    ##------------calculate min cumulated survival outcome and censored time-----------
    u[,1]    = pmin(y[,1], C)
    if (nstage == 2){
      u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1]
    }else if (nstage > 2){
      for (j in 3:nstage){
        u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
      }
    }
    for (j in 1:nstage){
      u[,j]    = u[,j] + 0.01
    }
    
    ##-----------------------censored indicator for each interval------------------------
    delta_mat[,1] = y[,1] < C
    if (nstage > 1){
      for (kkk in 2:nstage){
        delta_mat[,kkk] = (y[,1] < C)
      }
    }
    
    model0=fit.subgroup(X[,tt,], y = Surv(u[,tt], delta_mat[,tt]), as.vector(A1), propensity.func = prop.func, loss   = "cox_loss_lasso", method="a_learning",
                        nfolds = 5)
    model1=fit.subgroup(X[,tt,], y = Surv(u[,tt], delta_mat[,tt]), as.vector(A2), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
                        nfolds = 5)
    model2=fit.subgroup(X[,tt,], y = Surv(u[,tt], delta_mat[,tt]), as.vector(A3), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
                        nfolds = 5)
    
    
    set.seed(((rr+121)*555)+321)
    N=2000
    X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
    
    A=f=d=matrix(0,N,M)
    
    Pr=matrix(1/3,N,M)
    
    for (i in 1:N) {
      A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    
    f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
    f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
    
    d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
    d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
    f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
    f_21 = (f_21)^3
    f_22 = (f_22)^3
    
    d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
    d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    
    for (i in 1:N) {
      A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
    f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
    f_31 = (f_31)^3
    f_32 = (f_32)^3
    
    d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
    d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
    f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
    f_41 = (f_41)^3
    f_42 = (f_42)^3
    
    d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
    d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    for (i in 1:N) {
      A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
    }
    f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
    f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
    f_51 = (f_51)^3
    f_52 = (f_52)^3
    
    d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
    d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
    
    # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
    # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
    
    Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
    
    Y = apply(cbind(Tf,C),1,min)
    delta = (Tf<=C)
    censrate = c(censrate,mean(1-delta))
    ################################
    indexstep = c()
    for (j in 1:N){
      realsurv = min(Y[j],givent)
      if(realsurv>=mt[1]&realsurv<mt[2]){
        indexstep[j] = 1
      }
      if(realsurv>=mt[2]&realsurv<mt[3]){
        indexstep[j] = 2
      }
      if(realsurv>=mt[3]&realsurv<mt[4]){
        indexstep[j] = 3
      }
      if(realsurv>=mt[4]&realsurv<mt[5]){
        indexstep[j] = 4
      }
      if(realsurv>=mt[5]){
        indexstep[j] = 5
      }
    }
    #############################---get true survival curve---##########################
    
    S=timegrid(Y,delta,givent)
    grid=length(S)
    if(givent>=mt[1]&givent<mt[2]){
      phat=matrix(1/3,N,1)
      beta = beta1
    }
    if(givent>=mt[2]&givent<mt[3]){
      phat=matrix(1/3,N,2)
      beta = c(beta1,beta2)
    }
    if(givent>=mt[3]&givent<mt[4]){
      phat=matrix(1/3,N,3)
      beta = rep(0,9*dimco+4+3*dimco+2)
      #beta = c(beta1,beta2,beta3)
    }
    if(givent>=mt[4]&givent<mt[5]){
      phat=matrix(1/3,N,4)
      beta =c(beta1,beta2,beta3)
    }
    if(givent>=mt[5]){
      phat=matrix(1/3,N,5)
      beta = rep(0,25*dimco+16+5*dimco+4)
    }
    phat=as.matrix(phat)
    Nstep=length(phat[1,])
    
    
    a1=  predict(model0, newx = X[,tt,], type = "benefit.score")
    
    a2 = predict(model1, newx = X[,tt,], type = "benefit.score")
    
    a3 = predict(model2, newx = X[,tt,], type = "benefit.score")
    
    a = cbind(a1,a2,a3)
    
    dhat_sub[,tt] =  one_vs_all(a)
  }
  
  #
  set.seed(((rr+121)*555)+321)
  N=2000
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  ################################
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  }
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  
  KMtsub=KMcurve_wtc(grid,S,Y,delta,A,dhat_sub,phat,mt,indexstep)
  if (sum(KMtsub == "NaN") > 60){
    estsurv_sub=c(estsurv_sub,NA)
  } else {
    estsurv_sub=c(estsurv_sub,min(KMtsub,na.rm = T))
  }
  points(c(0,S),KMtsub,type = "l",col="orange")
  MCR_sub[rr] = mean((d[,1:Nstep]!=dhat_sub))
  
  mean(censrate)
  ########################---censored q learning---######################
  
  set.seed(((rr+121)*555)+123)
  N=500
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  ################################
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  } 
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  
  n = N
  nstage = Nstep
  y =  matrix(0, nrow=n, ncol=nstage)
  u =  matrix(0, nrow=n, ncol=nstage)
  delta_mat = matrix(0, nrow=n, ncol=nstage)
  d_qhat= matrix(0, nrow=2000, ncol=nstage)
  weight_mat = matrix(0, nrow = n, ncol = nstage)
  
  ##--------------------calculate interval survival outcome---------------------
  y[,1] = exp(0.4*((A[,1]==d[,1]))-3.5*abs(X[,1,3])+X[,2,2]) 
  
  y[,2] =  exp(0.4*apply((A[,1:2]==d[,1:2]),1,sum)
               -3.5*abs(X[,1,3])+X[,2,2])-  y[,1] 
  
  for (i in 3:nstage){
    y[,i] =  exp(0.4*apply((A[,1:i]==d[,1:i]),1,sum)
                 -3.5*abs(X[,1,3])+X[,2,2])- exp(0.4*apply((A[,1:(i-1)]==d[,1:(i-1)]),1,sum)
                                                 -3.5*abs(X[,1,3])+X[,2,2])
  }
  
  
  ##------------calculate min cumulated survival outcome and censored time-----------
  u[,1]    = pmin(y[,1], C)
  u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1]
  for (j in 3:nstage){
    u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
  }
  for (j in 1:nstage){
    u[,j]    = u[,j]
  }
  ##-----------------------censored indicator for each interval------------------------
  delta_mat[,1] = y[,1] < C
  if (nstage > 1){
    for (kkk in 2:nstage){
      # delta_mat[,kkk] = (apply(y[,1:kkk], 1, sum) < C)
      delta_mat[,kkk] = (apply(y[,1:kkk], 1, sum) < C)
    }
  }
  
  ##-----------------------generated histroy information matrix-----------------------
  
  
  d_qhat_pesudo= matrix(0, nrow=2000, ncol=nstage)
  
  hinfor   = list(nstage)
  for (j in 1:nstage){ 
    res <- model.matrix(~as.factor(A[,j]))
    hinfor[[j]] = cbind(1,X[,j,],res[,2:3])
  }
  ##---------------------caluculate conditional survival function--------------------
  
  for (i in 1:nstage){
    basisinfo =   (cbind(X[,i,]))
    if (i == 1){
      u_Vec  =  u[ , i]
    }
    if (i > 1){
      u_Vec  =  apply(u[, 1:i], 1, sum)
    }
    S_c_m  =  cal_survival_and_density_functions(u_Vec, basisinfo, (1 - delta_mat[, i]),u_Vec, basisinfo, fitting_method = "Cox")$S
    weight_mat[, i] =  delta_mat[ , i]/pmax(S_c_m, 0.05)
  }
  
  
  ######################################################################################################
  qlearningcoef = qlearning.censored(hinfor, u, delta_mat, weight_mat, nstage)
  ######################################################################################################
  
  set.seed(((rr+121)*555)+321)
  N=2000
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  ################################
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  } 
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  
  
  hinfor   = list(nstage)
  for (j in 1:nstage){ 
    res <- model.matrix(~as.factor(A[,j]))
    hinfor[[j]] = cbind(1,X[,j,],res[,2:3])
  }
  
  p0_tmp      =  ncol(hinfor[[1]])
  p_tmp       =  length(qlearningcoef[[1]])
  
  for (stg in 1:nstage){
    for (j in 1:nrow(hinfor[[1]])){
      
      tmpf_1 = hinfor[[stg]][j,1:(p0_tmp-2)]%*%qlearningcoef[[stg]][1:(p0_tmp-2)] 
      
      tmpf_2 = hinfor[[stg]][j,1:(p0_tmp-2)]%*%qlearningcoef[[stg]][1:(p0_tmp-2)] +
        hinfor[[stg]][j,p0_tmp-1]*qlearningcoef[[stg]][(p0_tmp-1)] +
        c((hinfor[[stg]][j, 2:(dim(hinfor[[stg]])[2]-2)])*1,
          (hinfor[[stg]][j, 2:(dim(hinfor[[stg]])[2]-2)])*0)%*%
        (qlearningcoef[[stg]][(p0_tmp+1):p_tmp])
      
      
      tmpf_3 = hinfor[[stg]][j,1:(p0_tmp-2)]%*%qlearningcoef[[stg]][1:(p0_tmp-2)] +
        hinfor[[stg]][j,p0_tmp]*qlearningcoef[[stg]][(p0_tmp)] +
        c((hinfor[[stg]][j, 2:(dim(hinfor[[stg]])[2]-2)])*0,
          (hinfor[[stg]][j, 2:(dim(hinfor[[stg]])[2]-2)])*1)%*%
        (qlearningcoef[[stg]][(p0_tmp+1):p_tmp])
      
      tmpf_all = c(tmpf_1,tmpf_2,tmpf_3)
      #print(tmpf_all)
      
      d_qhat[j,stg] = which.max(tmpf_all)
    }
  }
  
  
  #
  set.seed(((rr+121)*555)+321)
  N=2000
  X= array(runif(N*M*dimco,0,1),dim=c(N,M,dimco)); C = exp(0.5*abs(X[,1,2])+1*X[,1,1]+rnorm(N,0,1)-0.9); X[,1,2] = runif(N,0,1) ; X[,1,1] = runif(N,0,1)
  
  A=f=d=matrix(0,N,M)
  
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]; f_11 = (f_11)^3
  f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]; f_12 = (f_12)^3
  
  d[,1] = 1+apply(cbind(0,matrix(sign(f_11),N,1)),1,max) + apply(cbind(0,matrix(sign(f_12),N,1)),1,max)
  d[sample(1:N,0.05*N),1] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  f_21 = (f_21)^3
  f_22 = (f_22)^3
  
  d[,2] = 1+apply(cbind(0,matrix(sign(f_21),N,1)),1,max) + apply(cbind(0,matrix(sign(f_22),N,1)),1,max)
  d[sample(1:N,0.05*N),2] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  f_31 = (f_31)^3
  f_32 = (f_32)^3
  
  d[,3] = 1+apply(cbind(0,matrix(sign(f_31),N,1)),1,max) + apply(cbind(0,matrix(sign(f_32),N,1)),1,max)
  d[sample(1:N,0.05*N),3] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  f_41 = (f_41)^3
  f_42 = (f_42)^3
  
  d[,4] = 1+apply(cbind(0,matrix(sign(f_41),N,1)),1,max) + apply(cbind(0,matrix(sign(f_42),N,1)),1,max)
  d[sample(1:N,0.05*N),4] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(1,2,3),1,replace = T,prob = c(1/3,1/3,1/3))
  }
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  f_51 = (f_51)^3
  f_52 = (f_52)^3
  
  d[,5] = 1+apply(cbind(0,matrix(sign(f_51),N,1)),1,max) + apply(cbind(0,matrix(sign(f_52),N,1)),1,max)
  d[sample(1:N,0.05*N),5] = sample(x=c(1,2,3),0.05*N,replace = T,prob = c(1/3,1/3,1/3))
  
  # Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  # C = exp(0.5*abs(X[,1,2])+exp(X[,1,1])+rnorm(N,0,1))
  
  Tf = exp(0.4*apply((A==d),1,sum)-3.5*abs(X[,1,3])+X[,2,2])# # #
  
  Y = apply(cbind(Tf,C),1,min)
  delta = (Tf<=C)
  censrate = c(censrate,mean(1-delta))
  
  ################################
  indexstep = c()
  for (j in 1:N){
    realsurv = min(Y[j],givent)
    if(realsurv>=mt[1]&realsurv<mt[2]){
      indexstep[j] = 1
    }
    if(realsurv>=mt[2]&realsurv<mt[3]){
      indexstep[j] = 2
    }
    if(realsurv>=mt[3]&realsurv<mt[4]){
      indexstep[j] = 3
    }
    if(realsurv>=mt[4]&realsurv<mt[5]){
      indexstep[j] = 4
    }
    if(realsurv>=mt[5]){
      indexstep[j] = 5
    }
  }
  
  #############################---get true survival curve---##########################
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  if(givent>=mt[1]&givent<mt[2]){
    phat=matrix(1/3,N,1)
    beta = beta1
  }
  if(givent>=mt[2]&givent<mt[3]){
    phat=matrix(1/3,N,2)
    beta = c(beta1,beta2)
  }
  if(givent>=mt[3]&givent<mt[4]){
    phat=matrix(1/3,N,3)
    beta = rep(0,9*dimco+4+3*dimco+2)
    #beta = c(beta1,beta2,beta3)
  } 
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta =c(beta1,beta2,beta3)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = rep(0,25*dimco+16+5*dimco+4)
  }
  phat=as.matrix(phat)
  Nstep=length(phat[1,])
  
  KMtq=KMcurve_wtc(grid,S,Y,delta,A,d_qhat,phat,mt,indexstep)
  if (sum(KMtq == "NaN") > 60){
    estsurv_q=c(estsurv_q, NA)
  } else {
    estsurv_q=c(estsurv_q,min(KMtq,na.rm = T))
  }
  points(c(0,S),KMtq,type = "l",col="green4")
  MCR_q[rr] = mean((d[,1:Nstep]!=d_qhat))
  print(rr)
  
}


#####################################################################################
print(mean(MCR))
sd(MCR)
print(mean(MCR_q))
sd(MCR_q)
print(mean(MCR_sub))
sd(MCR_sub)
#####################################################################################
print(mean(truesurv,na.rm=T))
print(sd(truesurv,na.rm = T))
print(mean(estsurv,na.rm=T))
print(sd(estsurv,na.rm = T))
print(mean(estsurv_q,na.rm=T))
print(sd(estsurv_q,na.rm = T))
print(mean(estsurv_sub,na.rm=T))
print(sd(estsurv_sub,na.rm = T))
#####################################################################################

