library(personalized)
library(survival)
#source("utility_func.R")
mean(censrate)
Run =100
N = 500;
M = 5;
k=3
dimco = 25;
givent = 1.5;
#mt = c(0,0.5,1,1.5,2)
mt = c(0.5,1,1.5,2,2.5)
censrate=NULL
truesurv=NULL
estsurv=NULL
estsurv_q=NULL
estsurv_sub=NULL
MCR=matrix(0,Run,1)
MCR_q=matrix(0,Run,1)
MCR_sub=matrix(0,Run,1)
V=XI.gen(k)
beta_set1 = c(1,1,1,rep(0,dimco-3))
beta_set2 = c(1,-1,-1,rep(0,dimco-3))
beta1 = c(beta_set1,beta_set2) 
beta2 = c(0.5*beta_set1,-1.5*beta_set1,0.1,0.5*beta_set2,-1.5*beta_set2,0.1)
beta3 = c(0.25*beta_set1,-0.5*beta_set1,1*beta_set1,0.05,0.1,0.25*beta_set2,-0.5*beta_set2,1*beta_set2,0.05,-0.1)
beta4 = c(0.1*beta_set1,-0.25*beta_set1,0.5*beta_set1,1*beta_set1,0.05,-0.05,0.1,0.1*beta_set2,-0.25*beta_set2,0.5*beta_set2,1*beta_set2,0.05,-0.05,0.1)
beta5 = c(0.05*beta_set1,-0.1*beta_set1,0.25*beta_set1,0.5*beta_set1,1*beta_set1,0.05,-0.05,0.05,-0.1,0.05*beta_set2,-0.1*beta_set2,0.25*beta_set2,0.5*beta_set2,1*beta_set2,0.05,-0.05,0.05,-0.1)

probeta1 = c(0,0,0.25,rep(0,dimco-3))
probeta2 = c(0,-0.25,0.5,rep(0,dimco-3))

for (rr in 1:Run) {
  set.seed((rr+203)*44+123)
  N = 500;
  #######################---DATA---############################
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X; 
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  #############################---get true survival curve---##########################
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N,train)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  #########################----proposed method---###############################
  #lambda = 0.001   givent = 1.5
  #lambda = 0.0001  givent = 2.5
  
  # ii = 0
  # lambdaset = c(0.00001,0.0001,0.001)
  # bset = c(1,5,10)
  # lamvalue = c()
  # for (lambdatmp in lambdaset){
  #   for (btmp in bset){
  #     ii = ii+1
  #     lamvalue[ii] = crossvalidate(train,2,lambdatmp,btmp)
  #   }
  # }
  # b = bset[(which.max(lamvalue) - length(bset)*floor((which.max(lamvalue)-1)/length(bset)))]
  # lambda = lambdaset[ceiling(which.max(lamvalue)/length(bset))]
  
  
  lambda = 0.00001
  a=1; b=1;
  grid=length(S)
  l0_qt =matrix(0,grid,1);
  for (ss in 1:grid) {
    l0_qt[ss]=1-sum((Y==S[ss])*delta)/sum((Y>=S[ss]))
  }
  f_norm = -sum(log(l0_qt))/lambda - 1
  V_norm = norm(V)
  c = -(f_norm +  V_norm)
  
  betatmp = beta; betatmp[1:3] = 0.5
  
  betahat = optim(betatmp,Qfunc,lambda=lambda,X=X,V=V,A=A,a=a,b=b,c=c,
                  phat=phat,TT=S,Y=Y,delta=delta,dimco=dimco,mt=mt,indexstep=indexstep,method="BFGS",
                  control = list(maxit = 5, trace = TRUE, REPORT = 500))$par
  
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X 
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N,train)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  KMtc=KMcurve_wtc(grid,S,Y,delta,A,d,phat,mt,indexstep)
  if (sum(KMtc == "NaN") > 60){
    truesurv=c(truesurv,NA)
  } else {
    truesurv=c(truesurv,min(KMtc,na.rm = T))
  }
  if(rr==1){
    plot(c(0,S),KMtc,main="KM",type = "l",xlab = "t",ylim=c(0.2,1),ylab = "Survival Proportion",col="red")
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
  MCR[rr] = mean(apply((d[,1:Nstep]==dhat),1,sum) == Nstep)
  
  ########################---subgroup identification--######################
  
  dhat_sub = matrix(0,N,Nstep)
  
  for (tt in 1:Nstep) {
    
    set.seed((rr+203)*44+123)
    N=500
    train <- trainfunc(N); A = train$A; d = train$d; X=train$X
    Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
    Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
    
    censrate = c(censrate,mean(1-delta))
    
    indexstep = indexstepfunc(Y,N,mt,givent)
    
    #############################---get true survival curve---##########################
    S=timegrid(Y,delta,givent)
    grid=length(S)
    
    # giventmp = giventimefn(givent,mt,N,train)
    # phat = giventmp$phat
    # beta = giventmp$betatmp
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
    
    adjust_Time = adjust_func(N,Nstep,Ttime,C)
    
    model0=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A1), propensity.func = prop.func, loss   = "cox_loss_lasso", method="a_learning",
                        nfolds = 10)
    model1=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A2), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
                        nfolds = 10)
    model2=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A3), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
                        nfolds = 10)
    
    set.seed((rr+188)*44+321)
    N=2000
    
    train <- trainfunc(N); A = train$A; d = train$d; X=train$X
    Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
    Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
    
    censrate = c(censrate,mean(1-delta))
    
    indexstep = indexstepfunc(Y,N,mt,givent)
    
    #############################---get true survival curve---##########################
    S=timegrid(Y,delta,givent)
    grid=length(S)
    
    # giventmp = giventimefn(givent,mt,N,train)
    # phat = giventmp$phat
    # beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    
    a1=  predict(model0, newx = X[,tt,], type = "benefit.score")
    
    a2 = predict(model1, newx = X[,tt,], type = "benefit.score")
    
    a3 = predict(model2, newx = X[,tt,], type = "benefit.score")
    
    a = cbind(a1,a2,a3)
    
    dhat_sub[,tt] =  one_vs_all(a)
  }
  
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  # giventmp = giventimefn(givent,mt,N,train)
  # phat = giventmp$phat
  # beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  KMtsub=KMcurve_wtc(grid,S,Y,delta,A,dhat_sub,phat,mt,indexstep)
  if (sum(KMtsub == "NaN") > 60){
    estsurv_sub=c(estsurv_sub,NA)
  } else {
    estsurv_sub=c(estsurv_sub,min(KMtsub,na.rm = T))
  }
  points(c(0,S),KMtsub,type = "l",col="orange")
  MCR_sub[rr] = mean(apply((d[,1:Nstep]==dhat_sub),1,sum) ==   Nstep)
  
  
  ########################-------------------------######################
  ########################---censored q learning---######################
  ########################-------------------------######################
  
  set.seed((rr+203)*44+123)
  N=500
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N,train)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  n = N
  nstage = Nstep
  y =  matrix(0, nrow=n, ncol=nstage)
  u =  matrix(0, nrow=n, ncol=nstage)
  delta_mat = matrix(0, nrow=n, ncol=nstage)
  d_qhat= matrix(0, nrow=2000, ncol=nstage)
  weight_mat = matrix(0, nrow = n, ncol = nstage)
  
  ##--------------------calculate interval survival outcome---------------------
  y[,1] = Ttime$Tfy[,1]
  
  if (nstage != 1){
    for (i in 2:nstage){
      y[,i] =  Ttime$Tfy[,i]
    }
  }
  ##------------calculate min cumulated survival outcome and censored time-----------
  # u[,1]    = pmin(y[,1], C)
  # if (nstage == 2){
  #   u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1]
  # }else if (nstage > 2){
  #   for (j in 3:nstage){
  #     u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
  #   }
  # }
  u[,1]    = pmin(y[,1], C)
  u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1]
  for (j in 3:nstage){
    u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
  }
  
  ##-----------------------censored indicator for each interval------------------------
  delta_mat[,1] = y[,1] < C
  if (nstage > 1){
    for (kkk in 2:nstage){
      delta_mat[,kkk] = (apply(y[,1:kkk], 1, sum) < C)
    }
  }
  
  ##-----------------------generated histroy information matrix-----------------------
  
  d_qhat_pesudo= matrix(0, nrow=2000, ncol=nstage)
  set.seed((rr+203)*44+123)
  N=500
  train <- trainfunc(N); X=train$X
  
  hinfor   = list(nstage)
  for (j in 1:nstage){
    res <- model.matrix(~as.factor(A[,j]))
    hinfor[[j]] = cbind(1,X[,j,],res[,2:3])
  }
  
  ##---------------------caluculate conditional survival function--------------------
  basisinfo =   (cbind(X[,1,]))
  for (i in 1:nstage){
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
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X
  
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  #indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N,train)
  phat = giventmp$phat
  beta = giventmp$betatmp
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
  
  
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N,train)
  phat = giventmp$phat
  #beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  KMtq=KMcurve_wtc(grid,S,Y,delta,A,d_qhat,phat,mt,indexstep)
  if (sum(KMtq == "NaN") > 60){
    estsurv_q=c(estsurv_q, NA)
  } else {
    estsurv_q=c(estsurv_q,min(KMtq,na.rm = T))
  }
  points(c(0,S),KMtq,type = "l",col="green4")
  MCR_q[rr] = mean(apply((d[,1:Nstep] == d_qhat),1,sum) == Nstep)
  print(rr)
}
#####################################################################################
#####################################################################################
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

# 
# truesurv1 = truesurv
# estsurv1 = estsurv
# estsurv_q1 = estsurv_q
# estsurv_sub1 = estsurv_sub
# 
# mean(censrate)
