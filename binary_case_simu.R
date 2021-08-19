library(personalized)
library(survival)
#source("utility_func.R")
library(DTRreg)

Run =50
N = 500;
M = 5;
k=2
dimco = 25;
givent = 1.5;
#mt = c(0,0.5,1,1.5,2)
mt = c(0.5,1,1.5,2,2.5)
censrate=NULL
truesurv=NULL
estsurv=NULL
estsurv_q=NULL
estsurv_sub=NULL
estsurv_DWsurv=NULL
MCR=matrix(0,Run,1)
MCR_q=matrix(0,Run,1)
MCR_sub=matrix(0,Run,1)
MCR_DWsurv =matrix(0,Run,1)
V=XI.gen(k)
beta_set1 = c(1,1,1,rep(0,dimco-3))
#beta_set2 = c(1,-1,-1,rep(0,dimco-3))
beta1 = c(beta_set1) 
beta2 = c(0.5*beta_set1,-1.5*beta_set1,0.1)
beta3 = c(0.25*beta_set1,-0.5*beta_set1,1*beta_set1,0.05,0.1)
beta4 = c(0.1*beta_set1,-0.25*beta_set1,0.5*beta_set1,1*beta_set1,0.05,-0.05,0.1)
beta5 = c(0.05*beta_set1,-0.1*beta_set1,0.25*beta_set1,0.5*beta_set1,1*beta_set1,0.05,-0.05,0.05,-0.1)


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
  
  giventmp = giventimefn(givent,mt,N)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  #########################----proposed method---###############################
  #lambda = 0.001   givent = 1.5
  #lambda = 0.0001  givent = 2.5
  #lambda = 0.00001


  #lambda = 0.000001   5stage
  lambda = 0.000001
  a=1; b=1;
  grid=length(S)
  l0_qt =matrix(0,grid,1);
  for (ss in 1:grid) {
    l0_qt[ss]=1-sum((Y==S[ss])*delta)/sum((Y>=S[ss]))
  }
  f_norm = -sum(log(l0_qt))/lambda - 1
  V_norm = norm(V)
  c = -(f_norm +  V_norm)
  
  #betatmp = beta; betatmp[1:3] = 0.5
  
  betatmp = beta; betatmp[1:3] = 0.5
  
  betahat = optim(betatmp,Qfunc,lambda=lambda,X=X,V=V,A=A,a=a,b=b,c=c,
                  phat=phat,TT=S,Y=Y,delta=delta,dimco=dimco,mt=mt,indexstep=indexstep,method="BFGS",
                  control = list(maxit = 20, trace = TRUE, REPORT = 500))$par

  
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X 
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N)
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
  #fhat2=decifunc2(phat,X,A,betahat,dimco)

  dhat = matrix(0,N,Nstep)
  for (tt in 1:Nstep){
    inner.matrix=matrix(0,N,k)
    for (ii in 1:k){
      inner.matrix[,ii] = apply(V[,ii]*cbind(fhat1[,tt]),1,sum)
    }
    dhat[,tt]=apply(inner.matrix,1,pred)
  }
  dhat = transd(dhat,N,Nstep)
  
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
    
    giventmp = giventimefn(givent,mt,N)
    phat = giventmp$phat
    beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    # A1  = A[,tt]
    # A1[A1==3] = 0
    # A1[A1==2] = 0
    # 
    # A2   = A[,tt]
    # A2[A2==1] = 0
    # A2[A2==3] = 0
    # A2[A2==2] = 1
    # 
    # 
    # A3   = A[,tt]
    # A3[A3==1] = 0
    # A3[A3==2] = 0
    # A3[A3==3] = 1
    
    adjust_Time = adjust_func(N,Nstep,Ttime,C)

    model0=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A[,tt]), propensity.func = prop.func, loss   = "cox_loss_lasso", method="weighting",
                        nfolds = 10)
    # model1=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A2), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
    #                     nfolds = 10)
    # model2=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A3), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
    #                     nfolds = 10)
    
    # model0=fit.subgroup(X[,tt,], y = Surv(Y, delta), as.vector(A1), propensity.func = prop.func, loss   = "cox_loss_lasso", method="a_learning",
    #                     nfolds = 10)
    # model1=fit.subgroup(X[,tt,], y = Surv(Y, delta), as.vector(A2), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
    #                     nfolds = 10)
    # model2=fit.subgroup(X[,tt,], y = Surv(Y, delta), as.vector(A3), propensity.func = prop.func, loss   = "cox_loss_lasso",method="a_learning",
    #                     nfolds = 10)

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
    
    giventmp = giventimefn(givent,mt,N)
    phat = giventmp$phat
    beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    
    a1=  predict(model0, newx = X[,tt,], type = "trt.group")
    
    # a2 = predict(model1, newx = X[,tt,], type = "benefit.score")
    # 
    # a3 = predict(model2, newx = X[,tt,], type = "benefit.score")
    # 
    # a = cbind(a1,a2,a3)
    # 
    dhat_sub[,tt] =  as.vector(a1)
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
  
  giventmp = giventimefn(givent,mt,N)
  phat = giventmp$phat
  beta = giventmp$betatmp
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
  
  giventmp = giventimefn(givent,mt,N)
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

  
  pesudo_trt <- function(A_star){
    
    d_qhat_pesudo= matrix(0, nrow=2000, ncol=nstage)
    set.seed((rr+203)*44+123)
    N=500
    train <- trainfunc(N); X=train$X
    
    hinfor   = list(nstage)
    hinfor[[1]] = cbind(1,X[,1,])
    for (j in 2:nstage){ 
      hinfor[[j]] = cbind(1,X[,j,],A_star[,j-1])
    }
    
    ##---------------------caluculate conditional survival function--------------------
    basisinfo =   hinfor[[1]][,-1]
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
    qlearningcoef = qlearning.censored(hinfor, hinfor, a = A_star, u, delta_mat, weight_mat, nstage)
    
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
    
    giventmp = giventimefn(givent,mt,N)
    phat = giventmp$phat
    beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    
    hinfor   = list(nstage)
    hinfor[[1]] = cbind(1,X[,1,])
    for (j in 2:nstage){ 
      hinfor[[j]] = cbind(1,X[,j,],A_star[,j-1])
    }
    
    if (nstage == 5){
      f1hat = hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
      
      f2hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))]) 
               + hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))]))
      
      f3hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
               +hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))])
               +hinfor[[3]][,2:(ncol(hinfor[[3]]))]%*%as.matrix(qlearningcoef[[3]][(ncol(hinfor[[3]])+2):(2*ncol(hinfor[[3]]))]))
      
      f4hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
               +hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))])
               +hinfor[[3]][,2:(ncol(hinfor[[3]]))]%*%as.matrix(qlearningcoef[[3]][(ncol(hinfor[[3]])+2):(2*ncol(hinfor[[3]]))])
               +hinfor[[4]][,2:(ncol(hinfor[[4]]))]%*%as.matrix(qlearningcoef[[4]][(ncol(hinfor[[4]])+2):(2*ncol(hinfor[[4]]))]))
      
      f5hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
               +hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))])
               +hinfor[[3]][,2:(ncol(hinfor[[3]]))]%*%as.matrix(qlearningcoef[[3]][(ncol(hinfor[[3]])+2):(2*ncol(hinfor[[3]]))])
               +hinfor[[4]][,2:(ncol(hinfor[[4]]))]%*%as.matrix(qlearningcoef[[4]][(ncol(hinfor[[4]])+2):(2*ncol(hinfor[[4]]))])
               +hinfor[[5]][,2:(ncol(hinfor[[5]]))]%*%as.matrix(qlearningcoef[[5]][(ncol(hinfor[[5]])+2):(2*ncol(hinfor[[5]]))]))
      
      fhat = cbind(f1hat,f2hat,f3hat,f4hat,f5hat)
      
    }else if (nstage == 3){
      f1hat = hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
      
      f2hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))]) 
               + hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))]))
      f3hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
               +hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))])
               +hinfor[[3]][,2:(ncol(hinfor[[3]]))]%*%as.matrix(qlearningcoef[[3]][(ncol(hinfor[[3]])+2):(2*ncol(hinfor[[3]]))]))
      fhat = cbind(f1hat,f2hat,f3hat)
    }else if (nstage == 2){
      f1hat = hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))])
      f2hat = (hinfor[[1]][,2:(ncol(hinfor[[1]]))]%*%as.matrix(qlearningcoef[[1]][(ncol(hinfor[[1]])+2):(2*ncol(hinfor[[1]]))]) 
               + hinfor[[2]][,2:(ncol(hinfor[[2]]))]%*%as.matrix(qlearningcoef[[2]][(ncol(hinfor[[2]])+2):(2*ncol(hinfor[[2]]))]))
      fhat = cbind(f1hat,f2hat)
      
    }
    
    for (stg in 1:nstage){
      #  d_qhat_pesudo[,stg] = sign(fhat[,stg])
      d_qhat_pesudo[,stg] = sign(fhat[,stg])
    }
    return(d_qhat_pesudo)
  }
  
  d_qhat = pesudo_trt(A)
  

  # for (stg in 1:nstage){
  #   d_qhat[,stg] = one_vs_all(cbind(a5[,stg],a6[,stg],a7[,stg]))
  # }
  # 
  
  set.seed((rr+188)*44+321)
  N=2000
  train <- trainfunc(N); A = train$A; d = train$d; X=train$X
  Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
  Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
  
  censrate = c(censrate,mean(1-delta))
  
  indexstep = indexstepfunc(Y,N,mt,givent)
  
  S=timegrid(Y,delta,givent)
  grid=length(S)
  
  giventmp = giventimefn(givent,mt,N)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  KMtq=KMcurve_wtc(grid,S,Y,delta,A,d_qhat,phat,mt,indexstep)
  if (sum(KMtq == "NaN") > 60){
    estsurv_q=c(estsurv_q, NA)
  } else {
    estsurv_q=c(estsurv_q,min(KMtq,na.rm = T))
  }
  points(c(0,S),KMtq,type = "l",col="green4")
  MCR_q[rr] = mean(apply((d[,1:Nstep] == d_qhat),1,sum) == Nstep)
  
######################## DWSurv #######################
  
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
  
  giventmp = giventimefn(givent,mt,N)
  phat = giventmp$phat
  beta = giventmp$betatmp
  Nstep=length(phat[1,])
  
  adjust_Time = adjust_func(N,Nstep,Ttime,C)
  
  # model0=fit.subgroup(X[,tt,], y = Surv(adjust_Time$u[,tt], adjust_Time$delta_mat[,tt]), as.vector(A[,tt]), propensity.func = prop.func, loss   = "cox_loss_lasso", method="a_learning",
  #                     nfolds = 10)

  if (Nstep ==3){
u1 = adjust_Time$u[,1]
u2 = adjust_Time$u[,2]
u3 = adjust_Time$u[,3]

X1 = X[,1,]
X2 = X[,2,]
X3 = X[,3,]

A1 = A[,1]
A2 = A[,2]
A3 = A[,3]

X11 = X[,1,1]
X12 = X[,1,2]
X13 = X[,1,3]
X14 = X[,1,4]
X15 = X[,1,5]
X16 = X[,1,6]
X17 = X[,1,7]
X18 = X[,1,8]
X19 = X[,1,9]
X110= X[,1,10]

X21 = X[,2,1]
X22 = X[,2,2]
X23 = X[,2,3]
X24 = X[,2,4]
X25 = X[,2,5]
X26 = X[,2,6]
X27 = X[,2,7]
X28 = X[,2,8]
X29 = X[,2,9]
X210 = X[,2,10]

X31 = X[,3,1]
X32 = X[,3,2]
X33 = X[,3,3]
X34 = X[,3,4]
X35 = X[,3,5]
X36 = X[,3,6]
X37 = X[,3,7]
X38 = X[,3,8]
X39 = X[,3,9]
X310 = X[,3,10]


delta1 = adjust_Time$delta_mat[,1]
delta2 = adjust_Time$delta_mat[,2]
delta3 = adjust_Time$delta_mat[,3]

set.seed((rr+188)*44+321)
N=2000
train <- trainfunc(N); A = train$A; d = train$d; X=train$X


X11test = X[,1,1]
X12test = X[,1,2]
X13test = X[,1,3]
X14test = X[,1,4]
X15test = X[,1,5]
X16test = X[,1,6]
X17test = X[,1,7]
X18test = X[,1,8]
X19test = X[,1,9]
X110test = X[,1,10]


X21test = X[,2,1]
X22test = X[,2,2]
X23test = X[,2,3]
X24test = X[,2,4]
X25test = X[,2,5]
X26test = X[,2,6]
X27test = X[,2,7]
X28test = X[,2,8]
X29test = X[,2,9]
X210test = X[,2,10]

X31test = X[,3,1]
X32test = X[,3,2]
X33test = X[,3,3]
X34test = X[,3,4]
X35test = X[,3,5]
X36test = X[,3,6]
X37test = X[,3,7]
X38test = X[,3,8]
X39test = X[,3,9]
X310test = X[,3,10]

A1test = A[,1]
A2test = A[,2]
A3test = A[,3]


#mydata = data.frame(u1,u2,u3,X1,X2,X3,A1,A2,A3,X1test,X2test,X3test,X11,X13,X21,X23,X31,X33,delta1,delta2,delta3)
# mydata1 = data.frame(u1,u2,u3,
#                      X11,X12,X13,X14,X15,X16,X17,X18,X19,X110,
#                      X21,X22,X23,X24,X25,X26,X27,X28,X29,X210,
#                      X31,X32,X33,X34,X35,X36,X37,X38,X39,X310,
#                      X11test,X12test,X13test,X14test,X15test,X16test,X17test,X18test,X19test,X110test,
#                      X21test,X22test,X23test,X24test,X25test,X26test,X27test,X28test,X29test,X210test,
#                      X31test,X32test,X33test,X34test,X35test,X36test,X37test,X38test,X39test,X310test,
#                      A1,A2,A3,delta1,delta2,delta3)

mydata1 = data.frame(u1,u2,u3,
                     X11,X12,X13,
                     X21,X22,X23,
                     X31,X32,X33,
                     X11test,X12test,X13test,
                     X21test,X22test,X23test,
                     X31test,X32test,X33test,
                     A1,A2,A3,
                     A1test,A2test,A3test,
                     delta1,delta2,delta3)

mod <- DWSurv(time = list(~u1, ~u2,~u3),
              blip.mod = list(~1+X11+X12+X13,
                              ~1+X21+X22+X23+A1,
                              ~1+X31+X32+X33+A1+A2),
              treat.mod = list(A1~1,
                               A2~1,
                               A3~1),
              tf.mod =  list(~1+X11+X12+X13,
                             ~1+X21+X22+X23,
                             ~1+X31+X32+X33),
              cens.mod = list(delta1~1, 
                              delta2~1,
                              delta3~1),
              data=mydata1)


dhat_DWsurvtmp = matrix(NA,N,Nstep)

#summary(mod)

for (i in 1:2000){
  dhat_DWsurvtmp[i,1] = as.numeric(I(sum(mod$psi[[1]] * cbind(1,X11test[i],X12test[i],X13test[i])) >0))
}

for (i in 1:2000){
  dhat_DWsurvtmp[i,2] = as.numeric(I(sum(mod$psi[[2]] * cbind(1,X21test[i],X22test[i],X23test[i],A1test[i])) >0))
}

for (i in 1:2000){
  dhat_DWsurvtmp[i,3] = as.numeric(I(sum(mod$psi[[3]] * cbind(1,X31test[i],X32test[i],X33test[i],A1test[i],A2test[i])) >0))
}

opt.treat = list()
opt.treat[[1]] = dhat_DWsurvtmp[,1]
opt.treat[[2]] = dhat_DWsurvtmp[,2]
opt.treat[[3]] = dhat_DWsurvtmp[,3]
# mod <- DWSurv(time = list(~u1, ~u2,~u3), blip.mod = list(~X1, ~X2,~X3),
#               treat.mod = list(A1~1, A2~1,A3~1), tf.mod = list(~X11+X13, ~X21+X23,~X31+X33),
#               cens.mod = list(delta1~1, delta2~1,delta3~1), var.estim = "asymptotic",data=mydata)
# 
# mod <- DWSurv(time = list(~u1, ~u2,~u3), blip.mod = list(~X1+A1, ~X2+A2,~X3+A3),
#               treat.mod = list(A1~1, A2~1,A3~1), tf.mod = list(~X1test, ~X2test,~X3test),
#               cens.mod = list(delta1~1, delta2~1,delta3~1), var.estim = "none",data=mydata)
# # 
# 
# mod <- DWSurv(time = list(~u1, ~u2,~u3),
#               blip.mod = list(~1+X11+X12+X13+X14+X15+X16+X17+X18+X19+X110,
#                               ~1+X21+X22+X23+X24+X25+X26+X27+X28+X29+X210,
#                               ~1+X31+X32+X33+X34+X35+X36+X37+X38+X39+X310),
#               treat.mod = list(A1~X11+X12+X13+X14+X15+X16+X17+X18+X19+X110,
#                                A2~X21+X22+X23+X24+X25+X26+X27+X28+X29+X210,
#                                A3~X31+X32+X33+X34+X35+X36+X37+X38+X39+X310),
#               tf.mod = list(~X11+X12+X13+X14+X15+X16+X17+X18+X19+X110,
#                             ~X21+X22+X23+X24+X25+X26+X27+X28+X29+X210+A1,
#                             ~X31+X32+X33+X34+X35+X36+X37+X38+X39+X310+A1+A2),
#                cens.mod = list(delta1~1, delta2~1,delta3~1),data=mydata1)

# mod <- DWSurv(time = list(~u1, ~u2,~u3),
#               blip.mod = list(~X11+X12+X13,
#                               ~X21+X22+X23,
#                               ~X31+X32+X33),
#               treat.mod = list(A1~1+X11+X12+X13,
#                                A2~1+X21+X22+X23,
#                                A3~1+X31+X32+X33),
#               tf.mod =  list(~1+X11+X12+X13,
#                              ~1+X21+X22+X23+A1,
#                              ~1+X31+X32+X33+A1+A2),
#               cens.mod = list(delta1~1, delta2~1,delta3~1),data=mydata1)
# 

# mod <- DWSurv(time = list(~u1, ~u2,~u3),
#               blip.mod = list(~1+X11+X12+X13+X14+X15+X16+X17+X18+X19+X110,
#                               ~1+X21+X22+X23+X24+X25+X26+X27+X28+X29+X210+A1,
#                               ~1+X31+X32+X33+X34+X35+X36+X37+X38+X39+X310+A1+A2),
#               treat.mod = list(A1~1+X11test+X12test+X13test+X14test+X15test+X16test+X17test+X18test+X19test+X110test,
#                                A2~1+X21test+X22test+X23test+X24test+X25test+X26test+X27test+X28test+X29test+X210test,
#                                A3~1+X31test+X32test+X33test+X34test+X35test+X36test+X37test+X38test+X39test+X310test),
#               tf.mod = list(~1+X11test+X12test+X13test+X14test+X15test+X16test+X17test+X18test+X19test+X110test,
#                             ~1+X21test+X22test+X23test+X24test+X25test+X26test+X27test+X28test+X29test+X210test+A1,
#                             ~1+X31test+X32test+X33test+X34test+X35test+X36test+X37test+X38test+X39test+X310test+A1+A2),
#                cens.mod = list(delta1~1, delta2~1,delta3~1),data=mydata1)


dhat_DWsurv <-  tranDW(opt.treat,N,Nstep)

Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
censrate = c(censrate,mean(1-delta))

indexstep = indexstepfunc(Y,N,mt,givent)

S=timegrid(Y,delta,givent)
grid=length(S)

phat = giventmp$phat
Nstep=length(phat[1,])

KMtDWsurv=KMcurve_wtc(grid,S,Y,delta,A,dhat_DWsurv,phat,mt,indexstep)
if (sum(KMtq == "NaN") > 60){
  estsurv_DWsurv=c(estsurv_DWsurv, NA)
} else {
  estsurv_DWsurv=c(estsurv_DWsurv,min(KMtDWsurv,na.rm = T))
}
points(c(0,S),KMtDWsurv,type = "l",col="cyan")
MCR_DWsurv[rr] = mean(apply((d[,1:Nstep] == dhat_DWsurv),1,sum) == Nstep)

  }
  
  
  if (Nstep ==5){
    u1 = adjust_Time$u[,1]
    u2 = adjust_Time$u[,2]
    u3 = adjust_Time$u[,3]
    u4 = adjust_Time$u[,4]
    u5 = adjust_Time$u[,5]
    
    X1 = X[,1,]
    X2 = X[,2,]
    X3 = X[,3,]
    X4 = X[,4,]
    X5 = X[,5,]
    
    A1 = A[,1]
    A2 = A[,2]
    A3 = A[,3]
    A4 = A[,4]
    A5 = A[,5]
    
    X11 = X[,1,1]
    X12 = X[,1,2]
    X13 = X[,1,3]

    X21 = X[,2,1]
    X22 = X[,2,2]
    X23 = X[,2,3]

    X31 = X[,3,1]
    X32 = X[,3,2]
    X33 = X[,3,3]

    X41 = X[,4,1]
    X42 = X[,4,2]
    X43 = X[,4,3]
    
    X51 = X[,5,1]
    X52 = X[,5,2]
    X53 = X[,5,3]
    
    delta1 = adjust_Time$delta_mat[,1]
    delta2 = adjust_Time$delta_mat[,2]
    delta3 = adjust_Time$delta_mat[,3]
    delta4 = adjust_Time$delta_mat[,4]
    delta5 = adjust_Time$delta_mat[,5]
    
    
    set.seed((rr+188)*44+321)
    N=2000
    train <- trainfunc(N); A = train$A; d = train$d; X=train$X
    
    X11test = X[,1,1]
    X12test = X[,1,2]
    X13test = X[,1,3]
    
    X21test = X[,2,1]
    X22test = X[,2,2]
    X23test = X[,2,3]
    
    X31test = X[,3,1]
    X32test = X[,3,2]
    X33test = X[,3,3]

    X41test = X[,4,1]
    X42test = X[,4,2]
    X43test = X[,4,3]
    
    X51test = X[,5,1]
    X52test = X[,5,2]
    X53test = X[,5,3]
    
    
    A1test = A[,1]
    A2test = A[,2]
    A3test = A[,3]
    A4test = A[,4]
    A5test = A[,5]
    
   # mydata = data.frame(u1,u2,u3,u4,u5,X1,X2,X3,X4,X5,A1,A2,A3,A4,A5,X1test,X2test,X3test,X4test,X5test,delta1,delta2,delta3,delta4,delta5)
    
    mydata2 = data.frame(u1,u2,u3,u4,u5,
                         X11,X12,X13,
                         X21,X22,X23,
                         X31,X32,X33,
                         X41,X42,X43,
                         X51,X52,X53,
                         X11test,X12test,X13test,
                         X21test,X22test,X23test,
                         X31test,X32test,X33test,
                         X41test,X42test,X43test,
                         X51test,X52test,X53test,
                         A1,A2,A3,A4,A5,
                         A1test,A2test,A3test,A4test,A5test,
                         delta1,delta2,delta3,delta4,delta5)
    
    
    # mod <- DWSurv(time = list(~u1, ~u2,~u3), blip.mod = list(~X1, ~X2,~X3),
    #               treat.mod = list(A1~1, A2~1,A3~1), tf.mod = list(~X11+X13, ~X21+X23,~X31+X33),
    #               cens.mod = list(delta1~1, delta2~1,delta3~1), var.estim = "asymptotic",data=mydata)
    # mod <- DWSurv(time = list(~u1, ~u2,~u3),
    #               blip.mod = list(~1+X11+X12+X13,
    #                               ~1+X21+X22+X23+A1,
    #                               ~1+X31+X32+X33+A1+A2,
    #                               ~1+X41+X42+X43+A1+A2+A3,
    #                               ~1+X51+X52+X53+A1+A2+A3+A4),
    #               treat.mod = list(A1~1+X11+X12+X13,
    #                                A2~1+X21+X22+X23,
    #                                A3~1+X31+X32+X33,
    #                                A4~1+X41+X42+X43,
    #                                A5~1+X51+X52+X53),
    #               tf.mod =  list(~1+X11+X12+X13,
    #                              ~1+X21+X22+X23+A1,
    #                              ~1+X31+X32+X33+A1+A2,
    #                              ~1+X41+X42+X43+A1+A2+A3,
    #                              ~1+X51+X52+X53+A1+A2+A3+A4),
    #               cens.mod = list(delta1~1, delta2~1,delta3~1,delta4~1,delta5~1),data=mydata2)
    # 
    
    mod <- DWSurv(time = list(~u1, ~u2,~u3,~u4,~u5),
                  blip.mod = list(~1+X11+X12+X13,
                                  ~1+X21+X22+X23+A1,
                                  ~1+X31+X32+X33+A1+A2,
                                  ~1+X41+X42+X43+A1+A2+A3,
                                  ~1+X51+X52+X53+A1+A2+A3+A4),
                  treat.mod = list(A1~1,
                                   A2~1,
                                   A3~1,
                                   A4~1,
                                   A5~1),
                  tf.mod =  list(~1+X11+X12+X13,
                                 ~1+X21+X22+X23,
                                 ~1+X31+X32+X33,
                                 ~1+X41+X42+X43,
                                 ~1+X51+X52+X53),
                  cens.mod = list(delta1~1, delta2~1,delta3~1,delta4~1,delta5~1),data=mydata2)
    
    
    dhat_DWsurvtmp = matrix(NA,N,Nstep)
    
    #summary(mod)
    
    for (i in 1:2000){
      dhat_DWsurvtmp[i,1] = as.numeric(I(sum(mod$psi[[1]] * cbind(1,X11test[i],X12test[i],X13test[i])) >0))
    }
    
    for (i in 1:2000){
      dhat_DWsurvtmp[i,2] = as.numeric(I(sum(mod$psi[[2]] * cbind(1,X21test[i],X22test[i],X23test[i],A1test[i])) >0))
    }
    
    for (i in 1:2000){
      dhat_DWsurvtmp[i,3] = as.numeric(I(sum(mod$psi[[3]] * cbind(1,X31test[i],X32test[i],X33test[i],A1test[i],A2test[i])) >0))
    }
    
    for (i in 1:2000){
      dhat_DWsurvtmp[i,4] = as.numeric(I(sum(mod$psi[[4]] * cbind(1,X41test[i],X42test[i],X43test[i],A1test[i],A2test[i],A3test[i])) >0))
    }
    
    for (i in 1:2000){
      dhat_DWsurvtmp[i,5] = as.numeric(I(sum(mod$psi[[5]] * cbind(1,X51test[i],X52test[i],X53test[i],A1test[i],A2test[i],A3test[i],A4test[i])) >0))
    }
    
    opt.treat = list()
    opt.treat[[1]] = dhat_DWsurvtmp[,1]
    opt.treat[[2]] = dhat_DWsurvtmp[,2]
    opt.treat[[3]] = dhat_DWsurvtmp[,3]
    opt.treat[[4]] = dhat_DWsurvtmp[,4]
    opt.treat[[5]] = dhat_DWsurvtmp[,5]
    
    # mod <- DWSurv(time = list(~u1, ~u2,~u3,~u4,~u5), blip.mod = list(~X1, ~X2+A1,~X3+A1+A2,~X4+A1+A2+A3,~X5+A1+A2+A3+A4),
    #               treat.mod = list(A1~1, A2~1,A3~1, A4~1,A5~1), tf.mod = list(~X1test, ~X2test,~X3test ,~X4test,~X5test),
    #               cens.mod = list(delta1~1, delta2~1,delta3~1,delta4~1,delta5~1), var.estim = "none",data=mydata)
    # 
    dhat_DWsurv <-  tranDW(opt.treat,N,Nstep)
    
    Ttime = Tffunc(N,M,A,d,X); Tf = Ttime$Tf; C = Ctfunc(mt,N)
    Y = apply(cbind(Tf,C),1,min); delta = (Tf<=C)
    censrate = c(censrate,mean(1-delta))
    
    indexstep = indexstepfunc(Y,N,mt,givent)
    
    S=timegrid(Y,delta,givent)
    grid=length(S)
    
    phat = giventmp$phat
    Nstep=length(phat[1,])
    
    KMtDWsurv=KMcurve_wtc(grid,S,Y,delta,A,dhat_DWsurv,phat,mt,indexstep)
    if (sum(KMtq == "NaN") > 60){
      estsurv_DWsurv=c(estsurv_DWsurv, NA)
    } else {
      estsurv_DWsurv=c(estsurv_DWsurv,min(KMtDWsurv,na.rm = T))
    }
    points(c(0,S),KMtDWsurv,type = "l",col="cyan")
    MCR_DWsurv[rr] = mean(apply((d[,1:Nstep] == dhat_DWsurv),1,sum) == Nstep)
    
  }

  # mod <- DWSurv(time = list(~adjust_Time$u[,1], ~adjust_Time$u[,2],~adjust_Time$u[,3]), blip.mod = list(~X[,1,], ~X[,2,],~X[,3,]), 
  #               treat.mod = list(A[,1]~1, A[,2]~1,A[,3]~1), tf.mod = list(~X[,1,1]+X[,1,3], ~X[,2,1]+X[,2,3],~X[,3,1]+X[,3,3]), 
  #               cens.mod = list(adjust_Time$delta_mat[,1]~1, adjust_Time$delta_mat[,2]~1,adjust_Time$delta_mat[,3]~1), var.estim = "asymptotic")
  # 
  # 
  
  
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
print(mean(estsurv_DWsurv,na.rm=T))
print(sd(estsurv_DWsurv,na.rm = T))
#####################################################################################

# 
# truesurv1 = truesurv
# estsurv1 = estsurv
# estsurv_q1 = estsurv_q
# estsurv_sub1 = estsurv_sub
# 
mean(censrate)
