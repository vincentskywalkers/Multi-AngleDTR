
adjust_func <- function(N,Nstep,Ttime,C){
  
  n = N
  nstage = Nstep
  y =  matrix(0, nrow=n, ncol=nstage)
  u =  matrix(0, nrow=n, ncol=nstage)
  delta_mat = matrix(0, nrow=n, ncol=nstage)
  
  ##--------------------calculate interval survival outcome---------------------
  y[,1] = Ttime$Tfy[,1]
  
  if (nstage != 1){
    for (i in 2:nstage){
      y[,i] =  Ttime$Tfy[,i]
    }
  }
  
  ##------------calculate min cumulated survival outcome and censored time-----------
  u[,1]    = pmin(y[,1], C) 
  u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1] ; 
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
      delta_mat[,kkk] = (apply(y[,1:kkk], 1, sum) < C)
    }
  }
  
  datainfo = list(u=u,y=y,delta_mat=delta_mat)
  return(datainfo)
}

##----------------------------------------------------
XI.gen = function(k){
  XI = matrix(0,k-1,k)
  XI[,1]=rep((k-1)^(-1/2),k-1)
  for (ii in 2:k){
    XI[,ii]=rep( -(1+sqrt(k))/((k-1)^(1.5)), k-1)
    XI[ii-1,ii]=XI[ii-1,ii]+sqrt(k/(k-1))
  }
  return(XI)
}
##----------------------------------------------------

Y.matrix.gen=function(k,nobs,y.train ){
  Y.matrix = matrix(0,nobs,k-1)
  XI=XI.gen(k)
  for (ii in 1:nobs){
    Y.matrix[ii,] = XI[,y.train[ii]]
  }
  return(Y.matrix)
}
##----------------------------------------------------

pred=function(f){
  y=min(which(f==max(f)))
  return(y)
}

##----------------------------------------------------
timegrid <- function(Y,delta,givent){
  SS=sort(unique(Y*delta))
  SS=sort(unique(c(SS[which(SS<=givent)],givent)))
  gridss=length(SS)
  SS=SS[2:gridss]
  return(SS)
}
##-----------------------------------------------------

multiclassfier <- function(a){
  b = c()
  for (i in 1:N){
    
    if (a[i,3] == 1) 
    { 
      b[i]= 3
    }
    else if (a[i,3] == 0 & a[i,2] == 1 ){
      b[i]= 2
    } else {
      b[i]= 1
    }
    
  }
  return(b)
}

##-----------------------------------------------------

multiclassfier2 <- function(a){
  b = c()
  for (i in 1:N){
    
    if (a[i,3] == 1) 
    { 
      b[i]= 3
    }
    else if (a[i,3] == -1 & a[i,1] == 1 ){
      b[i]= 1
    } else {
      b[i]= 2
    }
    
  }
  return(b)
}
##-----------------------------------------------------

KMcurve_wtc <- function(grid,TT,Y,delta,A,d,phat,mt,indexstep){
  Wsi=wei(A,d,phat,indexstep)
  qt=matrix(0,grid,1)
  #qt[1]=1-sum((Y==TT[1])*delta*rep(1,N))/sum((Y>=TT[1])*rep(1,N))
  for (ss in 1:grid) {
    if(TT[ss]>=mt[1]&TT[ss]<mt[2]){
      ww=Wsi[1,]
    }
    if(TT[ss]>=mt[2]&TT[ss]<mt[3]){
      ww=Wsi[2,]
    }
    if(TT[ss]>=mt[3]&TT[ss]<mt[4]){
      ww=Wsi[3,]
    }
    if(TT[ss]>=mt[4]&TT[ss]<mt[5]){
      ww=Wsi[4,]
    }
    if(TT[ss]>=mt[5]){
      ww=Wsi[5,]
    }
    qt[ss]=1-sum((Y==TT[ss])*delta*ww)/sum((Y>=TT[ss])*ww)
  }
  KMt=matrix(0,grid+1,1)
  KMt[1]=1
  for (ii in 1:grid) {
    KMt[ii+1]=prod(qt[1:ii])
  }
  return(KMt)
}

################LUM function#####################################
orig <-function(u,a,c){
  z=1-u
  index <- which( u>=(c/(1+c)) )
  z[index]= 1/(1+c) * (a/( (1+c)*u[index]-c+a ))^a
  return(z)
}
################derivation of LUM function#####################################
dorig <-function(u,a,c){
  z=matrix(-1,length(u),1)
  index <- which( u>=(c/(1+c)) )
  z[index]= -(a/( (1+c)*u[index]-c+a ))^(a+1)
  return(z)
}
############hinge loss function######
hing <- function(u){
  z = apply(cbind(0,1+u),1,max)
  return(z)
}
############logistic surrogate function######
logistic <- function(a,b,c,u){
  z = a/(1+exp(-b*(u-c)))
  return(z)
}



##---------------------------------
decifunc<-function(phat,X,A,beta,dimco){
  N=length(phat[,1])
  Nstep=length(phat[1,])
  fg=matrix(0,N,Nstep)
  for (ss in 1:Nstep) {
    if(ss==1){
      fg[,ss]=X[,1,]%*%beta[1:dimco]
    }
    if(ss==2){
      fg[,ss]=X[,1,]%*%beta[(2*dimco+1):(3*dimco)]+X[,2,]%*%beta[(3*dimco+1):(4*dimco)]+A[,1]*beta[4*dimco+1]
    }
    if(ss==3){
      fg[,ss]=X[,1,]%*%beta[(4*dimco+2+2*dimco+1):(4*dimco+2+2*dimco+dimco)]+X[,2,]%*%beta[(4*dimco+2+2*dimco+dimco+1):(4*dimco+2+2*dimco+2*dimco)]
      +X[,3,]%*%beta[(4*dimco+2+2*dimco+2*dimco+1):(4*dimco+2+2*dimco+3*dimco)]+A[,1]*beta[4*dimco+2+2*dimco+3*dimco+1]+A[,2]*beta[4*dimco+2+2*dimco+3*dimco+2]
    }
    if(ss==4){
      fg[,ss]=X[,1,]%*%beta[(12*dimco+6+1):(12*dimco+6+dimco)]+X[,2,]%*%beta[(12*dimco+6+dimco+1):(12*dimco+6+2*dimco)]
      +X[,3,]%*%beta[(12*dimco+6+2*dimco+1):(12*dimco+6+3*dimco)]+X[,4,]%*%beta[(12*dimco+6+3*dimco+1):(12*dimco+6+4*dimco)]
      +A[,1]*beta[12*dimco+6+4*dimco+1]+A[,2]*beta[12*dimco+6+4*dimco+2]+A[,3]*beta[12*dimco+6+4*dimco+3]
    }
    if(ss==5){
      fg[,ss]=X[,1,]%*%beta[(20*dimco+12+1):(20*dimco+12+dimco)]+X[,2,]%*%beta[(20*dimco+12+dimco+1):(20*dimco+12+2*dimco)]
      +X[,3,]%*%beta[(20*dimco+12+2*dimco+1):(20*dimco+12+3*dimco)]+X[,4,]%*%beta[(20*dimco+12+3*dimco+1):(20*dimco+12+4*dimco)]
      +X[,5,]%*%beta[(20*dimco+12+4*dimco+1):(20*dimco+12+5*dimco)]
      +A[,1]*beta[20*dimco+12+5*dimco+1]+A[,2]*beta[20*dimco+12+5*dimco+2]+A[,3]*beta[20*dimco+12+5*dimco+3]+A[,4]*beta[20*dimco+12+5*dimco+4]
    }
    
  }
  return(fg)
}

decifunc2<-function(phat,X,A,beta,dimco){
  N=length(phat[,1])
  Nstep=length(phat[1,])
  fg=matrix(0,N,Nstep)
  for (ss in 1:Nstep) {
    if(ss==1){
      fg[,ss]=X[,1,]%*%beta[(dimco+1):(2*dimco)]
    }
    if(ss==2){
      fg[,ss]=X[,1,]%*%beta[(4*dimco+1+1):(4*dimco+1+dimco)]+X[,2,]%*%beta[(4*dimco+1+dimco+1):((4*dimco+1+2*dimco))]+A[,1]*beta[4*dimco+2+2*dimco]
    }
    if(ss==3){
      fg[,ss]=X[,1,]%*%beta[(9*dimco+4+1):(9*dimco+4+dimco)]+X[,2,]%*%beta[(9*dimco+4+dimco+1):(9*dimco+4+2*dimco)]
      +X[,3,]%*%beta[(9*dimco+4+2*dimco+1):(9*dimco+4+3*dimco)]+A[,1]*beta[9*dimco+4+3*dimco+1]+A[,2]*beta[9*dimco+4+3*dimco+2]
    }
    if(ss==4){
      fg[,ss]=X[,1,]%*%beta[(12*dimco+6+4*dimco+3+1):(12*dimco+6+4*dimco+3+dimco)]+X[,2,]%*%beta[(12*dimco+6+4*dimco+3+dimco+1):(12*dimco+6+4*dimco+3+2*dimco)]
      +X[,3,]%*%beta[(12*dimco+6+4*dimco+3+2*dimco+1):(12*dimco+6+4*dimco+3+3*dimco)]+X[,4,]%*%beta[(12*dimco+6+4*dimco+3+3*dimco+1):(12*dimco+6+4*dimco+3+4*dimco)]
      +A[,1]*beta[12*dimco+6+4*dimco+3+4*dimco+1]+A[,2]*beta[12*dimco+6+4*dimco+3+4*dimco+2]+A[,3]*beta[12*dimco+6+4*dimco+3+4*dimco+3]
    }
    if(ss==5){
      fg[,ss]=X[,1,]%*%beta[(25*dimco+16+1):(25*dimco+16+dimco)]+X[,2,]%*%beta[(25*dimco+16+dimco+1):(25*dimco+16+2*dimco)]
      +X[,3,]%*%beta[(25*dimco+16+2*dimco+1):(25*dimco+16+3*dimco)]+X[,4,]%*%beta[(25*dimco+16+3*dimco+1):(25*dimco+16+4*dimco)]
      +X[,5,]%*%beta[(25*dimco+16+4*dimco+1):(25*dimco+16+5*dimco)]
      +A[,1]*beta[25*dimco+16+5*dimco+1]+A[,2]*beta[25*dimco+16+5*dimco+2]+A[,3]*beta[25*dimco+16+5*dimco+3]+A[,4]*beta[25*dimco+16+5*dimco+4]
    }
    
  }
  return(fg)
}

## censored q-learning ##
qlearning.censored <- function(hinfor,y_mat, AtRisk, weight_mat,nstage){
  # hinfor, hinfor, a = A_star, u, delta_mat, weight_mat, nstage
  # h0 = hinfor; h1=hinfor; y_mat = u; AtRisk = delta_mat
  yTilde    =  y_mat[ , nstage]
  stageCoef       =  list()
  weight_mat[y_mat[ ,nstage] == 0, nstage] = 0
  
  
  for (stg in nstage:2){
    
    a1 = (hinfor[[stg]][, (dim(hinfor[[stg]])[2]-1)])[AtRisk[, stg] == 1]
    a2 = (hinfor[[stg]][, (dim(hinfor[[stg]])[2])])[AtRisk[, stg] == 1]
    X_at = (hinfor[[stg]][, 2:(dim(hinfor[[stg]])[2]-2)])[AtRisk[, stg] == 1,]
    
    stageModel  =  lm(yTilde[AtRisk[, stg]  == 1] ~ hinfor[[stg]][AtRisk[, stg] == 1, ] +
                        X_at:a1 + X_at:a2 - 1,
                      weights = weight_mat[AtRisk[, stg] == 1, stg])
    
    bTemp       =  coef(stageModel)
    
    bTemp[is.na(bTemp)] = 0
    stageCoef[[stg]] = bTemp
    ## create psuedo outcome for the first stage
    p0Temp      =  ncol(hinfor[[stg]])
    pTemp       =  length(bTemp)
    weight_mat[y_mat[, stg] == 0, stg] = 0
    
    
    a1_tmp = (hinfor[[stg]][, (dim(hinfor[[stg]])[2]-1)])
    a2_tmp = (hinfor[[stg]][, (dim(hinfor[[stg]])[2])])
    X_tmp = (hinfor[[stg]][, 2:(dim(hinfor[[stg]])[2]-2)])
    
    # yTilde      =  (hinfor[[stg]]%*%bTemp[1:p0Temp] + abs(cbind(X_tmp,
    #                                                             X_tmp)%*%bTemp[(p0Temp+1):pTemp])) * 
    #   
    yTilde      =  (hinfor[[stg]]%*%bTemp[1:p0Temp] + (cbind(X_tmp,
                                                             X_tmp)%*%bTemp[(p0Temp+1):pTemp])) * 
      
      (AtRisk[, stg] == 1) + y_mat[, stg - 1];
  }
  # yTilde      =  (hinfor[[stg]]%*%bTemp[1:p0Temp] + abs(cbind(X_tmp,
  #                    X_tmp)%*%bTemp[(p0Temp+1):pTemp]))*(AtRisk[, stg] == 1) + y_mat[, stg - 1];
  # stg  =  1
  for (stg in 1:nstage){
    a1_tmpfit = (hinfor[[stg]][, (dim(hinfor[[stg]])[2]-1)])
    a2_tmpfit = (hinfor[[stg]][, (dim(hinfor[[stg]])[2])])
    
    stage1Model      =  lm(yTilde ~ hinfor[[stg]] +
                             hinfor[[stg]][, 2:(dim(hinfor[[stg]])[2]-2)]:a1_tmpfit +
                             hinfor[[stg]][, 2:(dim(hinfor[[stg]])[2]-2)]:a2_tmpfit - 1,
                           weights = weight_mat[, stg])
    bTemp            =  coef(stage1Model)
    bTemp[is.na(bTemp)] = 0
    stageCoef[[stg]] = bTemp
  }
  return(stageCoef)
}

####weight######

wei<-function(A,d,phat,indexstep){
  N=length(phat[,1])
  Nstep=length(phat[1,])
  W=matrix(0,Nstep,N)
  for (ii in 1:N) {
    if (indexstep[ii] == Nstep){
      for (ss in 1:indexstep[ii]) {
        W[ss,ii]=prod((A[ii,1:ss]==d[ii,1:ss]))/prod(phat[ii,1:ss])  #min(lum[1:ss,ii])
      } 
    } else if (indexstep[ii] < Nstep){
      for (ss in 1:indexstep[ii]) {
        W[ss,ii]=prod((A[ii,1:ss]==d[ii,1:ss]))/prod(phat[ii,1:ss])  #min(lum[1:ss,ii])
      } 
    }
  }
  return(W)
}

####weight######
weight=function(V,A,f1,f2,a,b,c,phat,TT,indexstep){
  N=length(phat[,1])
  grid=length(TT)
  W=matrix(0,Nstep,N)
  U=matrix(0,Nstep,N)
  logit=matrix(0,Nstep,N)
  for (ss in 1:Nstep) {
    Y.matrix = Y.matrix.gen(k,N,A[,ss])
    U[ss,] = apply(Y.matrix*cbind(f1[,ss],f2[,ss]),1,sum) 
    logit[ss,]=logistic(a,b,c,U[ss,])
  }
  for (ii in 1:N) {
    if (indexstep[ii] == Nstep){
      for (ss in 1:indexstep[ii]) {
        W[ss,ii]= prod(logit[1:ss,ii])/prod(phat[ii,1:ss])
      } 
    } else if (indexstep[ii] < Nstep){
      for (ss in 1:indexstep[ii]) {
        logit[(indexstep[ii]),ii] = 1
        W[ss,ii]= prod(logit[1:ss,ii])/prod(phat[ii,1:ss])
      } 
    }
  }
  return(W)
}

###-----------Q function--------
Qfunc <- function(beta,lambda,X,V,A,a,b,c,phat,TT,Y,delta,dimco,mt,indexstep){
  grid=length(TT)
  N=length(Y)
  Nstep=length(phat[1,])
  
  f1=decifunc(phat,X,A,beta,dimco)
  f2=decifunc2(phat,X,A,beta,dimco)
  
  Wsi=weight(V,A,f1,f2,a,b,c,phat,TT,indexstep);
  qt=matrix(0,grid,1);
  # qt[1]=1-sum((Y==TT[1])*delta*rep(1,N))/sum((Y>=TT[1])*rep(1,N))
  for (ss in 1:grid) {
    if(TT[ss]>=mt[1]&TT[ss]<mt[2]){
      ww=Wsi[1,]
    }
    if(TT[ss]>=mt[2]&TT[ss]<mt[3]){
      ww=Wsi[2,]
    }
    if(TT[ss]>=mt[3]&TT[ss]<mt[4]){
      ww=Wsi[3,]
    }
    if(TT[ss]>=mt[4]&TT[ss]<mt[5]){
      ww=Wsi[4,]
    }
    if(TT[ss]>=mt[5]){
      ww=Wsi[5,]
    }
    qt[ss]=1-sum((Y==TT[ss])*delta*ww)/sum((Y>=TT[ss])*ww)
  }
  Sbeta=sum(log(qt))-lambda*sum(beta^2)  ##prod(qt) ##prod(qt) ##prod(qt)
  return(-Sbeta)
}

one_vs_all <- function(class_mat){
  n_mat = nrow(class_mat)
  recommend.trt = c()
  for (i in 1:n_mat){
    recommend.trt[i] = which.max(class_mat[i,])
  }
  return(recommend.trt)
}

KMcurve_wto <- function(grid,TT,Y,delta){
  qt=matrix(0,grid,1)
  for (ii in 1:grid) {
    qt[ii]=1-sum((Y==TT[ii])*delta)/sum((Y>=TT[ii]))
  }
  KMt=matrix(0,grid+1,1)
  KMt[1]=1
  for (ii in 1:grid) {
    KMt[ii+1]=prod(qt[1:ii])
  }
  return(KMt)
}


prop.func <- function(x, trt)
{
  # fit propensity score model
  propens.model <- cv.glmnet(y = trt,
                             x = x, family = "binomial")
  pi.x <- predict(propens.model, s = "lambda.min",
                  newx = x, type = "response")[,1]
  pi.x
}


cal_survival_and_density_functions    <- function(X, Z, status, X_new, Z_new, fitting_method = "Cox"){
  
  if (fitting_method == "Cox"){
    
    theta       <- coxph(Surv(X, status) ~ Z)$coefficients
    X_l         <- X[status == 1]
    Z_l         <- Z[status == 1, ]
    if (is.vector(Z_l)){
      Z_l = matrix(Z_l, ncol = 1)
    }
    lambda_l    <- 1/(sum.I(X_l, "<", X, exp(Z %*% theta)) + exp(Z_l %*% theta))
    Lambda_new      <- sum.I(X_new, ">=", X_l, lambda_l)
    S_new           <- exp( - exp(Z_new %*% theta)* Lambda_new)
  } else {
    print("Wrong fitting method")
  }
  return( list( S = S_new, Lambda = Lambda_new, theta = theta) )
}


sum.I   <- function (yy, FUN, Yi, Vi = NULL)
{
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  pos <- rank(c(yy, Yi), ties.method = "f")[1:length(yy)] -
    rank(yy, ties.method = "f")
  if (substring(FUN, 2, 2) == "=")
    pos <- length(Yi) - pos
  if (!is.null(Vi)) {
    if (substring(FUN, 2, 2) == "=")
      tmpind <- order( - Yi )
    else tmpind <- order(Yi) 
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ]) 
  }
  else return(pos)
}

