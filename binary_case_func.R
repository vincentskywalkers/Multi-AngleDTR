##----------------------------------------------------
Ctfunc <- function(mt,N){
  C = sample(c(mt,3),N,replace = T,prob=c(5,5,5,1,1,1))
 #C = sample(c(mt,3),N,replace = T,prob=c(0.1,0.25,0.25,0.25,1,1))
 # C = sample(c(mt,3),N,replace = T,prob=c(0.5,0.5,0.5,1,1,1))
  return(C)
}
##----------------------------------------------------


Tffunc <- function(N,M,A,d,X){
  Tfy = matrix(NA,N,M); 
  Tf = c()
  for (i in 1:N){
    Tfy[i,1] = I(d[i,1]==A[i,1])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,1,1])^3  + X[i,1,3] >0),0.15+0.05*I((X[i,1,1])^3  + X[i,1,3] <=0))) +
      I(d[i,1]!=A[i,1])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,1,5])^2-(X[i,1,6]) < 0.1),0.25+0.025*I((X[i,1,5])^2-(X[i,1,6]) >= 0.1),0.15))
    for (ttt in 2:(M-2)){
      Tmp =  I(d[i,ttt]==A[i,ttt])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] >0),0.15+0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] <=0))) + 
        I(d[i,ttt]!=A[i,ttt])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) < 0.1),0.25+0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) >= 0.1),0.15))
      Tfy[i,ttt] = I(Tfy[i,(ttt-1)]>0.25)*Tmp
    }
    for (ttt in (M-1):(M)){
      Tmp =  I(d[i,ttt]==A[i,ttt])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] >0),0.15+0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] <=0))) + 
        I(d[i,ttt]!=A[i,ttt])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) < 0.1),0.25+0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) >= 0.1),0.15))
      Tfy[i,ttt] = I(Tfy[i,(ttt-1)]>0.25)*Tmp
    }
    Tf[i] = sum(Tfy[i,]) + 0.5 
  }
  datainfo = list(Tf = Tf, Tfy = Tfy)
  return(datainfo)
}


# Tffunc <- function(N,M,A,d,X){
#   Tfy = matrix(NA,N,M); 
#   Tf = c()
#   for (i in 1:N){
#     Tfy[i,1] = I(d[i,1]==A[i,1])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,1,1])^3  + X[i,1,3] >0),0.15+0.05*I((X[i,1,1])^3  + X[i,1,3] <=0))) +
#     I(d[i,1]!=A[i,1])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,1,5])^2-(X[i,1,6]) < 0.1),0.25+0.025*I((X[i,1,5])^2-(X[i,1,6]) >= 0.1),0.15))
#     for (ttt in 2:(M-2)){
#       Tmp =  I(d[i,ttt]==A[i,ttt])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] >0),0.15+0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] <=0))) + 
#         I(d[i,ttt]!=A[i,ttt])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) < 0.1),0.25+0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) >= 0.1),0.15))
#       Tfy[i,ttt] = I(Tfy[i,(ttt-1)]>0.25)*Tmp
#     }
#     for (ttt in (M-1):(M)){
#       Tmp =  I(d[i,ttt]==A[i,ttt])*sample(c(0.5,0.25),size=1,prob=c(0.85-0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] >0),0.15+0.05*I((X[i,ttt,1])^3  + X[i,ttt,3] <=0))) + 
#         I(d[i,ttt]!=A[i,ttt])*sample(c(0.5,0.25,0),size=1,prob=c(0.6-0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) < 0.1),0.25+0.025*I((X[i,ttt,5])^2-(X[i,ttt,6]) >= 0.1),0.15))
#       Tfy[i,ttt] = I(Tfy[i,(ttt-1)]>0.25)*Tmp
#     }
#     Tf[i] = sum(Tfy[i,]) + 0.5 
#   }
#   datainfo = list(Tf = Tf, Tfy = Tfy)
#   return(datainfo)
# }


##----------------------------------------------------
trainfunc <- function(N){
  
  X=array(rnorm(N*M*dimco,0,1),dim=c(N,M,dimco));
  A=f=d=matrix(0,N,M);
  Pr=matrix(1/3,N,M)
  
  for (i in 1:N) {
    A[i,1]=sample(x=c(-1,1),1,replace = T,prob = c(1/2,1/2))
  }
  
  f_11=X[,1,]%*%beta1[1:dimco]
  #f_12=X[,1,]%*%beta1[(dimco+1):(2*dimco)]
  
  for (i in 1:N) {
    A[i,2]=sample(x=c(-1,1),1,replace = T,prob = c(1/2,1/2))
  }
  
  f_21=X[,1,]%*%beta2[1:dimco] + X[,2,]%*%beta2[(dimco+1):(2*dimco)] + A[,1]*beta2[2*dimco+1]
  #f_22=X[,1,]%*%beta2[(2*dimco+2):(3*dimco+1)] + X[,2,]%*%beta2[(3*dimco+2):(4*dimco+1)] + A[,1]*beta2[4*dimco+2]
  
  
  for (i in 1:N) {
    A[i,3]=sample(x=c(-1,1),1,replace = T,prob = c(1/2,1/2))
  }
  
  f_31=X[,1,]%*%beta3[1:dimco] + X[,2,]%*%beta3[(dimco+1):(2*dimco)] +X[,3,]%*%beta3[(2*dimco+1):(3*dimco)]+A[,1]*beta3[(3*dimco)+1]+A[,2]*beta3[(3*dimco)+2]
  #f_32=X[,1,]%*%beta3[(3*dimco+3):(4*dimco+2)] + X[,2,]%*%beta3[(4*dimco+3):(5*dimco+2)] +X[,3,]%*%beta3[(5*dimco+3):(6*dimco+2)]+A[,1]*beta3[6*dimco+3]+A[,2]*beta3[6*dimco+4]
  
  for (i in 1:N) {
    A[i,4]=sample(x=c(-1,1),1,replace = T,prob = c(1/2,1/2))
  }
  
  
  f_41=X[,1,]%*%beta4[1:dimco] + X[,2,]%*%beta4[(dimco+1):(2*dimco)] +X[,3,]%*%beta4[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta4[(3*dimco+1):(4*dimco)]+A[,1]*beta4[4*dimco+1]+A[,2]*beta4[4*dimco+2] +A[,3]*beta4[4*dimco+3]
  #f_42=X[,1,]%*%beta4[(4*dimco+4):(5*dimco+3)] + X[,2,]%*%beta4[(5*dimco+4):(6*dimco+3)] +X[,3,]%*%beta4[(6*dimco+4):(7*dimco+3)]+X[,4,]%*%beta4[(7*dimco+4):(8*dimco+3)]+A[,1]*beta4[8*dimco+4]+A[,2]*beta4[8*dimco+5] +A[,3]*beta4[8*dimco+6]
  
  for (i in 1:N) {
    A[i,5]=sample(x=c(-1,1),1,replace = T,prob = c(1/2,1/2))
  }
  
  f_51=X[,1,]%*%beta5[1:dimco] + X[,2,]%*%beta5[(dimco+1):(2*dimco)] +X[,3,]%*%beta5[(2*dimco+1):(3*dimco)]+X[,4,]%*%beta5[(3*dimco+1):(4*dimco)]+X[,5,]%*%beta5[(4*dimco+1):(5*dimco)]+A[,1]*beta5[5*dimco+1]+A[,2]*beta5[5*dimco+2] +A[,3]*beta5[5*dimco+3] +A[,4]*beta5[5*dimco+4]
  #f_52=X[,1,]%*%beta5[(5*dimco+5):(6*dimco+4)] + X[,2,]%*%beta5[(6*dimco+5):(7*dimco+4)] +X[,3,]%*%beta5[(7*dimco+5):(8*dimco+4)]+X[,4,]%*%beta5[(8*dimco+5):(9*dimco+4)]+X[,5,]%*%beta5[(9*dimco+5):(10*dimco+4)]+A[,1]*beta5[10*dimco+5]+A[,2]*beta5[10*dimco+6]+A[,3]*beta5[10*dimco+7] +A[,4]*beta5[10*dimco+8]
  
  f1=cbind(f_11,f_21,f_31,f_41,f_51)
  #f2=cbind(f_12,f_22,f_32,f_42,f_52)
  
  for (tt in 1:M) {
    #print(tt)
    d[,tt] = sign(f1[,tt])
    # inner.matrix=matrix(0,N,3)
    # for (ii in 1:3){
    #   inner.matrix[,ii] = apply(V[,ii]*cbind(f1[,tt],f2[,tt]),1,sum) 
    # }
    # d[,tt]=apply(inner.matrix,1,pred)
  }
  
  datainfo = list(X=X,d=d,A=A,f1=f1)
  return(datainfo)
  
}

##----------------------------------------------------

indexstepfunc <- function(Y,N,mt,givent){
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
  return(indexstep)
}

##----------------------------------------------------

giventimefn <- function(givent,mt,N){
  
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
    beta = c(beta1,beta2,beta3)
  }
  if(givent>=mt[4]&givent<mt[5]){
    phat=matrix(1/3,N,4)
    beta = c(beta1,beta2,beta3,beta4)
  }
  if(givent>=mt[5]){
    phat=matrix(1/3,N,5)
    beta = c(beta1,beta2,beta3,beta4,beta5)
  }
  phat=as.matrix(phat)
  betatmp = beta
  datainfo = list(phat=phat,betatmp=betatmp)
  return(datainfo)
}


#######################################################################

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

transd <- function(dtrt,N,Nstep){
for (i in 1:N){
  for (j in 1:Nstep){
    if (dtrt[i,j] ==2){
      dtrt[i,j] = -1
    }
  }
}
  return(dtrt)
}
##----------------------------------------------------

tranDW <- function(dtrt,N,Nstep){
  dopt=matrix(0,N,Nstep);
  for (i in 1:N){
    for (j in 1:Nstep){
      if (dtrt[[j]][i] ==0){
        dtrt[[j]][i] = -1
      }
    }
  }
  for (j in 1:Nstep){
  dopt[,j] = dtrt[[j]]
  }
  return(dopt)
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

one_vs_all <- function(class_mat){
  n_mat = nrow(class_mat)
  recommend.trt = c()
  for (i in 1:n_mat){
    recommend.trt[i] = which.max(class_mat[i,])
  }
  return(recommend.trt)
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
# dorig <-function(u,a,c){
#   z=matrix(-1,length(u),1)
#   index <- which( u>=(c/(1+c)) )
#   z[index]= -(a/( (1+c)*u[index]-c+a ))^(a+1)
#   return(z)
# }
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
      fg[,ss]=X[,1,]%*%beta[(dimco+1):(2*dimco)]+X[,2,]%*%beta[(2*dimco+1):(3*dimco)]+A[,1]*beta[3*dimco+1]
    }
    if(ss==3){
      fg[,ss]=X[,1,]%*%beta[(3*dimco+1+1):(3*dimco+1+dimco)]+X[,2,]%*%beta[(3*dimco+1+dimco+1):(3*dimco+1+2*dimco)]
      +X[,3,]%*%beta[(3*dimco+1+2*dimco+1):(3*dimco+1+3*dimco)]+A[,1]*beta[3*dimco+1+3*dimco+1]+A[,2]*beta[3*dimco+1+3*dimco+2]
    }
    if(ss==4){
      fg[,ss]=X[,1,]%*%beta[(3*dimco+1+3*dimco+2+1):(3*dimco+1+3*dimco+2+dimco)]+X[,2,]%*%beta[(3*dimco+1+3*dimco+2+dimco+1):(3*dimco+1+3*dimco+2+dimco+dimco)]
      +X[,3,]%*%beta[(3*dimco+1+3*dimco+2+dimco+dimco+1):(3*dimco+1+3*dimco+2+dimco+dimco+dimco)]+A[,1]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+1]+A[,2]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+2]
      +A[,3]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+3]
    }
    if(ss==5){
      fg[,ss]=X[,1,]%*%beta[(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+1):(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco)]+X[,2,]%*%beta[(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+1):(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco)]
      +X[,3,]%*%beta[(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+1):(3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+dimco)]+A[,1]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+dimco+1]+A[,2]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+dimco+2]
      +A[,3]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+dimco+3]+A[,4]*beta[3*dimco+1+3*dimco+2+dimco+dimco+dimco+3+dimco+dimco+dimco+4]
    }
    if(ss==6){
      fg[,ss]=X[,1,]%*%beta[(12*dimco+1+2+3+4+1):(12*dimco+1+2+3+4+dimco)]+X[,2,]%*%beta[(12*dimco+1+2+3+4+dimco+1):(12*dimco+1+2+3+4+dimco+dimco)]
      +X[,3,]%*%beta[(12*dimco+1+2+3+4+dimco+dimco+1):(12*dimco+1+2+3+4+dimco+dimco+dimco)]+A[,1]*beta[12*dimco+1+2+3+4+dimco+dimco+dimco+1]+A[,2]*beta[12*dimco+1+2+3+4+dimco+dimco+dimco+2]
      +A[,3]*beta[12*dimco+1+2+3+4+dimco+dimco+dimco+3]+A[,4]*beta[12*dimco+1+2+3+4+dimco+dimco+dimco+4]+A[,5]*beta[12*dimco+1+2+3+4+dimco+dimco+dimco+5]
    }
  }
  return(fg)
}

# decifunc<-function(phat,X,A,beta,dimco){
#   N=length(phat[,1])
#   Nstep=length(phat[1,])
#   fg=matrix(0,N,Nstep)
#   for (ss in 1:Nstep) {
#     if(ss==1){
#       fg[,ss]=X[,1,]%*%beta[1:dimco]
#     }
#     if(ss==2){
#       fg[,ss]=X[,1,]%*%beta[(2*dimco+1):(3*dimco)]+X[,2,]%*%beta[(3*dimco+1):(4*dimco)]+A[,1]*beta[4*dimco+1]
#     }
#     if(ss==3){
#       fg[,ss]=X[,1,]%*%beta[(4*dimco+2+2*dimco+1):(4*dimco+2+2*dimco+dimco)]+X[,2,]%*%beta[(4*dimco+2+2*dimco+dimco+1):(4*dimco+2+2*dimco+2*dimco)]
#       +X[,3,]%*%beta[(4*dimco+2+2*dimco+2*dimco+1):(4*dimco+2+2*dimco+3*dimco)]+A[,1]*beta[4*dimco+2+2*dimco+3*dimco+1]+A[,2]*beta[4*dimco+2+2*dimco+3*dimco+2]
#     }
#     if(ss==4){
#       fg[,ss]=X[,1,]%*%beta[(12*dimco+6+1):(12*dimco+6+dimco)]+X[,2,]%*%beta[(12*dimco+6+dimco+1):(12*dimco+6+2*dimco)]
#       +X[,3,]%*%beta[(12*dimco+6+2*dimco+1):(12*dimco+6+3*dimco)]+X[,4,]%*%beta[(12*dimco+6+3*dimco+1):(12*dimco+6+4*dimco)]
#       +A[,1]*beta[12*dimco+6+4*dimco+1]+A[,2]*beta[12*dimco+6+4*dimco+2]+A[,3]*beta[12*dimco+6+4*dimco+3]
#     }
#     if(ss==5){
#       fg[,ss]=X[,1,]%*%beta[(20*dimco+12+1):(20*dimco+12+dimco)]+X[,2,]%*%beta[(20*dimco+12+dimco+1):(20*dimco+12+2*dimco)]
#       +X[,3,]%*%beta[(20*dimco+12+2*dimco+1):(20*dimco+12+3*dimco)]+X[,4,]%*%beta[(20*dimco+12+3*dimco+1):(20*dimco+12+4*dimco)]
#       +X[,5,]%*%beta[(20*dimco+12+4*dimco+1):(20*dimco+12+5*dimco)]
#       +A[,1]*beta[20*dimco+12+5*dimco+1]+A[,2]*beta[20*dimco+12+5*dimco+2]+A[,3]*beta[20*dimco+12+5*dimco+3]+A[,4]*beta[20*dimco+12+5*dimco+4]
#     }
#     
#   }
#   return(fg)
# }
# 
# decifunc2<-function(phat,X,A,beta,dimco){
#   N=length(phat[,1])
#   Nstep=length(phat[1,])
#   fg=matrix(0,N,Nstep)
#   for (ss in 1:Nstep) {
#     if(ss==1){
#       fg[,ss]=X[,1,]%*%beta[(dimco+1):(2*dimco)]
#     }
#     if(ss==2){
#       fg[,ss]=X[,1,]%*%beta[(4*dimco+1+1):(4*dimco+1+dimco)]+X[,2,]%*%beta[(4*dimco+1+dimco+1):((4*dimco+1+2*dimco))]+A[,1]*beta[4*dimco+2+2*dimco]
#     }
#     if(ss==3){
#       fg[,ss]=X[,1,]%*%beta[(9*dimco+4+1):(9*dimco+4+dimco)]+X[,2,]%*%beta[(9*dimco+4+dimco+1):(9*dimco+4+2*dimco)]
#       +X[,3,]%*%beta[(9*dimco+4+2*dimco+1):(9*dimco+4+3*dimco)]+A[,1]*beta[9*dimco+4+3*dimco+1]+A[,2]*beta[9*dimco+4+3*dimco+2]
#     }
#     if(ss==4){
#       fg[,ss]=X[,1,]%*%beta[(12*dimco+6+4*dimco+3+1):(12*dimco+6+4*dimco+3+dimco)]+X[,2,]%*%beta[(12*dimco+6+4*dimco+3+dimco+1):(12*dimco+6+4*dimco+3+2*dimco)]
#       +X[,3,]%*%beta[(12*dimco+6+4*dimco+3+2*dimco+1):(12*dimco+6+4*dimco+3+3*dimco)]+X[,4,]%*%beta[(12*dimco+6+4*dimco+3+3*dimco+1):(12*dimco+6+4*dimco+3+4*dimco)]
#       +A[,1]*beta[12*dimco+6+4*dimco+3+4*dimco+1]+A[,2]*beta[12*dimco+6+4*dimco+3+4*dimco+2]+A[,3]*beta[12*dimco+6+4*dimco+3+4*dimco+3]
#     }
#     if(ss==5){
#       fg[,ss]=X[,1,]%*%beta[(25*dimco+16+1):(25*dimco+16+dimco)]+X[,2,]%*%beta[(25*dimco+16+dimco+1):(25*dimco+16+2*dimco)]
#       +X[,3,]%*%beta[(25*dimco+16+2*dimco+1):(25*dimco+16+3*dimco)]+X[,4,]%*%beta[(25*dimco+16+3*dimco+1):(25*dimco+16+4*dimco)]
#       +X[,5,]%*%beta[(25*dimco+16+4*dimco+1):(25*dimco+16+5*dimco)]
#       +A[,1]*beta[25*dimco+16+5*dimco+1]+A[,2]*beta[25*dimco+16+5*dimco+2]+A[,3]*beta[25*dimco+16+5*dimco+3]+A[,4]*beta[25*dimco+16+5*dimco+4]
#     }
#     
#   }
#   return(fg)
# }

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
  u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1] ;  pos=0.01;
  for (j in 3:nstage){
    u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
  }
  for (j in 1:nstage){
    u[,j]    = u[,j] + pos;
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

## censored q-learning ##
qlearning.censored <- function(h0,h1,a,y_mat, AtRisk, weight_mat,nstage){
  yTilde          =  y_mat[ , nstage]
  stageCoef       =  list()
  weight_mat[y_mat[ ,nstage] == 0, nstage] = 0
  for (stg in nstage:2){
    stageModel  =  lm(yTilde[AtRisk[, stg]  == 1] ~ h0[[stg]][AtRisk[, stg] == 1, ] + 
                        h1[[stg]][AtRisk[, stg] == 1, ]:a[AtRisk[, stg] == 1, stg] - 1, 
                      weights = weight_mat[AtRisk[, stg] == 1, stg])
    bTemp       =  coef(stageModel)
    bTemp[is.na(bTemp)] = 0
    stageCoef[[stg]] = bTemp
    ## create psuedo outcome for the first stage
    p0Temp      =  ncol(h0[[stg]])
    pTemp       =  length(bTemp)
    weight_mat[y_mat[, stg] == 0, stg] = 0
    yTilde      =  (h0[[stg]]%*%bTemp[1:p0Temp] + abs(h1[[stg]]%*%bTemp[(p0Temp+1):pTemp])) * 
      (AtRisk[, stg] == 1) + y_mat[, stg - 1];
  }
  stg              =  1
  stage1Model      =  lm(yTilde ~ h0[[stg]] + h1[[stg]]:a[, stg] - 1, 
                         weights = weight_mat[, stg])
  bTemp            =  coef(stage1Model)
  bTemp[is.na(bTemp)] = 0
  stageCoef[[stg]] = bTemp
  return(stageCoef)
}


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

weight=function(V,A,f1,a,b,c,phat,TT,indexstep){
  N=length(phat[,1])
  grid=length(TT)
  W=matrix(0,Nstep,N)
  U=matrix(0,Nstep,N)
  logit=matrix(0,Nstep,N)
  for (ss in 1:Nstep) {
    Y.matrix = Y.matrix.gen(k,N,A[,ss])
    U[ss,] = apply(Y.matrix*cbind(f1[,ss]),1,sum) 
    logit[ss,]=logistic(a,b,c,U[ss,])
  }
  for (ii in 1:N) {
    if (indexstep[ii] == Nstep){
      for (ss in 1:indexstep[ii]) {
        W[ss,ii]= prod(logit[1:ss,ii])/prod(phat[ii,1:ss])
        #W[ss,ii]=prod((A[ii,1:min((ss+1),Nstep)]==d[ii,1:min((ss+1),Nstep)]))/prod(phat[ii,1:min((ss+1),Nstep)])  #min(lum[1:ss,ii])
      } 
    } else if (indexstep[ii] < Nstep){
      for (ss in 1:indexstep[ii]) {
        logit[(indexstep[ii]+1),ii] = 1
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
  #f2=decifunc2(phat,X,A,beta,dimco)
  
  Wsi=weight(V,A,f1,a,b,c,phat,TT,indexstep);
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


cal_survival_and_density_functions <- function(X, Z, status, X_new, Z_new, fitting_method = "Cox"){
  
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


######################################################################

######################################################################
#######################################################################
#######################################################################

crossvalidate <- function(train,Kfold,lambda,b){
  estsurv_value = c()
  for (kkk in 1:Kfold){
    corss_N = (1-1/Kfold)*N
    index = sample(1:N,corss_N,replace = F)
    corss_X <- train$X[index,,] 
    corss_A <- train$A[index,] 
    corss_Y <- Y[index]
    corss_delta <- delta[index]
    indexstep = indexstepfunc(corss_Y ,corss_N,mt,givent)
    S=timegrid(corss_Y, corss_delta ,givent)
    grid=length(S)
    
    giventmp = giventimefn(givent,mt,corss_N)
    phat = giventmp$phat
    beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    lambda = lambda
    a=1; b=b;
    grid=length(S)
    l0_qt =matrix(0,grid,1);
    for (ss in 1:grid) {
      l0_qt[ss]=1-sum((corss_Y==S[ss])*corss_delta)/sum((corss_Y>=S[ss]))
    }
    f_norm = -sum(log(l0_qt))/lambda - 1
    V_norm = norm(V)
    c = -(f_norm +  V_norm)
    
    betatmp = beta; betatmp[1:6] = 0.5
    
    betahat = optim(betatmp,Qfunc,lambda=lambda,X=corss_X,V=V,A=corss_A,a=a,b=b,c=c,
                    phat=phat,TT=S,Y=corss_Y,delta=corss_delta,dimco=dimco,mt=mt,indexstep=indexstep,method="BFGS",
                    control = list(maxit = 5, trace = TRUE, REPORT = 500))$par
    
    validate_N = (1/Kfold)*N
    validate_X <- train$X[-index,,] 
    validate_A <- train$A[-index,] 
    validate_Y <- Y[-index]
    validate_delta <- delta[-index]
    validate_indexstep = indexstepfunc(validate_Y ,validate_N ,mt,givent)
    
    S=timegrid(validate_Y, validate_delta,givent)
    grid=length(S)
    
    giventmp = giventimefn(givent,mt,validate_N)
    phat = giventmp$phat
    Nstep=length(phat[1,])
    
    fhat1=decifunc(phat,validate_X,validate_A,betahat,dimco)
    #fhat2=decifunc2(phat,validate_X,validate_A,betahat,dimco)
    
    # dhat = matrix(0,validate_N,Nstep)
    # for (tt in 1:Nstep){
    #   inner.matrix=matrix(0,validate_N,3)
    #   for (ii in 1:k){
    #     inner.matrix[,ii] = apply(V[,ii]*cbind(fhat1[,tt],fhat2[,tt]),1,sum) 
    #   }
    #   dhat[,tt]=apply(inner.matrix,1,pred)
    # }
    
    dhat = matrix(0,N,Nstep)
    for (tt in 1:Nstep){
      inner.matrix=matrix(0,N,k)
      for (ii in 1:k){
        inner.matrix[,ii] = apply(V[,ii]*cbind(fhat1[,tt]),1,sum)
      }
      dhat[,tt]=apply(inner.matrix,1,pred)
    }
    dhat = transd(dhat,N,Nstep)
    
    
    KMtce=KMcurve_wtc(grid,S,validate_Y,validate_delta,validate_A,
                      dhat,phat,mt,validate_indexstep)
    
    estsurv_value[kkk] = min(KMtce,na.rm = T)
  }
  return(mean(estsurv_value))
}


#######################################################################

crossvalidate_org <- function(train,Kfold,lambda){
  estsurv_value = c()
  for (kkk in 1:Kfold){
    corss_N = (1-1/Kfold)*N
    index = sample(1:N,corss_N,replace = F)
    corss_X <- train$X[index,,] 
    corss_A <- train$A[index,] 
    corss_Y <- Y[index]
    corss_delta <- delta[index]
    indexstep = indexstepfunc(corss_Y ,corss_N,mt,givent)
    S=timegrid(corss_Y, corss_delta ,givent)
    grid=length(S)
    
    giventmp = giventimefn(givent,mt,corss_N)
    phat = giventmp$phat
    beta = giventmp$betatmp
    Nstep=length(phat[1,])
    
    lambda = lambda
    a=1; b=1;
    grid=length(S)
    l0_qt =matrix(0,grid,1);
    for (ss in 1:grid) {
      l0_qt[ss]=1-sum((corss_Y==S[ss])*corss_delta)/sum((corss_Y>=S[ss]))
    }
    f_norm = -sum(log(l0_qt))/lambda - 1
    V_norm = norm(V)
    c = -(f_norm +  V_norm)
    
    betahat = optim(beta,Qfunc,lambda=lambda,X=corss_X,V=V,A=corss_A,a=a,b=b,c=c,
                    phat=phat,TT=S,Y=corss_Y,delta=corss_delta,dimco=dimco,mt=mt,indexstep=indexstep,method="BFGS",
                    control = list(maxit = 5, trace = TRUE, REPORT = 500))$par
    
    validate_N = (1/Kfold)*N
    validate_X <- train$X[-index,,] 
    validate_A <- train$A[-index,] 
    validate_Y <- Y[-index]
    validate_delta <- delta[-index]
    validate_indexstep = indexstepfunc(validate_Y ,validate_N ,mt,givent)
    
    S=timegrid(validate_Y, validate_delta,givent)
    grid=length(S)
    
    giventmp = giventimefn(givent,mt,validate_N)
    phat = giventmp$phat
    Nstep=length(phat[1,])
    
    # fhat1=decifunc(phat,validate_X,validate_A,betahat,dimco)
    # fhat2=decifunc2(phat,validate_X,validate_A,betahat,dimco)
    # 
    # dhat = matrix(0,validate_N,Nstep)
    # for (tt in 1:Nstep){
    #   inner.matrix=matrix(0,validate_N,3)
    #   for (ii in 1:k){
    #     inner.matrix[,ii] = apply(V[,ii]*cbind(fhat1[,tt],fhat2[,tt]),1,sum) 
    #   }
    #   dhat[,tt]=apply(inner.matrix,1,pred)
    # }
    # 
    value = Qfunc(betahat,lambda,validate_X,V,validate_A,a,b,c,phat,S,validate_Y,validate_delta,dimco,mt,validate_indexstep)
    estsurv_value[kkk] = value 
  }
  return(mean(estsurv_value))
}



#######################################################################

# 
# multiclassfier <- function(a){
#   b = c()
#   for (i in 1:N){
#     
#     if (a[i,3] == 1) 
#     { 
#       b[i]= 3
#     }
#     else if (a[i,3] == 0 & a[i,2] == 1 ){
#       b[i]= 2
#     } else {
#       b[i]= 1
#     }
#     
#   }
#   return(b)
# }
# 
# 
# multiclassfier2 <- function(a){
#   b = c()
#   for (i in 1:N){
#     
#     if (a[i,3] == 1) 
#     { 
#       b[i]= 3
#     }
#     else if (a[i,3] == -1 & a[i,2] == 1 ){
#       b[i]= 2
#     } else {
#       b[i]= 1
#     }
#     
#   }
#   return(b)
# }

# adjust_func_DW <- function(N,Nstep,Ttime,C){
#   
#   n = N
#   nstage = Nstep
#   y =  matrix(0, nrow=n, ncol=nstage)
#   u =  matrix(0, nrow=n, ncol=nstage)
#   delta_mat = matrix(0, nrow=n, ncol=nstage)
#   
#   ##--------------------calculate interval survival outcome---------------------
#   y[,1] = Ttime$Tfy[,1]
#   
#   if (nstage != 1){
#     for (i in 2:nstage){
#       y[,i] =  Ttime$Tfy[,i]
#     }
#   }
#   ##------------calculate min cumulated survival outcome and censored time-----------
#   u[,1]    = pmin(y[,1], C) 
#   u[,2]    = pmin(apply(y[,1:2], 1, sum), C) -  u[,1] ; pos = 0.05
#   for (j in 3:nstage){
#     u[,j]    = pmin(apply(y[,1:j], 1, sum), C) - apply(u[,1:(j-1)], 1, sum)
#   }
#   for (j in 1:nstage){
#     u[,j]    = u[,j]
#   }
#   
#   ##-----------------------censored indicator for each interval------------------------
#   delta_mat[,1] = y[,1] < C
#   if (nstage > 1){
#     for (kkk in 2:nstage){
#       delta_mat[,kkk] = (apply(y[,1:kkk], 1, sum) < C)
#     }
#   }
#   
#   datainfo = list(u=u,y=y,delta_mat=delta_mat)
#   return(datainfo)
# }
# 
