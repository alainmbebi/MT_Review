library(openxlsx)
library(data.table)
library(matrixcalc)

rm(list=ls()) 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the joint estimation for B and Omega 
cmor<-function(lambda1, lambda2, lambda3, lambda4, X, Y, eps=1e-6, zeta=1e-8, tol.out=1e-6){
  #----------------------
  
  # t=0 initialization
  
  Sigma_0=diag(ncol(Y))                           
  Omega_0=diag(ncol(Y))                          
  invSigma0=solve(Sigma_0)
  invOmega0=solve(Omega_0)
  P0=t(chol(lambda1*Omega_0 +  lambda2*Sigma_0))  
  
  #----------------------
  SVD_P0=svd(P0)
  U2_0=SVD_P0$u
  Sigma2_0=diag(SVD_P0$d)
  Sigma2vect_0=c(SVD_P0$d)
  tV2_0=SVD_P0$v
  
  #----------------------
  SVD_X=svd(X)
  U1=SVD_X$u
  Sigma1=diag(SVD_X$d)
  Sigma1vect=c(SVD_X$d)
  tV1=t(SVD_X$v)
  
  #----------------------
  Btilde0= matrix(0, nrow=nrow(X), ncol=ncol(Y))
  S0=tV1%*%t(X)%*%Y%*%U2_0
  diagSigma1=(Sigma1vect)^2 #length p
  diagSigma2_0=(Sigma2vect_0)^2 #length s
  for(i in 1:(nrow(X))){
    for(j in 1:(ncol(Y))){
      Btilde0[i,j]<-S0[i,j]/(diagSigma1[i] + diagSigma2_0[j])
    }
  }
  B0=t(tV1)%*%Btilde0%*%t(U2_0)
  
  residual.cov_0=crossprod(Y-X%*%B0)/nrow(Y)
  obj_0=matrix.trace(residual.cov_0%*%invOmega0)-nrow(Y)*determinant(invOmega0, logarithm=TRUE)$mod[1] +
    lambda1*matrix.trace(t(B0)%*%B0) + lambda2*matrix.trace(invSigma0%*%t(B0)%*%B0) -ncol(X)*determinant(invSigma0, logarithm=TRUE)$mod[1] +
    lambda3*matrix.trace(invOmega0) +  lambda4*matrix.trace(invSigma0)
  #----------------------
  Sigma_1=(lambda2*t(B0)%*%B0 + lambda4*diag(ncol(Y)))/nrow(X)
  Omega_1=(t(Y-X%*%B0)%*%(Y-X%*%B0) + lambda3*diag(ncol(Y)))/nrow(Y)
  invSigma1=solve(Sigma_1)
  invOmega1=solve(Omega_1)
  P1=t(chol(lambda1*Omega_1 +  lambda2*Sigma_1)) 
  
  #----------------------
  SVD_P1=svd(P1)
  U2_1=SVD_P1$u
  Sigma2_1=diag(SVD_P1$d)
  Sigma2vect_1=c(SVD_P1$d)
  tV2_1=SVD_P1$v
  
  #----------------------
  Btilde1= matrix(0, nrow=nrow(X), ncol=ncol(Y))
  S1=tV1%*%t(X)%*%Y%*%U2_1
  diagSigma1=(Sigma1vect)^2 #length p
  diagSigma2_1=(Sigma2vect_1)^2 #length s
  for(i in 1:(nrow(X))){
    for(j in 1:(ncol(Y))){
      Btilde1[i,j]<-S1[i,j]/(diagSigma1[i] + diagSigma2_1[j])
    }
  }
  B1=t(tV1)%*%Btilde1%*%t(U2_1)
  
  residual.cov_1=crossprod(Y-X%*%B1)/nrow(Y)
  obj_1=matrix.trace(residual.cov_1%*%invOmega1)-nrow(Y)*determinant(invOmega1, logarithm=TRUE)$mod[1] +
    lambda1*matrix.trace(t(B1)%*%B1) + lambda2*matrix.trace(invSigma1%*%t(B1)%*%B1) -ncol(X)*determinant(invSigma1, logarithm=TRUE)$mod[1] +
    lambda3*matrix.trace(invOmega1) +  lambda4*matrix.trace(invSigma1)
  #----------------------
  
  tcontjmor=0
  conv_crit=tol.out*sum(diag(crossprod(Y))/nrow(Y))
  #while((sum(abs(B1))<=sum(abs(B0)))==TRUE  && tcontjmor<=itermax){ 
  while(((obj_0-obj_1)>conv_crit)==TRUE  && (min(diag(residual.cov_0)) > eps) && tcontjmor<=itermax){ 
    
    Sigma_0=Sigma_1
    Omega_0=Omega_1  
    invSigma0=invSigma1                                         
    invOmega0=invOmega1   
    P0=P1
    #----------------------
    #----------------------
    SVD_P0=SVD_P1
    U2_0=U2_1
    Sigma2vect_0=Sigma2vect_1
    tV2_0=tV2_1
    
    #----------------------
    Btilde0= Btilde1
    S0=S1
    B0=B1
    
    #----------------------
    residual.cov_0=residual.cov_1
    obj_0=obj_1
    
    #----------------------
    #----------------------
    
    SVD_P0=svd(P0)
    U2_0=SVD_P0$u
    Sigma2_0=diag(SVD_P0$d)
    Sigma2vect_0=c(SVD_P0$d)
    tV2_0=SVD_P0$v
    
    #---------------------- 
    Sigma_1=(lambda2*t(B0)%*%B0 + lambda4*diag(ncol(Y)))/nrow(X)
    Omega_1=(t(Y-X%*%B0)%*%(Y-X%*%B0) + lambda3*diag(ncol(Y)))/nrow(Y)
    invSigma1=solve(Sigma_1)
    invOmega1=solve(Omega_1)
    P1=t(chol(lambda1*Omega_1 +  lambda2*Sigma_1)) 
    
    #----------------------
    SVD_P1=svd(P1)
    U2_1=SVD_P1$u
    Sigma2_1=diag(SVD_P1$d)
    Sigma2vect_1=c(SVD_P1$d)
    tV2_1=SVD_P1$v
    
    #----------------------
    Btilde1= matrix(0, nrow=nrow(X), ncol=ncol(Y))
    S1=tV1%*%t(X)%*%Y%*%U2_1
    diagSigma1=(Sigma1vect)^2 #length p
    diagSigma2_1=(Sigma2vect_1)^2 #length s
    for(i in 1:(nrow(X))){
      for(j in 1:(ncol(Y))){
        Btilde1[i,j]<-S1[i,j]/(diagSigma1[i] + diagSigma2_1[j])
      }
    }
    B1=t(tV1)%*%Btilde1%*%t(U2_1)
    
    residual.cov_1=crossprod(Y-X%*%B1)/nrow(Y)
    obj_1=matrix.trace(residual.cov_1%*%invOmega1)-nrow(Y)*determinant(invOmega1, logarithm=TRUE)$mod[1] +
      lambda1*matrix.trace(t(B1)%*%B1) + lambda2*matrix.trace(invSigma1%*%t(B1)%*%B1) -ncol(X)*determinant(invSigma1, logarithm=TRUE)$mod[1] +
      lambda3*matrix.trace(invOmega1) +  lambda4*matrix.trace(invSigma1)
    #----------------------
    
    #----------------------
    tcontjmor <- sum(tcontjmor, 1)
    
    print(tcontjmor)
    
  }
  return(list(B1, invOmega1, invSigma1, tcontjmor))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in the MOR
CV_cmor<-function(lambda1, lambda2, lambda3, lambda4, X, Y, kfold=5) {
  lambda1 = sort(lambda1)
  lambda2 = sort(lambda2)
  lambda3 = sort(lambda3)
  lambda4 = sort(lambda4)
  n = nrow(Y)
  listt=list(lambda1, lambda2, lambda3,lambda4)
  dimll=length(lambda1)*length(lambda2)*length(lambda3)*length(lambda4)
  listindex=rep(list(rep(0,length(listt))),dimll )
  
  #CV_errors = array(0, c(dimll, kfold))
  CV_errors=rep(0,length(listindex))
  ind=sample(n)
  for (k in 1:kfold){
    
    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    
    # validation set
    Y.valid = Y[leave.out,,drop=FALSE ]
    X.valid = X[leave.out,,drop=FALSE ]
    
    # loop over all tuning parameters
    
    for (m in 1:dimll) {
      for(i in 1:length(lambda1)) {
        for(j in 1:length(lambda2)) {
          for (t in 1:length(lambda3)) {
            for (l in 1:length(lambda4)) {
              # compute the joint penalized regression matrix estimate with the training set
              Estimes.train = cmor(lambda1[i], lambda2[j], lambda3[t], lambda4[l], X.train, Y.train) 
              B_hat.train = Estimes.train[[1]]
              invOmega_hat_train =  Estimes.train[[2]]
              invSigma_hat_train =  Estimes.train[[3]]
              listindex[[m]]=c(lambda1[i], lambda2[j], lambda3[t], lambda4[l])
              CV_errors[m]= CV_errors[m] + mean((Y.valid-X.valid%*%B_hat.train)^2)
            }
          }
        }
      }
    }   
  }
  
  # determine optimal tuning parameters
  AVG = mean(CV_errors) 
  bestindexAVG=which.min(CV_errors)
  error = min(CV_errors)
  opt = listindex[[bestindexAVG]]
  opt.lam1 = opt[1]
  opt.lam2 = opt[2]
  opt.lam3 = opt[3]
  opt.lam4 = opt[4]
  
  
  # return best lambdaO, lambdaB and other meaningful values
  return(list(lambda1=opt.lam1, lambda2=opt.lam2, lambda3=opt.lam3, lambda4=opt.lam4, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Global parameters
itermax=50
eps=10^-5 

### Loading trait matrix
Y <- read.xlsx("../Supplementary Table 3_Trait_Matrix.xlsx", rowNames = TRUE)

### CV A Japonica-Japonica
Y <- Y[Y$Population == "Indica",  -c(1)]
target_accessions <- rownames(Y)
for(i in 1:length(target_accessions)){
  target_accessions[i] <- paste(cbind(target_accessions[i], target_accessions[i]), collapse = "_")
}
n_traits <- ncol(Y)

# Reading the markers
X <- fread("../Supplementary File 1_markers.txt")
colnames(X)[1] <- "Marker"
X <- X[X$Marker %in% target_accessions, -c(1)]

### Setting up Cross Validation approach 10
### Define folds
k_folds <- 5

# Create folds of equal size
set.seed(1)
folds <- cut(seq(1, nrow(Y)), breaks=k_folds, labels=FALSE)

# Retaining true test values to calculate prediction accuracy
true_fold_test_list <- list()

# List to allocate predictions for the test sets of each fold
prediction_list <- list()

for(i in 1:k_folds){
  
  sprintf("Current fold: %s", i)
  
  ### Getting the test and data partitions for this particular fold
  testIndexes <- which(folds==i,arr.ind=TRUE)
  
  X.TRN <- X[-testIndexes,]
  X.TRN <- as.matrix(X.TRN)
  Y.TRN <- Y[-testIndexes,]
  Y.TRN <- as.matrix(Y.TRN)
  
  Y.TST <- Y[testIndexes, ]
  Y.TST <- as.matrix(Y.TST)
  X.TST <- X[testIndexes]
  X.TST <- as.matrix(X.TST)
  
  true_fold_test_list[[i]] <- Y.TST
  
  # Hypeparameter tuning for this fold --> feeding the training partition
  L.cmor.opt=CV_cmor(lambda1 = 2, lambda2 = 2, lambda3 = 2,
                     lambda4 = 2, X.TRN, Y.TRN, kfold=5)
  lam1.opt.cmor=L.cmor.opt[[1]]
  lam2.opt.cmor=L.cmor.opt[[2]]
  lam3.opt.cmor=L.cmor.opt[[3]]
  lam4.opt.cmor=L.cmor.opt[[4]]
  
  Final_cmor_est_CV=cmor(lam1.opt.cmor, lam2.opt.cmor, lam3.opt.cmor, lam4.opt.cmor, X.TRN, Y.TRN)
  Bhat_cmor=Final_cmor_est_CV[[1]]
  
  Yhat_cmor=X.TST%*%Bhat_cmor
  
  prediction_list[[i]] <- Yhat_cmor
}
