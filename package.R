##########################################################
# The function to conduct gene-based association test
# 
# currently, only allow 2 covariates and 3 categories
#
# X: covariates matrix
# G: genotype matrix
# y: response vector
# smallAdjust: indecator to specify if resampling is applied due to small sample size
# per: permutation times for small sample size adjustment
##########################################################

function(X, G, y, smallAdjust=F, per=1000){
  
  # based on the threshold to generate ordinal response
  # index latent variable based on cut-off
  idxLatent1= which(y==1)
  idxLatent2= which(y==2)
  idxLatent3= which(y==3)

  # count number for each category
  n1=length(idxLatent1)
  n2=length(idxLatent2)
  n3=length(idxLatent3)
  n=n1+n2+n3
  
  dat= cbind(y, X)
  colnames(dat)= c("y", "x1", "x2")
  
  # estimate the parameters under the null
  fit= vglm( y ~ x1+x2 , family=cumulative(parallel=TRUE), data=dat)
  
  # summary(fit)
  
  ahat1=as.numeric(coef(fit)[1])
  ahat2=as.numeric(coef(fit)[2])
  bhat1=as.numeric(coef(fit)[3])
  bhat2=as.numeric(coef(fit)[4])
  
  beta= as.vector(-c(bhat1, bhat2))  # should be careful about the sign of beta 
  ########################################
  ####################
  # score statistic
  ####################
  
  psi1= 1 - (1/(1+exp(-ahat1 + X[idxLatent1,]%*%beta)))
  psi2= 1 - (1/(1+exp(-ahat2 + X[idxLatent2,]%*%beta))) - (1/(1+exp(-ahat1  + X[idxLatent2,]%*%beta)))
  psi3=   - (1/(1+exp(-ahat2 + X[idxLatent3,]%*%beta)))
  
  Q= apply(G,1,sum) # aggregate multiple SNPs, here assume weight equal to 1
  
  score= as.numeric( -( t(Q[idxLatent1])%*%psi1 + t(Q[idxLatent2])%*%psi2 + t(Q[idxLatent3])%*%psi3 ) ) # original score 
  
  
  if (smallAdjust==T) {
  
  ###########################
  # Fisher information matrix
  ###########################
  
  #########
  # J=1
  #########
  Fisher1=matrix(0, 4, 4)
  for (i in idxLatent1) {
    f1= function(x) {
      log(1/( exp(x[1]*Q[i]+x[2]*X[i,1]+x[3]*X[i,2]-x[4])+1 ) )
    }
    Fisher1=Fisher1 + hessian(func=f1, x=c(0, -bhat1, -bhat2, ahat1))
  }
  
  # add extra zero row and column to exisiting matrix
  F1=matrix(0, 5, 5)
  upper1=c( Fisher1[upper.tri(Fisher1, diag = T)], rep(0,5) )
  
  F1[upper.tri(F1, diag=TRUE)] <- upper1
  
  lower1=c( Fisher1[lower.tri(Fisher1, diag = F)] )
  
  F1[1:4,1:4][lower.tri(F1[1:4,1:4], diag=F)] <- lower1
  
  #########
  # J=2
  #########
  Fisher2=matrix(0, 5, 5)
  for (i in idxLatent2) {
    f1= function(x) {
      log(1/( exp(x[1]*Q[i]+x[2]*X[i,1]+x[3]*X[i,2]-x[5])+1 ) - 1/( exp(x[1]*Q[i]+x[2]*X[i,1]+x[3]*X[i,2]-x[4])+1 )   )
    }
    Fisher2=Fisher2 + hessian(func=f1, x=c(0, -bhat1, -bhat2, ahat1, ahat2))
  }
  
  #########
  # J=3
  #########
  Fisher3=matrix(0, 4, 4)
  for (i in idxLatent3) {
    f1= function(x) {
      log(1 - 1/( exp(x[1]*Q[i]+x[2]*X[i,1]+x[3]*X[i,2]-x[4])+1 )   )
    }
    Fisher3=Fisher3 + hessian(func=f1, x=c(0, -bhat1, -bhat2, ahat2))
  }
  
  
  F3=matrix(0, 5, 5)
  upper3=c( Fisher3[1:3,1:3][upper.tri(Fisher3[1:3,1:3], diag = T)], rep(0,4) )
  
  F3[1:4,1:4][upper.tri(F3[1:4,1:4], diag=T)] <- upper3
  
  lower3=c( Fisher3[1:3,1:3][lower.tri(Fisher3[1:3,1:3], diag = F)] )
  F3[1:3,1:3][lower.tri(F3[1:3,1:3], diag=F)] <- lower3
  
  F3[1:3,5]= Fisher3[1:3,4]
  F3[5,5]= Fisher3[4,4]
  F3[5,1:3]= Fisher3[1:3,4]
  ##################################
  # sum together
  
  Fisher= -( F1+Fisher2+F3 )
  
  InvFisher= solve(Fisher)
  
  
  testStat= score^2*InvFisher[1,1]
  # testStat
  
  pvalue = 1-pchisq(testStat, 1)

  } else {
  
  #####################
  # permutation
  #####################  
  
  permutation= function(m){
    
    idx= sample(seq(1:n))
    Qper= Q[idx]
    
    scorePer= -( -( t(Qper[idxLatent1])%*%psi1 + t(Qper[idxLatent2])%*%psi2 + t(Qper[idxLatent3])%*%psi3 ) )
    
    return(scorePer)
  }
  
  scoretest= sapply(1:per, permutation)
  
  pvalue= mean(abs(scoretest)>=abs(score))
  }
   
  return(pvalue)
  
}