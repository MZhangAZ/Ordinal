####################################
# This is a script for type I error comparison
# by using simulated haplotype in SKAT (add only rare variants)
# I consider covariate-adjusting (one continuous variable and one binary variable)
# I choose SCORE-Seq for comparison
# I merge two categories into one
# I use multiple arrary job and for loop in R
#
# For multiple arrary, I set #PBS -J 1-1000
# For R forloop, I set for (k in 1:1000)
# total replicates are 10e06
# I write single p-value set, denoted as "seedINV=seedNum", and also write 1000-set p-values, denoted as "seed=NAI"
#
# 
# date: 10/12/2016 Miao Zhang
####################################

################
# read the haplotype data and load packages
################
library(SKAT)

data(SKAT.haplotypes)
names(SKAT.haplotypes)

attach(SKAT.haplotypes)

library(HardyWeinberg)

library(VGAM)
library(numDeriv)
##################################
# function to compute observed MAF
# 
# genoCol: a vector of additive coding
##################################

MAF<- function(genoCol){
  count0= sum(genoCol==0)
  count1= sum(genoCol==1)
  count2= sum(genoCol==2)
  MAFtable= c(count0, count1, count2)
  p= maf(MAFtable)
  return(p)
}

################
# create additive genotype data
# randomly draw 2 haplotype and add them together
################


################
# function to draw genotype data 
# rare=T: only draw rare variants
# rare=F: common + rare variants
################


#####################################################
# generate genotype additive coding
# based on SKAT haplotype data base
# specify subregion to generate genotype data 
# For commone + rare case, remove those sample SNPs with observed MAF=0
# For rare only case, remove those sample SNPs with observed MAF=0 & MAF>=0.05
# 
# n: sample size
# length: length of subregion
# rare: T only rare variants; F common + rare variants
# update date: 07/23/2016  Miao Zhang
#####################################################


genotype= function(n, length=5000, rare=T){
  
  startPoint= sample(1:60000, 1) # draw start SNP position
  BPidx= which(SKAT.haplotypes$SNPInfo$CHROM_POS>=startPoint & SKAT.haplotypes$SNPInfo$CHROM_POS<=(startPoint+length) )
  
  sampleSNP=SKAT.haplotypes$SNPInfo[BPidx,]
  
  snpIdx= which(SKAT.haplotypes$SNPInfo$CHROM_POS %in% sampleSNP$CHROM_POS) # denote snp index (column) in haplotype file
  
  G= matrix(0, n, length(snpIdx)) # create empty genotype matrix, which snp number is not fixed
  
  G= SKAT.haplotypes$Haplotype[sample(1:10000, n), snpIdx]+SKAT.haplotypes$Haplotype[sample(1:10000, n), snpIdx]
  
  if (rare==T) {
    sMAF= apply(G, 2, MAF)    # sample MAF estimated from sample
    rareIdx= which(sMAF>0 & sMAF<0.05) # keep observed MAF>0 and <0.05 as rare variants
    G= G[, rareIdx]
    return(G)
    
  } else {
    sMAF= apply(G, 2, MAF)    # sample MAF estimated from sample
    rareIdx= which(sMAF>0) # keep observed MAF>0 as common + rare variants
    G= G[, rareIdx]
    return(G)
  }
  
}

################################################
# specify parameter values
# 
# n: sample size
# length: length of subregion in SKAT haplotype data
# rep: simulation times
#################################################

# NAI is index for multiple array job
AI <- Sys.getenv("PBS_ARRAY_INDEX")
NAI <- as.numeric(AI)

index01=1000
index02=2000

n= index01 
length= index02

########################
# simulation
########################

pvalue=matrix(0, 1000, 12)

for (k in 1:1000) {

# set seed
seedNum= (NAI-1)*1000+k # index is index for multiple R scripts
set.seed(seedNum)   
  
  
# generate latent variable
e=rlogis(n)  # draw error term from standard logistic distribution
geno=genotype(n, length, rare=T) # generate genotype data
x1= rnorm(n) # generate continuous covariate
x2= sample(c(0,1), n, replace=T) #generate binary covariate
b1= 1.5 # coef for continuous variable
b2= 0.6 # coef for binary variable
yl= x1*b1 + x2*b2 + e

# set threshold
orderIdx= order(yl)

a1= (yl[orderIdx])[ceiling((n/3)*1)]
a2= (yl[orderIdx])[ceiling((n/3)*2)]


# based on the threshold to generate ordinal response
# index latent variable based on cut-off
y=numeric(n)
idxLatent1= which(yl<=a1)
idxLatent2= which(yl>a1 & yl<=a2)
idxLatent3= which(yl>a2)
y[idxLatent1]= 1 
y[idxLatent2]= 2 
y[idxLatent3]= 3 

# count number for each category
n1=length(idxLatent1)
n2=length(idxLatent2)
n3=length(idxLatent3)


data=data.frame(y)

# estimate the parameters under the null
fit= vglm( y ~x1+x2 , family=cumulative(parallel=TRUE), data=data)

# summary(fit)

ahat1=as.numeric(coef(fit)[1])
ahat2=as.numeric(coef(fit)[2])
bhat1=as.numeric(coef(fit)[3])
bhat2=as.numeric(coef(fit)[4])

beta= as.vector(-c(bhat1, bhat2))  # should be careful about the sign of beta 
########################################

# score statistic

X= as.matrix( cbind(x1, x2) )

psi1= 1 - (1/(1+exp(-ahat1 + X[idxLatent1,]%*%beta)))
psi2= 1 - (1/(1+exp(-ahat2 + X[idxLatent2,]%*%beta))) - (1/(1+exp(-ahat1  + X[idxLatent2,]%*%beta)))
psi3=   - (1/(1+exp(-ahat2 + X[idxLatent3,]%*%beta)))

Q= apply(geno,1,sum) # aggregate multiple SNPs, here assume weight equal to 1


#    idx1=which(Q[,2]==1)
#    idx2=which(Q[,2]==2)
#    idx3=which(Q[,2]==3)

# S1=sum(Q[idxLatent1])
# S2=sum(Q[idxLatent2])
# S3=sum(Q[idxLatent3])

# SS1=sum( (Q[idxLatent1])^2 )
# SS2=sum( (Q[idxLatent2])^2 )
# SS3=sum( (Q[idxLatent3])^2 )

score= as.numeric( -( t(Q[idxLatent1])%*%psi1 + t(Q[idxLatent2])%*%psi2 + t(Q[idxLatent3])%*%psi3 ) ) # original score 

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

OrdinalP = 1-pchisq(testStat, 1)

#######################################
# SCORE-Seq
#
#######################################

# merge two categories into one
# group1 (y=1) is control group
# group3 (y=3) is case group
# group2 (y=2) is case group
bidata=numeric(n)
bidata[idxLatent1]=0
bidata[idxLatent3]=1

bidata[idxLatent2]=1
bidata= as.vector(bidata)

subid= seq(1,n) # create subject id
snpid= seq(1, dim(geno)[2]) # create snp id

# create variables for score-seq input
ScoreGeno= cbind(snpid, t(geno))
ScoreMap= cbind("Gene", snpid)
ScorePheno= cbind(subid, bidata, x1, x2)

phenoFile= paste("/12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqPheno-seed=", seedNum, "-sample=", index01, "-length=", index02,".txt", sep="")
genoFile= paste("/12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqGeno-seed=", seedNum, "-sample=", index01, "-length=", index02,".txt", sep="")
mapFile= paste("/12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqMap-seed=", seedNum, "-sample=", index01, "-length=", index02,".txt", sep="")

write.table(ScorePheno, file=phenoFile,  row.names = F, col.names = F, sep = "\t")
write.table(ScoreGeno, file=genoFile,  row.names = F, col.names = F, sep = "\t")
write.table(ScoreMap, file=mapFile,  quote = F, row.names = F, col.names = F, sep = "\t")

command= paste("./SCORE-Seq -pfile /12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqPheno-seed=", seedNum, "-sample=", index01, "-length=", index02, ".txt -gfile /12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqGeno-seed=", seedNum, "-sample=", index01, "-length=", index02, ".txt -mfile /12kx/gsfs2/home/u11/miaozhang/EREC/txt/Score-seqMap-seed=", seedNum, "-sample=", index01, "-length=", index02, ".txt -ofile /12kx/gsfs2/home/u11/miaozhang/EREC/txt/rareP-seed=", seedNum, "-sample=", index01, "-length=", index02, ".txt -MAF 0.05 -MAC 5 -resample -1", sep="")

system(command)

seqout= paste("/12kx/gsfs2/home/u11/miaozhang/EREC/txt/rareP-seed=", seedNum, "-sample=", index01, "-length=", index02, ".txt", sep="")

result= as.numeric(strsplit(readLines(seqout), "\t")[[2]][2:12])
result=  c(OrdinalP, result) 
file=paste("/12kx/gsfs2/home/u11/miaozhang/result/TypeI-seedINV=", seedNum, "-sample=", index01, "-length=", index02, ".Rdata", sep="")
save(result, file=file)

pvalue[k,]= result

#############################
# delete Score-seq txt files
#############################
file.remove(phenoFile)
file.remove(genoFile)
file.remove(mapFile)
file.remove(seqout)

}

file=paste("/12kx/gsfs2/home/u11/miaozhang/result/TypeI-seed=", NAI, "-sample=", index01, "-length=", index02, ".Rdata", sep="")
save(pvalue, file=file)

