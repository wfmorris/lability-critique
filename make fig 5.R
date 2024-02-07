# make panels of Fig 5 for Morris & Doak manuscript entitled 
# "What is adaptive demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()

# *************************************************************
# EXPLORE EFFECT OF INCREASING ENVTAL VARIATION CENTERED ON DIFFERENT 
# PORTIONS OF THE UNIMODAL (quadratic logistic) SURVIVAL FUNCTION FOR 
# THE SLOW LIFE HISTORY: BUT ADJUST OTHER 
# VITAL RATES TO GET LAM.1=1 IN A CONSTANT ENVIRONMENT

# mean vital rates for the slow life history at mean(z)=0
sj=.1
sa=.95
f=0.5
a.sa=log(sa/(1-sa)) # adult survival intercept at mean(z)=0

b=-0.5  # curvature parameter for quadratic logistic

Mz=fs=matrix(0,1,7)

# code to find values of fertility that produce lam.1=1 at the mean environment
# as the mean environment changes but the survival function is fixed

Mz[1]=3.5
(sa.low=1/(1+exp(-a.sa-b*(Mz[1])^2))) # unimodal survival function
fs[1]=9.601; A=matrix(c(0,sj,fs[1],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[2]=3
(sa.low=1/(1+exp(-a.sa-b*(Mz[2])^2))) 
fs[2]=8.257; A=matrix(c(0,sj,fs[2],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[3]=2.5
(sa.low=1/(1+exp(-a.sa-b*(Mz[3])^2))) # sa unimodal
fs[3]=5.45; A=matrix(c(0,sj,fs[3],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[4]=2
(sa.low=1/(1+exp(-a.sa-b*(Mz[4])^2))) # sa unimodal
fs[4]=2.8; A=matrix(c(0,sj,fs[4],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[5]=1.5
(sa.low=1/(1+exp(-a.sa-b*(Mz[5])^2))) # sa unimodal
fs[5]=1.395; A=matrix(c(0,sj,fs[5],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[6]=1
(sa.low=1/(1+exp(-a.sa-b*(Mz[6])^2))) # sa unimodal
fs[6]=0.798; A=matrix(c(0,sj,fs[6],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[7]=0.5
(sa.low=1/(1+exp(-a.sa-b*(Mz[7])^2))) # sa unimodal
fs[7]=.563; A=matrix(c(0,sj,fs[7],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[8]=0
fs[8]=f

# get lambda.g by simulation for each choice of mean(z)

tmax=5E5
n0=matrix(0.5,2,1)
sds=seq(0,.5,0.01)  # sd(z) values to use
lam.s=sa.mean=matrix(0,length(sds),8)

for(m in 1:8){
  
  A=matrix(c(0,sj,fs[m],0),2,2); # fill matrix with constant vr's
  
  for(i in 1:length(sds)){
    
    zs=rnorm(tmax,mean=Mz[m],sd=sds[i])  # generate random envts. with current mean and sd
    
    sas=1/(1+exp(-a.sa-b*(zs)^2)) # sa unimodal - get all adult survival rates
    sa.mean[i,m]=mean(sas)  # save arithmetic mean adult survival
    
    n=n0
    lambdas=matrix(0,1,tmax)
    for(t in 1:tmax){
      A[2,2]=sas[t]   # put variable adult survival into matrix
      n=A %*% n
      lambdas[t]=sum(n)
      n=n/sum(n)
    }
    
    lam.s[i,m]=exp(mean(log(lambdas)))
  }
}

# fig 5C
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
windows()
matplot(sds,sa.mean,type='l',xlab='SD(environment)',ylab='Arithmetic mean adult survival',
        lwd=4,col=cbbPalette,cex.axis=2,lty='solid')

# fig 5D
windows()
matplot(sds,lam.s,type='l',xlab='SD(environment)',ylab='Stochastic lambda',
        lwd=4,col=cbbPalette,cex.axis=2,lty='solid')

# fig 5A
z=seq(-5,5,.1)
sa.z=1/(1+exp(-a.sa-b*z^2))

windows()
matplot(z,cbind(sa.z),type="l",
        ylab="Adult survival",xlab='Environment, z',
        lwd=4,col=c('black'),cex.axis=2,lty='solid')

# fig 5B
M1=0
M2=-1.5
M3=-3
z=seq(-5,5,.0005)
zd1a=dnorm(z,mean=M1,sd=.1)
zd1b=dnorm(z,mean=M1,sd=.5)
zd2a=dnorm(z,mean=M2,sd=.1)
zd2b=dnorm(z,mean=M2,sd=.5)
zd3a=dnorm(z,mean=M3,sd=.1)
zd3b=dnorm(z,mean=M3,sd=.5)
windows()
matplot(z,cbind(zd1a,zd1b,zd2a,zd2b,zd3a,zd3b),type='l',
        xlab='Environment, z',ylab='Density',lwd=4,
        col=c('black'),cex.axis=2,lty='solid')


# PUSH THE EXAMPLE ABOVE WITH STRONGEST ADAPTIVE LABILITY WITH LOW SD(Z)
# TO HIGHER SD(Z)

tmax=5E5
n0=matrix(0.5,2,1)
sds=seq(0,10,0.1)  # wider range of SD(z)
lam.s=sa.mean=matrix(0,length(sds))

m=2 # this is the case where M(z)=-3 and f=8.257

A=matrix(c(0,sj,fs[m],0),2,2); # fill matrix with constant vr's

for(i in 1:length(sds)){
  
  zs=rnorm(tmax,mean=Mz[m],sd=sds[i])  # generate random envts. with current mean and sd
  
  sas=1/(1+exp(-a.sa-b*(zs)^2)) # sa unimodal - get all adult survival rates
  sa.mean[i]=mean(sas)  # save arithmetic mean adult survival
  
  n=n0
  lambdas=matrix(0,1,tmax)
  for(t in 1:tmax){
    A[2,2]=sas[t]   # but variable adult survival into matrix
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lam.s[i]=exp(mean(log(lambdas)))
}

# fig 5E
windows()
plot(sds,lam.s,type='l',
     xlab='SD(environment)',ylab='Stochastic lambda',
     lwd=4,
     col="#E69F00",cex.axis=2,lty='solid')
lines(c(0,10),c(1,1),lty='dashed',lwd=4,col='black')