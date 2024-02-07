# make panels of Figs S2, S3, and S4 for a manuscript entitled 
# "What is adaptive demographic lability and when might we expect to see it?"

graphics.off()

# fig S2 A

r0=0.5
a=0
b=.4
z=seq(-12,12,.1)
r=1/(1+exp(-a-b*z))
windows()
plot(z,r,xlab="Environment, z",ylab="Adult survival",
     type="l",lwd=4,cex.axis=2)
lines(c(-5,-5),c(0,.12),lty='dashed',lwd=4,col='red')
lines(c(-12,-5),c(.12,.12),lty='dashed',lwd=4,col='red')
lines(c(6,6),c(0,.92),lty='dashed',lwd=4,col='black')
lines(c(-12,6),c(.92,.92),lty='dashed',lwd=4,col='black')

# fig S2 B

windows()
zf=seq(-7,-3,.1)
rf=1/(1+exp(-a-b*zf))
zs=seq(4,8,.1)
rs=1/(1+exp(-a-b*zs))
x=seq(-2,2,.1)
matplot(x,cbind(rf,rs),xlab="Environment, z",ylab="Adult survival",
        type="l",lwd=4,cex.axis=2,lty='solid',col=c('red','black'))

# fig S2 C

a=log(2/3)
b=.4
z=seq(-12,12,.1)
r=1/(1+exp(-a-b*z))
r1=r
r2=2*r
r3=3*r
windows()
matplot(z,cbind(r1,r3),xlab="Environment, z",ylab="Fertility",
        type="l",lwd=4,cex.axis=2,col=c('black','red'),lty='solid')
abline(v=0,lty='dashed',lwd=3)

# fig S2 D

z=seq(-2,2,.1)
r=1/(1+exp(-a-b*z))
windows()
matplot(z,cbind(r),xlab="Environment, z",ylab="Fertility",
        type="l",lwd=4,cex.axis=2)
abline(v=0,lty='dashed',lwd=3)


# fig S3

# slow life history; adult survival varies (quadratic function) with beta environmental stochasticity

sj=.1
f=.5
sa=.95

# plot quadratic functions
b=.06 # slope parameter
cs=c(-.06,-.03,0,.03,.06)  # cuvature: convex if c>0, concave if c<0
z=seq(-.5,.5,.01)
z2=z^2

sas=matrix(0,length(z),length(cs))
for(i in 1:length(cs)){
  sas[,i]= sa+b*z+cs[i]*z2  # quadratic function
}

# fig S3 A

windows()
matplot(z,sas,type="l",ylim=c(0.9,1.0),ylab="s.a",xlab="Environment, z",lty="solid",col="black",lwd=3,cex.axis=2)

# simulate log lambda.g for different values of b, beta environment

sj=.1
f=.5
sa=.95

m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

# fig S3 B

windows()
z=seq(-.5,.5,.01)
plot(z,dbeta(z+0.5,aa,bb),xlab="Environment, z",ylab="Density",type="l",lwd=4,cex.axis=2)


tmax=50000
zs=rbeta(tmax,aa,bb)-0.5  # z's centered on zero
zs2=zs^2
n0=matrix(c(0.5,0.5),2,1)

# sa varies, f and sj constant

# quadratic parameters
b=.06 # slope parameter
cs=seq(-.06,.06,0.002)  # curvature: convex if c>0, concave if c<0

LLg.c=array(0,length(cs))

for(i in 1:length(cs)){
  
  lambdas=array(0,tmax)
  n=n0
  
  sas = sa+b*zs+cs[i]*zs2  # quadratic function
  
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  LLg.c[i]=mean(log(lambdas))
  
}

# fig S3 C

windows()
plot(cs,LLg.c,xlab="Curvature parameter, c",ylab="log(geometric mean)",
     type="l",lwd=3,cex.axis=2)
abline(0,0,lty='dashed',lwd=3)
abline(v=0,lty='dotted',lwd=3)

# Fig S4 A
# plot power functions
cs=c(.75,1,2,4,8)  # convex if c>1, concave if c<1
as=((2^cs)*sa-1)/((2^cs)-1)
z=seq(0,1,.01)

sas=matrix(0,length(z),length(cs))
for(i in 1:length(cs)){
  sas[,i]= as[i]+(1-as[i])*(z^cs[i])  # power function
}
windows()
matplot(z,sas,type="l",ylab="s.a",xlab="Environment, z",lty="solid",col="black",lwd=3,cex.axis=2)


# simulate log lambda.g for different values of c, Beta environment

sj=.1
f=.5
sa=.95

tmax=50000

m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

# Fig S4 B

# plot beta pdf
windows()
plot(z,dbeta(z,aa,bb),xlab="Environment, z",ylab="Density",type="l",lwd=3,cex.axis=2)

zs=rbeta(tmax,aa,bb)
n0=matrix(c(0.5,0.5),2,1)

# sa varies, f and sj constant

cs=seq(.8,6,.2)  # convex if c>1, concave if c<1
as=((2^cs)*sa-1)/((2^cs)-1)

LLg.c=array(0,length(cs))

for(i in 1:length(cs)){
  
  lambdas=array(0,tmax)
  n=n0
  sas = as[i]+(1-as[i])*(zs^cs[i])  # power function
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  LLg.c[i]=mean(log(lambdas))
  
}

# Fig S4 C

windows()
plot(cs,LLg.c,xlab="Power function exponent, c",ylab="log(geometric mean)",
     type="l",lwd=3,cex.axis=2)
abline(0,0,lty='dashed',lwd=3)
abline(v=1,lty='dotted',lwd=3)
