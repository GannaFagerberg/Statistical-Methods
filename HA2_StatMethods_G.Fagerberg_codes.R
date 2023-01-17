
#########################
##      HOMEWORK 2     ##
##         in          ##
## Statistical Methods ##
#########################

#######################
##  Ganna Fagerberg  ##
#######################

library(ICAOD)
library(knitr)
library(plotly)
library(nloptr)


#######################
###    PROBLEM 1    ###
#######################

## 1A. Calculate the standardized infromation matrix, (SIM).

## Compute SIM

# Function to compute SIM

Inf_linearreg<-function(w,x){
  Mi<-matrix(0,2,2)
  for(i in 1:length(w)){
    z<-c(1,x[i]^2)
    z1<-z%*%t(z)
    Mi<-Mi+w[i]*z1
  }
  (Mi)
}

x<-c(0,3,6) # design points
w<-c(1/3,1/3,1/3) # weights at design points

SIM=Inf_linearreg(w,x) # compute SIM

################################################

## 1B. Plot the standardized prediction variance for design 1

# Find the inverse of SIM.

SIM_inv=solve(SIM)

# Plot the prediction function.
(0.004*x^4-0.128*x*2+1.962)
x<-seq(0,6,length=100) # 
plot(x,(0.004*x^4-0.128*x*2+1.962),type="l", ylim=c(1,6), 
main="Figure 1. Prediction function for design 1")
abline(h=2)

test=function(x){(0.004*x^4-0.128*x*2+1.962)}
test(0); test(6)


################################################

## 1C. Plot the standardized prediction variance for the optiml design

# Function to compute SIM

x<-c(0,3,6) # design points
w<-c(1/2,0,1/2) # weights at design points

SIM_opt=Inf_linearreg(w,x) # compute SIM
SIM_inv_opt=solve(SIM_opt) # find the inverse

# Plot the function

x<-seq(0,6,length=100) 
plot(x,(0.00308642*x^4-0.11111112*x^2+2),type="l", ylim=c(1,3), 
     main="Figure 2. Prediction function for design 2")
abline(h=2)

test=function(x){(0.00308642*x^4-0.11111112*x^2+2)}
test(6); test(0)


################################################

## 1D. Compute relative D-efficiency of design 1

RE=sqrt(det(SIM_inv_opt)/det(SIM_inv))
mult=1/RE
extra=round((mult-1)*100) # how many more subjects needed (in procent)

################################################

## 1E. Obtain A-optimal design by modifying R-code provided in lecture-2

# Criterion function for A-optimality

crit<-function(w,x){
  Mi<-matrix(0,2,2)
  for(i in 1:length(w)){
    z<-c(1,x[i]^2)
    z1<-z%*%t(z)
    Mi<-Mi+w[i]*z1
  }
  
  # for other optimality criterion function change here 
  sum(diag(solve(Mi)))
  
}


# Optimum value of w while x values are fixed

critw<-function(w,x){
  g<-w^2/sum(w^2)  # why g is based on w2?
  crit(g,x) 
}


# Optimum values of x and w simulatneously

critwx<-function(wx){
  r<-length(wx)/2
  c<-length(wx)/r
  mat<-matrix(wx,r,c,byrow =FALSE)
  x<-mat[,1]
  w<-mat[,2]
  w<-w^2/sum(w^2)
  crit(w,x) 
}


# Removing values with low weight and average the value wich
# has less difference than 0.001 and sum-up the weight


Consolidate<-function(x,w,ex,ew){
  de<-data.frame(cbind(x,w))
  de<-de[which(de$w>ew),]
  x<-de[,1]
  
  k<-dim(de)[1]
  g <-rep(0,k)
  for(i in 1:k){
    
    if(g[i]==0){
      g[i]<-max(g)+1
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)
          {g[j]<-g[i]} 
      }}}
  de$g<-g
  f1<-aggregate(x~g,de,mean)[,2]
  f2<-aggregate(w~g,de,sum)[,2]
  newx<-data.frame(x=f1,w=f2) 
  newx
}


# Use initial design now and minimize varience for both x and w 

final<-function(xop,wop,iter,accuracy){         
  i=1;improvement=1            
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last=cbind(xop,wop)
    n<-length(xop)
    
    wxinit<-c(xop,sqrt(wop))
    
    opdes<-neldermead(x0=wxinit,fn=critwx, 
                      lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r<-length(opdes )/2
    c<-length(opdes )/r
    
    mat<-matrix(opdes,r,c,byrow =FALSE)
    x<-mat[,1]
    w<-mat[,2]
    w<-w^2/sum(w^2)  
    
    #again simplify the design
    sdes<-Consolidate(x,w,0.01,0.001)
    xop<-sdes$x
    wop<-sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr<-data.frame(x=xop,w=wop,trace=crit(wop,xop),improvement=improvement)
  
  rr[order(rr$x),]
}


# Use all the above functions

library(nloptr)
xl=0; xu=6

# Starting values of design points with weight
x=runif(10,xl,xu);                       
winit = rep(1,length(x)) / length(x)
winit = sqrt(winit)


# Fixed value of x; find optimal value for w by maximizing  critw
wop<-neldermead(x0=winit,fn=critw,x=x)$par

# Covert wop between 0 and 1 so that that the sum is 1
wop= (wop^2)/sum(wop^2);
initialdes<-data.frame(x,w=wop,trace=crit(wop,x))

# Simplify the design
sdes<-Consolidate(x,wop,0,0.0001)
xop<-sdes$x; wop<-sdes$w
siminitialdes<-data.frame(x=xop,w=wop,trace=crit(wop,xop))

# Final result
final.res<- final(xop,wop,iter=100,accuracy=0.000000001)




#######################
###    PROBLEM 2    ###
#######################

## 2A. Estimate the parameters in the logistic model 

# Beetle data
beetle.dat = cbind(c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113,1.8369, 1.861, 1.8839), 
             c(59, 60, 62, 56, 63, 59, 62, 60), 
             c(6, 13, 18, 28, 52, 53, 61, 60))
colnames(beetle.dat) <- c("Dose", "N","killed") # rename the columns

# Display the data
kable(beetle.dat, digits = 4)

ns = length(beetle.dat[,1]) # number of responses
ni <- beetle.dat[,2] # numbers of beetles
yi <- beetle.dat[,3] # numbers of killed 
xi <- beetle.dat[,1] # dose
wi <- ni/sum(ni) #weights
N  <- sum(ni)

# Plot the proportion killed a function of dose:
plot(xi, yi/ni,pch=16,col='blue',
     ylim=c(0,1),xlab="log(dose)",ylab="Proportion killed", )

# Fit logistic regression:
log.regr <- glm((cbind(yi,ni-yi) ~ xi + I(xi^2)), family=binomial(logit), 
                 data=data.frame(beetle.dat))
summary(log.regr)


# Plot the fitted line:
x.grid = seq(1.65, 1.95, 0.01)
pi.hat.logistic = predict(log.regr, newdata=data.frame(xi=x.grid), type="response")
lines(x.grid,pi.hat.logistic,lwd=2,col='red')

################################################

## 2B. D-optimal design by modyfying codes given in lecture 2

# Criterion function for D-optimality

crit<-function(w,x, b0, b1, b2){
  p1<- 1 #the numerator
  p2<- 1+exp(-(b0+b1*x+b2*x^2)) # the denominator
  e<-p1/p2  # the function
  v<-e*(1-e) # the variance
  Mi<-matrix(0,3,3)
  for(i in 1:length(x)){
    z<-c(1, x[i], x[i]^2)
    z1<-z%*%t(z)
    Mi<-Mi+w[i]*v[i]*z1
  }
  
  -log(det(Mi))
  
}


# Optimum value of w while x values are fixed

critw<-function(w,x,b0,b1, b2){
  g<-w^2/sum(w^2)  # why g is based on w2?
  crit(g,x,b0,b1, b2) 
}


# Optimum values of x and w simulatneously

critwx<-function(wx,b0,b1, b2){
  r<-length(wx)/2
  c<-length(wx)/r
  mat<-matrix(wx,r,c,byrow =FALSE)
  x<-mat[,1]
  w<-mat[,2]
  w<-w^2/sum(w^2)
  crit(w,x,b0,b1, b2) 
}

# Simplify design by removing values with low weight and average the value which
# has less difference than 0.001; sum-up their weights

Consolidate<-function(x,w,ex,ew){                   
  
  de<-data.frame(cbind(x,w))
  de<-de[which(de$w>ew),]                 
  
  x<-de[,1]
  
  k<-dim(de)[1] #number of obs in de
  
  g <-rep(0,k)
  
  for(i in 1:k){
    
    if(g[i]==0){
      g[i]<- max(g)+1
      
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)             #ex
            
          {g[j]<-g[i]} 
      }}}
  
  de$g<-g
  
  f1<-aggregate(x~g,de,mean)[,2]  # check aggregate
  f2<-aggregate(w~g,de,sum)[,2]
  
  newx<-data.frame(x=f1,w=f2) 
  newx
}


# Minimize variance for both x and w 

final<-function(xop,wop,b0,b1, b2,iter,accuracy){         
  
  i=1;improvement=1
  
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last=cbind(xop,wop)
    n<-length(xop)
    
    wxinit<-c(xop,sqrt(wop))
    
    opdes<-neldermead(x0=wxinit,fn=critwx,b0=b0,b1=b1, b2=b2, 
                      lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r<-length(opdes )/2
    c<-length(opdes )/r
    
    mat<-matrix(opdes,r,c,byrow =FALSE)
    x<-mat[,1]
    w<-mat[,2]
    w<-w^2/sum(w^2)  
    
    #again simplify the design
    sdes<-Consolidate(x,w,0.01,0.001)
    xop<-sdes$x
    wop<-sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr<-data.frame(x=xop,w=wop,var=crit(wop,xop,b0,b1, b2),improvement=improvement)
  
  rr[order(rr$x),]
}


# Get final results
library(nloptr)
xl=1.6907; xu=1.8839
b0=431.11; b1=-520.62; b2=156.41

# Starting values of design point with weight
x=runif(41,xl,xu);                 
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit);


# Find optimal value for w (for fixed x) by maximizing  critw
wop<-neldermead(x0=winit,fn=critw,x=x,b0=431.11,b1=-520.62, b2=156.41)$par

# Convert wop so that it lies between 0 and 1 and its sum is 1
wop= (wop^2)/sum(wop^2);
initialdes<-data.frame(x,w=wop,var=crit(wop,x,b0, b1,b2))

# Simplify the design
sdes<-Consolidate(x,wop,0,0.0001)
xop<-sdes$x; wop<-sdes$w
siminitialdes<-data.frame(x=xop,w=wop,var=crit(wop,xop,b0=431.11,b1=-520.62, b2=156.41))

final.res<- final(xop,wop,b0,b1,b2,iter=100,accuracy=0.000000001)


# Compare the result with that from the ICAOD package

library(ICAOD)
res <-locally(formula = ~(1/(1+exp(-(b0+b1*x+b2*x^2)))),
               predvars=c("x") , parvars = c("b0", "b1","b2"),
               family=binomial(), lx=c(1.6907),ux=c(1.8839),
               inipars=c(431.11,-520.62, 156.41), iter=30, k=3)

################################################

## 2C. The sensitivity analysis


final.res1<- final(xop,wop,b0=300,b1=-400,b2=100,iter=100,accuracy=0.000000001)
final.res2<- final(xop,wop,b0=250,b1=-200,b2=100,iter=100,accuracy=0.000000001)
final.res3<- final(xop,wop,b0=500,b1=-600,b2=250,iter=100,accuracy=0.000000001)
final.res4<- final(xop,wop,b0=400,b1=-450,b2=130,iter=100,accuracy=0.000000001)



################################################

## 2E. E-optimal design by modyfying codes from lecture 2


crit<-function(w,x, b0, b1, b2){
  p1<- 1 
  p2<- 1+exp(-(b0+b1*x+b2*x^2)) 
  e<-p1/p2  # the function
  v<-e*(1-e) # the variance
  Mi<-matrix(0,3,3)
  for(i in 1:length(x)){
    z<-c(1, x[i], x[i]^2)
    z1<-z%*%t(z)
    Mi<-Mi+w[i]*v[i]*z1
  }
  
  max(eigen(solve(Mi))$values)
  
}



# Optimum value of w while x values are fixed

critw<-function(w,x,b0,b1, b2){
  g<-w^2/sum(w^2)  # why g is based on w2?
  crit(g,x,b0,b1, b2) 
}


# Optimum values of x and w simulatneously

critwx<-function(wx,b0,b1, b2){
  r<-length(wx)/2
  c<-length(wx)/r
  mat<-matrix(wx,r,c,byrow =FALSE)
  x<-mat[,1]
  w<-mat[,2]
  w<-w^2/sum(w^2)
  crit(w,x,b0,b1, b2) 
}


# Simplify design by removing values with low weight; average the values which
# have less difference than 0.001 and sum-up their weights

Consolidate<-function(x,w,ex,ew){                
  
  de<-data.frame(cbind(x,w))
  de<-de[which(de$w>ew),]                 
  
  x<-de[,1]
  
  k<-dim(de)[1] #number of obs in de
  
  g <-rep(0,k)
  
  for(i in 1:k){
    
    if(g[i]==0){
      g[i]<- max(g)+1
      
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)             #ex
            
          {g[j]<-g[i]} 
      }}}
  
  de$g<-g
  
  f1<-aggregate(x~g,de,mean)[,2]  # check aggregate
  f2<-aggregate(w~g,de,sum)[,2]
  
  newx<-data.frame(x=f1,w=f2) 
  newx
}

# Minimize variance for both x and w 

final<-function(xop,wop,b0,b1, b2,iter,accuracy){       
  
  i=1;improvement=1
  
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last=cbind(xop,wop)
    n<-length(xop)
    
    wxinit<-c(xop,sqrt(wop))
    
    opdes<-neldermead(x0=wxinit,fn=critwx,b0=b0,b1=b1, b2=b2, 
                      lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r<-length(opdes )/2
    c<-length(opdes )/r
    
    mat<-matrix(opdes,r,c,byrow =FALSE)
    x<-mat[,1]
    w<-mat[,2]
    w<-w^2/sum(w^2)  
    
    #again simplify the design
    sdes<-Consolidate(x,w,0.01,0.001)
    xop<-sdes$x
    wop<-sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr<-data.frame(x=xop,w=wop,var=crit(wop,xop,b0,b1, b2),improvement=improvement)
  
  rr[order(rr$x),]
}


# Get final results

library(nloptr)
xl=1.6907; xu=1.8839
b0=431.11; b1=-520.62; b2=156.41

# Starting values of design point with weight
x=runif(41,xl,xu);                       # does it matter that 41 is chosen?
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit)

#F optimal value for w by maximizing critw (x is fixed)
wop<-neldermead(x0=winit,fn=critw,x=x,b0=431.11,b1=-520.62, b2=156.41)$par

# Convert wop to lie between 0 and 1 so thatits sum is 1
wop= (wop^2)/sum(wop^2);
initialdes<-data.frame(x,w=wop,var=crit(wop,x,b0, b1,b2))

# Simplify the design
sdes<-Consolidate(x,wop,0,0.0001)
xop<-sdes$x
wop<-sdes$w
siminitialdes<-data.frame(x=xop,w=wop,var=crit(wop,xop,b0=431.11,b1=-520.62, b2=156.41))


final.res<- final(xop,wop,b0,b1,b2,iter=100,accuracy=0.000000001)


# Verify the result 

library(matrixcalc)
library(ICAOD)

Eopt <-function(x,w,b0,b1,b2,fimfunc){
  e1<-eigen(svd.inverse(fimfunc(x = x, w = w, b0=b0,b1=b1,b2=b2)))$values
  max(e1)
}

res<-locally(formula=~(exp(b0+b1*x+b2*x^2)/(1+exp(b0+b1*x+b2*x^2))), predvars = "x", 
              parvars = c("b0", "b1","b2"),
              lx =1.6907, ux =1.8839,crtfunc=Eopt, inipars =c(431.11,-520.62,156.41),
              k=3,iter=100)



#######################
###    PROBLEM 3    ###
#######################

## 3A.  The Emax model and D-optimality. P. 289.

# Design region X=[0,100]
# Theta=(E0, EMAX, ED50)=(22, 11.2, 70)
# This is the best guess (maybe Bayesian analysis?)

library(ICAOD)

res3A <-locally(formula=~E0+((EMAX*x)/(ED50+x)), predvars = "x", 
              parvars = c("E0", "EMAX","ED50"),
              lx =0, ux =100, inipars = c(22, 11.2, 70),k=3,iter=100)


################################################

## 3B. The D-optimal design for the 3PL model by modyfying the codes in lecture 2

# Criterion function for D-optimality

crit<-function(w,x,a,b,d){
  #p1 <-1+c*exp(-a*(x-b))
  #p2<-1+exp(-a*(x-b))
  
  e<-d+((1-d)/(1+exp(-a*(x-b))))
  
  v<-e*(1-e)
  
  Mi<-matrix(0,3,3)
  for(i in 1:length(x)){
    
    z<-c(x[i]-b, -a, (1+exp(-a*(x[i]-b)))/(1-d))
    
    scalar<-1/(1+d*exp(-a*(x[i]-b)))
    
    z1<- scalar*(z%*%t(z))
    
    Mi<-Mi+w[i]*v[i]*z1
  }
  

  -log(det(Mi))
}


# Optimum value of w while x values are fixed

critw<-function(w,x,a,b,d){
  g<-w^2/sum(w^2)  
  crit(g,x,a,b, d) 
}


# Optimum values of x and w simulatneously

critwx<-function(wx,a,b,d){
  r<-length(wx)/2
  c<-length(wx)/r
  mat<-matrix(wx,r,c,byrow =FALSE)
  x<-mat[,1]
  w<-mat[,2]
  w<-w^2/sum(w^2)
  crit(w,x,a,b,d) 
}


# Simplify the design by removing values with low weights; average the values which
# have less difference than 0.001 and sum-up their weights

Consolidate<-function(x,w,ex,ew){
  de<-data.frame(cbind(x,w))
  de<-de[which(de$w>ew),]
  x<-de[,1]
  
  k<-dim(de)[1]
  g <-rep(0,k)
  for(i in 1:k){
    
    if(g[i]==0){
      g[i]<-max(g)+1
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)
          {g[j]<-g[i]} 
      }}}
  de$g<-g
  f1<-aggregate(x~g,de,mean)[,2]
  f2<-aggregate(w~g,de,sum)[,2]
  newx<-data.frame(x=f1,w=f2) 
  newx
}

# Minimize the variances for both x and w 

final<-function(xop,wop,a,b,d, iter,accuracy){
  i=1;improvement=1
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last=cbind(xop,wop)
    n<-length(xop)
    
    wxinit<-c(xop,sqrt(wop))
    
    opdes<-neldermead(x0=wxinit,fn=critwx,a=a,b=b,d=d,
                      lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r<-length(opdes )/2
    c<-length(opdes )/r
    mat<-matrix(opdes,r,c,byrow =FALSE)
    x<-mat[,1]
    w<-mat[,2]
    w<-w^2/sum(w^2)  
    
    #again simplify the sedign
    sdes<-Consolidate(x,w,0.01,0.001)
    xop<-sdes$x
    wop<-sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr<-data.frame(x=xop,w=wop,var=crit(wop,xop,a,b,d),improvement=improvement)
  rr[order(rr$x),]
}


# Get final results

library(nloptr)
xl=-3; xu=3
a=0.5; b=1; d=0.05

# Starting values of design points with weight
x=runif(10,xl,xu);
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit)

# Optimal values for w by maximizing critw (x are fixed)
wop<-neldermead(x0=winit,fn=critw,x=x,a=0.5,b=1, d=0.05)$par

#convert wop so that it lies between 0 and 1 and sums up to 1
wop= (wop^2)/sum(wop^2);
initialdes<-data.frame(x,w=wop,var=crit(wop,x,a,b,d))

# Simplify the design
sdes<-Consolidate(x,wop,0,0.0001)
xop<-sdes$x; wop<-sdes$w
siminitialdes<-data.frame(x=xop,w=wop,var=crit(wop,xop,a,b,d))


res3B<-final(xop,wop,a,b,d,iter=100,accuracy=0.000000001)


# Verify the result with ICAOD package

library(ICAOD)
res3B_alt <-locally(formula = ~d+((1-d)/(1+exp(-a*(x-b)))),
               predvars=c("x") , parvars = c("a", "b", "d"),
               family=binomial(), lx=c(-3),ux=c(3),
               inipars=c(0.5,1,0.05), iter=100, k=3)









