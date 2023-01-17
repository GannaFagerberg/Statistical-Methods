###############################################
## HOME ASSIGNEMENT 1 in Statistical methods###
###############################################

## GANNA FAGERBERG ##
 


################
### Question 1##
################

#Packages used
library(ggplot2)
library(emdbook)
library(spatstat)
library(coda)
library(mcmcplots)
library(MASS) 
library(MCMCpack)
library(bayestestR)
library(dplyr)
library(HDInterval)
library(ggpubr)

################

## Reproduction of Figure 1##

# Generate data
M=500; k=10 # set number of iterations
n=16; a=2; b=4; # params of the beta-binomial distr
X=Y=array(0,dim=c(k,M)) # initialization vector

# Gibbs sampling method
for (j in 1:M){ 
  
  # seed for reproducibility
  set.seed(2012*j) 
  
  Y[1,j]=rbeta(1,a,b)  #set starting values
  X[1,j]=rbinom(1,n,Y[1,j])  
  
  for (i in 2:k) { #sampling loop for k 
    
    X[i,j]=rbinom(1,n,Y[i-1,j]) 
    Y[i,j]=rbeta(1,a+X[i,j],n-X[i,j]+b) 
  }
  
}

# Extract X and formate the data
count.X_gibbs=data.frame(table(X[10,]))
gibbs <- rep("gibbs", each=nrow(count.X_gibbs))

count.X_gibbs<-data.frame(count.X_gibbs, gibbs)
names(count.X_gibbs)[names(count.X_gibbs) == "gibbs"] <- "type"
names(count.X_gibbs)[names(count.X_gibbs) == "Var1"] <- "X"

## Direct sampling from the Beta-Binomial
set.seed(20102009)
prob.Beta <- rbeta(n,a,b) 
bin.out <- rbinom(M,n,prob.Beta) #bin outcomes based on prob

# Format the data
count.X_dir <-data.frame(table(bin.out))
direct <- rep("direct", each=nrow(count.X_dir))
count.X_dir<-data.frame(count.X_dir, direct)
names(count.X_dir)[names(count.X_dir) == "direct"] <- "type"
names(count.X_dir)[names(count.X_dir) == "bin.out"] <- "X"

# Reproduce Figure 1.
# Combine the two datasets obtaine from previous simulations
df.merged<-rbind(count.X_gibbs, count.X_dir)
pX.merged<-ggplot(data=df.merged, aes(x=X, y=Freq,  fill=type)) +
  geom_bar(stat="identity", color="black", position ="dodge", width=0.8)+ theme_classic()+
  labs(caption="Figure 1. Comparison of Two Histograms of Samples of Size 
       m = 500 From the Beta-Binomial Distribution With n = 16, a = 2, 
       and b = 4. The black histogram sample was obtained using Gibbs 
       sampling with k = 10. The white histogram sample was generated 
       directly from the beta-binomial distribution.")+
  scale_y_continuous(name=" ", limits=c(0, 70), breaks=seq(0,70,10))+
  scale_x_discrete(name=" ")+ theme(legend.position="none")+ 
  scale_fill_manual(values=c('black','white')); pX.merged

############################################

## Reproduction of Figure 3##

# Set initial values
Theta<-Y[10, ]
x=c(0:16) #range of x outcomes

prob.X<- sapply(x, function(x){mean(dbinom(x,n,Theta))})

# Format the data
count.X_theta<-data.frame(x, prob.X)
names(count.X_theta)[names(count.X_theta) == "x"] <- "X"
names(count.X_theta)[names(count.X_theta) == "prob.X"] <- "perc"
gibbs <- rep("gibbs", each=nrow(count.X_theta))
count.X_theta<-data.frame(count.X_theta, gibbs)
names(count.X_theta)[names(count.X_theta) == "gibbs"] <- "type"

# Exact probabilities
Gam_fun<- function (x, n, a, b) {
  
  result <- choose(n,x)*(gamma(a+b)/(gamma(a)*gamma(b)))*
    ((gamma(x+a)*gamma(n-x+b)/gamma(a+b+n)))
  
}
exact_res<- Gam_fun(x,n,a,b)

# Extract the probabilites and format the data  
count.X_exact<-data.frame(x, exact_res)
names(count.X_exact)[names(count.X_exact) == "x"] <- "X"
names(count.X_exact)[names(count.X_exact) == "exact_res"] <- "perc"
exact <- rep("exact", each=nrow(count.X_exact))
count.X_exact<-data.frame(count.X_exact, exact)
names(count.X_exact)[names(count.X_exact) == "exact"] <- "type"


# Reproduce Figure 3.
df.merged_theta<-rbind(count.X_theta, count.X_exact)

pX.merged_theta<-ggplot(data=df.merged_theta, aes(x=X, y=perc,  fill=type)) +
  geom_bar(stat="identity", color="black", position ="dodge", width=0.8)+ theme_classic()+
  labs(caption="Figure 3. Comparison of Two Probability Histograms of the Beta-Binomial 
                Distribution With n = 16, ct = 2, and f3 = 4. The black histogram represents 
                estimates of the marginal distribution of using Equation (2.11), 
                based on a sample of Size m = 500 from 
                the pair of conditional distributions in (2.6). The Gibbs sequence 
                had length k = 10. The white histogram represents the
                exact beta-binomial probabilities")+
  scale_y_continuous(name=" ", limits=c(0, 0.12), breaks=seq(0,0.12,0.02))+
  scale_x_discrete(name=" ")+ theme(legend.position="none")+ 
  scale_fill_manual(values=c('black','white')); pX.merged_theta



############################################

## Reproduction of Figure 5## 

# Gibbs sampler
lambda=16 # parameter of a poisson
X=Y=N=array(0,dim=c(k,M)) #matricies to hold the sim values

for (j in 1:M){
  
  set.seed(2012*j)
  #set initial values
  N[1,j] <- rpois(1,lambda)
  Y[1,j] <- rbeta(1,a,b) 
  X[1,j] <- rbinom(1, N[1,j], Y[1,j])
  
  
  for (i in 2:k){ #sampling loop 
    X[i,j]=rbinom(1,N[i-1,j],Y[i-1,j]) 
    Y[i,j]=rbeta(1,a+X[i,j],N[i-1,j]-X[i,j]+b) 
    N[i,j]=X[i,j]+rpois(1, lambda=(lambda*(1-Y[i,j])))
  }
}

# Extract and format the data
Y.vec <- Y[10,]; N.vec<-N[10,]
x=seq(0, 20, by=1) #new range for x

prob.X_pois = as.data.frame(cbind(x,sapply(x,function(x){mean(dbinom(x,N.vec,Y.vec))})))

# the plot
pX.pois<-ggplot(data=prob.X_pois, aes(x=x, y=V2)) +
  geom_bar(stat="identity", width=0.5, fill="black")+ theme_classic()+
  labs(caption="Figure 5. Estimates of Probabilities of the Marginal Distribution of 
                          X Using Equation (2. 11), Based on a Sample of Size m = 500 
                          From the Three Conditional Distributions in (4.9) With A = 16, a= 
                          2, and:1 = 4. The Gibbs sequences had length k =10.")+
  scale_y_continuous(name=" ", breaks=seq(0,0.13,0.02))+
  scale_x_continuous(name=" ", limits=c(0, 19), breaks=seq(0,20, 4)); pX.pois

#########################################

## An alternative Gibbs sampling method ## 

M = 1000 #take larger number of iterations
burn.in= 1:500 # set burn.in M/2
X=Y=array(0,dim=c(M,1)) # initialization
set.seed(2009)
# Perform Gibbs iterations 
for (i in 2:M) { 
  Y[1]=rbeta(1,a,b)
  X[1]=rbinom(1, n, T[1]) 
  X[i] =rbinom(1,size=n,prob=Y[i-1]) 
  Y[i] = rbeta(1,a+X[i],b+n-X[i]) 
}


# Discard burn-in
X.out = X[-burn.in]; Y.out = Y[-burn.in]
# Prepare the data for plotting
df.X_out<-data.frame(table(X.out))
names(df.X_out)[names(df.X_out) == "X.out"] <- "X"
alt <- rep("alt", each=nrow(df.X_out))
df.X_out<-data.frame(df.X_out, alt)
names(df.X_out)[names(df.X_out) == "alt"] <- "type"

# Generate directly from beta-binomial
X.out_direct<-rbetabinom(500, size=16, shape1=2, shape2=4)

df.X_direct<-data.frame(table(X.out_direct))
names(df.X_direct)[names(df.X_direct) == "X.out_direct"] <- "X"
direct<- rep("direct", each=nrow(df.X_direct))
df.X_direct<-data.frame(df.X_direct, direct)
names(df.X_direct)[names(df.X_direct) == "direct"] <- "type"

# Merge two datasets
df.merged_alt<-rbind(df.X_direct, df.X_out)

# Plot with Gibbs sampling
p.X_alt<-ggplot(data=df.merged_alt, aes(x=X, y=Freq,  fill=type)) +
  geom_bar(stat="identity", color="black", position ="dodge", width=0.8)+ theme_classic()+
  labs(caption="Comparison of Two Histograms of the Beta- 
       Binomial Distribution With n = 16, a = 2, and b = 14 based on a sample of 500.
       The black histogram represents estimates of the marginal distribution based on 
       a Gibbs sample of length 1000. The white histogram represents direct sampling
       from the beta-binomial dostribution")+
  scale_y_continuous(name=" ", limits=c(0, 70), breaks=seq(0,70,10))+
  scale_x_discrete(name=" ")+ theme(legend.position="none")+ 
  scale_fill_manual(values=c('black','white')); p.X_alt


###################################
## The diagnostics for convergence##


# KS test
#ks.test(X.out, X.out_direct)

# Cdfs
#plot(ecdf(X.out), xlim=range(c(X, bin.out)), col="dodgerblue")
#plot(ecdf(X.out_direct), add=TRUE, lty="dashed", col="purple")


## Trace plots and running average plots

# Running average
# Teoretical mean X
n=16;a=2;b=4; mean_X=n*a/(a+b)
cummean.Xout <- function(X.out) cumsum(X.out) / seq_along(X.out)
cummean.Xout <-cummean.Xout(X.out)

# Teoretical mean Y
mean_Y=1/(1+b/a)
cummean.Yout <- function(Y.out) cumsum(Y.out) / seq_along(Y.out)
cummean.Yout <-cummean.Yout(Y.out)

#Plots
index=c(1:500)
par(mar = c(2, 2, 2, 2)); par(mfcol=c(2,2))
trX=plot(index, X.out, type="l", main="Trace plot for X") # Trace plots
trY=plot(index, Y.out, type="l", main="Trace plot for Y")
cumX=plot(cummean.Xout, type="l", ylim=c(0, 7), xlab="X", ylab="Cumulative mean, X", 
          panel.first = abline(h=5.3333, col="red", las=1), 
          main="Running average plot for X")
cumY=plot(cummean.Yout, type="l", ylim=c(0.3, 0.5), xlab="Y", ylab="Cumulative mean, Y", 
          panel.first = abline(h=0.3333, col="red", las=1), 
          main="Running average plot for Y")

# ACF plots
par(mfcol=c(1,2))
acf(X.out)
acf(Y.out)


## Gelman diagnostics and plots##
## Produce four chains at different starting values

gibbs.diag<-function(seed){
  M = 1000 #take larger number of iterations
  a=2
  b=4
  X=Y=array(0,dim=c(M,1))
  
  # Perform Gibbs iterations
  for (i in 1:M){
    
    set.seed(seed)
    if (seed > 100){
      Y[1]=rbeta(1,a,b)
      X[1]=rbinom(1, n, Y[1]) 
      
    } else{
      Y[1]=runif(1)
      X[1]=rbinom(1, n, Y[1]) 
    }
    
    for (i in 2:M) { 
      
      X[i] =rbinom(1,size=n,prob=Y[i-1]) 
      Y[i] = rbeta(1,a+X[i],b+n-X[i]) 
    }
  }
  df=data.frame(X,Y)
  return(df)
}




ch1<- as.mcmc(gibbs.diag(25))
ch2<- as.mcmc(gibbs.diag(50))
ch3<- as.mcmc(gibbs.diag(110))
ch4<- as.mcmc(gibbs.diag(120))

chains.diag<-mcmc.list(ch1,ch2,ch3, ch4) 

# Gelman diagnostics
gelm.shrink_factor<- gelman.diag(chains.diag, autoburnin = FALSE, multivariate=TRUE)
gelm.shrink_factor 

# Gelman plots
gelman.plot<-gelman.plot(chains.diag, autoburnin = TRUE); gelman.plot


#####################################################################

#####################################################################

################
### Question 2##
################

# Create a function to reproduce a plot with varying priors and posteriors

gamma_fun<-function(shape, rate){
  # set seed for reproducibility
  set.seed(1979)
  # generate data from the poisson as according the observed distribution
  y = rpois(576,0.9373) 
  # set a grid for a range of theta values
  grid = seq(0.01,1.3,0.01) 
  #compute densities of the prior given grid
  d.prior=dgamma(grid, shape, rate) 
  # compute densities of the conjugate posterior
  d.post=dgamma(grid, shape+537, rate+576)
  # compute normalising constant
  const=max(d.post)/max(d.prior)
  # update the prior
  d.prior=const*d.prior
  # parameters of the graphs
  par(mar = c(2, 2, 2, 2))  
  # posterior density plot 
  p=plot(grid, d.post,  type="l", lty = 1, lwd=2,  xlim=c(0.01, 1.3), yaxt='n',
         xlab="", ylab="", col="blue") 
  # graph plots in the same box 
  par(new=TRUE) 
  # prior density plot
  p=plot(grid, d.prior, type = "l", pch=15,  xlab="", ylab="", xlim=c(0.01, 1.3), 
         axes=FALSE, col="red", main=substitute(paste(alpha,"=",shape,", " , beta,"=",rate)))
  
}

# Reproduce the graph for nine combinations of the parameter
par(mfcol=c(3,3)) #set the parameters of the output
gamma_fun(0.1, 0.1)
gamma_fun(10, 0.1)
gamma_fun(100, 0.1)
gamma_fun(0.1, 10)
gamma_fun(10, 10)
gamma_fun(100, 10)
gamma_fun(0.1, 100)
gamma_fun(10, 100)
gamma_fun(100, 100)
legend("topleft",legend=c("Posterior","Prior"),
       text.col=c("black","red"),lty=c(1,1),lwd=c(2,1), col=c("black","red"))

########################################

## Credible intervals ##

## Approximate 95% interval
## Simulate 1000 observations from the posterior
set.seed(1979); ps = rgamma(100000, 537+0.5, 576+0.5); ps
mean(ps); median(ps)

up.appr=round(mean(ps)+1.96*sd(ps),4)
low.appr=round(mean(ps)-1.96*sd(ps),4)
app_ci=c(low.appr, up.appr)

## Exact95% interval
excat_ci=round(qgamma(c(0.025, 0.975), 537, 576),4)

## HDI
#library(HDInterval)
set.seed(1979); ps = rgamma(100000, 537, 576); 
hdi_ci<-round(hdi(ps, credMass = 0.95),4)

## Manual alternative to HPI
# dx<-density(ps); dn <- cumsum(dx$y)/sum(dx$y); li <- which(dn>=0.025)[1]; 
#ui <- which(dn>=0.975)[1]; dx$x[c(li,ui)]

# Plot credible intervals
par(mfcol=c(1,1))
plot(density(ps), main="Credible intervals, 95%")
abline(v=hdi_ci, col="red", lty = 4, lwd=3)
abline(v=excat_ci, col="green", lty = 2, lwd=2)
abline(v=app_ci, col="blue", lty = 3, lwd=3)
legend("topright",legend=c("Appr","Exact", "HPI"),
       lty=c(3,2,4),lwd=c(3,2,3), col=c("blue","green","red"))


#####################################################################

#####################################################################

################
### Question 3##
################

# Generate bivariate data
set.seed(2009); N=1000 
cov_m = matrix(c(1,0.4,0.4,1),2,2,byrow=TRUE)
df = as.data.frame(mvrnorm(1000, mu = c(0,0), Sigma = cov_m))
x=df[1];y=df[2]
#cor(df); x=df[,1]; y=df[,2]; plot(x,y)

#################################
# Metropolis-Hastings algorithm
MH.fun = function(M, half.int,  verbose = TRUE) { 
  
  # create the function that evaluates the log of the post density
  logpost = function(rho){
    if(rho>-1 & rho <1) return(-3./2*log(1.-rho**2) - N*log((1.-rho**2)**(1/2)) -                  
                                 sum(1./(2.*(1.-rho**2))*(x**2-2.*rho*x*y+y**2)))
    else return(-Inf)
  }
  
  # start counter for the number of accepted samples
  accepted_number=0 
  
  # vector to store the samples
  chain_rho = vector("numeric", M) 
  
  rho=0 # starting value
  
  #start the loop
  for (i in 1:M) { 
    
    # draw a sample from the proposal distr (Equation 7)
    rho_candidate=rho + runif(1,-half.int, half.int)  
    #rho_candidate=rnorm(1,0,1) # ADD TERMS TO ACCEPT
    
    # Compute the acceptance probability (Eq. 8 and Eq. 6) (in log domain) 
    num=logpost(rho_candidate)
    den=logpost(rho)
    
    accept=num-den
    
    accept=min(0,accept) # 0 as we are operating in log domain
    
    accept=exp(accept) # transform to original scale
    
    # # Accept rho_candidate with probability accept
    if (runif(1,0,1) < accept) {
      
      rho=rho_candidate 
      accepted_number=accepted_number+1
    }
    
    chain_rho[i]=rho  
  }
  
  if(verbose) cat("acceptance rate:", ((accepted_number/M)), "\n")  
  
  list(chain_rho, acceptance_rate=accepted_number/M)
} 

res<- MH.fun(M=10000, half.int=0.07)

#mean(unlist(res[1]), names=FALSE); sd(unlist(res[1])); range(res[1])

####################################################################
# Reproduce Figure 1 ##

# Prepare the data 
rho.sim <- unlist(res[1], use.names=FALSE)
index<-c(1:10000)

# Scatterplot
scatter.xy<-ggplot(data=df, aes(x=V1, y=V2)) + 
  geom_point(size=1.5, shape=21, color="black", fill="blue")+
  scale_y_continuous(name=" ", limits=c(-4, 4), breaks=seq(-4,4,1))+
  scale_x_continuous(name=" ", limits=c(-4, 4), breaks=seq(-4,4,1))+
  theme_bw()

# Traceplot
df.rho<-data.frame(index, rho.sim) #for traceplot
trace.rho<-ggplot(df.rho, aes(x=index, y=rho.sim))+ geom_line(color="blue")+
  scale_y_continuous(name="rho", limits=c(0, 0.6), breaks=seq(0,.6,.1))+theme_bw()+
  scale_x_continuous(name="", limits=c(0, 10000), breaks=seq(0, 10000,2000))


# Histogram, rho #CHANGE binwidth
rho.sim_df<-data.frame(rho.sim) #for histogram
hist.rho <-ggplot(data=rho.sim_df, aes(x=rho.sim)) +
  geom_histogram(color="black", fill="blue", bins=55)+
  scale_x_continuous(name="rho", limits=c(0, .6), breaks=seq(0, .6,.1))+theme_bw()+
  scale_y_continuous(name="rho", limits=c(0, 2000), breaks=seq(0,2000,200))

# Arrange three plots 
ggarrange(scatter.xy,trace.rho, hist.rho, ncol = 1, nrow = 3)


#############################################################

## Modify the length of the interval to get different levels of acceptance ##

MH.fun1 = function(M, rho, half.int,  verbose = TRUE) { 
  
  set.seed(2009)
  # create the function that evaluates the log of the post density
  logpost = function(rho){
    if(rho>-1 & rho <1) return(-3./2*log(1.-rho**2) - N*log((1.-rho**2)**(1/2)) -                  
                                 sum(1./(2.*(1.-rho**2))*(x**2-2.*rho*x*y+y**2)))
    else return(-Inf)
  }
  
  # start counter for the number of accepted samples
  accepted_number=0 
  
  # vector to store the samples
  chain_rho = vector("numeric", M) 
  
  #start the loop
  for (i in 1:M) { 
    
    # draw a sample from the proposal distr (Equation 7)
    rho_candidate=rho + runif(1,-half.int, half.int)  
    #rho_candidate=rnorm(1,0,1) # ADD TERMS TO ACCEPT
    
    # Compute the acceptance probability (Eq. 8 and Eq. 6) (in log domain) 
    num=logpost(rho_candidate)
    den=logpost(rho)
    
    accept=num-den
    
    accept=min(0,accept) # 0 as we are operating in log domain
    
    accept=exp(accept) # transform to original scale
    
    # # Accept rho_candidate with probability accept
    if (runif(1,0,1) < accept) {
      
      rho=rho_candidate 
      accepted_number=accepted_number+1
    }
    
    chain_rho[i]=rho  
    
    
  }
  
  if(verbose) cat("acceptance rate:", ((accepted_number/M)), "\n")  
  
  return(c(mean=mean(chain_rho), sd=sd(chain_rho), acceptance_rate=accepted_number/M))
  
} 


grid<-seq(0.05, 0.5, 0.02) # set values for half.width
res.accep<-sapply(grid, function(grid){MH.fun1(M, rho, half.int=grid, verbose = FALSE)})
res.accep_df<-as.data.frame(t(res.accep))

# Plot acceptance rate
res.accep_df<-data.frame(grid, res.accep_df)
p.accept_rates<- ggplot(data=res.accep_df, aes(x=grid, y=acceptance_rate)) +
  geom_line()+ theme_classic()+
  labs(caption=" MH acceptance rates depending on 
                the half-length parameter")+
  scale_y_continuous(name="Acceptance rate", limits=c(0, 0.7), breaks=seq(0,0.7,0.1))+
  scale_x_continuous(name="Half intervals", limits=c(0.05, 0.5), breaks=seq(0.05, 0.5, 0.05))
p.accept_rates


##############################
## Normal symmetric proposal##

MH.fun_nm = function(M, rho, sigma,  verbose = TRUE) { 
  
  set.seed(2009)
  # create the function that evaluates the log of the post density
  logpost = function(rho){
    if(rho>-1 & rho <1) return(-3./2*log(1.-rho**2) - N*log((1.-rho**2)**(1/2)) -                  
                                 sum(1./(2.*(1.-rho**2))*(x**2-2.*rho*x*y+y**2)))
    else return(-Inf)
  }
  
  
  accepted_number=0 # start counter
  
  chain_rho = vector("numeric", M) # store the samples
  rho=0
  
  for (i in 1:M) { #start the loop
    
    # draw a sample from the proposal distr, Equation 7
    rho_candidate=rho + rnorm(1,0, sigma)  
    #rho_candidate=rnorm(1,0,1) # ADD TERMS TO ACCEPT
    
    num=logpost(rho_candidate)
    den=logpost(rho)
    
    accept=num-den
    
    accept=min(0,accept) # 0 as we are operating in log domain
    
    accept=exp(accept) # transform to original scale
    
    # # Accept rho_candidate with probability accept
    if (runif(1,0,1) < accept) {
      
      rho=rho_candidate 
      accepted_number=accepted_number+1
    }
    
    chain_rho[i]=rho  
    
    
  }
  
  if(verbose) cat("acceptance rate:", ((accepted_number/M)), "\n")  
  
  return(c(mean=mean(chain_rho), sd=sd(chain_rho), acceptance_rate=accepted_number/M))
  #chain_rho
  
} 

sigma=seq(0.01, 0.5, 0.02)

res_nm<-sapply(sigma, function(sigma){MH.fun_nm(M, rho, sigma=sigma, verbose = FALSE)})
res_nm<-as.data.frame(t(res_nm))

# Plot acceptance rate
df.accept_nm<-data.frame(sigma, res_nm)

p.accept_nm<- ggplot(data=df.accept_nm, aes(x=sigma, y=acceptance_rate)) +
  geom_line()+ theme_classic()+
  labs(caption=" Acceptance rates depending on 
                the standard errors of the normal proposal")+
  scale_y_continuous(name="Acceptance rate", limits=c(0, 0.9), breaks=seq(0,0.9,0.1))+
  scale_x_continuous(name="Sigmas", limits=c(0.01, 0.5), breaks=seq(0.01, 0.5, 0.05))
p.accept_nm


##################################
## Non-symmetric normal proposal##

set.seed(1979)
MH.fun_asy = function(M, sigma, verbose=TRUE) { 
  
  f.ln <- function(rho){    
    if(rho>-1 & rho <1) 
      return(-3/2*log(1-rho**2) - N*log((1-rho**2)**(1/2)) -                  
               sum(1/(2*(1-rho**2))*(x**2-2.*rho*x*y+y**2)))  
    else return(-Inf)} 
  
  q.ln <- function(rho){    
    if(rho>-1 & rho <1) 
      return(dnorm(rho, 0, sd=sigma, log=TRUE))  
    else return(-Inf)} 
  
  accepted_number=0
  
  chain_rho = vector("numeric", M)
  
  rho = 0
  #rnorm(1,0,sd=sigma)
  
  for (i in 1:M) { #start the loop
    
    
    rho_candidate=rnorm(1,0, sd=sigma)  # candidate value
    
    # Compute the acceptance probability, Eq. 8 and Eq. 6 (in log domain) 
    num=f.ln(rho_candidate)+q.ln(rho)
    den=f.ln(rho)+q.ln(rho_candidate)
    
    accept=num-den
    accept=min(0,accept) # 0 as we are operating in log domain
    
    accept=exp(accept) # transform to original scale
    
    # # Accept rho_candidate with probability accept
    if (runif(1,0,1) < accept) {
      
      rho=rho_candidate 
      accepted_number=accepted_number+1
    }
    
    chain_rho[i]=rho  
    
    
  }
  
  if(verbose) cat("acceptance rate:", ((accepted_number/M)), "\n")  
  
  return(c(mean=mean(chain_rho), sd=sd(chain_rho), acceptance_rate=accepted_number/M))
  #chain_rho
  
} 

res_asy<-MH.fun_asy(M=10000, sigma=0.2)
sigma=seq(0.05, 0.2, 0.02)
res_asy_nm<-sapply(sigma, function(sigma){MH.fun_asy(M=10000, sigma=sigma, verbose = FALSE)})
res_asy_nm<-as.data.frame(t(res_asy_nm))
res_asy_nm<-data.frame(sigma, res_asy_nm)

p.accept_nm<- ggplot(data=df.accept_nm, aes(x=sigma, y=acceptance_rate)) +
  geom_line()+ theme_classic()+
  labs(caption="Figure 7. Acceptance rates depending on 
                the standard errors of the proposed normal distribution")+
  scale_y_continuous(name="Acceptance rate", limits=c(0, 0.9), breaks=seq(0,0.9,0.1))+
  scale_x_continuous(name="Sigmas", limits=c(0.01, 0.5), breaks=seq(0.01, 0.5, 0.05))
p.accept_nm

###########################################################
## Extend the analysis to the three streams of normal obs##

N=1000 
sigma = matrix(c(1,0.4,0.4,0.4,1,0.4,0.4,0.4,1),3,3,byrow=TRUE) 
df = as.data.frame(mvrnorm(1000, mu = c(0,0,0), Sigma = sigma)) 
x=df[,1] 
y=df[,2] 
z=df[,3] 

# Estimated correlation matrix
#sigma_emp=cor(df)

# Metropolis-Hastings algorithm

MH.fun_ml = function(M, rho, half.int,  verbose = TRUE) { 
  
  
  logpost = function(rho){
    if(rho>-1 & rho <1 & rho !=- 0.5) 
      return(-log((2*rho+1)**2)-log((rho-1)**4)-1/2*log((2*rho+1)**N)-N*log(abs(rho-1))-
               sum(1/(2*(2*rho+1)*(1-rho))*((1+rho)*(x**2+y**2+z**2)-2*rho*(x*y+x*z+y*z))))
    else return(-Inf)
  }
  
  accepted_number=0 # start counter
  
  chain_rho = vector("numeric", M) # store the samples
  
  for (i in 1:M) { #start the loop
    
    # draw a sample from the proposal distr, Equation 7
    rho_candidate=rho + runif(1,-half.int, half.int)  
    #rho_candidate=rnorm(1,0,1) # ADD TERMS TO ACCEPT
    
    # Compute the acceptance probability, Eq. 8 and Eq. 6 (in log domain) 
    num=logpost(rho_candidate)
    den=logpost(rho)
    
    accept=num-den
    
    accept=min(0,accept) # 0 as we are operating in log domain
    
    accept=exp(accept) # transform to original scale
    
    # # Accept rho_candidate with probability accept
    if (runif(1,0,1) < accept) {
      
      rho=rho_candidate 
      accepted_number=accepted_number+1
    }
    
    chain_rho[i]=rho  
    
    
  }
  
  if(verbose) cat("acceptance rate:", ((accepted_number/M)), "\n")  
  
  #list(chain_rho, acceptance_rate=accepted_number/M)
  return(c(mean=mean(chain_rho), sd=sd(chain_rho), acceptance_rate=accepted_number/M))
  
} 

res_ml<- MH.fun_ml(M=10000, rho=0, half.int=0.07)

df_ml<-unlist(res_ml[1])
mean(df_ml); sd(df_ml)

#Trace and ACF plot
index=c(1:10000)
par(mfcol=c(1,2))
plot(index, df_ml, type="l" )
acf(df_ml)

## Diffreent values of grid

grid<-seq(0.05, 0.3, 0.02) # set values for half.width
res.ml<-sapply(grid, function(grid){MH.fun_ml(M=10000, rho=0, half.int=grid, verbose = FALSE)})
res.ml_df<-as.data.frame(t(res.ml))
res.ml_df<-data.frame(grid, res.ml_df)


# Plot acceptance rate
res.accep_df<-data.frame(grid, res.ml_df)
p.accept_rates<- ggplot(data=res.accep_df, aes(x=grid, y=acceptance_rate)) +
  geom_line()+ theme_classic()+
  labs(caption="Figure 7. MH acceptance rates depending on 
                the half-length parameter")+
  scale_y_continuous(name="Acceptance rate", limits=c(0, 0.7), breaks=seq(0,0.7,0.1))+
  scale_x_continuous(name="Half intervals", limits=c(0.05, 0.5), breaks=seq(0.05, 0.5, 0.05))
p.accept_rates

























