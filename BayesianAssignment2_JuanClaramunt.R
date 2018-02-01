##########################################
#Outline:
#
#1)Set prior parameters (line 46)
#2)Gibbs sampler with MH step (line 79)
#3)Convergence (line 216)
#4)Estimates and credible intervals (line 283)
#5)DIC (line 308)
#6)Bayes Factor (line 470)
#7)Generate new datasets (line 601)
#8)Posterior predictive check (line 620)
#9)Frequentist approach (line 697)
#
##########################################




#Open needed packages
require(MASS)
#Read the data
data<-read.csv2('dataNBAR.csv',header = TRUE)

#Transform the variables to the desired format
y<-as.numeric(as.character(data[1:474,3]))
x1<-as.numeric(as.character(data[1:474,1]))
x2<-as.numeric(as.character(data[1:474,2]))

#compute the length of the data
l.data<-length(y)
N<-l.data


#set the nuber of desired iterations and burn-in period
n.iter<-10000
burn.period<-1000


#set the nuber of desired chains and their initial values
n.chains<-2
init.values<-matrix(rep(0,n.chains*4),nrow = n.chains, ncol=4)
init.values[1,]<-c(400,1,3,100)#initial values chain 1
init.values[2,]<-c(500,0,1,30)#initial values chain 2


#parameters of the priors

#Mean and standard deviation of the a0 prior (normal distribution)
mu00<-466.676#mean
tau00<-1.922#standard deviation

#b1 prior distribution is a t-distribution   
nu0<-1 #degrees of freedom
m0<-0.971 #mean
s0<-0.041 #standard deviation

#Mean and standard deviation of the b2 prior (normal distribution)
mu20<-1.022#mean
tau20<-0.11#standard deviation

#Sigma prior is an uninformative gamma function 
#with the following hyperparameters
a.prior<-0.001
b.prior<-0.001

#Set the variables
post.a0<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
post.b1<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)  
post.b2<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
post.sigma<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
accept.rate<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)

p.a0<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
p.b1<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
p.b2<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
p.sigma<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))


#Gibbs sampler with MH step: Sampling from the conditional posterior distributions 

for (j in 1:n.chains) {
  #read the initial values of each chain
  a0<-init.values[j,1]
  b1<-init.values[j,2]
  b2<-init.values[j,3]
  sigma<-init.values[j,4]
  
  #start sampling
  for(i in 1:n.iter){
    #compute the mean and standard deviation of the posterior distribution of a0 (normal)
    sum1<-0
    for (k in 1:N) {
      sum1<-sum1+y[k]-b1*x1[k]-b2*x2[k]
    }
    #mean
    mu01<-(sum1/(sigma^2)+mu00/(tau00^2))/(N/(sigma^2)+1/(tau00^2))
    #standard deviation
    tau01<-1/(N/(sigma^2)+1/(tau00^2))
    #sample from the conditional posterior distribution of a0
    a0<-rnorm(1,mean=mu01,sd=tau01)
    #Save the value 
    post.a0[j,i]<-a0
    
    # Metropolis-Hastings step
    
    #proposal density
    beta.ast<-rnorm(1,mean=b1,sd=0.1)
    
    #sample a random value from a uniform(0,1)
    u<-runif(1,min=0,max=1)
    
    #conditional posterior computations
    sum.x1sq<-0
    
    for (k in 1:N) {
      sum.x1sq<-sum.x1sq+x1[k]^2
    }
    
    sum2<-0
    for (k in 1:N) {
      sum2<-sum2+x1[k]*(y[k]-a0-b2*x2[k])
    }
    
    
    #Compute the acceptance rate. We do it directly to avoid computational problems when the denominator tends to 0
    accept.rate[j,i]<-exp(((-(beta.ast^2)*sum.x1sq/(2*sigma^2))+beta.ast*sum2/sigma^2)-(-(b1^2)*sum.x1sq/(2*sigma^2)+b1*sum2/sigma^2))*((1+((b1-m0)^2)/(nu0*s0^2))^(-(nu0+1)/2))/(((1+((b1-m0)^2)/(nu0*s0^2)))^(-(nu0+1)/2))
    
    #keep the previous value if the computed acceptance 
    #rate is below the random value from the uniform distribution
    #or save the sampled value from the proposed distribution if 
    #the acceptance rate is above the random value from the uniform distribution
    if (accept.rate[j,i]>u){
      b1<-beta.ast
    }else{
      if(i==1){
        b1<-b1
      }else{
        b1<-post.b1[j,i-1]
      }
    }  
    #Save the value
    post.b1[j,i]<-b1
    
    #compute the mean and standard deviation of the posterior distribution of b2 (normal)
    sum3<-0
    for (k in 1:N) {
      sum3<-sum3+(y[k]-a0-b1*x1[k])*x2[k]
    }
    #mean
    mu21<-(sum3/(sigma^2)+mu20/(tau20^2))/(sum(x2^2)/(sigma^2)+1/(tau20^2))
    #standard deviation
    tau21<-1/(sum(x2^2)/(sigma^2)+1/(tau20^2))  
    #sample from the conditional posterior distribution of b2
    b2<-rnorm(1,mean=mu21,sd=tau21)
    #save the value
    post.b2[j,i]<-b2
    
    
    
    #compute the parameters of the conditional posterior distribution of sigma (inverse gamma function)
    a.post<-N/2+a.prior
    
    S<-0
    for (k in 1:N) {
      S<-S+(y[k]-a0-b1*x1[k]-b2*x2[k])^2
    }
    b.post<-S/2+b.prior
    #sample from the conditional posterior distribution of sigma
    sigma<-sqrt(1/rgamma(1,shape=a.post,rate=b.post))
    #save the values
    post.sigma[j,i]<-sigma
  }
  
  #Delete the burn-in period
  for (i in (burn.period+1):n.iter) {
    p.a0[j,(i-burn.period)]<-post.a0[j,i]
    p.b1[j,(i-burn.period)]<-post.b1[j,i]
    p.b2[j,(i-burn.period)]<-post.b2[j,i]
    p.sigma[j,(i-burn.period)]<-post.sigma[j,i]
  }
  
}#end sampling

#compute the acceptance rate for each chain.
a.rate<-rep(0,2)
for(j in 1:n.chains){
  a.rate[j]<-mean(accept.rate[j,burn.period:n.iter])
}

#compute the length of the sampling
l.sample<-dim(p.a0)[2]

#Combine the values of the different chains 
f.a0<-rep(0,l.sample)
f.b1<-rep(0,l.sample)
f.b2<-rep(0,l.sample)
f.sigma<-rep(0,l.sample)
for (i in 1:l.sample) {
  f.a0[i]<-mean(p.a0[,i])
  f.b1[i]<-mean(p.b1[,i])
  f.b2[i]<-mean(p.b2[,i])
  f.sigma[i]<-mean(p.sigma[,i])
}

#compute the length of the new vectors which contain the values of each parameter
l.post<-length(f.a0)

#histograms of posterior distributions which gives us an aproximation of the distributions

post.dist.a0<-hist(f.a0,main = "Posterior distribution of a0")
post.dist.b1<-hist(f.b1,main = "Posterior distribution of b1")
post.dist.b2<-hist(f.b2,main = "Posterior distribution of b2")
post.dist.sigma<-hist(f.sigma,main = "Posterior distribution of sigma")


#Assess convergence

#History plots

#History plot of a0
plot(post.a0[1,],type="l",main = "History plot of a0", xlab = "iteration", ylab="a0",col="red")
for (i in 2:n.chains) {
  lines(post.a0[i,])
}

#History plot of b1
plot(post.b1[1,],type="l",main = "History plot of b1", xlab = "iteration", ylab="b1",col="red")
for (i in 2:n.chains) {
  lines(post.b1[i,])
}

#History plot of b2
plot(post.b2[1,],type="l",main = "History plot of b2", xlab = "iteration", ylab="b2",col="red")
for (i in 2:n.chains) {
  lines(post.b2[i,])
}

#History plot of sigma
plot(post.sigma[1,],type="l",main = "History plot of sigma", xlab = "iteration", ylab="St.deviation",col="red")
for (i in 2:n.chains) {
  lines(post.sigma[i,])
}


#Autocorrelation

#set the vectors
a0.autocor<-rep(0,l.post)
b1.autocor<-rep(0,l.post)
b2.autocor<-rep(0,l.post)
sigma.autocor<-rep(0,l.post)

#set the number of desired lags
n.lag<-100
#set the first correlations to 1
a0.autocor[1]<-1
b1.autocor[1]<-1
b2.autocor[1]<-1
sigma.autocor[1]<-1

for(i in 2:n.lag){
  #compute the lagged vectors
  a0.lag<-c(rep(NA,i),f.a0[i:l.post-i])#add as many NA as desired lags
  b1.lag<-c(rep(NA,i),f.b1[i:l.post-i])#add as many NA as desired lags
  b2.lag<-c(rep(NA,i),f.b2[i:l.post-i])#add as many NA as desired lags
  sigma.lag<-c(rep(NA,i),f.sigma[i:l.post-i])#add as many NA as desired lags
  #compute the correlations between the original vectors containing the parameters
  #and the lagged vectors
  a0.autocor[i]<-cor(f.a0,a0.lag,use = "complete.obs")
  b1.autocor[i]<-cor(f.b1,b1.lag,use = "complete.obs")
  b2.autocor[i]<-cor(f.b2,b2.lag,use = "complete.obs")
  sigma.autocor[i]<-cor(f.sigma,sigma.lag,use = "complete.obs")
}

#Plot the autocorrelations previously computed
barplot(a0.autocor, main="Autocorrelation a0", xlab="lag", ylim = c(-0.3,1), col="blue", xlim = c(0,n.lag))
barplot(b1.autocor, main="Autocorrelation b1", xlab="lag",ylim = c(-0.3,1), col="blue", xlim = c(0,n.lag))
barplot(b2.autocor, main="Autocorrelation b2", xlab="lag",ylim = c(-0.3,1), col="blue", xlim = c(0,n.lag))
barplot(sigma.autocor, main="Autocorrelation Standard deviation", xlab="lag",ylim = c(-0.3,1), col="blue", xlim = c(0,n.lag))



#Estimates and credible intervals

#Mean and 95%CI of parameter a0
a0.mean<-mean(f.a0)#mean
a0.ord<-sort(f.a0)#order the values
CIa0<-quantile(x=a0.ord,probs = c(0.025,0.975))# 95% credible interval

#Mean and 95%CI of parameter b1
b1.mean<-mean(f.b1)#mean
b1.ord<-sort(f.b1)#order the values
CIb1<-quantile(x=b1.ord,probs = c(0.025,0.975))# 95% credible interval

#Mean and 95%CI of parameter b2
b2.mean<-mean(f.b2)#mean
b2.ord<-sort(f.b2)#order the values
CIb2<-quantile(x=b2.ord,probs = c(0.025,0.975))# 95% credible interval

#Mean and 95%CI of parameter sigma
sigma.mean<-mean(f.sigma)#mean
sigma.ord<-sort(f.sigma)#order the values
CIsigma<-quantile(x=sigma.ord,probs = c(0.025,0.975))# 95% credible interval




#DIC compuation

#First, for the model with 2 predictors

#Dhat computation

#set parameter
dhat<-0
#compute the length of the vectors containing the parameters
Q<-length(f.a0)

#compute the dhat/(-2)
for (i in 1:N) {
  dhat<-dhat+((-1*(y[i]-(a0.mean+b1.mean*x1[i]+b2.mean*x2[i]))^2)/(2*sigma.mean^2))-log(sqrt(2*pi)*sigma.mean)
}

#multiply the previous result by -2 to obtain the Dhat
Dhat<-(-2*dhat)


#Dbar

#set the vectors
dbar<-rep(0,Q)
dbar2<-rep(0,Q)
#Previos computation for the dbar
for (q in 1:Q) {
  for (i in 1:N) {
    dbar[q]<-dbar[q]+((-1*(y[i]-(f.a0[q]+f.b1[q]*x1[i]+f.b2[q]*x2[i]))^2)/(2*f.sigma[q]^2))-log(sqrt(2*pi)*f.sigma[q])
    dbar2[q]<-(-2)*dbar[q]
    }
}
#compute the Dbar
Dbar<-sum(dbar2)/Q
#Compute pD
pD<-(Dbar-Dhat)

#Compute the DIC
DIC<-Dhat+2*pD


#DIC one predictors model model (b1 is removed):

#Using the same priors for a0,b2 and b3 we performed again the Gibbs sampler

#Set matrices for the parameters
dic.a0<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
dic.b2<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
dic.sigma<-matrix(rep(0,n.chains*n.iter),nrow = n.chains,ncol=n.iter)
d.a0<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
d.b2<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
d.sigma<-matrix(rep(0,n.chains*(n.iter-burn.period)),nrow = n.chains,ncol=(n.iter-burn.period))
#sample from the posteriors using a Gibbs sampler without the MH step as b1 is not part of the second model
for (j in 1:n.chains) {
  #read initial values
  a0<-init.values[j,1]
  b2<-init.values[j,3]
  sigma<-init.values[j,4]
  #sample a0
  for(i in 1:n.iter){
    #compute the mean and standard deviation of the posterior distribution of a0 
    sum1<-0
    for (k in 1:N) {
      sum1<-sum1+y[k]-b2*x2[k]
    }
    #mean
    mu01<-(sum1/(sigma^2)+mu00/(tau00^2))/(N/(sigma^2)+1/(tau00^2))
    tau01<-1/(N/(sigma^2)+1/(tau00^2)) #standard deviation
    a0<-rnorm(1,mean=mu01,sd=tau01) #sample from the conditional posterior distribution of a0
    dic.a0[j,i]<-a0#Save the value
    
    
    #compute the mean and standard deviation of the posterior distribution of b2
    sum3<-0
    for (k in 1:N) {
      sum3<-sum3+(y[k]-a0)*x2[k]
    }
    #mean
    mu21<-(sum3/(sigma^2)+mu20/(tau20^2))/(sum(x2^2)/(sigma^2)+1/(tau20^2))
    tau21<-1/(sum(x2^2)/(sigma^2)+1/(tau20^2))#standard deviation  
    #sample from the conditional posterior distribution of b2
    b2<-rnorm(1,mean=mu21,sd=tau21)
    dic.b2[j,i]<-b2#save the value
    
    
    #compute the parameters of the posterior distribution of sigma (inverse gamma function)
    a.post<-N/2+a.prior
    
    S<-0
    for (k in 1:N) {
      S<-S+(y[k]-a0-b2*x2[k])^2
    }
    b.post<-S/2+b.prior
    #sample from the conditional posterior distribution of sigma
    sigma<-sqrt(1/rgamma(1,shape=a.post,rate=b.post))
    dic.sigma[j,i]<-sigma #save the value
  }
  
  #Remove the burn-in period.
  for (i in (burn.period+1):n.iter) {
    d.a0[j,(i-burn.period)]<-dic.a0[j,i]
    d.b2[j,(i-burn.period)]<-dic.b2[j,i]
    d.sigma[j,(i-burn.period)]<-dic.sigma[j,i]
  }
  
}#end sampling

l.sample2<-dim(d.a0)[2]

#Combine the values of the different chains 
f.d.a0<-rep(0,l.sample2)
f.d.b1<-rep(0,l.sample2)
f.d.b2<-rep(0,l.sample2)
f.d.sigma<-rep(0,l.sample2)
for (i in 1:l.sample2) {
  f.d.a0[i]<-mean(d.a0[,i])
  f.d.b2[i]<-mean(d.b2[,i])
  f.d.sigma[i]<-mean(d.sigma[,i])
}

#compute the mean values for the dhat
da0.mean<-mean(f.d.a0)
db2.mean<-mean(f.d.b2)
dsigma.mean<-mean(f.d.sigma)

#compute DIC


#Compute Dhat as previously

dhat2<-0
Q2<-length(d.a0)
#Previous computation for the Dhat
for (i in 1:N) {
  dhat2<-dhat2+((-1*(y[i]-(da0.mean+db2.mean*x2[i]))^2)/(2*dsigma.mean^2))-log(sqrt(2*pi)*dsigma.mean)
}

#Resulting Dhat
Dhat2<-(-2*dhat2)

#Compute Dbar as previously

dbar.2<-rep(0,Q2)
dbar2.2<-rep(0,Q2)
for (q in 1:Q2) {
  for (i in 1:N) {
    dbar.2[q]<-dbar.2[q]+((-1*(y[i]-(d.a0[q]+d.b2[q]*x2[i]))^2)/(2*d.sigma[q]^2))-log(sqrt(2*pi)*d.sigma[q])
    dbar2.2[q]<-(-2)*dbar.2[q]
  }
}
#Resulting Dbar
Dbar2<-sum(dbar2.2)/Q2

#Compute the new pD
pD2<-(Dbar2-Dhat2)

#Compute the DIC for the second model
DIC2<-Dhat2+2*pD2




#Bayes factor

# hypothesis: b2>b1 

#First we define new priors for b1 and b2 with same means and variances
#to compute the bayes factor. These new priors are quite similar to the original ones.

#Define the new means
means<-c(1,1)
#Define the variance-covariance matrix of b1 and b2with the new variances
covar<-matrix(c(0.0064,0,0,0.0064),2,2)

#sample from prior to obtain complexity
pre.f1<-mvrnorm(n = n.iter, mu=means, Sigma=covar, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

#compute the complexity by computing the proportion the amount of times the values sampled from the prior agree with our hypothesis 
cont.c1<-0
for (i in 1:n.iter) {
  if(pre.f1[i,2]>pre.f1[i,1]){
    cont.c1<-cont.c1+1
  }
}

#compute the proportion
c1<-cont.c1/n.iter


#Sample from the posteriors using the Gibbs sampler with the new priors

#Define the parameters of the new priors
mu10.f1<-1#mean
tau10.f1<-0.08#variance
mu20.f1<-1#mean
tau20.f1<-0.08#variance


#Initial values. Chosen to converge fast, although it should not be necessary.
a0.f1<-a0.mean
b1.f1<-1
b2.f1<-1
sigma.f1<-sigma.mean
  
#Set the matrices
post.f1.a0<-matrix(rep(0,n.iter),nrow = 1,ncol=n.iter)
post.f1.b1<-matrix(rep(0,n.iter),nrow = 1,ncol=n.iter)  
post.f1.b2<-matrix(rep(0,n.iter),nrow = 1,ncol=n.iter)
post.f1.sigma<-matrix(rep(0,n.iter),nrow = 1,ncol=n.iter)
  
  
p.f1.a0<-matrix(rep(0,(n.iter-burn.period)),nrow = 1,ncol=(n.iter-burn.period))
p.f1.b1<-matrix(rep(0,(n.iter-burn.period)),nrow = 1,ncol=(n.iter-burn.period))
p.f1.b2<-matrix(rep(0,(n.iter-burn.period)),nrow = 1,ncol=(n.iter-burn.period))
p.f1.sigma<-matrix(rep(0,(n.iter-burn.period)),nrow = 1,ncol=(n.iter-burn.period))
  
#Gibbs sampler.Similar to the previous ones.  
for(i in 1:n.iter){
  #compute the mean and standard deviation of the posterior distribution of a0 
    sum1<-0
    for (k in 1:N) {
      sum1<-sum1+y[k]-b1.f1*x1[k]-b2.f1*x2[k]
    }
    #mean
    mu01<-(sum1/(sigma.f1^2)+mu00/(tau00^2))/(N/(sigma.f1^2)+1/(tau00^2))
    tau01<-1/(N/(sigma.f1^2)+1/(tau00^2))#standard deviation
    a0.f1<-rnorm(1,mean=mu01,sd=tau01)#sample from the conditional posterior distribution of a0
    post.f1.a0[i]<-a0.f1#save the value
    
    #compute the mean and standard deviation of the posterior distribution of b1 
    sum2<-0
    for (k in 1:N) {
      sum2<-sum2+(y[k]-a0.f1-b2.f1*x2[k])*x1[k]
    }
    #mean
    mu11.f1<-(sum2/(sigma.f1^2)+mu10.f1/(tau10.f1^2))/(sum(x1^2)/(sigma.f1^2)+1/(tau10.f1^2))
    tau11.f1<-1/(sum(x1^2)/(sigma.f1^2)+1/(tau10.f1^2)) #standard deviation 
    b1.f1<-rnorm(1,mean=mu11.f1,sd=tau11.f1)#sample from the conditional posterior distribution of b1
    post.f1.b1[i]<-b1.f1#save the value
    
    #compute the mean and standard deviation of the posterior distribution of b2
    sum3<-0
    for (k in 1:N) {
      sum3<-sum3+(y[k]-a0.f1-b1.f1*x1[k])*x2[k]
    }
    #mean
    mu21.f1<-(sum3/(sigma.f1^2)+mu20.f1/(tau20.f1^2))/(sum(x2^2)/(sigma.f1^2)+1/(tau20.f1^2))
    tau21.f1<-1/(sum(x2^2)/(sigma.f1^2)+1/(tau20.f1^2))#sd  
    b2<-rnorm(1,mean=mu21.f1,sd=tau21.f1)#sample from the conditional posterior distribution of b2
    post.f1.b2[i]<-b2.f1
    
    #compute the parameters of the posterior distribution of sigma (inverse gamma function)
    a.post<-N/2+a.prior
    S<-0
    for (k in 1:N) {
      S<-S+(y[k]-a0-b1.f1*x1[k]-b2.f1*x2[k])^2
    }
    b.post<-S/2+b.prior
    #sample from the conditional posterior distribution of sigma
    sigma.f1<-sqrt(1/rgamma(1,shape=a.post,rate=b.post))
    post.f1.sigma[i]<-sigma.f1#save the value
}#end sampling
  
  #Remove burn-in period
  
  for (i in (burn.period+1):n.iter) {
    p.f1.a0[(i-burn.period)]<-post.f1.a0[i]
    p.f1.b1[(i-burn.period)]<-post.f1.b1[i]
    p.f1.b2[(i-burn.period)]<-post.f1.b2[i]
    p.f1.sigma[(i-burn.period)]<-post.f1.sigma[i]
  }
  
#compute the length of the sample
l.f1<-length(p.f1.b1)  

#compute the fit by computing the proportion the amount of times the values sampled from the posterior agree with our hypothesis 
cont.f1<-0
for (i in 1:l.f1) {
  if(p.f1.b2[i]>p.f1.b1[i]){
    cont.f1<-cont.f1+1
  }
}
#compute the proportion
f1<-cont.f1/l.f1

#Compute the bayes factor
BF1=f1/c1 #where f1 is the fit and c1 is the complexity

#Compute the posterior model probabilities
PMP1<-BF1/(BF1+1)
PMP2<-1/(BF1+1)


#Generate new datasets under H0

#Set the vectors
y.pred<-rep(0,l.data)
pred.y<-matrix(rep(0,l.post*l.data),nrow = l.data, ncol = l.post)
#Generate the data for the outcome variable
for (n in 1:l.post) {
  for (i in 1:l.data) {
    #the mean of the new y's are the expected values given the previous parameters of the regression
    pred.mean<-f.a0[n]+f.b1[n]*x1[i]+f.b2[n]*x2[i]
    #standard deviation is the one obtained by sampling from the posterior of sigma
    pred.sig<-f.sigma[n]
    #we sample the new y's from a normal distribution with the previous mean and standard deviation
    y.pred[i]<-rnorm(1,pred.mean,pred.sig)
  }
  #save the values
  pred.y[i,n]<-y.pred[i]  
}

#Posterior predictive check

#Apply test to original dataset

#set vector of residuals
res<-rep(0,l.data)
for (i in 1:l.data) {
  #compute the residuals
  res[i]<-y[i]-(f.a0[i]+f.b1[i]*x1[i]+f.b2[i]*x2[i])
}
#compute the standard deviation of the residuals
sd.res<-sd(res)
#set vector of pretest computation, t.res.
t.res<-rep(0,l.data)
#compute the residuals minus 3 times their standard deviation
for (i in 1:l.data) {
  t.res[i]<-res[i]-3*sd.res
}

#Select the maximum value
test1<-max(t.res)




#Apply the test to the generated datasets. 

#set the vectors we are going to use
res2<-matrix(rep(0,l.data*l.post),nrow=l.data, ncol=l.post)
t.res2<-matrix(rep(0,l.data*l.post),nrow=l.data, ncol=l.post)
pred.test<-rep(0,l.post)
sd.res2<-rep(0,l.post)
pred.test1<-rep(0,l.post)
#The following computations may take some time (around 30min).
for (n in 1:l.post) {
  for (i in 1:l.data) {
    #compute the residuals
    res2[i,n]<-pred.y[i,n]-(f.a0[i]+f.b1[i]*x1[i]+f.b2[i]*x2[i])
  }
  #compute the standard deviations of the residuals
  sd.res2[n]<-sd(res2[1:l.data,])
  
  #compute the residuals minus 3 times the standard deviations
  for (j in 1:l.data) {
    t.res2[j,n]<-res2[j,n]-3*sd.res2[n]
  }
  #Select the maximum value of the previous computation
  pred.test1[n]<-max(t.res2[1:l.data,])
}

l.pval<-length(pred.test1)
#posterior predictive p-value
pval<-rep(0,l.pval)
#compute the proportion of times in which the test 
#applied to the generated datasets is greater than the 
#test applied to the original dataset
for (i in 1:l.pval) {
  if (test1<pred.test1[i]){
    pval[i]<-1
  }else{
    pval[i]<-0
  }
  
}

#compute the p-value by obtaining the mean of the previous computation
pvalue<-mean(pval)
#compute the 95%credible interval of the p-value by ordering
#their distribution and selecting the 2.5 and 97.5 percentile
p.ord<-sort(pval)
CIpval<-quantile(x=p.ord,probs = c(0.025,0.975))


#Check potential outliers by checking those values in the original 
#dataset with a higher test value than the highest in the generated datasets
i.pos<-t.res>max(pred.test1)
pos.outliers<-t.res[i.pos]

#compare results with frequentist approach

freq.res<-lm(y~x1+x2)
#confidence intervals
freqCI<-confint(freq.res, level=0.95)
