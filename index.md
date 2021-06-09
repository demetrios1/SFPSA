<img src="{{site.baseurl}}/going_vs_bank_reg-1.png" width="200" height="200" />
![going_vs_bank_reg-1.png]({{site.baseurl}}/going_vs_bank_reg-1.png =100x20)
## Introduction and usage

This code implements the methods from the paper **Do forecasts of bankruptcy cause bankruptcy? A machine learning sensitivity analysis.**, which can be found [here](https://arxiv.org/pdf/2106.04503.pdf).  The abstract of the paper is reproduced here:

> It is widely speculated that auditors’ public forecasts of bankruptcy are,at least in part, self-fulfilling prophecies in the sense that they might actuallycause bankruptcies that would not have otherwise occurred. This conjectureis hard to prove, however, because the strong association between bankrupt-cies and bankruptcy forecasts could simply indicate that auditors are skillfulforecasters with unique access to highly predictive covariates. In this paper,we investigate the causal effect of bankruptcy forecasts on bankruptcy usingnonparametric sensitivity analysis. We contrast our analysis with two alterna-tive approaches: a linear bivariate probit model with an endogenous regressor,and a recently developed bound on risk ratios called E-values. Additionally,our machine learning approach incorporates a monotonicity constraint corre-sponding to the assumption that bankruptcy forecasts do not make bankrupt-cies less likely. Finally, a tree-based posterior summary of the treatment effectestimates allows us to explore which observable firm characteristics moderatethe inducement effect.


While the paper is focused on whether or not going concern opinions affect probability of bankruptcy, these methods are widely applicable to any setting with a binary treatment variable and a binary response variable (though the methods could also be extended to multinomial response as well).  


### Monotone Bart
```markdown
library(foreach)
library(bcf)
library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)
#Set a seed and generate data from the bivariate probit
library(MASS)
set.seed(0)
N <- 500 # Number of random samples
a=1

x1=runif(N, -a,a)
x2=runif(N, -a,a)
x3=runif(N,-a,a)
x4=runif(N,-a,a)
x5=runif(N, -a,a)
beta1= -0.2
alpha1= 0.7
beta0= -0
alpha0= -0.5
mu1 <- beta0+beta1*(x1+x2+x3+x4+x5)
mu2 <- alpha0+alpha1*(x1+x2+x3+x4+x5)
mu<-matrix(c(mu1, mu2), nrow=N, ncol=2)
rho=.5
gamma=1
B1.true=pnorm(mu2+gamma)
B0.true=pnorm(mu2)
sigma <- matrix(c(1, rho,rho,1),
               2) # Covariance matrix
sim_data=t(sapply(1:N, function(i)mvrnorm(1, mu = mu[i,], Sigma = sigma )))
#generate the binary treatments
G=sapply(1:N, function(i)ifelse(sim_data[i,1]>=0, 1,0))
#generate the binary outcomes
B=sapply(1:N, function(i)ifelse(sim_data[i,2]>=-1*gamma*G[i], 1,0))
print(table(G,B))
covariates=data.frame(x1,x2,x3,x4,x5,B, G)
vars=c('x1','x2','x3','x4','x5')
intframe=BARTpred(covariates, treat='G', Outcome='B',vars, mono=T)
```

### Integration
This is the code for the integration.  Right now, you can use bart, monotone bart, or random forest, although more models will likely be added.  We also have a balanced cross-validation function for the monotone bart and bart model.  

This integration takes a mean of the posterior draws and uses that as the passed probability for each obvservation in the integration optimization step.  One could also optimize the entire dataset for all the posterior draws of the joint probabilities of treatment/outcome for each observation.  However, doing so is quite the computational strain, so we do not currently include this feature. 

```markdown
#pick a function to integrate over
#here the standard deviation is chosen to match the generated data
sd_chosen = sqrt(rho/(1-rho))
f=function(u){
dnorm(u, mean=0, sd = sd_chosen) #for the time being, this has to be done as a number not a variable
}
treat_frame=integrate_function(intframe, constraint=T, f=f, n_cores=1, lambda=0)
library(rpart)
#merge the irds and covariates
dataset=data.frame(covariates, diff=treat_frame[, 'diff'])
tree_fit<-rpart(diff~x1+x2+x3+x4+x5,data=dataset, minbucket=100)
rpart.plot::rpart.plot(tree_fit)
####Optional, if you want to make a quick histogram of IRD's
library(tidyverse)
library(dplyr)
dataset%>%
ggplot(aes(x=diff))+
geom_histogram(aes(y=..count../sum(..count..)),color='white',
               fill='black')+
               ggtitle('All observations')+
               ylab('Density')+xlab('Individual Risk Differences')+
               xlim(0,.50)+theme_minimal(base_size = 16)+
               theme(plot.title = element_text(hjust = 0.5,size=16))


```



### Support or Contact

Additional detail can be found [here](https://github.com/demetrios1/bankruptcy_sensitivity), where more codes and a supplemental file is available.
