## Introduction and usage
This code implements the methods from the paper Do forecasts of bankruptcy cause bankruptcy? A machine learning sensitivity analysis.
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
#pick a function to integrate over
#here the standard deviation is chosen to match the generated data
sd_chosen = sqrt(rho/(1-rho))
f=function(u){
dnorm(u, mean=0, sd = 1) #for the time being, this has to be done as a number not a variable
}
treat_frame=integrate_function(intframe, constraint=T, f=f, n_cores=1, lambda=0)
library(rpart)
#merge the irds and covariates
dataset=data.frame(covariates, diff=treat_frame[, 'diff'])
tree_fit<-rpart(diff~x1+x2+x3+x4+x5,data=dataset, minucket=100)
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

[Link](url) and ![Image](src)
For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Integration
This is the code for the integration.  Right now, you can use bart, monotone bart, or random forest, although more models will likely be added.  We also have a balanced cross-validation function for the monotone bart and bart model.  

This integration takes a mean of the posterior draws and uses that as the passed probability for each obvservation in the integration optimization step.  One could also optimize the entire dataset for all the posterior draws of the joint probabilities of treatment/outcome for each observation.  However, doing so is quite the computational strain, so we do not currently include this feature. 
Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/demetrios1/Causallysensitive/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
