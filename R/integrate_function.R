#' The main integration function
#'
#' This function does the integration from the paper
#' Here, we pass the mean of the posterior estimates of observed joint probabilities
#' rather than the full posterior estimates.
#' @name integrate_function
#' @param intframe the dataframe of interest
#' @param constraint boolean, whether or not you want to implement b1>=b0 in integration
#' @param f the function for the distribution of f(u)
#' @param n_cores the number of cores (note detect_cores only gives total number of cores, not necessarily available)
#' @param lambda regularization term in constraint, multiplied by exp(b1)+b0-b0)^2
#' @keywords cats
#' @export
#' @examples
#' #Set a seed and generate data from the bivariate probit
#'set.seed(0)
#'
#' N <- 500 # Number of random samples
# set.seed(123) #for variation in the x's
# Target parameters for univariate normal distributions
#' a=1
#'
#' x1=runif(N, -a,a)
#' x2=runif(N, -a,a)
#' x3=runif(N,-a,a)
#' x4=runif(N,-a,a)
#' x5=runif(N, -a,a)

#' beta1= -0.2
#' alpha1= 0.7
#' beta0= -0
#' alpha0= -0.5
#' mu1 <- beta0+beta1*(x1+x2+x3+x4+x5)
#' mu2 <- alpha0+alpha1*(x1+x2+x3+x4+x5)
#' mu<-matrix(c(mu1, mu2), nrow=N, ncol=2)
#' rho=.5
#' gamma=1
#' B1.true=pnorm(mu2+gamma)
#' B0.true=pnorm(mu2)
#' sigma <- matrix(c(1, rho,rho,1),
#'                2) # Covariance matrix


#' sim_data=t(sapply(1:N, function(i)MASS::mvrnorm(1, mu = mu[i,], Sigma = sigma )))
#' #generate the binary treatments
#' G=sapply(1:N, function(i)ifelse(sim_data[i,1]>=0, 1,0))

#' #generate the binary outcomes
#' B=sapply(1:N, function(i)ifelse(sim_data[i,2]>=-1*gamma*G[i], 1,0))
#' print(table(G,B))
#' covariates=data.frame(x1,x2,x3,x4,x5,B, G)

#' vars=c('x1','x2','x3','x4','x5')
#' intframe=BARTpred(covariates, treat='G', Outcome='B',vars, mono=T)
#' #pick a function to integrate over
#' #here the standard deviation is chosen to match the generated data
#' f=function(u){
#'  dnorm(u, mean=0, sd=sqrt(rho/(1-rho)))
#' }
#' treat_frame=integrate_function(intframe, constraint=T, f=f, n_cores=1, lambda=0)
#' library(rpart)
#' #merge the ites and covariates
#' dataset=data.frame(covariates, tau=treat_frame[, 'tau'])
#' tree_fit<-rpart(tau~x1+x2+x3+x4+x5,data=dataset)
#' rpart.plot::rpart.plot(tree_fit)
#' ####Optional, if you want to make a quick histogram of ITE's
#' library(tidyverse)
#' library(dplyr)
#' dataset%>%
#' ggplot(aes(x=tau))+
#' geom_histogram(aes(y=..count../sum(..count..)),color='white',
#' fill='#1d2951')+
#' ggtitle('All observations')+
#' ylab('Density')+xlab('Individual Treatment Effects')+
#' xlim(0,.25)+theme_minimal(base_size = 16)+
#' theme(plot.title = element_text(hjust = 0.5,size=16))


integrate_function=function(intframe, constraint=T, f=f, n_cores=n_cores, lambda=0){
  n_cores=n_cores
  registerDoParallel(cores = n_cores)
  set.seed(12296)
  mBG_out<-c()


  newoutcome=as.matrix(intframe)

  ptm = proc.time()
  mBG_out = foreach( i = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
    #better global definition.
    optimdist=function(vals){
      b1 = vals[1]
      b0 = vals[2]
      g = vals[3]
      a=integrate( function(u) pnorm(b1+u)*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
      b=integrate(function(u) pnorm(b0+u)*(1-pnorm(g+u))*f(u), lower = -Inf, upper = Inf )$value
      c=integrate(  function(u) (1-pnorm(b1+u))*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
      #return(c(a,b,c))
      #vProbBG is lhs, global variable
      return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
    }
    optimdist_constraint=
      function(vals){
        b1 = vals[1]
        b0 = vals[2]
        g = vals[3]
        #replace b1 with b1+b0
        a=integrate( function(u) pnorm((exp(b1)+b0)+u)*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
        b=integrate(function(u) pnorm(b0+u)*(1-pnorm(g+u))*f(u), lower = -Inf, upper = Inf )$value
        c=integrate(  function(u) (1-pnorm((exp(b1)+b0)+u))*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
        #return(c(a,b,c))
        #vProbBG is lhs, global variable
        return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2)+lambda*((exp(b1)+b0) -b0)^2)
      }
    #calculate the treatment
    treatment_Compute = function( vb){
      integrate( function(u) (pnorm(vb[1]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value

    }
    treatment_Compute_constraint = function( vb){
      integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value

    }
    counterfac_Compute = function(vb){
      b1 =(integrate( function(u) (pnorm(vb[1]+u))*f(u), lower = -Inf, upper = Inf )$value)
      b0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)

      return(c(b1, b0))

    }

    counterfac_Compute_constraint = function(vb){
      b1 =(integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
      b0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)

      return(c(b1, b0))

    }


    # Probability to fit.
    vProbBG     = newoutcome[i,c( "ProbB1G1", "ProbB1G0", "ProbB0G1" )]

    # Starting value

    start0 = c( newoutcome[i,"ProbB1G1"], newoutcome[i,"ProbB1G0"],newoutcome[i,"ProbB0G1"] )

    names(start0) = c("b1","b0","g")




    if(constraint==T){
      vBG_out     = try(optim( qnorm(start0), optimdist_constraint,  control = list(maxit = 500)
      )
      , silent=TRUE)
      if ( class( vBG_out ) == "try-error" ){
        out        = as.numeric(newoutcome[i,c("indexes")])
        out        = c( out, rep(NA,5) )
      } else {
        out        = as.numeric(newoutcome[i,c("indexes")] )
        vBG        = vBG_out$par
        # tau      = treatment_Compute( vb = vBG[c(1,2)])

        #  counter=counterfac_Compute(vb = vBG[c(1,2)])

        tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])

        counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])

        out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter)
      }
      return( out)
    }else{
      vBG_out     = try(optim( qnorm(start0), optimdist,  control = list(maxit = 500)
      )
      , silent=TRUE)
      if ( class( vBG_out ) == "try-error" ){
        out        = as.numeric(newoutcome[i,c("indexes")])
        out        = c( out, rep(NA,5) )
      } else {
        out        = as.numeric(newoutcome[i,c("indexes")] )
        vBG        = vBG_out$par
        tau      = treatment_Compute( vb = vBG[c(1,2)])

        counter=counterfac_Compute(vb = vBG[c(1,2)])

        # tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])

        #counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])

        out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter)
      }
      return( out)
    }
  }

  print( proc.time() - ptm )

  colnames(mBG_out) = c( "b1", "b0", "g", "tau", "convergence", "fnvalue" ,'B1','B0')

  return(mBG_out)
}
