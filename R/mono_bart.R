#' @import foreach doParallel dbarts bcf monbart
#' @import dplyr ggplot2 caret
NULL



#' Monotone Bart
#'
#' Implements bart with a monotonicity constraint, i.e. $b_1(x)> b_0(x)$
#' @name monotone_bart_function
#' @param y observed outputs
#' @param z 1 if y==0, 0 o.w.
#' @param x Dataframe of covariates
#' @param xpred Dataframe of covariates we predict on (the test set)
#' @param nskip Number of burn in draws
#' @param ndpost Number of posterior draws we keep
#' @param m Number of trees
#' @keywords BART
#' @export
#' @examples
#' monotone_bart_function()

set.seed(12296)


monotone_bart_function = function(y, z, x, xpred, nskip=5000, ndpost=5000, m = 50) {

  sort_ix = order(z, y)
  x = x[sort_ix,]
  z = z[sort_ix]
  y = y[sort_ix]
  n = nrow(x)
  n0 = sum(z==0)
  n00 = sum(z==0 & y==0)
  yobs = y
  yobs0 = rbinom(n, 1, 0.5)
  yobs[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
  yobs0[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
  yobs0[1:n00][yobs[1:n00]==1] = 0 # To statisfy DA constraints
  yobs0[z==0 & y==1] = 1
  yobs0 = yobs0[1:n0]

  offset =  0#qnorm(mean(y[z==1]))
  offset0 = 0# qnorm(mean(flu$grp==0 & flu$fluy2==1)/mean(flu$grp==1 & flu$fluy2==1)) #<- wtf

  zz = offset + 3*yobs - 3*(1-yobs)
  z0 = offset0 + 3*yobs0 - 3*(1-yobs0)

  ################################################################################
  # MCMC
  ################################################################################
#the seed is set in the data generating process
 # set.seed(1022)


  xi = lapply(1:ncol(x), function(i) bcf:::.cp_quantile(x[,i]))
  fit.mono = bartRcppMono(yobs, zz, t(as.matrix(x)), t(xpred),
                          yobs0, z0, t(as.matrix(x)),t(xpred),
                          n00,
                          xi,
                          nskip, ndpost, m, 3.0,
                          offset, offset0)

  #xpred.exp = rbind(data.frame(xpred, z=1), data.frame(xpred, z=0))
  #fit =  bart(cbind(x, z), y, xpred.exp,
  #               nskip=5000, ndpost=5000,
  #               ntree=50, usequants=T)


  # P(Y|Z=1, X)
  pr1 = pnorm(fit.mono$postpred)

  # P(Y|Z=0, X)
  pr0 = pr1*pnorm(fit.mono$postpred0)

  return(list(pr0 = pr0, pr1 = pr1))

}
#' A shark fin generator
#' Generating from the shark fin distribution
#' @name shark_fin
#' @param z vector of quantiles
#' @param q the q-parameter
#' @param varexp the variance you want
#' @keywords shark_fin
#' @export
#' @examples
#' n=1000
#' xgrid = seq(-4, 4, length.out=n)
#' q=.5
#' var=1
#' plot(xgrid,shark_fin(xgrid, q, var), type='l', lwd=2,col='firebrick4',
#'    main='Sharkfin plot', xlab='u', ylab='f(u)')
#' #draw from a sharkfin
#'   ind = rbinom(n,1,q)
#'   u = rep(NA,n)
#' a=(1-(dnorm(0)/(1-pnorm(0)))^2)
#' b=((dnorm(0)/(1-pnorm(0))))
#' sig=sqrt(var/((1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)))
#' u[ind==1] = truncnorm::rtruncnorm(sum(ind==1),a = -Inf, b = 0, sd = sig)
#' u[ind==0] = truncnorm::rtruncnorm(sum(ind==0),a = 0, b = +Inf,sd = sig*(1-q)/q)  #prob 1-q over 0, 1/s and s cancel
#'plot(density(u), lwd=2)
#'lines(xgrid, shark_fin(xgrid, q, var), type='l', lwd=2, col='firebrick4')

shark_fin = function(z, q, var){

  a=(1-(dnorm(0)/(1-pnorm(0)))^2)
  b=((dnorm(0)/(1-pnorm(0))))
  varexp=(1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)

  sig<-sqrt(var/varexp)

  val = rep(NA, length(z))
  val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
  val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
  return(val)

}


#' A helper function.  Allows you to fit observed joint probabilities with BART, or with a monotonicity constraint
#' @name BARTpred
#' @param df dataframe including the covariates plus the treatment and outcome columns
#' @param treat String of name of treatment column
#' @param Outcome string, name of outcome column
#' @param vars list of names of columns of covariates
#' @param mono boolean, whether or not to impose monotonicity constraint,
#' use if monbart is installed
#' @param nd_post number of posterior draws to keep
#' @param n_skip number of burn in samples
#' @keywords BART
#' @export
#' @examples
#'#' N <- 500 # Number of random samples
#set.seed(123) #for variation in the x's
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
#' CV_pred=BARTpred_CV(covariates, treat='G', Outcome='B', vars=vars,mono=T,
#'  nd_post=200, n_skip=500, nfold=5)
#'
BARTpred=function(df, treat='G', Outcome='B',vars, mono=T, nd_post=2000, n_skip=2000){

 covtreat=df[df[treat]==1,]
 covcontrol=df[df[treat]==0,]
 covtreat[treat]=as.factor(covtreat[treat])
  #case 1, the treatment
  if (mono==F){
    bart1=bart(as.matrix(covtreat[vars]),as.factor(unlist(covtreat[Outcome])),
               as.matrix(df[vars]),ndpost = nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)

    #case 2 control
    bart2=bart(as.matrix(covcontrol[vars]),as.factor(unlist(covcontrol[Outcome])),
               as.matrix(df[vars]),ndpost =nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)

    #case 3 propensity
    bart3=bart(as.matrix(df[vars]), as.factor(unlist(df[treat])),
               as.matrix(df[vars]),ndpost = nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
    #PB1G1
    pred1 = colMeans(pnorm(bart1$yhat.test))
    #PB1G0
    pred2=colMeans(pnorm(bart2$yhat.test))
    #PG1
    pred3=colMeans(pnorm(bart3$yhat.test))
  }
  else if(mono==T){
    xtest=df[vars]
    xtraincontrol=covcontrol[vars]
    xtraintreat=covtreat[vars]
    ytrain0    = as.factor( unlist(covcontrol[Outcome] ))
    ytrain1    =as.factor( unlist(covtreat[Outcome]))

    # mono fits

    bart_mono = monotone_bart_function(y = as.numeric(c(ytrain1, ytrain0)==1),
                              z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                              x = rbind(xtraintreat, xtraincontrol),
                              xpred = xtest, nskip = n_skip, ndpost = nd_post,m=100)


    bart3 = bart(x.train=as.matrix(df[vars]), y.train=as.factor(unlist(df[treat])),
                 x.test=as.matrix(df[vars]), ndpost = nd_post, nskip = n_skip, ntree=100, usequants=T,  verbose=F)


    #use this method for prediction on binary
    #PB1G1
    pred1= 1-colMeans(bart_mono$pr0)
    #PB1G0
    pred2= 1-colMeans(bart_mono$pr1)
    #PG1
    pred3=colMeans(pnorm(bart3$yhat.test))
  }

  else{
    print('mono argument is T or F Boolean')
  }




  expoutcomesfun	   =cbind( df[,], data.frame(  treatpred = pred1, notreatpred=pred2), propensity=pred3 )

  expoutcomesfun2    = rbind( NULL, expoutcomesfun )
  outcomenew=expoutcomesfun2[,c('propensity','notreatpred','treatpred')]


  ####need the joints####

  eps     = 1e-6

  outcomenew[,c("prop", "B_G0", "B_G1")] = 0.5*eps + (1-eps)*outcomenew[,c('propensity','notreatpred','treatpred')]

  outcomenew[,"ProbB1G1"]=outcomenew[,"prop"]*outcomenew[,"treatpred"]
  outcomenew[,"ProbB1G0"] = (1-outcomenew[,"prop"])*outcomenew[,"notreatpred"]
  outcomenew[,"ProbB0G1"] = outcomenew[,"prop"]*(1-outcomenew[,"treatpred"])

  #for error analysis, not super necessary
  indexes=seq(from=1, to=length(outcomenew$treatpred), by=1)
  outcomenew=as.data.frame(cbind(indexes, outcomenew))
  outcomenew$indexes=as.numeric(outcomenew$indexes)
  row.names(outcomenew)<-NULL
  newoutcome4 =as.matrix(outcomenew)
  return(outcomenew)
}


#' A helper function.  Allows you to fit observed joint probabilities with linear model, or with a monotonicity constraint
#' also have an option for random forest
#'
#' This function allows you to express your love of cats.
#' @name prediction_function
#' @param df The data you wanna pass.  Should have your covariates, treatment column, and outcome column
#' @param treat the treatment column name (pass as a string)
#' @param Outcome the outcome column name (pass as a string)
#' @param vars the column names of the covariates you wanna pass, should all be in the df dataframe
#' @param model options include random forest, logit
#' @keywords random forest
#' @export
#' @examples
#' set.seed(0)
#' N <- 500 # Number of random samples
#' a=1
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
#' sigma <- matrix(c(1, rho,rho,1), 2) # Covariance matrix
#' sim_data=t(sapply(1:N, function(i)MASS::mvrnorm(1, mu = mu[i,], Sigma = sigma )))
#' #generate the binary treatments
#' G=sapply(1:N, function(i)ifelse(sim_data[i,1]>=0, 1,0))
#' #generate the binary outcomes
#' B=sapply(1:N, function(i)ifelse(sim_data[i,2]>=-1*gamma*G[i], 1,0))
#' print(table(G,B))
#' covariates=data.frame(x1,x2,x3,x4,x5,B, G)
#' vars=c('x1','x2','x3','x4','x5')
#' prediction_function(covariates, 'G', 'B', vars, model='random forest')
prediction_function=function(df, treat='G', Outcome='B',vars, model='random forest'){
  covtreat=df[df[treat]==1,]
  covcontrol=df[df[treat]==0,]
  covtreat[treat]=as.factor(covtreat[treat])
  #case 1, the treatment
  if (model=='random forest'){

    mod1=randomForest::randomForest(covtreat[vars],as.factor(unlist(covtreat[Outcome])) )
    #print(mod1)
    #PB1G1
    pred1=predict(mod1, df[vars], type='prob')[,1]
    #case 2 control
    mod2=randomForest::randomForest(covcontrol[vars],as.factor(unlist(covcontrol[Outcome])) )
    #print(mod2)
    #PB1G0
    pred2=predict(mod2, df[vars], type='prob')[,1]
    #case 3 propensity
    mod3=randomForest::randomForest(df[vars],as.factor(unlist(df[treat])) )
    #print(mod3)
    #PG1


    pred3=predict(mod3, df[vars], type='prob')[,1]

  }
  else if(model=='logit'){
  print('do later')
  }

  else{
    print('This model is not yet available')
  }




  expoutcomesfun	   =cbind( df[,], data.frame(  treatpred = pred1, notreatpred=pred2), propensity=pred3 )

  expoutcomesfun2    = rbind( NULL, expoutcomesfun )
  outcomenew=expoutcomesfun2[,c('propensity','notreatpred','treatpred')]


  ####need the joints####

  eps     = 1e-6

  outcomenew[,c("prop", "B_G0", "B_G1")] = 0.5*eps + (1-eps)*outcomenew[,c('propensity','notreatpred','treatpred')]

  outcomenew[,"ProbB1G1"]=outcomenew[,"prop"]*outcomenew[,"B_G1"]
  outcomenew[,"ProbB1G0"] = (1-outcomenew[,"prop"])*outcomenew[,"B_G0"]
  outcomenew[,"ProbB0G1"] = outcomenew[,"prop"]*(1-outcomenew[,"B_G1"])

  #for error analysis, not super necessary
  indexes=seq(from=1, to=length(outcomenew$B_G1), by=1)
  outcomenew=as.data.frame(cbind(indexes, outcomenew))
  outcomenew$indexes=as.numeric(outcomenew$indexes)
  row.names(outcomenew)<-NULL
  newoutcome4 =as.matrix(outcomenew)
  return(outcomenew)
}




#' A helper function.  Allows you to fit observed joint probabilities with BART, or with a monotonicity constraint
#' @name BARTpred_CV
#' @param df dataframe including the covariates plus the treatment and outcome columns
#' @param treat String of name of treatment column
#' @param Outcome string, name of outcome column
#' @param vars list of names of columns of covariates
#' @param mono Boolean, whether or not to use monotonicity constraint
#' @param nfold how many cross validation folds you want, integer value
#' @param nd_post number of posterior draws to keep
#' @param n_skip number of burn in samples
#' @keywords BART
#' @export
#' @examples
#' N <- 500 # Number of random samples
#set.seed(123) #for variation in the x's
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
#' CV_pred=BARTpred_CV(covariates, treat='G', Outcome='B', vars=vars,mono=T,
#'  nd_post=200, n_skip=500, nfold=5)
#' output=do.call('rbind', CV_pred$imp_frame)
#' cc=ggrocs(list(PB1G0=roc(output$B[output$G==0],
#'                    output$BG0[output$G==0]),
#'               PB1G1=roc(output$B[output$G==1],
#'                         output$BG1[output$G==1]),
#'               PG1=roc(output$G[],
#'                         output$G1)),
#'          breaks = seq(0,1,0.1),
#'          tittle = "ROC perf. predicting our probs")
#' #pick a function to integrate over
#' #here the standard deviation is chosen to match the generated data
#' f=function(u){
#'   dnorm(u, mean=0, sd=sqrt(rho/(1-rho)))
#' }
#' treat_frame=integrate_function(as.matrix(output), constraint=T, f=f, n_cores=1)
BARTpred_CV=function(df, treat='G', Outcome='B',vars,mono=T, nd_post=20, n_skip=20, nfold=5){
  fold=nfold
  covtreat=df[df[treat]==1,]
  covcontrol=df[df[treat]==0,]
  covtreat[treat]=as.factor(covtreat[treat])
#  x_here=covcontrol
#  x_here_2=covtreat
 # y_here=as.factor(df[df[treat]==0, Outcome])
#  y_here_2=as.factor(df[df[treat]==1, Outcome])

  df$y_here=ifelse(df[Outcome]==1&df[treat]==0, 3,
                   ifelse(df[Outcome]==0&df[treat]==0, 2,
                      ifelse(df[Outcome]==1&df[treat]==1, 1,0)))
  test_index <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                           y_1=split(sample(which(df$y_here==3)), 1:fold),
                           y_0=split(sample(which(df$y_here==2)), 1:fold))

  test_index_2 <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                             y_1=split(sample(which(df$y_here==1)), 1:fold),
                             y_0=split(sample(which(df$y_here==0)), 1:fold))


  y_train <- list()
  y_test <- list()


  mbart_train_pred <- list()
  mbart_test_pred <- list()

  mbart_train_pred_prob <- list()
  mbart_test_pred_prob <- list()
  pred1<-list()
  pred2<-list()

  #### CV: filter with ROC in parallel way
  y_train<-c()
  y_test<-c()
  y_train0<-c()
  y_train1<-c()
  y_test0<-c()
  y_test1<-c()
  true_test<-c()
  xtest<-c()
  bart_train_pred_prob<-list()
  bart_test_pred_prob<-list()

  imp_frame<-list()
  for (cv in 1:fold) {
    train_x=t(df)[vars,-test_index[[cv]] ]
    test_x=t(df)[vars,test_index[[cv]] ]
    x_train=sapply(data.frame(t(train_x)), as.numeric)
    x_test=sapply(data.frame(t(test_x)), as.numeric)
    train_x_0 = t(df)[vars,-test_index[[cv]] ]
    test_x_0 = t(df)[vars,test_index[[cv]] ]

    train_y_0 = as.factor(df[, Outcome])[-test_index[[cv]]]
    test_y_0 = as.factor(df[, Outcome])[test_index[[cv]]]
    train_x_1 = t(df)[vars,-test_index_2[[cv]] ]
    test_x_1 = t(df)[vars,test_index_2[[cv]] ]
    train_y_1 =as.factor(df[, Outcome])[-test_index_2[[cv]]]
    test_y_1 =as.factor(df[, Outcome])[test_index_2[[cv]]]
    x_train_0=data.frame(t(train_x_0))
    x_test_0=data.frame(t(test_x_0))
    x_train_1=data.frame(t(train_x_1))
    x_test_1=data.frame(t(test_x_1))
    y_train0[[cv]] = train_y_0
    y_test0[[cv]] = test_y_0
    y_train1[[cv]] = train_y_1
    y_test1[[cv]] = test_y_1
    train_y =as.factor(df[, treat])[-test_index[[cv]]]
    test_y = as.factor(df[, treat])[test_index[[cv]]]
    y_train[[cv]] = train_y
    y_test[[cv]] = test_y
   # ytrain0    = as.factor( unlist(covcontrol[Outcome] ))
  #  ytrain1    = as.factor( unlist(covtreat[Outcome] ))

    xtest[[cv]]=df[c(test_index[[cv]],test_index_2[[cv]]), vars]
    x_train_0<-sapply(x_train_0, as.numeric)
    x_train_1<-sapply(x_train_1, as.numeric)
  if(mono==T){
    bart_mono = monotone_bart_function(y = as.numeric(c(y_train1[[cv]], y_train0[[cv]])==1),
                              z = 1-c(rep(1, length(y_train1[[cv]])),
                                      rep(0, length(y_train0[[cv]]))),
                              x = rbind(x_train_1, x_train_0),
                              xpred = xtest[[cv]], nskip = n_skip, ndpost = nd_post,m=100)
    bart_fit=bart(x_train, y_train[[cv]],xtest[[cv]],ndpost = nd_post, nskip = n_skip,ntree=100,
                  verbose=T,usequants = TRUE,numcut = 1000)


    bart_train_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.train))
    bart_test_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.test))

    pred1[[cv]]=1-colMeans(bart_mono$pr0)

    pred2[[cv]]= 1-colMeans(bart_mono$pr1)
  }else if(mono==F){
    bart_fit_1=bart(x_train_0, y_train0[[cv]],xtest[[cv]],ndpost = nd_post, nskip = n_skip,ntree=100,
                    verbose=T,usequants = TRUE,numcut = 1000)
    bart_fit_2=bart(x_train_1, y_train1[[cv]],xtest[[cv]],ndpost = nd_post, nskip = n_skip,ntree=100,
                    verbose=T,usequants = TRUE,numcut = 1000)
    bart_fit=bart(x_train, y_train[[cv]],xtest[[cv]],ndpost = nd_post, nskip = n_skip,ntree=100,
                  verbose=T,usequants = TRUE,numcut = 1000)


    bart_train_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.train))
    bart_test_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.test))

    pred1[[cv]]=colMeans(pnorm(bart_fit_1$yhat.test))

    pred2[[cv]]= colMeans(pnorm(bart_fit_2$yhat.test))

  }
    true_test[[cv]]=df[c(test_index[[cv]],test_index_2[[cv]]), c(Outcome, treat)]
    imp_frame_2<-data.frame(true_test[[cv]], BG1=pred1[[cv]],
                               BG0=pred2[[cv]], G1=bart_test_pred_prob[[cv]])
    expoutcomesfun2    = rbind( NULL, imp_frame_2 )

    outcomenew=expoutcomesfun2[,c('B', 'G', 'BG1','BG0','G1')]


    ####need the joints####

    eps     = 1e-6

    outcomenew[,c("prop", "B_G0", "B_G1")] = 0.5*eps +
      (1-eps)*outcomenew[,c('BG1','BG0','G1')]

    outcomenew[,"ProbB1G1"]=outcomenew[,"prop"]*outcomenew[,"B_G1"]
    outcomenew[,"ProbB1G0"] = (1-outcomenew[,"prop"])*outcomenew[,"B_G0"]
    outcomenew[,"ProbB0G1"] = outcomenew[,"prop"]*(1-outcomenew[,"B_G1"])

    #for error analysis, not super necessary
    indexes=seq(from=1, to=length(outcomenew$B_G1), by=1)
    outcomenew=as.data.frame(cbind(indexes, outcomenew))
    outcomenew$indexes=as.numeric(outcomenew$indexes)
    row.names(outcomenew)<-NULL
    #newoutcome4 =as.matrix(outcomenew)
    imp_frame[[cv]]<-outcomenew
  }


BG1_mbart=pROC::auc(imp_frame[[cv]][Outcome][imp_frame[[cv]][treat]==1],
              imp_frame[[cv]]$BG1[imp_frame[[cv]][treat]==1])
BG0_mbart=pROC::auc(imp_frame[[cv]][Outcome][imp_frame[[cv]][treat]==0],
              imp_frame[[cv]]$BG0[imp_frame[[cv]][treat]==0])
G1_mbart=pROC::auc(unlist(imp_frame[[cv]][treat]),
              imp_frame[[cv]]$G1)
cat(paste0('AUC for BG1 ',round(BG1_mbart[1],3), ' For cv run', cv ))
cat(paste0('AUC for BG0 ',round(BG0_mbart[1],3), ' For cv run', cv ))
cat(paste0('AUC for G1 ',round(G1_mbart[1],3), ' For cv run', cv ))


  print(cv)
  return(list(true_test = unlist(true_test),
              true_train_0=unlist(y_train0),
              true_test_0=unlist(y_test0),
              true_train_1=unlist(y_train1),
              true_test_1=unlist(y_test1),
              pred1=unlist(pred1),
              pred2=unlist(pred2),
              imp_frame=sapply(1:5, function(cv)list(imp_frame[[cv]]))
  ))
}






#' Functions plots multiple 'roc' objects into one plot
#' @name ggrocs
#' @param rocs
#'   A list of 'roc' objects. Every list item has a name.
#' @param breaks
#'   A vector of integers representing ticks on the x- and y-axis
#' @param legentTitel
#'   A string which is used as legend titel
#' @export
ggrocs <- function(rocs, breaks = seq(0,1,0.1), tittle = "Fit Model") {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {

    # Store all sensitivities and specifivities in a data frame
    # which an be used in ggplot
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(1-rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        auc = rep(sprintf("%.3f",rocs[[rocName]]$auc), length(rocs[[rocName]]$sensitivities)),
        trials = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringsAsFactors = FALSE
      )
    })


    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = trials)) +
      scale_color_manual(values=c('dodgerblue4', 'firebrick4', 'darkgreen'))+
      #scale_colour_viridis(discrete = TRUE)+
      geom_line(size = 1.5, alpha = 0.4) +
      geom_segment(aes(x = 0, y = 0, xend = 1,yend = 1), alpha = 0.4, colour = "gray") +
      geom_step() +
      scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,1), breaks = breaks) +
      scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), breaks = breaks) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))+

      ggtitle(tittle) +
      guides(fill = guide_legend(title=NULL)) +
      coord_equal()


    # auc_table <- unique(RocVals[, c("trials", "auc")])
    #  df.table <- gridExtra::tableGrob(auc_table, theme = gridExtra::ttheme_default(base_size = 8), rows = NULL)
    # Plot chart and table into one object
    return(gridExtra::grid.arrange(rocPlot,# df.table,
                        ncol=2,
                        as.table=TRUE,
                        widths=c(8,1)))
  }
}
