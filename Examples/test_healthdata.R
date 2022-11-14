library(devtools)
#devtools::install_github("demetrios1/SFPSA", ref="main")
library(SFPSA)
SFPSA::load.package()
#devtools::install_github("jaredsmurray/monbart", ref='main')
#library(monbart)

library(ISLR)

data("Wage")

# here we ask whether or not having health insurance leads to probability of being healthy
Wage2=Wage%>%
  mutate(healthy = ifelse(health%in% '2. >=Very Good', 1,0),
         health_ins = ifelse(health_ins %in% '1. Yes', 1,0),
         ID=seq(from=1, to=nrow(Wage)))%>%
  tidyr::spread(race, race)%>%
  tidyr::spread(jobclass, jobclass)%>%
  tidyr::spread(education, education)%>%
  tidyr::spread(year, year)%>%
  tidyr::spread(maritl, maritl)%>%
  dplyr::select(-c(region, health))
#remove the continuous covariates and add back later
Wage3 = Wage2 %>%
  dplyr::select(-c(age,health_ins, logwage, wage, healthy, ID))

#turn NA's into 1's
Wage3 <- +(!is.na(Wage3))
colnames(Wage3) <- c('White', 'Black', 'Asian', 'Other', 'Industrial',
                     'Information', 'LessHS', 'HS', 'SomeCollege', 'College', 'Advanced','2003',
                     '2004', '2005', '2006', '2007', '2008', '2009', 'NeverMarry', 'Married',
                     'Widowed', 'Divorced', 'Separated')

Wage4 <- cbind(Wage2%>%dplyr::select(c(age,health_ins, logwage, wage, healthy, ID)),
               Wage3)%>%
  dplyr::select(-c(ID, wage))


print(table(Wage4$health_ins,Wage4$healthy))

'%!in%' <- function(x,y)!('%in%'(x,y))
vars=colnames(Wage4)[colnames(Wage4)%!in% c('health_ins', 'healthy')]


intframe=SFPSA::BARTpred(Wage4, treat='health_ins', Outcome='healthy',
                  vars, mono=T)
#full train ROC
pROC::auc( Wage4$health_ins, intframe$propensity)
#intframe=prediction_function(Wage4, treat='health_ins', Outcome='healthy',vars)
#sd_chosen = sqrt(rho/(1-rho))
sd_chosen = .25
f=function(u){
  dnorm(u, mean=0, sd = 0.25) #for the time being, this has to be done as a number not a variable
}
treat_frame=integrate_function(intframe, constraint=T, f=f, n_cores=detectCores()-1, lambda=0)
head(treat_frame)
#library(rpart)

covariates <- Wage4[,vars]

#merge the irds and covariates
dataset=data.frame(covariates, outcome=treat_frame[, 'B1']-treat_frame[, 'B0'])
dataset
fit_tree(df=dataset, minbucket=500)

####Optional, if you want to make a quick histogram of IRD's
library(ggplot2)
dataset%>%
  ggplot(aes(x=outcome))+
  geom_histogram(aes(y=..count../sum(..count..))
                 ,color='white'
                 ,fill='#1d2951', bins=25)+
  ggtitle('All observations')+
  ylab('Density')+xlab('Individual Risk Differences')+theme_minimal(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5,size=14))

intframe


#compare to E-value...our nemesis...
#and the bivariate probit
library(EValue)
library(fmsb)
result_table <- table(Wage4$health_ins,Wage4$healthy)
RR_val <- (result_table[1]/(result_table[1]+result_table[3]))/(result_table[2]/(result_table[2]+result_table[4]))
RR_val
#because getting the confidence intervals are ANNOYING!
RR2=riskratio(result_table[1], result_table[2], result_table[1]+result_table[3],
              result_table[2]+result_table[4], p.calc.by.independence=FALSE)
RR=RR2$estimate
RR
RRlow=RR2$conf.int[1]
RRhigh=RR2$conf.int[2]

eval = EValue::evalues.RR(est=RR, RRlow, RRhigh)

bias_plot(RR,10)

evalue_nocontrol <- eval[2]


RR_control <- intframe$treatpred/intframe$notreatpred
#in their paper they turn RR<1 into 1/RR, we keep those as 1 with
#our monotonicity constraint
#RR_control <- ifelse(RR_control <1, 1/RR_control, RR_control)
Eval_control <- RR_control+sqrt(RR_control*(RR_control-1))
is.na(Eval_control)<- 1
data.frame(outcome=Eval_control)%>%
  ggplot(aes(x=outcome))+
  geom_histogram(aes(y=..count../sum(..count..))
                 ,color='white'
                 ,fill='#1d2951', bins=25)+
  geom_vline(xintercept=evalue_nocontrol)+
  geom_label(aes(label='No Control E-value', x=evalue_nocontrol+.35, y=.06))+
  geom_segment(aes(x = 2.0, y = 0.055, xend = evalue_nocontrol,
                   yend = 0.05),
     arrow = arrow(length=unit(0.20,"cm"),
                             ends="last", type = "closed"), lwd=.8)+
  ggtitle('Evalues')+
  ylab('Density')+xlab('Evalue Plots')+theme_minimal(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5,size=14))


data.frame(outcome=Eval_control, our = treat_frame[,'B1']/treat_frame[,'B0'])%>%
  ggplot(aes(x=outcome, y=our))+geom_point()
library(GJRM)
factors <- colnames(dataset)[1:(ncol(dataset)-1)]

formula_vars <- as.formula(paste("health_ins~", paste(factors, collapse="+")))
factors2 <- c('health_ins', factors)

formula_vars_2 <-  as.formula(paste("healthy~", paste(factors2, collapse="+")))

dataset$healthy <- Wage4$healthy
dataset$health_ins <- Wage4$health_ins
out<-gjrm(list(formula_vars,
               formula_vars_2),
          data=dataset,
          margins =c("probit", "probit"),
          Model = "B"
)
X.int <- as.matrix(out$X2)
ind.int<-(1:out$X2.d2)+out$X1.d2
coef.int <- as.numeric(out$coefficients[ind.int])
d0 <- d1 <- X.int
d0[, 'health_ins'] <- 0
d1[, 'health_ins'] <- 1
eti1 <- d1 %*% coef.int
eti0 <- d0 %*% coef.int

p.int1 <- probm(eti1, out$margins[2], min.dn = out$VC$min.dn,
                min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
p.int0 <- probm(eti0, out$margins[2],  min.dn = out$VC$min.dn,
                min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE)
n.sim=500
bs <- rMVN(n.sim, mean = out$coefficients, sigma = out$Vb)
eti1s <- d1 %*% t(bs[, ind.int])
eti0s <- d0 %*% t(bs[, ind.int])
#just get the coefficient for gamma for this!
gamma=qnorm(p.int1)-qnorm(p.int0)



rho_guess=summary(out)[]$theta
implied_U <- sqrt(rho_guess*(1-rho_guess))

data.frame(outcome=p.int1-p.int0)%>%
  ggplot(aes(x=outcome))+
  geom_histogram(aes(y=..count../sum(..count..))
                 ,color='white'
                 ,fill='#1d2951', bins=25)+
  ggtitle('All observations')+
  ylab('Density')+xlab('Individual Risk Differences')+theme_minimal(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5,size=14))



CV_pred=BARTpred_CV(Wage4, treat='health_ins', Outcome='healthy',
                    vars, mono=T,
                    nd_post=200, n_skip=500, nfold=5)
output=do.call('rbind', CV_pred$imp_frame)
output
cc=ggrocs(list(PB1G0=pROC::roc(output$healthy[output$health_ins==0],
                         output$BG0[output$health_ins==0]),
               PB1G1=pROC::roc(output$healthy[output$health_ins==1],
                         output$BG1[output$health_ins==1]),
               PG1=pROC::roc(output$health_ins[],
                       output$G1)),
          breaks = seq(0,1,0.1),
          tittle = "ROC perf. predicting our probs")


