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


intframe=BARTpred(Wage4, treat='health_ins', Outcome='healthy',
                  vars, mono=T)
intframe=prediction_function(Wage4, treat='health_ins', Outcome='healthy',vars)
#sd_chosen = sqrt(rho/(1-rho))
f=function(u){
  dnorm(u, mean=0, sd = .1) #for the time being, this has to be done as a number not a variable
}
treat_frame=integrate_function(intframe, constraint=F, f=f, n_cores=detectCores()-1, lambda=0)

library(rpart)

covariates <- Wage4[,vars]

#merge the irds and covariates
dataset=data.frame(covariates, outcome=treat_frame[, 'B1']-treat_frame[, 'B0'])

fit_tree(df=dataset, minbucket=500)

####Optional, if you want to make a quick histogram of IRD's
library(ggplot2)
dataset%>%
  ggplot(aes(x=outcome))+
  geom_histogram(aes(y=..count../sum(..count..))
                 ,color='white'
                 ,fill='#1d2951', bins=25)+
  ggtitle('All observations')+
  ylab('Density')+xlab('Individual Risk Differences')+theme_minimal(base_size = 16)+
  theme(plot.title = element_text(hjust = 0.5,size=16))


