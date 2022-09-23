#2022 
#statistics scripts by Statistics Services Ltd (Reading,UK)
#Mixed effects models for granule per chloroplast numbers
#Poisson and negative binomial models
#eventually we used negative binomial models

#NON-BACKCROSSED LINES
dat<-read.csv("2021_gpc_nonbackcrossed.csv")
dat$rep<-factor(dat$rep)
dat$group<-factor(dat$group)
dat$counts.f<-ordered(dat$counts)
colnames(dat)

# MAKE UNIQUE REP ID
dat$rep.id<-dat$rep:dat$group
length(levels(dat$rep.id))
head(dat)
# THE FITTING FUNCTION REQUIRES DISAGGREGATED DATA
dat.disag<-dat[rep(1:nrow(dat),dat[,2]),-2]
head(dat.disag)
aggregate(dat$freq~dat$rep.id,FUN=sum)
aggregate(dat.disag$counts~dat.disag$rep.id,FUN=length)

library(MASS)
nb0<-glm.nb(counts~1,data=dat.disag) #MUST BE FITTED BEFORE LOADING GLMMadaptive TO PREVENT MASKING

####################################################
# MIXED EFFECTS MODELS, FITTED BY MAXIMUM LIKELIHOOD
####################################################
install.packages("GLMMadaptive")
library(GLMMadaptive); library(emmeans)
####################
# POISSON REGRESSION
####################
poi0<-glm(counts~1,data=dat.disag,family=poisson)
poi1<-mixed_model(counts~1,~1|rep.id,data=dat.disag,family=poisson)
poi2<-mixed_model(counts~1+group,~1|rep.id,data=dat.disag,family=poisson)
anova(poi1,poi0) #TEST OF RANDOM EFFECT COMPONENT, RESULTING P-VALUE MUST BE HALVED
#there is an effect
anova(poi1,poi2) #TEST OF FIXED EFFECT ADJUSTED FOR RANDOM EFFECT
#there is an effect

# MUST REGRID, OR THE ARG type IS IGNORED
poi2r<-update(ref_grid(poi2),tran='log') 
emmeans(poi2r,~group,type='response') # PREDICTED MEAN COUNT
pairs(emmeans(poi2r,~group,type='response'),inf=T) #estimating performed on log scale, then antilog at the end
pairs(emmeans(poi2r,~group,type='response'),rev=T,inf=T)
pairs(emmeans(poi2r,~group),inf=T) #log scale

# ILLUSTRATE FITTED CONDITIONAL DISTRIBUTIONS
curve(dpois(x,5.21),0,21,22,lwd=2,ylim=c(0,.25),xlab='count',ylab='pmf')
curve(dpois(x,2.6),0,21,22,lwd=2,add=T,col=2)
curve(dpois(x,3.75),0,21,22,lwd=2,add=T,col=3)
curve(dpois(x,4.78),0,21,22,lwd=2,add=T,col=4)
legend('topright',levels(dat$group),lwd=2,col=1:4,title='group:')
title('Fitted conditional Poisson distribution')

# INFORMAL GOODNESS OF FIT CRITERION, e.g. FROM GBUR ET AL 2012, SECTION 5.4, PAGE 128
# SUM OF CONDITIONAL PEARSON RESIDUALS SQUARED / MODEL DF
# A RATIO >>1 SUGGESTS MORE OBSERVED VARIANCE THAN EXPECTED FROM A POISSON
dat.disag$pred.ss<-predict(poi2,dat.disag,type='sub')
dat.disag$pearson<-(dat.disag$counts-dat.disag$pred.ss)/sqrt(dat.disag$pred.ss)
pearchi2<-sum(dat.disag$pearson^2)
pearchi2/nobs(poi2) 
# MILD OVERDISPERSION, SO TRY NEGATIVE BINOMIAL

##############################
# NEGATIVE BINOMIAL REGRESSION
##############################

nb1<-mixed_model(counts~1,~1|rep.id,data=dat.disag,family=negative.binomial)
nb2<-mixed_model(counts~group,~1|rep.id,data=dat.disag,family=negative.binomial)
anova(nb1,nb0) #TEST OF RANDOM EFFECT COMPONENT, RESULTING P-VALUE MUST BE HALVED
#there is an effect

# IGNORE WARNING, IT'S CAUSED BY MASKING
anova(nb1,nb2) #TEST OF FIXED EFFECT ADJUSTED FOR RANDOM EFFECT
#there is an effect 

# MUST REGRID, OR THE ARG type IS IGNORED
nb2r<-update(ref_grid(nb2),tran='log')
emmeans(nb2r,~group,type='response') # PREDICTED MEAN COUNT
pairs(emmeans(nb2r,~group,type='response'),inf=T)
pairs(emmeans(nb2r,~group,type='response'),rev=T,inf=T)

# ILLUSTRATE FITTED CONDITIONAL DISTRIBUTIONS
theta<-exp(nb2$phis); theta
# NOTE ARG SIZE=THETA
curve(dnbinom(x,size=theta,mu=5.21),0,21,22,lwd=2,ylim=c(0,.25),xlab='count',ylab='pmf')
curve(dnbinom(x,size=theta,mu=2.61),0,21,22,lwd=2,add=T,col=2)
curve(dnbinom(x,size=theta,mu=3.75),0,21,22,lwd=2,add=T,col=3)
curve(dnbinom(x,size=theta,mu=4.77),0,21,22,lwd=2,add=T,col=4)
legend('topright',levels(dat$group),lwd=2,col=1:4,title='group:')
title('Fitted conditional Negative Binomial distribution')

# INFORMAL GOODNESS OF FIT CRITERION, e.g. FROM GBUR ET AL 2012, SECTION 5.4, PAGE 128
# SUM OF CONDITIONAL PEARSON RESIDUALS SQUARED / MODEL DF
# A RATIO >>1 SUGGESTS MORE OBSERVED VARIANCE THAN EXPECTED FROM A POISSON
dat.disag$pred.ss<-predict(nb2,dat.disag,type='sub')
dat.disag$pearson<-(dat.disag$counts-dat.disag$pred.ss)/sqrt(dat.disag$pred.ss+dat.disag$pred.ss^2/theta)
pearchi2<-sum(dat.disag$pearson^2)
pearchi2/nobs(nb2) 
# OVERDISPERSION NOW CURED

##################
# SAME ANALYSIS FOR BACKCROSSED LINES

dat2<-read.csv("2021_gpc_backcrossed.csv")
colnames(dat2)
dat2$rep<-factor(dat2$rep)
dat2$group<-factor(dat2$group)
dat2$counts.f<-ordered(dat2$counts)
dat2$counts<-as.integer(dat2$counts)
dat2$freq<-as.integer(dat2$freq)

sapply(dat2, class)
dat2$rep.id<-dat2$rep:dat2$group
length(levels(dat2$rep.id))
head(dat2)

# THE FITTING FUNCTION REQUIRES DISAGGREGATED DATA
dat2.disag<-dat2[rep(1:nrow(dat2),dat2[,2]),-2]
head(dat2.disag)
aggregate(dat2$freq~dat2$rep.id,FUN=sum)
aggregate(dat2.disag$counts~dat2.disag$rep.id,FUN=length)


##############################
# NEGATIVE BINOMIAL REGRESSION
##############################
library(MASS)
nb02<-glm.nb(counts~1,data=dat2.disag) #MUST BE FITTED BEFORE LOADING GLMMadaptive TO PREVENT MASKING

library(GLMMadaptive); library(emmeans)
nb12<-mixed_model(counts~1,~1|rep.id,data=dat2.disag,family=negative.binomial)
nb22<-mixed_model(counts~group,~1|rep.id,data=dat2.disag,family=negative.binomial)
anova(nb12,nb02) #TEST OF RANDOM EFFECT COMPONENT, RESULTING P-VALUE MUST BE HALVED
#there is an effect

# IGNORE WARNING, IT'S CAUSED BY MASKING
anova(nb12,nb22) #TEST OF FIXED EFFECT ADJUSTED FOR RANDOM EFFECT
#there is an effect

# MUST REGRID, OR THE ARG type IS IGNORED
nb2r2<-update(ref_grid(nb22),tran='log')
emmeans(nb2r2,~group,type='response') # PREDICTED MEAN COUNT
pairs(emmeans(nb2r2,~group,type='response'),inf=T)
pairs(emmeans(nb2r2,~group,type='response'),rev=T,inf=T)

# ILLUSTRATE FITTED CONDITIONAL DISTRIBUTIONS
theta2<-exp(nb22$phis); theta2

# INFORMAL GOODNESS OF FIT CRITERION, e.g. FROM GBUR ET AL 2012, SECTION 5.4, PAGE 128
# SUM OF CONDITIONAL PEARSON RESIDUALS SQUARED / MODEL DF
# A RATIO >>1 SUGGESTS MORE OBSERVED VARIANCE THAN EXPECTED FROM A MODEL
dat2.disag$pred.ss<-predict(nb22,dat2.disag,type='sub')
dat2.disag$pearson<-(dat2.disag$counts-dat2.disag$pred.ss)/sqrt(dat2.disag$pred.ss+dat2.disag$pred.ss^2/theta2)
pearchi22<-sum(dat2.disag$pearson^2)
pearchi22/nobs(nb22) 
