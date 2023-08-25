library(nlme)
library(graphics)
rockfish <- read.csv(file="Data/rockfish.csv",header=TRUE)
rockfish_spp <- c("chilippeper","bocaccio","widow","canary","black","yellowtail","perch_us","perch_goose","perch_ak","perch_ai","perch_ebs")
rockfish$species <- factor(rockfish_spp[rockfish$SP],levels=rockfish_spp)
StockRec <- groupedData(REC~SSB|species,data=as.data.frame(rockfish),
                        labels=list(x="SSB",y="REC"),units=list(x="(Biomass)",y="(log(Recruits*SBPR/SSB)"))

m1 = lme(log(REC*SBPR/SSB) ~ SSB:species,data=rockfish,random = ~1|species,method="REML")

boxplot(split(residuals(m1),rockfish$species),ylab="Residual",xlab="Rockfish Species",csi=0.2)
abline(h=0,col="red",lty=2,lwd=2)
ln_a_est <- m1$coefficients$fixed[[1]]
exp(ln_a_est)
print(ln_a_est+m1$coefficients$random[[1]])
m1$coefficients$random$species[,1]
spp=c(5,4,8,6,7,2,11,1,9,3,10)
plot(spp,ln_a_est+m1$coefficients$random[[1]],pch=16,cex=1.5,xlab="Rockfish Species",ylab="Most Likely Random Effect Conditional upon Fixed Effect")
abline(h=ln_a_est,lty=2,col="red",lwd=2)

norm.lik = function(ln_a_i,sig_e,sig_a,ln_a,group,beta)
{
  sub_dat=subset(data,data$group_id==group)
  nobs <- nrow(sub_dat) #number of observations per group
  lik.a=lik.fix=rep(NA,nobs)
  lik.random=1/sqrt(2*pi)*exp(-(ln_a_i^2)/2) # what does this statement imply?
  for(k in 1:nobs)
  {
    mu <- ln_a-beta*sub_dat$x[k]-ln_a_i*sig_a
    lik.fix[k]=1/sqrt(2*pi*sig_e^2)*exp(-((sub_dat$y[k]-mu)^2)/(2*sig_e^2))
  }
  lik.fix.a=prod(lik.fix)
  lik.a=lik.random*lik.fix.a
  return(lik.a)
}

simpson = function(lowerb, upperb, nbins,group,ln_a,beta, sig_e, sig_a) 
{
  step = (upperb-lowerb)/nbins
  ln_a_vec = seq(lowerb, upperb, by=step)
  frac_vec = c(1/3, rep(c(4/3, 2/3), times=((nbins+1)-3)/2), 4/3, 1/3) # divide the possible range of values into sections
  sum_comp = vector(length=length(ln_a_vec))
  for(i in 1:length(ln_a_vec)){
    sum_comp[i] = frac_vec[i]*norm.lik(ln_a_i=ln_a_vec[i], sig_e=sig_e, sig_a=sig_a,beta=beta, ln_a=ln_a,group=group)
  }
  full_like = sum(sum_comp)*step
  return(full_like)
}

loglike=function(l_sig_e,l_sig_a,ln_a,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,beta_9,beta_10,beta_11)
{
  sig_e=exp(l_sig_e)
  sig_a=exp(l_sig_a)
  beta <- c(beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,beta_9,beta_10,beta_11)
  lowerb=-5
  upperb=5
  nbins=100
  group=as.numeric(unique(data$group))
  like.group=vector(length=length(group))
  for (i in group)
  {
    like.group[i]=simpson(lowerb=lowerb, upperb=upperb, nbins=nbins, 
                        sig_e=sig_e, sig_a=sig_a, ln_a=ln_a,beta=beta[i], group=group[i])
  }
  nll=(-1)*sum(log(like.group))
  return(nll)
}
rockfish$y <- log(rockfish$REC*rockfish$SBPR/rockfish$SSB)
rockfish$group_id <- as.numeric(rockfish$species)
rockfish$x <- rockfish$SSB
data <- rockfish
beta <- as.vector(m1$coefficients$fixed[-1])
betas <- NULL
for(i in 1:length(rockfish_spp))
{
  assign(paste("beta",i,sep="_"),beta[i])
  #beta_obj = get(paste("beta",i,sep="_"))
}
betas <- list(beta_1=beta_1,beta_2=beta_2,beta_3=beta_3,beta_4=beta_4,
              beta_5=beta_5,beta_6=beta_6,beta_7=beta_7,beta_8=beta_8,
              beta_9=beta_9,beta_10=beta_10,beta_11=beta_11)
sapply(1:length(beta),function(x){assign(paste("beta", x, sep="_"),beta[x])})
fit1=mle(loglike,start=c(list("l_sig_e"=log(10),"l_sig_a"=log(10),"ln_a"=2),betas),nobs=nrow(data),method="BFGS")
summary(fit1)
l_sig_e=coef(fit1)[[1]]
l_sig_a=coef(fit1)[[2]]
sig_e=exp(coef(fit1)[[1]])
sig_b=exp(coef(fit1)[[2]])
ln_a=coef(fit1)[[3]]
betas = -1*coef(fit1)[4:length(coef(fit1))] # remember to get the correct direction here
plot(betas,beta,ylab="nlme",xlab="Simpson")
abline(b=1,a=0)

AIC(fit1)
(-1)*logLik(fit1)[1] # the -1 converts to the negative log likelihood


##
##
##
## BETA PROFILE
##
##
##
##

ln_a_profile=seq(from=ln_a-0.75*ln_a,to=ln_a+0.75*ln_a,length=1000)
theta=c(beta,sig_e,sig_b)
####################################
## Linear mixed model for streams ##
####################################
b.pro.fun = function (theta) #(theta)
{
  beta=b.pro
  sig_e=sig_e
  sig_b=sig_b
  lowerb=-5
  upperb=5
  nbins=100
  stream=c(1,2,3,4,5,6)
  like.str=vector(length=length(stream))
  for (i in stream)
  {
    like.str[i]=simpson(lowerb=lowerb, upperb=upperb, nbins=nbins, 
                        sig_e=sig_e, sig_b=sig_b, beta=beta, stream=stream[i])
  }
  nll=(-1)*sum(log(like.str))
  return(nll)
}

beta.loglik=rep(NA,length(b.profile))
for(i in 1:length(b.profile))
{
  b.pro=b.profile[i]
  fit.pro=optim(theta,b.pro.fun,method="BFGS",control=list(maxit=1000))
  beta.loglik[i]=fit.pro$value
}

v2=cbind(b.profile,beta.loglik)
plot(b.profile,beta.loglik,type="l",lty=1,lwd=2,col="black",ylab="Negative Log-Likelihood",xlab="Beta")
abline(v=v2[,1][which(v2[,2]==min(v2[,2]))],col="red",lty=1)
x=(-1)*beta.loglik>=logLik(fit1)[1]-1.92
print(x)
abline(v=v2[401,1],col="red",lty=2) #printing x and finding the minimum TRUE row number
abline(v=v2[600,1],col="red",lty=2) #printing x and finding the maximum TRUE row number
abline(v=65.49489,col="darkgreen",lty=3) #taken from the LME fits below
abline(v=91.76256,col="darkgreen",lty=3) #taken from the LME fits below
legend("topright",c("ML Beta Value","Profile 95% CI","LME 95% CI"),col=c("red","red","darkgreen"),lty=c(1,2,3))

lm3 <- lme(fixed = Density ~ 1,data=Streams,random = ~ 1 | Stream,method="ML") 
boxplot(split(lm3$residuals,Streams$Stream),ylab="Residual",xlab="Stream",csi=0.2)
abline(0,0,lwd=3)
print(summary(lm3))
coef(lm3)
lm3[[4]]$fixed[[1]]+2*6.202739
lm3[[4]]$fixed[[1]]-2*6.202739
intervals(lm3)