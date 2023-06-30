library(dplyr)
library(ggplot2)
library(mvtnorm)
library(stats4)
library(MASS)
library(coda)
library(ggplot2)
library(ggmcmc)

lice <- read.csv("Data/sealice_fish_data.csv")
sites <- read.csv("Data/sealice_site_data.csv")

lice <- lice %>%
  group_by(site_id,year,day,month,location,species) %>%
  mutate(n=n()) %>%
  ungroup()

lice_data <- merge(lice,sites,by=c("site_id","year","day","month","location"),all.x=TRUE)
lice_data <- lice_data %>%
  filter(species!="sockeye")
lice_data$p.total <- apply(lice_data[,c(11:25)], 1, sum, na.rm=TRUE)
lice_data$infested <- as.numeric(lice_data$p.total>0)
lice_data <- lice_data[lice_data$year>2010,]
#Keep the years for plotting
years <- unique(lice_data$year)
lice_data$year <- as.factor(lice_data$year)
lice_data <- lice_data[!is.na(lice_data$length),]
# length (from the tip of the snout to the fork in the tail) and height (the maximum body depth).
# height can be assumed to be a proxy for weight here
# Fit linear model with year as categorical variable

#Fit model with categorical interaction between length and infested
fit1<-lm(height ~ length * infested, data=lice_data)
summary(fit1)

#Compare te above interaction model to a model w/o interaction
fit1b <- lm(height ~ length, data=lice_data)
summary(fit1b)
#Use LRT to comare the models
#Report the loglikelihoods, the test statistic and the p-value
#This is the logLik, not the negative loglik
loglik1 <- logLik(fit1)[1]
loglik1b <- logLik(fit1b)[1]

#Simple model - complex model 
teststat <- (-2*loglik1b) - (-2*loglik1) # twice the negative log likelihood is the deviance
#Get pval from chi-square test 
#Pval is significant 
#The more complex model can be accepted
#The df here is 2 
pval <- pchisq(teststat, df = 2, lower.tail = FALSE)
pval

# implement the MCMC sampling algorithm to explore how many parasites, on average, are on chum and pink salmon given their body sizes
# P ~ Poisson(species * length)
m1 <- glm(p.total~length*species,data=lice_data,family = poisson(link = "log"))
coef(m1)

plot(lice_data$p.total[!is.na(lice_data$length)],predict(m1,type="r"))
Xvar <- model.matrix(p.total~length*species,data=lice_data)
coefs <- 1:ncol(Xvar)
names(coefs) <- names_coefs <- colnames(Xvar)
n_pars <- length(coefs)

mcmc_covar <- as.matrix(vcov(m1))
beta_init <- c(-0.4088346,0.0090470,0.2285305,-0.0082515)

data <- list(N=length(lice_data$p.total),p.total=lice_data$p.total,Xvar = Xvar,npars=n_pars_alt)

GetPosterior <- function(beta,data)
{
  Xvar <- data$Xvar
  prior.lambda <- -dnorm(beta[1],mean=-1e-4,sd=5,log=TRUE)
  prior.beta <- -dnorm(beta[2:ncol(Xvar)],mean=0,sd=5,log=TRUE)
  pred <- exp(Xvar %*% beta)
  
  nll <- -dpois(data$p.total,pred,log=TRUE)
  posterior <- sum(sum(nll),sum(prior.beta),sum(prior.lambda))
  ifelse(is.na(posterior),posterior<-1e6,posterior<-posterior)
  return(list(posterior=posterior,dev=2*sum(nll),nll=nll))
}

DoMCMC<-function(beta_init,data,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1,jump,scalar)
{
  pb <- txtProgressBar(min=1,max=Nsim,style=3)
  beta_curr <- beta_init
  Get1st <- GetPosterior(beta_curr,data)
  Fcurr <- -1*Get1st$posterior
  Deviance <- Get1st$dev
  Like <- -Get1st$nll
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+2)) # store the parameter coefficients and the posterior + deviance
  Likes <- matrix(0,nrow=(Nsim-Nburn),ncol=data$N) # store the likelihood for model selection
  Icnt <- 0
  covar.new <- covar
  accept <- 0
  rate <- c()
  for (Isim in 1:Nsim)
  {
    rate[Isim] <- accept/Isim
    if(Isim<=Nburn & Isim %% 50 == 0)
    {
      if(rate[Isim]>=0.30 & rate[Isim]<=0.40)
      {
        scalar <- scalar
      }else{
        if(rate[Isim]<0.30)
        {
          scalar <- max(0.01,min(1,scalar-0.01))
        }
        if(rate[Isim]>0.40)
          scalar <- min(1,max(0.01,scalar+0.01))
      }
    }
    if(Isim==Nburn){scalar.save <- scalar}
    covar.new <- covar*scalar
    beta_next <- as.numeric(rmvnorm(1, mean=beta_curr, sigma=as.matrix(covar.new)))
    Get2nd <- GetPosterior(beta_next,data)
    Fnext <- -1*Get2nd$posterior
    Dev2 <- Get2nd$dev
    Like2 <- -Get2nd$nll
    Rand1 <- log(runif(1,0,jump))
    if (Fnext > Fcurr+Rand1)
    {
      Fcurr <- Fnext
      beta_curr <- beta_next
      Deviance <- Dev2
      Like <- Like2
      accept <- accept + 1
    }   
    if (Isim > Nburn & Isim %% Nthin == 0)
    {
      Icnt <- Icnt + 1
      Outs[Icnt,] <- c(beta_curr,Fcurr,Deviance)
      Likes[Icnt,] <- Like
      cat("saving",Icnt,"\n")     
    }
    setTxtProgressBar(pb,Isim)
  } 
  xx <- seq(1:Icnt)
  par(mfrow=c(3,3),mar=c(4,4,1,1))
  for (II in 1:(Ndim+2))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- names(beta_init)[II]
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- (yy)
    plot(xx,yy,xlab="Cycle number",ylab=lab1,type='b',pch=16,cex=0.02)
  }
  par(mfrow=c(3,3),mar=c(4,4,1,1))
  for (II in 1:(Ndim+2))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- names(beta_init)[II] #paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- (yy)
    hist(yy,ylab=lab1,main="",col="steelblue")
  }
  return(list(results=Outs[1:Icnt,],accept_rate=rate,scalar.burn=scalar.save,like=Likes))
}


Nsim <- 1e4
burn_rate <- 0.5
Nthin <- 1
####
## 1st chain
####
jmp <- 1 
scale <- (2.4/sqrt(n_pars))^2#0.5 
a <- proc.time();
Outs <- DoMCMC(beta_init=as.numeric(beta_init),data=data,Ndim=n_pars,covar=mcmc_covar,Nsim=Nsim,Nburn=burn_rate*Nsim,Nthin=Nthin,jump=jmp,scalar=scale)
(proc.time() - a)/Nsim
par(mfrow=c(1,1))
plot(1:Nsim,Outs$accept_rate,type="l",lwd=2,lty=2,xlab="MCMC Chain",ylab="Acceptance Rate")

Output <- matrix(Outs$results,ncol=n_pars+2,byrow=F)  
MCMC1 <- mcmc(Output) #drop the posterior density, just keep the deviance
colnames(MCMC1) <- c(names(coefs),"Posterior","Deviance")

####
## 2nd chain
####

jmp <- 1
scale <- 0.5
a <- proc.time();
Outs2 <- DoMCMC(beta_init=as.numeric(beta_init),data=data,Ndim=n_pars,covar=mcmc_covar,Nsim=Nsim,Nburn=burn_rate*Nsim,Nthin=Nthin,jump=jmp,scalar=scale)
(proc.time() - a)/Nsim
par(mfrow=c(1,1))
plot(1:Nsim,Outs2$accept_rate,type="l",lwd=2,lty=2,xlab="MCMC Iteration",ylab="Acceptance Rate")
Output2 <- matrix(Outs2$results,ncol=n_pars+2,byrow=F)  
MCMC2 <- mcmc(Output2) #drop the posterior density, just keep the deviance
colnames(MCMC2) <- c(names(coefs),"Posterior","Deviance")

# Task 2: model selection

Xvar_alt <- model.matrix(p.total~length,data=lice_data)
coefs_alt <- 1:ncol(Xvar_alt)
names(coefs_alt) <- names_coefs_alt <- colnames(Xvar_alt)
n_pars_alt <- length(coefs_alt)

mcmc_covar_new <- mcmc_covar[1:2,1:2]
beta_init_alt <- c(-0.4088346,0.0090470)

data_alt <- list(N=length(lice_data$p.total),p.total=lice_data$p.total,Xvar = Xvar_alt,npars=n_pars_alt)

Outs_alt1 <- DoMCMC(beta_init=as.numeric(beta_init_alt),data=data_alt,Ndim=n_pars_alt,covar=mcmc_covar_new,Nsim=Nsim,Nburn=burn_rate*Nsim,Nthin=Nthin,jump=jmp,scalar=scale)
Outs_alt2 <- DoMCMC(beta_init=as.numeric(beta_init_alt),data=data_alt,Ndim=n_pars_alt,covar=mcmc_covar_new,Nsim=Nsim,Nburn=burn_rate*Nsim,Nthin=Nthin,jump=jmp,scalar=scale)

like_m1 <- rbind(Outs$like,Outs2$like)
like_m2 <- rbind(Outs_alt1$like,Outs_alt2$like)
loo::loo_compare(loo::loo(like_m1),loo::loo(like_m2))

# Task 3: histogram and convergence criteria for MCMC model
# merge the chains into an MCMC list object
MCMC_chains <- mcmc.list(MCMC1,MCMC2)
ggmcmc(ggs(MCMC_chains),file="Labs/Lab 3/mcmc_output.pdf")

print(gelman.diag(MCMC_chains))
gelman.plot(MCMC_chains)
ggs_Rhat(ggs(MCMC_chains))
print(heidel.diag(MCMC_chains))
# Graphical diagnostics
traceplot(MCMC_chains)

posterior <- as.matrix(MCMC_chains)

# plot histogram of results
par(mfrow=c(2,2),mar=c(5,4,1,1),omi=c(0.05,0.05,0.05,0.05))
for (i in 1:n_pars)
{
  hist(posterior[,i],main="",xlab=names_coefs[i],col="lightblue")
}

# Task 4: posterior predictive distribution for a pink salmon and chum salmon of length 60
names_coefs
pink_ppd <- rpois(1:nrow(posterior),exp(posterior[,1:n_pars]%*%c(1,60,1,60)))
hist(pink_ppd)
mean(pink_ppd);quantile(pink_ppd,probs=c(0.025,0.975))

chum_ppd <- rpois(1:nrow(posterior),exp(posterior[,1:n_pars]%*%c(1,60,0,0)))
hist(chum_ppd)
mean(chum_ppd);quantile(chum_ppd,probs=c(0.025,0.975))


# Task 5: posterior predictive distribution of chum salmon along their body sizes
body_size <- 25:150
chum_ppd <- sapply(body_size,function(x){rpois(1:nrow(posterior),exp(posterior[,1:n_pars]%*%c(1,x,0,0)))})
chum_mn <- colMeans(chum_ppd)
chum_ci <- apply(chum_ppd,2,quantile,probs=c(0.025,0.975))
pink_ppd <- sapply(body_size,function(x){rpois(1:nrow(posterior),exp(posterior[,1:n_pars]%*%c(1,x,1,x)))})
pink_mn <- colMeans(pink_ppd)
pink_ci <- apply(pink_ppd,2,quantile,probs=c(0.025,0.975))
layout(matrix(1:2,nrow=2))
plot(body_size,chum_mn,ylim=range(chum_ci,pink_ci),col=0,xlab="Body size (mm)",ylab="Parasites")
polygon(c(body_size,rev(body_size)),c(chum_ci[1,],rev(chum_ci[2,])),border=NA,col=adjustcolor("dodgerblue",0.5))
lines(body_size,chum_mn,lwd=2)

plot(body_size,pink_mn,ylim=range(chum_ci,pink_ci),col=0,xlab="Body size (mm)",ylab="Parasites")
polygon(c(body_size,rev(body_size)),c(pink_ci[1,],rev(pink_ci[2,])),border=NA,col=adjustcolor("dodgerblue",0.5))
lines(body_size,pink_mn,lwd=2)
points(lice_data$length[lice_data$species=="pink"],lice_data$p.total[lice_data$species=="pink"])
