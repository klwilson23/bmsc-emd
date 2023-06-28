# read in data for west coast vancouver island fisheries dataset from Fisheries & Oceans Canada
wcvi <- read.csv("Data/WCVI_biology.csv")
wcvi$length <- ifelse(!is.na(wcvi$Fork.length..mm.),wcvi$Fork.length..mm.,
                      ifelse(!is.na(wcvi$Total.length..mm.),wcvi$Total.length..mm.,wcvi$Standard.length..mm.))
wcvi$length_method <- ifelse(!is.na(wcvi$Fork.length..mm.),"fork_length",
                      ifelse(!is.na(wcvi$Total.length..mm.),"total_length","standard_length"))

wcvi <- data.frame("year"=wcvi$Survey.Year,"trip_id"=wcvi$Trip.identifier,"tow"=wcvi$Tow.number,"sample"=wcvi$Sample.identifier,"species"=gsub(" ","_",tolower(wcvi$English.common.name)),"length"=wcvi$length,"weight"=wcvi$Weight..g.,"age"=wcvi$Age)
wcvi <- wcvi[complete.cases(wcvi$weight,wcvi$length),]
wcvi$length <- wcvi$length/10 # convert from mm to cm
plot(weight~length,wcvi,xlab="Length (cm)",ylab="Weight (g)",xlim=c(0,200),ylim=c(0,80000))

# lets focus on a particular species, yellowtail rockfish, which is an important species for fisheries
df <- wcvi[wcvi$species=="yellowtail_rockfish",]
sigma <- round(sd(log(df$weight))/mean(log(df$weight)),2) # coefficient of variation in length-weight relationship
set.seed(10072023) # set the seed to the start-date of the course: 10th of July 2023
df <- df[sample(1:nrow(df),50,replace=FALSE),] # let's reduce the sample size down to something manageable
plot(log(weight)~log(length),df,ylim=c(-4.5,8),xlim=c(0,6.4))
abline(a=log(10^(-1.25)),b=3.04)
log10_a_prior <- c(-1.95,0.173) # prior on log10 a, from Froese et al. 2012: mean = -1.95, sd = 0.173
b_prior <- c(3.04,0.0857) # prior on b, from Froese et al. 2012: mean = 3.04, sd = 0.0857

# https://fishbase.us/popdyn/BayesianAnalysis.php
# build a function in R that estimates the relationship between length and weight of yellowtail under the following assumptions: 
# (1) assume that sigma is 0.07
# (2) the error around the mean of the log(weight) is normally distributed
# (3) the prior on alpha is based on the log10, not the natural log - be warned!
# (4) the prior

# generate a function that returns the expected body size of a fish for a given weight
# preference for the log-linear form so that it looks like y = b0 + b1*x
lwr <- function(log_alpha,beta,length,type="link")
{
  if(type=="link")
  {
    pred <- log_alpha+beta*log(length)
  }
  if(type=="response")
  {
    pred <- exp(log_alpha)*length^(beta)
  }
  return(pred)
}
# generate a function that calls the length-weight fx, and returns the joint log-likelihood of the data given the model
nll <- function(theta,sigma,sigma_est="fixed",data)
{
  # we pass some input parameters and objects to this function: theta is a the 'set' or 'vector' of parameters we are interested in
  log_alpha <- theta[1]
  beta <- theta[2]
  if(sigma_est=="estimated")
  {
    sigma <- theta[3]
  }else{
    sigma <- sigma
  }
  log_pred <- lwr(log_alpha=log_alpha,beta=beta,length=data$length)
  ss <- (log_pred-log(data$weight))^2 # this is the sums of squared errors - what does it seem similar to?
  like <- (1/(sigma*sqrt(2*pi)))*exp(-0.5*((log_pred-log(data$weight))/sigma)^2) # the normal probability density fx
  log_like <- log(like)
  ll <- sum(log_like) # we sum the log-likelihood here
  return(ll)
}

# generate some starting parameters
theta <- c(-4,3,sigma)
# call the function 'optim'
# optim is a fn that minimizes a target function 'nll', but we can flip it to maximize using the control argument
fit <- optim(par=theta,fn=nll,method="BFGS",hessian=TRUE,data=df,sigma=sigma,sigma_est="estimated",control=list(fnscale=-1))

# generate 95% confidence intervals from optim by using the Hessian matrix
fisher_info <- solve(-fit$hessian) # solve the Hessian matrix
prop_sigma <- sqrt(diag(fisher_info)) # take the square roots of the main diagonals of the variance-covariance matrix (squared root of variance is the std. deviation)
upper <- fit$par+1.96*prop_sigma # asymptotic confidence intervals is mean +/- 1.96*sigma
lower <- fit$par-1.96*prop_sigma
# store the maximum likelihood estimates and 95% CI in a data frame
approx_interval <- data.frame(metric=c("ln(a)","b","sigma"),method="optim (BFGS)",value=fit$par,sd=prop_sigma, upper=upper, lower=lower)


# make a function that returns the posterior distribution and marginal distributions 
# this fn should input a fixed sigma, a parameter that controls the size of the grid, the data, and the Froese et al. 2012 priors for log10(a) and b
# this fn should set up a grid, and loop through each element of the grid to calculate the prior, likelihood, and the posterior for each paramater combination
# remember that the posterior = prior * likelihood
# use the log-likelihoods to sum to avoid overload of many small numbers
# transform to likelihood
# normalize the posterior to be 0-1 by dividing by the sum of the posterior
# search the ln(a) from -5.5 to -3.5
# search the b from 2.9 to 3.2
# calculate the marginal distributions for each parameter by summing the posterior across all other parameters
make_grid <- function(grid.size,sigma,data,pr_alpha,pr_beta)
{
  log_alpha <- seq(-5.5,-3.5,length=grid.size)
  beta <- seq(2.9,3.2,length=grid.size)
  prior <- like <- posterior <- matrix(NA,nrow=grid.size,ncol=grid.size,dimnames=list(log_alpha,beta))
  for(i in 1:grid.size)
  {
    for(j in 1:grid.size)
    {
      pr.alpha <- dnorm(log10(exp(log_alpha[i])),mean=pr_alpha[1],sd=pr_alpha[2],log=FALSE)
      pr.beta <- dnorm(beta[j],mean=pr_beta[1],sd=pr_beta[2],log=FALSE)
      prior[i,j] <- (pr.alpha*pr.beta) # what is the prior probability
      ln_pred <- lwr(log_alpha=log_alpha[i],beta=beta[j],length=data$length,type="link") # get the predictions
      ll <- log(1/(sigma*sqrt(2*pi))*exp(-(log(data$weight)-ln_pred)^2/(2*sigma^2))) # what is the log likelihood
      like[i,j] <- exp(sum(ll)) # convert into normal space
      posterior[i,j] <- prior[i,j]*like[i,j] # calculate the posterior
    }
  }
  norm_post <- posterior/sum(posterior) # normalize
  marg.alpha <- rowSums(norm_post) # derive the marginal likelihood
  marg.beta <- colSums(norm_post) # derive the marginal likelihood
  return(list("norm_post"=norm_post,"log_alpha"=log_alpha,"beta"=beta,"marg.alpha"=marg.alpha,"marg.beta"=marg.beta))
}

results <- make_grid(grid.size=200,sigma=sigma,data=df,pr_alpha=c(-1.95,0.173),pr_beta=c(3.04,0.0857)) # Change grid.size here to assess sensitivity

persp(x=results$log_alpha,y=results$beta,
      z=results$norm_post,
      xlab="\u03B1",ylab="\u03B2",zlab="Posterior Density",theta=40,phi=20)

results_2 <- make_grid(grid.size=100,sigma=sigma,data=df,pr_alpha=c(-1.95,0.05),pr_beta=c(3.04,0.01)) # Change grid.size here to assess sensitivity

pr_alpha <- dnorm(log10(exp(seq(-8,-1,length=100))),mean=-1.95,sd=0.173,log=FALSE)
plot(names(results$marg.alpha),results$marg.alpha,ylab="Posterior Density",xlab=expression(alpha),type="l",lwd=2,xlim=c(-7,-2),ylim=c(0,0.1))
lines(names(results_2$marg.alpha),results_2$marg.alpha,lwd=2,col="red")
lines(seq(-8,-1,length=100),pr_alpha/sum(pr_alpha),lwd=2,col="blue")

plot(names(results$marg.beta),results$marg.beta,ylab="Posterior Density",xlab=expression(beta),type="l",lwd=2,xlim=c(2.9,3.2),ylim=c(0,0.1))
lines(names(results_2$marg.beta),results_2$marg.beta,lwd=2,col="red")
pr_beta <- dnorm(seq(2.9,3.2,length=100),mean=3.04,sd=0.0857,log=FALSE)
lines(seq(2.9,3.2,length=100),pr_beta/sum(pr_beta),lwd=2,col="blue")

# lets make uninformed priors on ln(a) and b, but fix sigma at MLE values
sigma <- approx_interval$value[approx_interval$metric=="sigma"]
results <- make_grid(grid.size=200,sigma=sigma,data=df,pr_alpha=c(-1.95,100),pr_beta=c(3.04,100)) # Change grid.size here to assess sensitivity

MAP <- which(results$norm_post==max(results$norm_post),arr.ind=TRUE) # what is the maximum a posteriori
ln_alpha_MAP <- results$log_alpha[MAP[1]] # the mode
beta_MAP <- results$beta[MAP[2]] # the mode

# norm_post <- as.data.frame(as.table(results$norm_post))
# colnames(norm_post) <- c("log_alpha","beta","posterior")
# norm_post$log_alpha <- as.numeric(as.character(norm_post$log_alpha))
# norm_post$beta <- as.numeric(as.character(norm_post$beta))

# calculate the MAP, MEAN, 95% CI, and SD for ln(a) and b
ln_a_mn <- sum(results$log_alpha*results$marg.alpha)
#posterior_samps_lna <- sample(x=norm_post$log_alpha,size=1e4,prob=norm_post$posterior,replace=TRUE)
posterior_samps_lna <- sample(x=results$log_alpha,size=1e4,prob=results$marg.alpha,replace=TRUE)
ln_a_ui <- quantile(posterior_samps_lna,probs=0.975)
ln_a_li <- quantile(posterior_samps_lna,probs=0.025)
ln_a_sd <- sqrt(sum(results$log_alpha^2*results$marg.alpha) - sum(results$log_alpha*results$marg.alpha)^2)

b_mn <- sum(results$beta*results$marg.beta) # calculate the mean
#posterior_samps_beta <- sample(x=norm_post$beta,size=1e4,prob=norm_post$posterior,replace=TRUE)
posterior_samps_beta <- sample(x=results$beta,size=1e4,prob=results$marg.beta,replace=TRUE)
b_ui <- quantile(posterior_samps_beta,probs=0.975)
b_li <- quantile(posterior_samps_beta,probs=0.025)
b_sd <- sqrt(sum(results$beta^2*results$marg.beta) - sum(results$beta*results$marg.beta)^2)

grid_interval <- data.frame(metric=c("ln(a)","b","sigma"),method="Bayesian grid",value=c(ln_a_mn,b_mn,sigma),sd=c(ln_a_sd,b_sd,NA), upper=c(ln_a_ui,b_ui,NA), lower=c(ln_a_li,b_li,NA))
intervals <- merge(grid_interval,approx_interval,all=TRUE)

# based on the mean and sigma from either of the above approaches:  
# construct a 95% predictive interval (either frequentist or Bayesian) for a Yellowtail Rockfish measured at 60 cm

weight_ppd <- exp(rnorm(length(posterior_samps_lna),mean=posterior_samps_lna+posterior_samps_beta*log(60),sd=sigma))
mean(weight_ppd);quantile(weight_ppd,probs=c(0.025,0.975))
# how likely is a 60 cm to exceed 3000g?
sum(weight_ppd>3000)/length(weight_ppd)

# compare the results from each approach above to the lm() function in R
summary(lm(log(weight)~log(length),data=df))
