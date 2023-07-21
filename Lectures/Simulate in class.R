# hello
# let's simulate a GLMM with a Poisson distribution

ngroups <- 12 # could be sites
sigma <- 1 # the standard deviation for the random effects among groups
thetas <- rnorm(ngroups,mean=0,sd=sigma)
names(thetas) <- 1:ngroups

# specify some regression coefficients
betas <- c("int"=0.5,"b1"=0.3,"b2"=-0.1)
### everything above here is truth

### everything below here is simulation
do_sim <- function()
{
  nsamp <- pmax(3,rpois(ngroups,12))
  N <- sum(nsamp)
  group_id <- rep(1:ngroups,times=nsamp)
  
  x1 <- rnorm(N,mean=0,sd=1)
  x2 <- rbinom(N,1,0.5)
  
  sim_data <- data.frame("y"=0,"x1"=x1,"x2"=x2,"group_id"=factor(group_id))
  sim_data$thetas <- thetas[sim_data$group_id]
  x_var <- model.matrix(y~x1+x2,data=sim_data)
  # log(lambda) = betas * X + thetas
  # lambda <- exp(betas * X + thetas)
  sim_data$lambda <- exp(x_var %*% betas + sim_data$thetas)
  sim_data$y <- rpois(N,lambda=sim_data$lambda)
  
  ## everything above here simulates
  
  ## everything here and below is estimating
  library(lme4)
  try(warnings(GLMER)!="TRUE")
  m1 <- lme4::glmer(y~x1+x2+(1|group_id),family="poisson",data=sim_data)
  summary(m1)
  coef(m1)
  lme4::fixef(m1)
  lme4::ranef(m1)
  plot(unlist(lme4::ranef(m1)$group_id),thetas,xlab="Random effects (lme4)",ylab="True effects")
  abline(b=1,a=0)
  
  ci <- confint(m1)
  true_pars <- c(sigma,betas)
  coverage <- as.numeric(data.table::between(true_pars,ci[,1],ci[,2]))
  bias <- 100*((true_pars-c(attr(summary(m1)$varcor$group_id,"stddev"),fixef(m1)))/true_pars)
  return(list("coverage"=coverage,"bias"=bias))
}

boot <- replicate(n=100,do_sim())
coverage <- do.call("rbind",boot[which(dimnames(boot)[[1]]=="coverage"),])
bias <- do.call("rbind",boot[which(dimnames(boot)[[1]]=="bias"),])
colnames(coverage) <- colnames(bias) <- c("sigma",names(betas))

colMeans(bias)
colMeans(coverage)





