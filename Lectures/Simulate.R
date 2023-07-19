# simulation
# let's simulate a glmm
ngroups <- 12
sigma <- 1
betas <- c("int"=0.5,"b1"=0.3,"b2"=-0.1)
thetas <- rnorm(ngroups,mean=0,sd=sigma)
names(thetas) <- 1:ngroups
do_sim <- function()
{
  nsamp <- rpois(ngroups,12)
  N <- sum(nsamp)
  group_id <- rep(1:ngroups,times=nsamp)
  x1 <- rnorm(N,0,1)
  x2 <- rbinom(N,1,0.5)
  
  sim_data <- data.frame("y"=0,"x1"=x1,"x2"=x2,"group_id"=factor(group_id))
  sim_data$thetas <- thetas[sim_data$group_id]
  x_var <- model.matrix(y~x1+x2,data=sim_data)
  sim_data$lambda <- exp(x_var %*% betas + thetas[sim_data$group_id])
  sim_data$y <- rpois(N,lambda=sim_data$lambda)
  
  m1 <- lme4::glmer(y~x1+x2+(1|group_id),family="poisson",data=sim_data)
  summary(m1)
  coef(m1)
  fixef(m1)
  ranef(m1)
  plot(unlist(ranef(m1)$group_id),as.numeric(thetas),xlab="Random effects (lme4)",ylab="True effects")
  abline(b=1,a=0)
  
  ci <- confint(m1)
  pars <- c(sigma,betas)
  coverage <- as.numeric(data.table::between(pars,ci[,1],ci[,2]))
  bias <- 100*((pars-c(attr(summary(m1)$varcor$group_id, "stddev"),fixef(m1)))/pars)
  return(list("coverage"=coverage,"bias"=bias))
}
boot <- replicate(100,do_sim())
coverage <- do.call("rbind", boot[which(dimnames(boot)[[1]]=="coverage"),])
bias <- do.call("rbind", boot[which(dimnames(boot)[[1]]=="bias"),])
colnames(coverage) <- colnames(bias) <- c("sigma",names(betas))
colMeans(bias)
colMeans(coverage)
