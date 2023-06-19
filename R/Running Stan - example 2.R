library("rstan") # observe startup messages
parallel::detectCores()

options(mc.cores = 4)

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



y <- as.matrix(read.table('https://raw.github.com/wiki/stan-dev/rstan/rats.txt', header = TRUE))
x <- c(8, 15, 22, 29, 36)
xbar <- mean(x)
N <- nrow(y)
T <- ncol(y)
rats_fit <- stan(file="Week 1/rats.stan", data = list(N=N, T=T, y=y, x=x, xbar=xbar))

print(rats_fit)
plot(rats_fit)
pairs(rats_fit, pars = c("mu_alpha", "mu_beta", "sigma_y"))

la <- extract(rats_fit, permuted = TRUE) # return a list of arrays 
mu <- la$alpha 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(rats_fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(rats_fit)
m <- as.matrix(rats_fit)
d <- as.data.frame(rats_fit)


