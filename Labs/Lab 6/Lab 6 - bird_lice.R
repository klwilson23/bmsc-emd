library(lme4)
library(tidyverse)
library(rstan)
library(bayesplot)
library(bayesAB)
library(loo)
library(DHARMa)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
birds <- read.csv("Data/allopreening_data.csv")
birds$scale_allo <- as.vector(scale(birds$perc_allo))
birds$scale_self <- as.vector(scale(birds$perc_self))
m1 <- glmer(louse~scale_allo+scale_self+(1|aviary/breed),data=birds,family="poisson")
summary(m1)
ranef(m1)
DHARMa::simulateResiduals(m1,plot=T)

Xvar <- model.matrix(louse~scale_allo+scale_self,data=birds)
burb_ave <- as.numeric(as.factor(birds$breed[match(unique(birds$aviary),birds$aviary)]))
x_predict <- expand.grid(1,scale_allo=seq(-2,3,by=0.1),scale_self=seq(-2,3,by=0.1))
stan_data <- list(N=nrow(Xvar),
                  Ncovar = ncol(Xvar),
                  Xvar= Xvar,
                  lice = birds$louse,
                  N_aviary = length(unique(birds$aviary)),
                  N_burbs = length(unique(birds$breed)),
                  aviary_id = as.numeric(as.factor(birds$aviary)),
                  burb_id = as.numeric(as.factor(birds$breed)),
                  burb_ave = burb_ave,
                  N_predict = nrow(x_predict),
                  x_predict = x_predict)

init_fx <- function(chain_id)
{
  list("beta" = rep(0,ncol(Xvar)),
       "phi" = 1,
       "sigma_aviary"=1,
       "sigma_burb"=1,
       "burb_eff" = rep(0,stan_data$N_burbs),
       "aviary_eff" = rep(0,stan_data$N_burbs))
}
#model_string <- cmdstanr::cmdstan_model(stan_file="Labs/Lab 6/bird_lice.stan")

fit <- rstan::stan(file = "Labs/Lab 6/bird_lice.stan",
                   data = stan_data,
                   iter = 2000,
                   chains = 4,
                   init = init_fx,
                   control=list("max_treedepth"=10,"adapt_delta"=0.95))
post_summ <- summary(fit,pars=c("betas","sigma_aviary","sigma_burb","phi","burb_eff","aviary_eff"),probs=c(0.025,0.975))$summary
summary(fit,pars=c("mu"))$summary

## let's check this in brms too?
library(brms)
brm_fit <- brm(louse~scale_allo+scale_self+(1|breed/aviary),data=birds,family=negbinomial())
summary(brm_fit)
brms::ranef(brm_fit)

## back to our Stan fit
# pp_checks
y_ppd <- rstan::extract(fit)$y_ppd
ci <- apply(y_ppd,2,function(x){c(mean(x),quantile(x,probs=c(0.025,0.975)))})

plot(stan_data$lice,ci[1,],ylim=range(ci))
segments(x0=stan_data$lice,y0=ci[2,],y1=ci[3,])
points(stan_data$lice,ci[1,])

## follow along with Hannah's lab 4
post_summ %>% 
  as.data.frame() %>% 
  mutate(rhat = round(Rhat,2)) %>% 
  filter(rhat != 1)

posterior <- as.array(fit)
bayesplot::mcmc_areas(posterior, pars = c("betas"))
rstan::traceplot(fit, pars = c('betas', 'sigma_aviary', 'sigma_burb','phi'))

# do a quick pp_check
b1 <- bayesplot::ppc_dens_overlay(stan_data$lice, y_ppd[sample(1:nrow(y_ppd),250),])
b2 <- bayesplot::ppc_stat_2d(stan_data$lice, y_ppd, stat = c("mean", "sd"))

# Dharma residuals - but I just use all 4000 samples?
dharma_resids <- function(model, iter, response, y_ppd = "y_ppd"){
  #extract model predictions based on our generated quantities section (the y_ppd variable)
  predictions <- as.data.frame(rstan::extract(model,
                                              pars = y_ppd)[[1]]) %>%
    #get things in the right format
    rownames_to_column() %>%
    pivot_longer(!rowname,
                 names_to = "obs",
                 values_to = "y_ppd") %>%
    mutate(rowname = as.numeric(rowname))
  # randomly select 250 draws
  sim_response_matrix <-
    data.frame(rowname = sample(1:iter, 4000, replace = FALSE)) %>%
    #left_join will grab only the rows in predictions that line up with the row
    #names we randomly generated above
    left_join(predictions, by = "rowname") %>%
    #remove obs so we can pivot wider without pissing it off
    dplyr::select(-obs) %>%
    group_by(rowname) %>%
    mutate(row = row_number()) %>%
    #get everything back in the right format
    pivot_wider(id_cols = everything(), names_from = 'rowname', values_from = 'y_ppd') %>%
    dplyr::select(-row) %>%
    as.matrix()
  #calculate the residuals using the simulated responses and the actual raw data
  resids <- DHARMa::createDHARMa(simulatedResponse = sim_response_matrix,
                                 observedResponse = response,
                                 fittedPredictedResponse = apply(t(sim_response_matrix), 2, mean),
                                 integerResponse = FALSE)
  return(resids)
}

#run our custom function on our model and with our raw data
resids <- dharma_resids(fit, iter=nrow(y_ppd), stan_data$lice)
plot(resids)

# plot the mean fits to the data
mean_lice <- extract(fit)$mu
mean_mat <- data.frame(mean_fit=NA,l.95.ci=NA,u.95.ci=NA)
#then, since we have one column for each observation in our original dataset,
#we'll run a for loop to do each calculation for each column
for(i in 1:ncol(mean_lice)){
  #take the median within each column
  mean_mat[i,1]=mean(mean_lice[,i])
  #the the upper and lower CIs
  mean_mat[i,2]=quantile(mean_lice[,i],0.025)
  mean_mat[i,3]=quantile(mean_lice[,i],0.975)
}

# mean posterior predictives to the data
mean_pred_mat <- data.frame(mean_pred=NA,l.95.pi=NA,u.95.pi=NA)
for(i in 1:ncol(mean_lice)){
  mean_pred_mat[i,1]=mean(y_ppd[,i])
  mean_pred_mat[i,2]=quantile(y_ppd[,i],0.025)
  mean_pred_mat[i,3]=quantile(y_ppd[,i],0.975)
}

#connect these to our original dataframe to make plotting easier
df_plot2 <- birds %>% 
  bind_cols(as.data.frame(mean_mat), as.data.frame(mean_pred_mat))

ggplot(data = df_plot2, aes(x = scale_allo, y = louse)) +
  geom_ribbon(aes(x = scale_allo, y = mean_pred, ymin = l.95.pi, ymax = u.95.pi),fill='darkgreen', alpha = 0.2) +
  geom_line(aes(x=scale_allo, y = mean_pred), lty=5,lwd=0.8, col = "darkgreen") +
  geom_ribbon(aes(x = scale_allo, y = mean_fit, ymin = l.95.ci, ymax = u.95.ci),fill='darkcyan', alpha = 0.2) +
  geom_line(aes(x=scale_allo, y = mean_fit), lwd=0.8, col = "darkcyan") +
  geom_point() +
  theme_classic() +
  labs(x = "% time allopreening", y = "No. lice")


# plot the smoothed relationships from our generated quantities
mean_lice <- extract(fit)$new_mu
mean_mat <- data.frame(mean_fit=NA,l.95.ci=NA,u.95.ci=NA)
#then, since we have one column for each observation in our original dataset,
#we'll run a for loop to do each calculation for each column
for(i in 1:ncol(mean_lice)){
  #take the median within each column
  mean_mat[i,1]=mean(mean_lice[,i])
  #the the upper and lower CIs
  mean_mat[i,2]=quantile(mean_lice[,i],0.025)
  mean_mat[i,3]=quantile(mean_lice[,i],0.975)
}
# make a plot of each preening variable when the other equals 0 (their average)
pred_plot <- x_predict[x_predict$scale_self==0,] %>% 
  bind_cols(as.data.frame(mean_mat[x_predict$scale_self==0,]))

g1 <- ggplot(data = birds, aes(x = scale_allo, y = louse)) +
  geom_ribbon(data = pred_plot,aes(x = scale_allo, y = mean_fit, ymin = l.95.ci, ymax = u.95.ci),fill='darkcyan', alpha = 0.2) +
  geom_line(data = pred_plot,aes(x=scale_allo, y = mean_fit), lwd=0.8, col = "darkcyan") +
  geom_point(data = birds, aes(x = scale_allo, y = louse),fill="grey",pch=21) +
  theme_classic() +
  labs(x = "% time allopreening", y = "No. lice")

# repeate, but for the other variable
pred_plot <- x_predict[x_predict$scale_allo==0,] %>% 
  bind_cols(as.data.frame(mean_mat[x_predict$scale_allo==0,]))

g2 <- ggplot(data = birds, aes(x = scale_self, y = louse)) +
  geom_ribbon(data = pred_plot,aes(x = scale_self, y = mean_fit, ymin = l.95.ci, ymax = u.95.ci),fill='darkcyan', alpha = 0.2) +
  geom_line(data = pred_plot,aes(x=scale_self, y = mean_fit), lwd=0.8, col = "darkcyan") +
  geom_point(data = birds, aes(x = scale_self, y = louse),fill="grey",pch=21) +
  theme_classic() +
  labs(x = "% time self-preening", y = "No. lice")

bayes_diag <- ggpubr::ggarrange(g2,g1, ncol = 1,nrow=2)
print(bayes_diag)
