norm.lik = function(bi,sig_e,sig_b,beta,stream,dat)
{
  nobs=3 #number of observations per stream
  lik.a=lik.fix=rep(NA,nobs)
  sub_dat=subset(dat,dat$pop==stream)
  lik.random=1/sqrt(2*pi)*exp(-(bi^2)/2)
  for(k in 1:nobs)
  {
    lik.fix[k]=1/sqrt(2*pi*sig_e^2)*exp(-((sub_dat$x[k]-beta-bi*sig_b)^2)/(2*sig_e^2))
  }
  lik.fix.a=prod(lik.fix)
  lik.a=lik.random*lik.fix.a
  return(lik.a)
}

simpson = function(lowerb, upperb, nbins, stream,beta, sig_e, sig_b) 
{
  step = (upperb-lowerb)/nbins
  b_vec = seq(lowerb, upperb, by=step)
  frac_vec = c(1/3, rep(c(4/3, 2/3), times=((nbins+1)-3)/2), 4/3, 1/3)
  sum_comp = vector(length=length(b_vec))
  for(i in 1:length(b_vec)){
    sum_comp[i] = frac_vec[i]*norm.lik(bi=b_vec[i], sig_e=sig_e, sig_b=sig_b, beta=beta, stream=stream)
  }
  full_like = sum(sum_comp)*step
  return(full_like)
}

loglike=function(l_sig_e,l_sig_b,beta)
{
  sig_e=exp(l_sig_e)
  sig_b=exp(l_sig_b)
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