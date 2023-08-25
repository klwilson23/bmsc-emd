library(nlme)
library(graphics)
rockfish <- read.csv(file="Data/rockfish.csv",header=TRUE)
rockfish_spp <- c("chilippeper","bocaccio","widow","canary","black","yellowtail","perch_us","perch_goose","perch_ak","perch_ai","perch_ebs")
rockfish$species <- factor(rockfish_spp[rockfish$SP],levels=rockfish_spp)

m1 = lme(log(REC*SBPR/SSB) ~ SSB:species,data=rockfish,random = ~1|species,method="REML")

boxplot(split(residuals(m1),rockfish$species),ylab="Residual",xlab="Rockfish Species",csi=0.2)
abline(h=0,col="red",lty=2,lwd=2)
ln_a_est <- m1$coefficients$fixed[[1]]
print(ln_a_est+m1$coefficients$random[[1]])
m1$coefficients$random$species[,1]
spp=c(5,4,8,6,7,2,11,1,9,3,10)
plot(spp,ln_a_est+m1$coefficients$random[[1]],pch=16,cex=1.5,xlab="Rockfish Species",ylab="Most Likely Random Effect Conditional upon Fixed Effect")
abline(h=ln_a_est,lty=2,col="red",lwd=2)
# ricker: CR * S * e(-log(CR)/K * S)
log(rockfish_traits$alpha[1] / rockfish_traits$SBPR[1])/rockfish_traits$beta[1]
beta <- as.vector(m1$coefficients$fixed[-1])

# finding MSY
rockfish_traits <- data.frame("alpha"=exp(as.numeric(ln_a_est+m1$coefficients$random[[1]])),"beta"=beta,"species"=rockfish_spp)
rockfish_traits$SBPR <- rockfish$SBPR[match(rockfish_traits$species,rockfish$species)]
rockfish_traits$K <- -(rockfish_traits$alpha/rockfish_traits$SBPR)/rockfish_traits$beta

yieldRoot<-function(alpha,beta,SBPR,Ucur)
{
  # population model is R ~ alpha / SSB * S*e(-beta * S); where K is carrying capacity equal to -(alpha/SPBR)/beta
  adults <- -((alpha/SBPR)/beta) * (1-Ucur) # K * (1-U)
  recruits <- alpha / SBPR * adults*exp(beta * adults)
  yield <- recruits*Ucur
  return(list(recruits=recruits,adults=adults,yield=yield))
}

fun_yield <- function(Ucur,alpha,beta,SBPR,delta,model)
{
  y1 <- yieldRoot(alpha=alpha,beta=beta,SBPR=SBPR,Ucur=Ucur-delta/2)$yield
  y2 <- yieldRoot(alpha=alpha,beta=beta,SBPR=SBPR,Ucur=Ucur+delta/2)$yield
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

###############################################
## uniroot the differentiation for Umsy #######
###############################################
findMSY <- function(data){
  U_msy <- NULL
  for(i in 1:nrow(data)){
    a_p <- data$alpha[i]
    b_p <- data$beta[i]
    SBPR_p <- data$SBPR[i]
    U_msy[i] <- uniroot(fun_yield, interval=c(0,1),extendInt="yes",alpha=a_p,beta=b_p,SBPR=SBPR_p, delta=0.0001)$root 
  }
  MSY <- sapply(1:nrow(data),function(x){yieldRoot(alpha=data$alpha[x],beta=data$beta[x],SBPR=data$SBPR[x],Ucur=U_msy[[x]][1])})
  return(list("U_msy"=U_msy,"MSY"=MSY))
}

MSY <- findMSY(data=rockfish_traits)
rockfish_traits$UMSY <- MSY$U_msy
rockfish_traits$R_MSY <- MSY$MSY["recruits",]
rockfish_traits$SSB_MSY <- MSY$MSY["adults",]  
rockfish_traits$MSY <- MSY$MSY["yield",]
