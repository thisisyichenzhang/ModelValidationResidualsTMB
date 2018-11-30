library(TMB)
library(ggplot2)# for plotting
compile("poisvalidation.cpp")
dyn.load(dynlib("poisvalidation"))

## Simulate the date 
set.seed(1357236)
N <- 200
r=0.8 # Growth rate 
K=10 # Carrying capacity
Q=0.1 # Process noise 
S=50 # Sample volume controlling measurement uncertainty
p <- 0.2 # probability of structural zero in zero-inflated Poisson model

X <- numeric(N)
Y <- numeric(N)
proc.err=rnorm(200,sd=sqrt(Q))

X[1]=3.0
for(t in 2:200){
  X[t] <- X[t-1] + r * ( 1 - (exp(X[t-1])/K) ) + proc.err[t]; 
}
# Simulate observations with zero-inflated poission 
library("VGAM") # zero-inflated possion simulation
Y <- rzipois(length(X), lambda = S*exp(X),  pstr0 = p) 




## TMB Data
data <- list(Y=Y)

## Initial guess on parameters for the fitted Ricker model
parameters0 <- list(
  X=X,
  logr=log(r),
  logK=log(K),
  logQ=log(Q),
  logS=log(S)
)

## Make TMB model 
fixed <- factor(NA)
obj <- MakeADFun(data,
                 parameters0,
                 random="X",
                 DLL="poisvalidation",
                 map=list(logS=fixed))

## Fit the model 
## par:  Parameters
## fn: The likelihood function
## gr: The gradient function 
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt
opt[c("par","objective","convergence")]

## Estimates and Standard Error
report <- sdreport(obj,opt$par)
sdreport <- summary(report,c("fixed","report"))
sdreport 

## Generate one step predictions
onestepresid <- oneStepPredict(obj,observation.name="Y",data.term.indicator="keep",
                               method="oneStepGeneric",discrete=TRUE,range=c(0,Inf),
                               conditional=1, ## Skip first residual
                               parallel=FALSE)

# K-S test for normality 
ks1<-ks.test(onestepresid$residual,pnorm,0,1) 
ks1$p.value

# Ljun-Box test
Box.test(onestepresid$residual, lag = 1, type = c( "Ljung-Box"), fitdf = 0)

### Plots

width <- 4
height <- 4
## Plot residual against previous observation:
pdf(file="Resid2.pdf",width=width,height=height)

plotdata<- data.frame(obs = head(Y,-1),
                      res = onestepresid$residual)
ggplot(plotdata,aes(obs,res))+
  geom_point(shape = 1,na.rm = T)+
  geom_smooth(na.rm = T)+
  geom_abline(slope = 0,intercept = 0)+
  labs(x=expression(y[i]),y=expression(z[i+1]))

dev.off()

with(plotdata,summary(lm(res ~ obs)))
## histogram 
hist(onestepresid$residual)
## qqplot
pdf(file="poisqq.pdf",width=width,height=height) 
qqnorm(onestepresid$residual)
abline(0,1)
dev.off()