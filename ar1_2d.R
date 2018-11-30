library(TMB)
compile("ar1_2d.cpp")
dyn.load(dynlib("ar1_2d"))

library(MASS)
simdata <- function(Dim1 = 10, Dim2 = 10, rho = 0.8, phi = 1, rho2 = 0, sds = rep(1,Dim1), sdObs = rep(1, Dim1), seed = 357246) {
  set.seed(seed)
  
  corrMat <- matrix(0, Dim1, Dim1)
  corrMat2 <- matrix(0, Dim1, Dim1)
  for (i in 1:Dim1) {
    for (j in 1:Dim1) {
      corrMat[i, j] <- rho^abs(i - j)
      corrMat2[i, j] <- rho2^abs(i - j)
    }
  }
  Sigma <- corrMat * (sds %o% sds)
  Sigma2 <- corrMat2 * (sdObs %o% sdObs)
  d <- matrix(NA, Dim2, Dim1)
  obs <- d
  ## init state
  d[1, ] <- rnorm(Dim1)
  i <- 1
  obs[i, ] <- d[i, ] + mvrnorm(1, rep(0, Dim1), Sigma = Sigma2)
  for (i in 2:Dim2) {
    d[i, ] <- phi * d[i - 1, ] + mvrnorm(1, rep(0, Dim1), Sigma = Sigma)
    obs[i, ] <- d[i, ] + mvrnorm(1, rep(0, Dim1), Sigma = Sigma2)
  }
  return(list(d = d, obs = obs, sds = sds, sdObs = sdObs))
}

Dim1 <- 10
Dim2 <- 10


## Simulate data
sim <- simdata(Dim1, Dim2)

d <- sim$d
obs <- sim$obs

data <- list(obs = t(obs))
parameters <- list(transf_rho = 0.5, transf_rho2 = 0, transf_phi = 0.5, logsds = sim$sds * 0, logsdObs = sim$sdObs * 
                     0, u = d)


fixed <- factor(NA)
obj <- MakeADFun(data, parameters,
                 random = "u", 
                 DLL = "ar1_2d",
                 map = list(logsds = rep(fixed,Dim1),
                            logsdObs = rep(fixed,Dim1)))
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
report <- sdreport(obj,opt$par)
sdreport <- summary(report,c("fixed","report"))
head(sdreport,6)

onestepresid <- oneStepPredict(obj, "obs", method = "fullGaussian")

res1 <- onestepresid$residual
# K-S test for normality 
ks1<-ks.test(res1,pnorm,0,1) 
ks1$p.value
# Ljun-Box test
Box.test(res1, lag = 1, type = c( "Ljung-Box"), fitdf = 0)

Box.test(as.numeric(t(matrix(res1,10))), lag = 1, type = c( "Ljung-Box"), fitdf = 0)

hist(res1)
qqnorm(res1)
abline(0,1)
plot(res1, obs)
abline(0,0)

########################################################################################
## Wrong, model 2, phi = 0 no dependence on the time t
########################################################################################

map <- list(transf_rho2 = fixed, 
            transf_phi = fixed, 
                       logsds = rep(fixed,Dim1),
                       logsdObs = rep(fixed,Dim1))
parameters$transf_rho = 0.8
parameters$transf_phi = 0

obj2 <- MakeADFun(data, parameters, random = "u", DLL = "ar1_2d", 
                 map = map)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
opt2
report2 <- sdreport(obj2,opt2$par)
sdreport2 <- summary(report2,c("fixed","report"))
head(sdreport2,4)

onestepresid2 <- oneStepPredict(obj2, "obs", method = "fullGaussian")

res2 <- onestepresid2$residual
# K-S test for normality 
ks2<-ks.test(res2,pnorm,0,1) 
ks2$p.value
# Ljun-Box test
Box.test(res2, lag = 1, type = c( "Ljung-Box"), fitdf = 0)

Box.test(as.numeric(t(matrix(res2,10))), lag = 1, type = c( "Ljung-Box"), fitdf = 0)

hist(res2)
qqnorm(res2)
abline(0,1)
plot(res2, obs)
abline(0,0)

########################################################################################
## Repeat with wrong model 2 repeated measurements on states are independent 
########################################################################################

map <- list(transf_rho2 = fixed,
            transf_rho = fixed, 
            logsds = rep(fixed,Dim1),
            logsdObs = rep(fixed,Dim1))
parameters$transf_rho = 0
parameters$transf_phi = 0.5


obj3 <- MakeADFun(data, parameters, random = "u", DLL = "ar1_2d", 
                  map = map)
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
opt3
report3 <- sdreport(obj3,opt3$par)
sdreport3 <- summary(report3,c("fixed","report"))
head(sdreport3,5)

onestepresid3 <- oneStepPredict(obj3, "obs", method = "fullGaussian")

res3 <- onestepresid3$residual
# K-S test for normality 
ks3<-ks.test(res3,pnorm,0,1) 
ks3$p.value
# Ljun-Box test
Box.test(res3, lag = 1, type = c( "Ljung-Box"), fitdf = 0)

Box.test(as.numeric(t(matrix(res3,10))), lag = 1, type = c( "Ljung-Box"), fitdf = 0)

hist(res3)
qqnorm(res3)
abline(0,1)
plot(res3, obs)
abline(0,0)


