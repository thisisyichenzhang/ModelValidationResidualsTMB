## This file serves as the R code for validation of beaver model
## with prediction quantile residuals

## Model and Data: 
## Reynolds, P. S. (1994)
## Time-series analyses of beaver body temperatures. 

library(datasets)
library(TMB)
library(ggplot2)
data(beavers)

#################### Plotting ################

## Beaver Temperature Time-series
with(beaver2 , {
  time <- strptime(paste(1990, day, time %/% 100, time %% 100),
                   "%Y %j %H %M")
  plot(time, temp, type = "l") # axis at 4-hour intervals.
  # now label every hour on the time axis
  plot(time, temp, type = "l",xaxt = "n", xlab = "Time", ylab = "Temperature", main="Time Series of Beaver Temperature")
  # add points colored by activity
  points(time, temp,col = factor(activ))
  legend("topleft", legend=c("Dormant","Active"), pch=16, col=unique(factor(activ)))
  r <- as.POSIXct(round(range(time), "hours"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = "hour"), format = "%H:00")
})

## Boxplot of Temperatures
with(beaver2, {
  activ<-ifelse(beaver2$activ==1,"Active","Dormant")
  ggplot2::qplot(x=activ,y=temp,geom="boxplot", colour=activ)+
    theme_bw(base_size = 12, base_family = "")+
    labs(title = "Boxplot of Temperatures",x = "Activity State",y="Temperature",colour="State")
})
################### END of Plotting #############


set.seed(1357121)
## Data wrangling
temp <- beaver2$temp + rnorm(length(beaver2$temp),0,0.1)

activity <- factor(beaver2$activ)
## Compile the C++ code
compile("beavervalidation.cpp")
## Load the dynamic library
dyn.load(dynlib("beavervalidation"))
## Create Optimization Object
obj <- MakeADFun(data = list(Y = temp,
                             act = activity),
                 parameters = list(logitPhi = 0,
                                   logSdState = 0,
                                   logSdObs = 0,
                                   mu = numeric(nlevels(activity)),
                                   X = numeric(length(temp))),
                 silent = TRUE, # make the output neat
                 random = "X",  # Specify the random effect/Use Laplace Approximation
                 DLL = "beavervalidation")

## par: Parameters
## fn: The likelihood function
## gr: The gradient fucntion 
opt <- nlminb(obj$par,obj$fn,obj$gr)
opt

## Estimates and Standard Error
report <- sdreport(obj,opt$par)
sdreport <- summary(report,c("fixed","report"))
sdreport # two mu parameters corresponds to mu of two states
 

## One-step Prediction Residual
onestepresid <- oneStepPredict(obj, observation.name = "Y",
                        data.term.indicator = "keep",
                        method = "oneStepGaussianOffMode")

plot(onestepresid$residual)
t.test(onestepresid$residual)

# Wrong model:  two mu are the same
obj2 <- MakeADFun(data = list(Y = temp,
                              act = activity),
                  parameters = list(logitPhi = 0,
                                    logSdState = 0,
                                    logSdObs = 0,
                                    mu = numeric(nlevels(activity)),
                                    X = numeric(length(temp))),
                  map = list(mu = factor(c(1,1))), # set mu[1] and mu[2] to be the same
                  silent = TRUE,
                  random = "X",
                  DLL = "beavervalidation")
opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr)
opt2
onestepresid2 <- oneStepPredict(obj2, observation.name = "Y",
                        data.term.indicator = "keep",
                        method = "oneStepGaussianOffMode")

onestepresid2
plot(onestepresid2$residual)
t.test(onestepresid2$residual)


#likelihood ratio test
testval <- 2 * (opt2$objective - opt$objective)

# p-value 
1 - pchisq(testval,length(opt$par) - length(opt2$par))

# plots
sumsdr <- summary(sdr)
Xest <- sumsdr[rownames(sumsdr) == "X",]
ciLow <- c(1,-2) %*% t(Xest)
ciHigh <- c(1,2) %*% t(Xest)

mu <- ss[rownames(ss) == "mu",]
muCiL <- c(1,-2) %*% t(mu)
muCiH <- c(1,2) %*% t(mu)
