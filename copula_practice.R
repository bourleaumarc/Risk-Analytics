# Tutorial adapted from https://datascienceplus.com/modelling-dependence-with-copulas/ 
library(copula)
library(psych) #for pairs.panels plots
library(VineCopula)
data(loss) #some insurance data from the copula package
# loss correspond to the loss amount and alae is the allocated loss adjustment expense

# visualise the data
df <- loss[,c(1,2)] # focus on the first two columns, we have ties (duplicated columns). 
plot(df) # difficult to assess the structure of dependence on this scale 
# we can't assess the independence / dependnece 
plot(df$loss) # kind of exponential but not too much info 
pairs.panels(df, method="kendall") #we have a distribution of each of the column. Local regression shows the positive association 
# kendall: measure of association.  hard to read as scales are different. We need to transform to a common scale.
# distribution is now uniform in the variables. 
# transformation to pseudo-observations
df.pobs <- pobs(df)
pairs.panels(df.pobs, method="kendall") # we have uniform margins. Now i see something helpful. In the right corner there seems to be more positive association. More concentration of points: positive asymptotic dependence. 
# copulas give you a way to represent points in a parametric fashion. 
# 
# We could retransform the data into any new marginals, with the same dependence structure
x1 <- qgamma(df.pobs[,1], shape=2, scale=1)
x2 <- qbeta(df.pobs[,2], shape1=2, shape2=2)
df2 <- cbind(x1, x2)
pairs.panels(df2)
pairs.panels(df2, method="kendall") # not uniform anymore. we have gamma and beta. we keep the same dependence structure. 

# Build copula: need to be fitted on uniform margina
cop.t <- tCopula(dim=2) # 2 dimensions because two variables. By default degrees of freedom = 4 
# we want to trick the parameters to fit the points 
fit.t <- fitCopula(cop.t, df.pobs, method="ml")
fit.t # rho.1 is the pearson correlation, df is the degrees of freedom as fitted 
coef(fit.t)
# to know if it's good: the lower the better 
aic.t <- 2*2-2*fit.t@loglik # 2*number of parameters - 2*fit.t
??copula
# Build copula 2, specifying an exact form 
# first obtain the correlation by the formula in the slides
rho <- sin(cor(df.pobs, method="kendall")*pi/2)[1,2]
# elliptical distribution: very arbitrary 
gen.t2 <- mvdc(copula=ellipCopula(family="t", param=rho, df=10), margins=c("gamma","gamma"), # shape of the margins
               paramMargins=list(list(shape=4,scale=3), list(shape=3,scale=2)))
gen.t3 <- mvdc(copula=ellipCopula(family="t", param=rho, df=10), margins=c("beta","beta"),
               paramMargins=list(list(shape1=2,shape2=1), list(shape1=3,shape2=2)))
gen.t4 <- mvdc(copula=ellipCopula(family="t", param=rho, df=10), margins=c("t","t"),
               paramMargins=list(list(df=10), list(df=10)))

fit.t2 <- fitCopula(gen.t2@copula, df.pobs, start=c(param=rho, df=10))
fit.t3 <- fitCopula(gen.t3@copula, df.pobs, start=c(param=rho, df=10))
fit.t4 <- fitCopula(gen.t4@copula, df.pobs, start=c(param=rho, df=10))

aic.t2 <- 2*2-2*attributes(fit.t2)$loglik
aic.t3 <- 2*2-2*attributes(fit.t3)$loglik
aic.t4 <- 2*2-2*attributes(fit.t4)$loglik
# all the aic is the same We just changed the margins but the dependence structure does not change. 
# fit a distribution to the margins fitdistr() from MASS this is better. 
# advantage of the second approach, even if no change in AIC, is the generator that mvdc provides

rdm.variates <- rMvdc(10000, gen.t2) # generated random observations from a model 
plot(rdm.variates)

# Automatic Bivariate copula selection: to see if the copula is good or not. 
select.copula <- BiCopSelect(df.pobs[,1], df.pobs[,2], familyset=NA)
select.copula # compute the aic to see if its good or not. 
select.copula$AIC
# check pertinence 
df.pobs2 <- pobs(rdm.variates)
select.copula2 <- BiCopSelect(df.pobs2[,1], df.pobs2[,2], familyset=NA) 
select.copula2 # 
# might not be the best for such model initially, here at n=1000, does not necessarily work
#################### TAIL DEPENDENCE #################### #################### 
# theoretical tail dependence coefficients
# first with copula package, easy with fitted model
lambda(fit.t@copula)  # theoretical tail dependence coefficient --> parametric way !! 
# lower tail = 0.06 and upper tail = 0.06 --> elliptical has the same coefficient for lower and upper 

# non-parametric estimator for the tail dependence 
# can differ from the above, as above there is a theoretical specification 
# non parametric way !! without assuming any model you can find the fit lambda 
fitLambda(df.pobs, lower.tail=F) # so we have the upper tail coefficient, much more than the t copula 
# we can try with the bb copula 
##### we choose: BiCopPar2TailDep 
# asymptotic dependence ? ? 

# compare to chi and chiplot 
library(evd)
chiplot(df)