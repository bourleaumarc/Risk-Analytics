## Simulate Gumbel distribution 
# Go through quantile function of the Gumbel distribution, using unif(0,1) variates
x <- runif(1000)
mu <- 1
beta <- 5
y <- mu-beta*log(-log(x))
plot(y, x) # cdf of a Gumbel distribution 
seqx <- seq(-50, 50, by=0.1)
z <- (seqx-mu)/beta
plot(seqx, 1/beta*exp(-z-exp(-z)), type="l")
# You can do pretty much the same for the other extreme-value distributions
# Anderson-Darling test
library(DescTools)
AndersonDarlingTest(y, null="pGumbel", loc=mu, scale=beta)

# Fitting a GEV distribution to Gumbel observations
library(ismev)
mod0 <- gev.fit(y)
# Does it fit? 
qqplot(rGumbel(length(y), loc=mu, scale=beta), y)
qqline(y, distribution=function(p) qGumbel(p, loc=mu, scale=beta))

# What about a distribution you might not know which domain of attraction it falls into? 
library(lubridate)
obs.date <- seq(from=as.Date("1993-01-01", format = "%Y-%m-%d"), to=as.Date("2022-12-31", format = "%Y-%m-%d"), by="day")
obs <- rt(length(obs.date), df=1.5)
df <- data.frame(obs.date)
df <- data.frame(df, obs)
max <- data.frame(0, 0, 0)
years <- unique(year(df$obs.date))
months <- unique(month(df$obs.date))
for(i in 1:length(years)){
  for(j in 1:length(months)){
    subset <- df[(month(df$obs.date)==months[j])&(year(df$obs.date)==years[i]),]
    max <- rbind(max, c(month=months[j], year=years[i], max(subset$obs)))
  }
}
max <- max[-1, ]
colnames(max) <- c("month", "year", "maximum")

mod1 <- gev.fit(max$maximum)
gev.diag(mod1) #diagnostics plot 

mod2 <- gev.fit(max$maximum, max, mul=2)
mod2
1-pchisq(-2*(mod2$nllh-mod1$nllh),1) #chi-square test 
# The trend is not significant, this is expected! 

## Data example
## Reproduced from great walkthrough by Hugo Winter 
## available at: https://younghydrologicsociety.files.wordpress.com/2018/04/eva_training_exercises_egu_2018.pdf
data(rain) #from ismev
years <- rep(1:48, rep(c(365,365,366,365), times = 12))[-17532]
rain.ann.max <- unlist(lapply(X = split(rain,years), FUN = max))
# annual maxima from this set of data in South-West England
# from extRemes package, can use fevd()
library(extRemes)
mod3 <- fevd(rain.ann.max, type="GEV", time.units="years")
plot(mod3)
mod3$results$par #gives parameters of the GEV
# suggest heavy-tailed model here, but only point estimates.
# What about building confidence intervals? 
ci.fevd(mod3, alpha=0.05, type="parameter")
# the CI includes 0, so not sure we're that heavy-tailed
gev.rl <- return.level(x = mod3, return.period = c(10,100,1000,10000),
                       do.ci = TRUE, alpha = 0.05)
gev.rl
#notice the very wide CI for the 10000-year return level. 
# We clearly do not have enough data and the uncertainty
# concerning the shape parameter drives this. 

## How to compute return levels (by hand) using GEV
# e.g. in the rain example, compute the 10-year return level
# understand here: the value exceeded one out of every 10 365-days block (1 year)
# approximated by
as.numeric(mod3$results$par[1]+mod3$results$par[2]*((-log(1-1/10))^(-mod3$results$par[3])-1)/mod3$results$par[3])
# or, as it is the 1-1/k = 1-1/10th quantile
qgev(1-1/10, location=mod3$results$par[1], scale=mod3$results$par[2], shape=mod3$results$par[3])

## Return period associated with level u
# recall: 1/(1-H(u))
as.numeric(1/(1-pgev(qgev(1-1/10, location=mod3$results$par[1], scale=mod3$results$par[2], shape=mod3$results$par[3]), location=mod3$results$par[1], scale=mod3$results$par[2], shape=mod3$results$par[3])))
# the above corresponds to 1/10 period as expected! 

unique.years <- unique(years)
plot(unique.years, rain.ann.max)
# For simplicity: suppose you want to assess the number of blocks (years here) needed to exceed level 70, i.e. {Mn > 50}
as.numeric(1/(1-pgev(70, location=mod3$results$par[1], scale=mod3$results$par[2], shape=mod3$results$par[3])))
# you expect an average of 14.01 blocks (years) to wait before exceeding 70.  

# making predictions with linear models (remainder)
mod4 <- lm(rain.ann.max~unique.years)
# making predictions for 5 years
predict.lm(mod4,newdata=data.frame("unique.years"=c(49:53)),se=T)