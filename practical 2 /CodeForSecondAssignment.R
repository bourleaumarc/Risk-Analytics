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
# the CI includes 0 for the shape parameter, so not sure we're that heavy-tailed
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




######################################################################################################################
library(xts)
library(extRemes)
library(RCurl)
# get data from eHYD

# adapted from example available at: https://www.gis-blog.com/eva-intro-3/

read.ehyd <- function(ehyd_url) {
  # separate the header, open the connection with correct encoding
  con <- url(ehyd_url, encoding = "latin1")
  header <- readLines(con, n=50)
  lines.header <- grep("Werte:", header, fixed = T)
  # read data, define time and values
  infile <- read.csv2(con, header = F, skip = lines.header,
                      col.names = c("time", "value"),
                      colClasses = c("character", "numeric"),
                      na.strings = "Lücke",
                      strip.white = TRUE, as.is = TRUE, fileEncoding = "latin1")
  infile$time <- as.POSIXct(infile$time, format = "%d.%m.%Y %H:%M:%S")
  # return time series object of class xts
  return(xts(infile$value, order.by = infile$time))
}

ehyd_url <- "http://ehyd.gv.at/eHYD/MessstellenExtraData/nlv?id=105700&file=2"
precipitation_xts <- read.ehyd(ehyd_url)
precipitation_xts <- precipitation_xts[-14611,]

# Visualise the data
plot(precipitation_xts)
hist(precipitation_xts, breaks = 30)
# seems to be right-skewed 

# MRL plot: 
mrlplot(precipitation_xts, main="Mean Residual Life Plot")
# Find the lowest threshold where the plot is somewhat linear

# fitting the GPD model over a range of thresholds
threshrange.plot(precipitation_xts, r = c(30, 45), nint = 16)
# nint is the number of threshold, r the range of the thresholds 
# ismev implementation is faster:
ismev::gpd.fitrange(precipitation_xts, umin=30, umax=45, nint = 16)
# set threshold
th <- 40

# Visualise the threshold 
plot(as.vector(precipitation_xts), type = 'l')
abline(h=th,col=2)

# maximum likelihood estimation
pot_mle <- fevd(as.vector(precipitation_xts), method = "MLE", type="GP", threshold=th)
# diagnostic plots
plot(pot_mle)
rl_mle <- return.level(pot_mle, conf = 0.05, return.period= c(2,5,10,20,50,100), do.ci=T)

# alternative using evd 
pot_mle_evd <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25)
pot_mle_evd2 <- fpot(as.vector(precipitation_xts), threshold = th)
par(mfrow = c(2,2))
plot(pot_mle_evd)
confint(profile(pot_mle_evd))

# return levels with evd, e.g. 50-year
rl_mle_evd <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25, mper=50)
plot(rl_mle_evd)
# return level plots
par(mfcol=c(1,1))
# return level plot w/ MLE
plot(pot_mle, type="rl",
     main="Return Level Plot for Oberwang w/ MLE",
     ylim=c(0,200), pch=16)
loc <- as.numeric(return.level(pot_mle, conf = 0.05, return.period=50))
segments(50, 0, 50, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,50, loc, col='midnightblue', lty=6)


# Check if there is clustering of extremes (Ferro-Segers estimate)
exi(precipitation_xts, u = 35)
# seems to indicate that there slighty is some level of clustering, but weak. 

pot_decl <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25, cmax = TRUE, ulow = th, mper = 50)
# some changes in the estimates because we have less data, if there is clustering of extremes (expected)
confint(profile(pot_decl))
# compare with
confint(profile(rl_mle_evd))

plot(pot_decl)

# 
# decluster and obtain plot, according to years for example
years <- c()
k <- 1
for(i in 1:nrow(precipitation_xts)){
  if(is.na(precipitation_xts[i])){
    next
  }else{
    years[k] <- year(precipitation_xts[i])
    k <- k+1
  }
}
years <- years-1999

decl1 <- decluster(as.vector(precipitation_xts), threshold=th, groups=years, na.action=na.omit)
decl1
plot(decl1)


## mini example of bootstrap 
