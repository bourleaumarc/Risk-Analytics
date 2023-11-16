library(xts)
library(extRemes)
library(RCurl)
library(evd)
## PEAK OVER THRESHOLD APPROACH 
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
plot(precipitation_xts) # xts allow us to plot time-series 
hist(precipitation_xts, breaks = 30) # there are a lot of data points close to 0. It looks long tail. 
# we can check with mean() and max() and quantile(precipitation, 0.95) this suggests we need extreme analysis 
# seems to be right-skewed 
########## FIRST THING: what is the threshold? mrlplot is  ########## 
# MRL plot: 
mrlplot(precipitation_xts, main="Mean Residual Life Plot") # mean excess (straight line) and dashed line (confidence intervals)
# We need to select a threshold that it should look linear. We should take 40-45 (very linear)
# Find the lowest threshold where the plot is somewhat linear

# fitting the GPD model over a range of thresholds
threshrange.plot(precipitation_xts, r = c(30, 45), nint = 16) # you see the parameters for a specific threshold (between 30 and 45). Tradeoff between 
# variance and bias. It looks like the confidence interval is limited to 40. 
# nint is the number of threshold, r the range of the thresholds 
# ismev implementation is faster:
ismev::gpd.fitrange(precipitation_xts, umin=30, umax=45, nint = 16) # we use this one for part 3 
# set threshold
th <- 40 # you can see what quantile does this threshold represent. corresponds to 99 quantile-> quantile(precipation, 0.99)
## main assumption of the data: independence. 
# Visualise the threshold 
plot(as.vector(precipitation_xts), type = 'l')
abline(h=th,col=2)

# maximum likelihood estimation
pot_mle <- fevd(as.vector(precipitation_xts), method = "MLE", type="GP", threshold=th) 
# diagnostic plots
plot(pot_mle) 
# fevd will always assume the data to be daily and the period to be yearly. 
rl_mle <- return.level(pot_mle, conf = 0.05, return.period= c(2,5,10,20,50,100), do.ci=T)
# CI increases with the period, which is logical. 
# alternative using evd 
#pot_mle_evd <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25)
#pot_mle_evd2 <- fpot(as.vector(precipitation_xts), threshold = th)
#par(mfrow = c(2,2))
#plot(pot_mle_evd)
#confint(profile(pot_mle_evd))

# return levels with evd, e.g. 50-year
#rl_mle_evd <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25, mper=50)
#plot(rl_mle_evd)
# return level plots
par(mfcol=c(1,1))
# return level plot w/ MLE: 
plot(pot_mle, type="rl",
     main="Return Level Plot for Oberwang w/ MLE",
    ylim=c(0,200), pch=16)
loc <- as.numeric(return.level(pot_mle, conf = 0.05, return.period=50))
segments(50, 0, 50, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,50, loc, col='midnightblue', lty=6)
# for fevd you need to decluster the data first. 
#### 
# Check if there is clustering of extremes (Ferro-Segers estimate)
exi(precipitation_xts, u = th) # extremalindex for fevd : 1/extremal index = approximate cluster size; only do clustering when the inverse is higher than 2. 
# seems to indicate that there slighty is some level of clustering, but weak. 
## npp: data you have in your data set to specific your period. ?fpot in this case 365.25 
pot_decl <- fpot(as.vector(precipitation_xts), threshold = th, npp = 365.25, cmax = TRUE, ulow = th, mper = 50)
# return level has changed. tells you the number of clusters. 
# some changes in the estimates because we have less data, if there is clustering of extremes (expected)
confint(profile(pot_decl))
# compare with
confint(profile(rl_mle_evd))

plot(pot_decl)


# We need to use this declustering function for FEVD. decluster and obtain plot, according to years for example 
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
years <- years-1980 # 1980 is the first year 

decl1 <- decluster(as.vector(precipitation_xts), threshold=th, groups=years, na.action=na.omit) # vector and groups need to have the same size 
decl1
plot(decl1) # shows in grey the points that are not retained in the selection. 
## mini example of bootstrap 
# 1st threshold 
# 2nd fit the model and mak
# 3rd 
