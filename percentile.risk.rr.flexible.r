# Function to get the percentile specific absolute risks
# betat: vector of the age specific betas
# nbins: number of bins to approximate the normal distribution. Use at least ~1000
# inc: average incidence
# Default values for "inc" are BRCA1 carrier breast cancer incidences from Kuchenbaecker et al 2017 JAMA
percentile.risk.rr.flexible <- function(betat=rep(0,80),nbins=1000,
inc=c(rep(0,20),rep(5.9,10),rep(23.5,10),rep(28.3,10),rep(25.7,10),rep(25.0,10),rep(16.5,10))/1000,
percentiles=c(5,10,50,90,95),min.max=FALSE) {
# BRCA2 carrier breast cancer incidences [Kuchenbaecker, Hopper, Barnes et al 2017 JAMA]
#   c(rep(0,20),rep(4.8,10),rep(10.8,10),rep(27.5,10),rep(30.6,10),rep(22.9,10),rep(21.9,10))/1000
# BRCA1 carrier ovarian cancer incidences [Kuchenbaecker, Hopper, Barnes et al 2017 JAMA]
#   c(rep(0,20),rep(0,10),rep(1.8,10),rep(7.0,10),rep(13.8,10),rep(29.4,10),rep(5.7,10))/1000
# BRCA2 carrier ovarian cancer incidences [Kuchenbaecker, Hopper, Barnes et al 2017 JAMA]
#   c(rep(0,20),rep(0,10),rep(0.3,10),rep(0,10),rep(6.5,10),rep(10.3,10),rep(2.3,10))/1000

stopifnot(class(betat)=="numeric")
stopifnot(class(nbins)=="numeric" | class(nbins)=="integer")
stopifnot(class(inc)=="numeric")

stopifnot(nbins>0)
stopifnot(length(betat)==length(inc))

stopifnot(class(percentiles)%in%c("numeric","integer"))
stopifnot(min(percentiles)>0)
stopifnot(max(percentiles)<100)

stopifnot(min(inc)>=0)


# log-RR
logrr.mat <- matrix(rep(NA,length(inc)*(nbins+1)),(nbins+1),length(inc))
for(i in 1:(nbins+1)){
	logrr.mat[i,] <- betat*((i-1)-nbins/2)/(sqrt(nbins/4))
}


# Percentile frequencies
vec1.freq <- c(rep(NA,(nbins+1)))
for(i in 1:(nbins+1)) {
	vec1.freq[i] <- dbinom(i-1,nbins,0.5)
}
cumfreq <- cumsum(vec1.freq)
combined.freq.mat <- data.frame(vec1.freq,cumfreq)
names(combined.freq.mat) <- c("freq","cumfreq")


# Matrix of survival functions ngen x length(inc)+1
surv.mat <- matrix(rep(NA,(nbins+1)*(1+length(inc))),(nbins+1),(1+length(inc)))
lambda0 <- rep(NA,length(inc))
lambda.mat <- matrix(rep(NA,(nbins+1)*length(inc)),(nbins+1),length(inc))


# First column is for age 0 so survival is 1
surv.mat[,1] <- rep(1,nbins+1)
surv.mat[,2] <- rep(1,nbins+1)
lambda.mat[,1] <- rep(0,nbins+1)
lambda0[1] <- 0
group <- c(1:(nbins+1))
for(i in 2:length(inc)) {
  lambda0[i] <- inc[i]*sum(combined.freq.mat$freq*surv.mat[,i])/(sum(combined.freq.mat$freq*exp(logrr.mat[,i])*surv.mat[,i]))
  lambda.mat[,i] <- lambda0[i]*exp(logrr.mat[,i])
  surv.mat[,i+1] <- exp(-rowSums(lambda.mat[,1:i]))
}


# Find the percentiles
per <- sort(unique(percentiles))/100
for(i in 1:length(per)) {
  r <- dim(combined.freq.mat[which(combined.freq.mat$cumfreq<=per[i]),])[1]
  if(i==1) { is <- r } else { is <- c(is,r) }
}
is[length(is)] <- is[length(is)]+1
is <- c(1,is,(nbins+1))

ri <- data.frame(t(1-surv.mat[is,]))
ri <- cbind(0:length(inc),ri)
colnames(ri)[1:2] <- c("age","min")
colnames(ri)[dim(ri)[2]] <- "max"
colnames(ri)[3:(dim(ri)[2]-1)] <- paste0(percentiles,"%")

if(min.max==FALSE) {
  ri <- ri[,setdiff(colnames(ri),c("min","max"))]
}

# Return results
return(ri)
}
