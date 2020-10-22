load("rerandN50sharp")
mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)

load("rerandN1000sharp")
mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)

load("rerandN50weak")
mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)

load("rerandN1000weak")
mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)


par(mfrow = c(2,1), mgp=c(1.8,0.5,0), mar=c(3,3,2,1), oma = c(0,0,2,0))
load("rerandN50sharp")
ecdfref = function(x){x}
ecdfstat = ecdf(1-refdistpre)
r = range(refdistpre) + c(-.01, + .01)
curve(ecdfref(x), from = r[1], to = r[2], main = "Sharp Null, N=50", lty = 2, lwd = 2, ylab = expression(paste("P(", "P-value"<= x, ")", sep ="")), xlab = "x", col = "#9ecae1")
curve(ecdfstat(x), from = r[1], to = r[2], col = "orange", add = T, lty = 1, lwd = 2)
legend("bottomright", c("True Distribution of P-Values", "Uniform(0,1)"), lty = c(1,2), col = c("orange", "#9ecae1"), bty = "n", lwd = c(2,2))

load("rerandN1000sharp")
ecdfref = function(x){x}
ecdfstat = ecdf(1-refdistpre)
r = range(refdistpre)+ c(-.01, + .01)
curve(ecdfref(x), from = r[1], to = r[2], main = "Sharp Null, N=1000", lty = 2, lwd = 2, ylab = expression(paste("P(", "P-value"<= x, ")", sep ="")), xlab = "x", col = "#9ecae1")
curve(ecdfstat(x), from = r[1], to = r[2], col = "orange", add = T, lty = 1, lwd = 2)
title("True Distribution of Large-Sample P-Values", outer = T)


mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)
load("rerandN1000weak")
par(mfrow = c(2,1), mgp=c(1.8,0.5,0), mar=c(3,3,2,1), oma = c(0,0,3,0))
  ecdfref = ecdf(refdiststu)
  ecdfstat = ecdf(statstu)
  r = range(refdiststu)+ c(-.01, + .01)
  curve(ecdfref(x), from = r[1], to = r[2], main = "Studentized, Weak Null, N=1000", lty = 2, lwd = 2, ylab = expression(paste("P(", X <= x, ")", sep ="")), xlab = "x", col = "#9ecae1")
  curve(ecdfstat(x), from = r[1], to = r[2], col = "orange", add = T, lty = 1, lwd = 2)
  legend("bottomright", c("True Distribution", "Reference Distribution"), lty = c(1,2), col = c("orange", "#9ecae1"), bty = "n", lwd = c(2,2))

  ecdfref = ecdf(refdistpre)
  ecdfstat = ecdf(statpre)
  r = range(refdistpre)+ c(-.01, + .01)
  curve(ecdfref(x), from = r[1], to = r[2], main = "Prepivoted, Weak Null, N=1000", lty = 2, lwd = 2, ylab = expression(paste("P(", X <= x, ")", sep ="")), xlab = "x", col = "#9ecae1")
  curve(ecdfstat(x), from = r[1], to = r[2], col = "orange", add = T, lty = 1, lwd = 2)

  title("True Distribution versus Reference Distribution", outer = T)



load("rerandN1000sharp")
mean(pvalabs <= .1)
mean(pvalstu <= .1)
mean(pvalpre <= .1)
mean(1-statpre <= .1)

hist(refdistpre)
hist(1-statpre)
plot(density(1-statpre,from=0, to = 1) , col = "black", lty = 4, xlab = "", ylab = "", lwd = 3)


load("rerandN1000weak")
hist(refdistpre)
hist(1-statpre)

hist(statstu)
plot(density(1-statpre, kernel = "cosine", from = 0, to = 1))
