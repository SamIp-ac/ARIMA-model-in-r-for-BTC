## Time series forecasting model
library(forecast)
library(TSA)
library(tseries)
## BTC.USD
BTC.USD <- read.csv('~/Desktop/Crypto_Data/BTC-USD.csv')
# training set
train<-BTC.USD[1:365,]
dBT<-train$Adj.Close
par(mfrow=c(2,2))
plot.ts(dBT)
acf(dBT)
plot.ts(log(dBT))
acf(log(dBT))
pacf(log(dBT))

# ACF calculation with lag 22
ACF_lag22 <- sum((dBT[1:(365-21)] - mean(dBT[1:(365-21)]))*(dBT[22:365] - mean(dBT[22:365])))/
  sqrt(sum((dBT[1:(365-21)] - mean(dBT[1:(365-21)]))^2)*sum((dBT[22:365] - mean(dBT[22:365]))^2))

# 
try <- rep(FALSE, length(diff(c(1,2,5))))
try[diff(c(1,2,5)) > 2] <- TRUE
#

lndBT<-log(dBT)
dlndBT<-diff(lndBT)

adf.test(lndBT)
adf.test(dlndBT)

par(mfrow=c(2,2))
plot.ts(dlndBT)
acf(dlndBT)
pacf(dlndBT,30)

d2lndBT<-diff(dlndBT)
par(mfrow=c(2,2))
plot.ts(d2lndBT)
acf(d2lndBT)
pacf(d2lndBT,30)

# ARIMA(0,2,1)model
d1<-arima(lndBT,order=c(0,2,1))
tsdiag(d1)
d1$aic

# ARIMA(0,2,4)model
d2<-arima(lndBT,order=c(0,2,4))
tsdiag(d2)
d2$aic

# ARIMA(2,2,0)model
d3<-arima(lndBT,order=c(2,2,0))
tsdiag(d3)
d3$aic

# ARIMA(1,1,1)model
d4<-arima(lndBT,order=c(1,1,1))
tsdiag(d4)
d4$aic

## ARIMA (1,1,1)has the lowest AIC
par(mfrow=c(2,2))
fit1<-arima(lndBT,order=c(1,1,1))
tfore1<-exp(predict(fit1,n.ahead=7)$pred[1:7])
test<-BTC.USD[366:372,]
testA<-test$Adj.Close

ts.plot(cbind(testA,tfore1),lty=1:2,col=c("Blue","Orange"),main="BTC-USD with ARIMA(1,1,1) model")
leg.names<-c("Ture","Prediction")
legend("bottomleft",leg.names,lty=1:2,col=c("Blue","Orange"),box.lty = 0)

##
library(stats)
par(mfrow=c(3,1))
fore
plot(forecast(fit1))
# when nexxt command is error
fit1$x <- ts(lndBT)
fore<-exp(as.vector(forecast(fit1, h=60)$mean))
lower=exp(forecast(fit1, h=60)$lower[,2]);upper=exp(forecast(fit1, h=60)$upper[,2])
ts.plot(c(dBT,fore));lines(c(dBT,lower),lty=2);lines(c(dBT,upper),lty=2)

sum((tfore1-testA)^2)/7
############3
model <- auto.arima(ts(lndBT), seasonal = TRUE, test='kpss', ic='aic')
summary(model)

exp(as.vector((forecast(model, h=6)$mean)))
fore <- exp(predict(model,n.ahead=7)$pred[1:7])

ts.plot(cbind(testA,fore),lty=1:2,col=c("Blue","Orange"),main="BTC-USD with ARIMA(1,1,1) model")
leg.names<-c("Ture","Prediction")
legend("bottomleft",leg.names,lty=1:2,col=c("Blue","Orange"),box.lty = 0)

lower=exp(forecast(model)$lower[,2]);upper=exp(forecast(model)$upper[,2])
ts.plot(c(dBT,fore));lines(c(dBT,lower),lty=2);lines(c(dBT,upper),lty=2)

sum((fore-testA)^2)/7
######

per = periodogram(lndBT)
require(data.table)
data.table(period=1/per$freq, spec=per$spec)[order(-spec)][1:5]

#Base Model
bestfit <- list(aic=model$aic, p=0, q=0, fit=model)
#Add fourior as regressor
for (i in 1:6){
  for (j in 1:4){
    #對於兩個超勢做FourierTransform，參數K為FourierTerms
    z1=fourier(ts(lndBT, frequency = 187), K=i)
    z2=fourier(ts(lndBT, frequency = 125), K=j)
    #重新建立新的ARIMA頻測模型
    fit=auto.arima(lndBT, xreg=cbind(z1, z2), seasonal=F)
    #如果新模型的AIC更小，就將p，q參數取代
    if (fit$aic<bestfit$aic)
      bestfit = list(aic=fit$aic, p=i, q=j, fit=fit)}}

bestfit

z1=fourier(ts(lndBT, frequency = 187), K=6)
z2=fourier(ts(lndBT, frequency = 125), K=4)

par(mfrow=c(3,1))
pred <- exp(predict(fit, newxreg = cbind(z1, z2))$pred[1:7])
ts.plot(cbind(testA,pred),lty=1:2,col=c("Blue","Orange"),main="BTC-USD with ARIMA(1,1,1) model")
leg.names<-c("Ture","Prediction")
legend("bottomleft",leg.names,lty=1:2,col=c("Blue","Orange"),box.lty = 0)

pred <- exp(predict(fit, newxreg = cbind(z1, z2))$pred)
ts.plot(c(dBT,pred))
lower=exp(forecast(fit, xreg = cbind(z1, z2))$lower[,2]);upper=exp(forecast(fit, xreg = cbind(z1, z2))$upper[,2])
lines(c(dBT,lower),lty=2);lines(c(dBT,upper),lty=2)
Box.test(fit$residuals)

###### Basic model
bestfit <- list(aic=model$aic, p=0, q=0, fit=model)
#Add fourier as regressor
for (i in 1:6){
  for (j in 1:4){
    #對於兩個超勢做Fourier Transform，參數K為Fourier Terms
    z1=fourier(ts(lndBT, frequency = 187), K=i)
    z2=fourier(ts(lndBT, frequency = 125), K=j)
    #重新建立新的ARIMA頻測模型
    fit2=Arima(lndBT, xreg=cbind(z1, z2), seasonal=F, include.drift = TRUE)
    #如果新模型的AIC更小，就將p，q參數取代
    if (fit2$aic<bestfit$aic)
      bestfit = list(aic=fit2$aic, p=i, q=j, fit=fit2)}}

bestfit
######
fit2 <- Arima(lndBT, seasonal=F)
pred2 <- exp(predict(fit2,n.ahead=7)$pred[1:7])
ts.plot(cbind(testA,pred2),lty=1:2,col=c("Blue","Orange"),main="BTC-USD with ARIMA(1,1,1) model")
leg.names<-c("Ture","Prediction")
legend("bottomleft",leg.names,lty=1:2,col=c("Blue","Orange"),box.lty = 0)

pred2 <- exp(predict(fit2,n.ahead=60)$pred[1:60])
ts.plot(c(dBT,pred2))
lower=exp(forecast(fit2)$lower[,2]);upper=exp(forecast(fit2)$upper[,2])
lines(c(dBT,lower),lty=2);lines(c(dBT,upper),lty=2)
Box.test(fit2$residuals)

rm(list = ls())
cat("/014")
######
dev.off()