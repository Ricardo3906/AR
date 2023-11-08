rm(list = ls())  # clean up all variables 
library(forecast)
library(tseries)
library(ggplot2)
library(readxl)
library(lubridate)
library(dplyr)
library(caret)



###load data
setwd("C:/Users/15059/Desktop/")
data = read.csv("cn_cpi.csv")
colnames(data) = c("date","cpi")
data$date = as.Date(data$date)
adf.test(data$cpi, alternative = "stationary")  ###Augmented Dickey-Fuller Test (检验单位根)
plot(x=data$date, y=data$cpi,type = 'l')    ###画图检查平稳性
cpi = data$cpi

##trans to y2y: Integrate
cpi_y2y = rep(NA,length(cpi))
for(i in 13:length(cpi)){
  cpi_y2y[i] = (cpi[i]-cpi[i-12])/cpi[i-12]*100
}
cpi_y2y = data.frame(data$date,cpi_y2y)
adf.test(cpi_y2y$cpi_y2y[13:dim(cpi_y2y)[1]], alternative = "stationary")   ###再次检验
colnames(cpi_y2y) = c("date","cpi")
plot(x=cpi_y2y$date, y=cpi_y2y$cpi,type = 'l')    ###画图检查平稳性


###season adjust
cpi_season = ts(cpi_y2y$cpi,frequency = 12)     ###将每年同一个月的数据构成一个时间序列，总共12个时间序列
cpi_dc = decompose(cpi_season)                  ###分解
plot(cpi_dc)
cpi_ad = cpi_season - cpi_dc$seasonal           ###remove seasonality
cpi_after_season = data.frame(data$date,cpi_ad)
names(cpi_after_season) = c("date","cpi")
cpi_after_season = na.omit(cpi_after_season)

######1-period-ahead: AR(1)
cpi_h1 = lead(cpi_after_season$cpi, n=1)        ###构造前置1期的数据，即pi_t+1
cpih1_dataset = data.frame(cbind(cpi_after_season$date,cpi_h1,cpi_after_season$cpi))   
cpih1_dataset[,1] = cpi_after_season$date
colnames(cpih1_dataset) = c("date","cpi_h1","cpi")
##Modeling
x = cpih1_dataset[-dim(cpih1_dataset)[1],3]     
y = cpih1_dataset[-dim(cpih1_dataset)[1],2]
ar1_model = lm(y~x)                             ###OLS线性回归

##predict the latest period:
x_test = cpih1_dataset[dim(cpih1_dataset)[1],3] ###提取最新一期的CPI数据
pred = c(1,x_test) %*% ar1_model$coefficients   ###预测
print(pred)



######Multistep-forecast: direct forecast AR(n)  
lead_h = 6                                      ###设定forecast step
cpi_h6 = lead(cpi_after_season$cpi, n = lead_h) ###构造前置lead_h期的数据，即pi_t+h
##获取滞后的pi作为解释变量
pi = cpi_after_season$cpi                   
lag_pi = matrix(NA, nrow = length(pi), ncol = 12)
colnames(lag_pi) = paste("lag", seq(1, 12, 1), sep = "")
for(l in 1:ncol(lag_pi)){
  lag_pi[,l] = lag(pi, n = as.integer(l))       ###滞后，从1到12期
}
cpih6_dataset = as.matrix(cbind(cpi_h6, pi, lag_pi))
dataset_With_lag = cpih6_dataset
cpih6_dataset = as.matrix(cpih6_dataset[1:(length(cpih6_dataset[,1]) - lead_h),])   ###去掉有空缺值的行,构成训练集
###AIC
cand_model = list()   ###候选的model存储
cand_aic = list()
#Modeling
for(ii in 0:12){
  y=cpih6_dataset[,1]
  x=cpih6_dataset[,2:(2+ii)]                    ###对不同数量的滞后分别训练模型
  cand_model[[ii+1]] = lm(y~x)
  cand_aic[[ii+1]] = AIC(cand_model[[ii+1]])
}
which(unlist(cand_aic) == min(unlist(cand_aic))) ###寻找最小的AIC值
model = cand_model[[which(unlist(cand_aic) == min(unlist(cand_aic)))]]   
lagnum = model$rank - 2                         ###减去一个常数项，和pi本身
##predict
test_predictors = dataset_With_lag[(dim(dataset_With_lag)[1] - lead_h + 1),2:model$rank]  ###获取pi_n及其滞后项，预测
pred=t(c(1,test_predictors))%*%model$coefficients
print(pred)


cand_model = list()   ###候选的model存储
cand_bic = list()
#Modeling
for(ii in 0:12){
  y=cpih6_dataset[,1]
  x=cpih6_dataset[,2:(2+ii)]
  cand_model[[ii+1]] = lm(y~x)
  cand_bic[[ii+1]] = BIC(cand_model[[ii+1]])
}
which(unlist(cand_bic) == min(unlist(cand_bic)))
model = cand_model[[which(unlist(cand_bic) == min(unlist(cand_bic)))]]   ###寻找最小的BIC
lagnum = model$rank - 2 ###减去一个常数项，和pi本身
##predict
test_predictors = dataset_With_lag[(dim(dataset_With_lag)[1] - lead_h + 1),2:model$rank]
pred=t(c(1,test_predictors))%*%model$coefficients
print(pred)




######Multistep-forecast: recursive forecast AR(n) 
lead_h = 6                                               ###设定forecast step
cand_model = list()   ###候选的model存储
cand_aic = list()
cpi_h12 = lead(cpi_after_season$cpi, n=1)                 ###构造前置1期的数据，即pi_t+1
pi = cpi_after_season$cpi
lag_pi = matrix(NA, nrow = length(pi), ncol = 12)
colnames(lag_pi) = paste("lag", seq(1, 12, 1), sep = "")
for(l in 1:ncol(lag_pi)){
  lag_pi[,l] = lag(pi, n = as.integer(l))
}
cpih12_dataset = as.matrix(cbind(cpi_h12, pi, lag_pi))
dataset_With_lag1 = cpih12_dataset
cpih12_dataset = as.matrix(cpih12_dataset[1:(length(cpih12_dataset[,1]) - lead_h),])   ###去掉有空缺值的行,构成训练集
for(ii in 0:12){
  y=cpih12_dataset[,1]
  x=cpih12_dataset[,2:(2+ii)]
  cand_model[[ii+1]] = lm(y~x)
  cand_aic[[ii+1]] = AIC(cand_model[[ii+1]])
}
which(unlist(cand_aic) == min(unlist(cand_aic)))
model = cand_model[[which(unlist(cand_aic) == min(unlist(cand_aic)))]]   ###寻找最小的AIC
lagnum = model$rank - 2 ###减去一个常数项，和pi本身
print(lagnum)
##predict
test_predictors = dataset_With_lag[(dim(dataset_With_lag)[1] - lead_h + 1),2:model$rank] ###获取pi_n及其滞后项，预测
pred=t(c(1,test_predictors))%*%model$coefficients   ### 得到^_pi_n+1

pi_for_pred = pi[1:(length(pi) - lead_h)]           ###保持pi_1 到 pi_n-1 不变，
for(i in 2:lead_h){
  pi_for_pred = append(pi_for_pred,pred)
  lag_pi_for_pred = matrix(NA, nrow = length(pi_for_pred), ncol = lagnum)
  colnames(lag_pi_for_pred) = paste("lag", seq(1, lagnum, 1), sep = "")
  for(l in 1:ncol(lag_pi_for_pred)){
    lag_pi_for_pred[,l] = lag(pi_for_pred, n = as.integer(l))
  }
  test_dataset = as.matrix(cbind(pi_for_pred, lag_pi_for_pred))
  test_predictors = test_dataset[dim(lag_pi_for_pred)[1],]
  pred=t(c(1,test_predictors))%*%model$coefficients   ###^_pi_n+i (i>=2)
}
print(pred)




