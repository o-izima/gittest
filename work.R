library(glmnet)
library(caret)
library(dplyr)
library(ggplot2)
library(leaps)
library(MASS)
library(data.table)

#Read in files
X = read.delim("periodicx_nc.dat", header = TRUE, sep = ",")
Y = read.delim("periodicy_nc.dat", header = TRUE, sep = ",")
TCP = read.delim("periodicload.dat", header = TRUE, sep = ",")
X0_cnsm = read.csv("X0_cnsm2015.csv", header = TRUE)
Y_cnsm = read.csv("Y_cnsm2015.csv", header = TRUE)
#Y_cnsm = Y_cnsm[,-1]

########## Histogram of TCPSCK
qplot(tcp, geom = "histogram")

ggplot(data=tcp, aes(tcp$X19)) + 
  geom_histogram(breaks=seq(19, 110, by=1), 
                 col="red", 
                 fill="green", 
                 alpha = .2) +
  ylim(c(0,2000))

### Combine x and y; rename the response variable
XY = cbind(Y,X)
names(XY)[1]<-"dispframe"

head(XY[, 1:5])
#head(XY, n = 10)
table(is.na(XY))

#### lets filter based on tcp values between 24 & 72 where counts are more than 500
tcp = filter(XY, X19 >= 24 & X19 <= 72)
t30 = filter(XY, X19 == 30)
t60 = filter(XY, X19 == 60)
t72 = filter(XY, X19 == 72)

min(tcp$X19)
max(tcp$X19)

attach(XY)

############ Data partition

set.seed(102)
id = sample(2, nrow(XY), replace = T, prob = c(0.6, 0.4))
train = XY[id == 1,]
test = XY[id == 2,]

x_train = model.matrix(dispframe~. -1, data = train)
x_test = model.matrix(dispframe~.-1, data = test)

y_train = train %>%
  select(dispframe) %>%
  unlist() %>%
  as.numeric()

y_test = test %>%
  select(dispframe) %>%
  unlist() %>%
  as.numeric()

######Grid values for lambda
grid = 10^seq(10, -2, length = 100)

#save files to PC
write.table(test,"test", sep = " ",row.names = FALSE)
b = read.table("X_train", header = TRUE, sep = " ")

################## Unadjusted Learning ##############################
################### Ridge
ridge_mod = glmnet(x_train, y_train, alpha = 0, lambda = grid)
par(mfrow = c(1,2))
plot(ridge_mod,xlab="L2 Norm", ylab="Standardized Coefficients")
plot(ridge_mod, xvar = "lambda", ylab="Standardized Coefficients", label = TRUE)

#Use CV to select a good lambda value using cv.glmnet
ridge_mod_cv = cv.glmnet(x_train, y_train, alpha = 0)
par(mfrow = c(1,1))
plot(ridge_mod_cv)
ridge_mod_cv$lambda.min
ridge_mod_cv$lambda.1se

#prediction
ridge_pred = predict(ridge_mod_cv, x_test)
rmse = mean((y_test - ridge_pred)^2)
rmse # 4004.043

plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(ridge_pred, type = "l", col = "red")
#######################
##################Lasso
lasso_mod = glmnet(x_train, y_train, alpha = 1, lambda = grid)
par(mfrow = c(1,2))
plot(lasso_mod, ylab="Standardized Coefficients")
plot(lasso_mod, xvar = "lambda", ylab="Standardized Coefficients", label = TRUE)

#Use CV to select a good lambda value using cv.glmnet
lasso_mod_cv =  cv.glmnet(x_train, y_train, alpha = 1)
par(mfrow = c(1,1))
plot(lasso_mod_cv)
lasso_mod_cv$lambda.min
lasso_mod_cv$lambda.1se

#prediction
lasso_pred = predict(lasso_mod_cv, x_test)
rmse = mean((y_test - lasso_pred)^2)
rmse # 774.8073

plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(lasso_pred, type = "l", col = "red")
########################
####################Elastic Net using caret


en = train(
  dispframe ~., data = train,
  method = "glmnet",
  trControl = en_cv
)
en
plot(en$finalModel, xvar = "lambda", ylab="Standardized Coefficients", label = T)
plot(en$finalModel, xvar = "dev", label = T)
en$bestTune
best = en$finalModel
coef(best, s = en$bestTune$lambda)
##################################################################
####### Unadjusted learning using EN for dispframes ##############
#custome control
en_cv = trainControl(method = "cv", number = 10)

en2 = train(dispframe ~ .,
           train,
           method = 'glmnet',
           tuneGrid = expand.grid(alpha = 0:1,
                                  lambda = 10^seq(10, -2, length = 100)),
           trControl = en_cv)

en2$bestTune
#  alpha    lambda
#113     1 0.2848036
#####################################################
############# Predict on the EN Model - Dispframes
#RMSE for train
p1 = predict(en2, train)
rmse_train_UA = sqrt(mean((p1-train$dispframe)^2))
rmse_train_UA   #28.30061

#RMSE for test
p2 = predict(en2, test)
rmse_test_UA = sqrt(mean((p2-test$dispframe)^2))
rmse_test_UA  #27.6736

par(mfrow = c(1,1))
plot(test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(p2, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)
str(p2)


# Finish the data.frame() call
my_solution <- data.frame(rtp = XY$dispframe, tcpsck = XY$X19)

# Use nrow() on my_solution
nrow(my_solution)
ncol(my_solution)
# Finish the write.csv() call
write.csv(my_solution, file = "rtp_tcpsck.csv", row.names = FALSE)

#predicted and real in a df
pred_real = data.frame(real_y = test$dispframe, pred_y = p2)
ncol(pred_real)
nrow(pred_real)
# Finish the write.csv() call
write.csv(pred_real, file = "real_pred.csv", row.names = FALSE)







###############################################
#Train and Test predictions
plot(errors, type = "l", lty = 1.8, col = "blue")
lines(p2, type = "l", col = "red")
legend("topright",legend=c("Train RMSE", "Test RMSE"), col=c("blue", "red"),lty=1.2,cex=0.6)

errors = test$dispframe-p2

en_t60 = train(
  dispframe ~., data = t60,
  method = "glmnet",
  trControl = en_cv
)
get_best_result(en_t60)

### prediction UA
en2_pred = predict(en2, test)
plot(test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(en2_pred, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)

############## Compare Models #################################
model_list = list(EN = en, EN2 = en2, la = lasso_mod)
res = resamples(model_list)
summary(res)
xyplot(res, metric = 'RMSE')


en2
plot(en)
plot(en2)
par(mfrow = c(1,2))
plot(en2$finalModel, xvar = "lambda", ylab="Standardized Coefficients", label = T)
plot(en2$finalModel, xvar = "dev", ylab="Standardized Coefficients",label = T)


get_best_result = function(caret_fit){
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}
get_best_result(en2)
#alpha    lambda     RMSE  Rsquared      MAE    RMSESD RsquaredSD
#1     1 0.2848036 28.43464 0.8036022 19.55824 0.8818206 0.01133255
#MAESD
#1 0.3637756

################################################################
##################### Linear Model ########################
lm_model = train(dispframe~.,
                 train,
                 method = 'lm',
                 trControl = en_cv)
summary(lm_model)

###################### Predictions on xy_train ##############
train$pred = predict(lm_model)

rmse_train <- RMSE(train$pred, train$dispframe)

rmse_train #25.42566 / 25.59475

###################### Predictions on xy_test ##############
test$pred = predict(lm_model, newdata = x_test)

rmse_test <- RMSE(test$pred, test$dispframe)
rmse_test  #2727107416






#######################################################################
##################################Load Adjusted
############ Data partition
set.seed(103)
id = sample(2, nrow(tcp), replace = T, prob = c(0.6, 0.4))
la_train = tcp[id == 1,]
la_test = tcp[id == 2,]

x_la_train = model.matrix(dispframe~. -1, data = la_train)
x_la_test = model.matrix(dispframe~.-1, data = la_test)

y_la_train = la_train %>%
  select(dispframe) %>%
  unlist() %>%
  as.numeric()

y_la_test = la_test %>%
  select(dispframe) %>%
  unlist() %>%
  as.numeric()
####################Ridge
#### LA Ridge
la_ridge_mod = glmnet(x_la_train, y_la_train, alpha = 0, lambda = grid)
par(mfrow = c(1,2))
plot(la_ridge_mod,xlab="L2 Norm",ylab="Standardized Coefficients")
plot(la_ridge_mod, xvar = "lambda",ylab="Standardized Coefficients", label = TRUE)

#Use CV to select a good lambda value using cv.glmnet
la_ridge_mod_cv = cv.glmnet(x_la_train, y_la_train, alpha = 0)
par(mfrow = c(1,1))
plot(la_ridge_mod_cv)

#prediction
la_ridge_pred = predict(la_ridge_mod_cv, x_la_test)
rmse = mean((la_ridge_pred - y_la_test)^2)
rmse # 4193.526
sqrt(la_ridge_mod_cv$cvm[la_ridge_mod_cv$lambda == la_ridge_mod_cv$lambda.min])

plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(la_ridge_pred, type = "l", col = "red")

################LA Lasso
la_lasso_mod = glmnet(x_la_train, y_la_train, alpha = 1, lambda = grid)
par(mfrow = c(1,2))
plot(la_lasso_mod,ylab="Standardized Coefficients")
plot(la_lasso_mod, xvar = "lambda", ylab="Standardized Coefficients",label = TRUE)

#Use CV to select a good lambda value using cv.glmnet
la_lasso_mod_cv =  cv.glmnet(x_la_train, y_la_train, alpha = 1)
par(mfrow = c(1,1))
plot(la_lasso_mod_cv)
sqrt(la_lasso_mod_cv$cvm[la_lasso_mod_cv$lambda == la_lasso_mod_cv$lambda.min])

#prediction
la_lasso_pred = predict(la_lasso_mod_cv, x_la_test)
rmse = mean((y_la_test - la_lasso_pred)^2)
rmse # 830.4624

par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(la_lasso_pred, type = "l", col = "red")

### EN la
la_en = train(dispframe ~ .,
            la_train,
            method = 'glmnet',
            tuneGrid = expand.grid(alpha = 0:1,
                                   lambda = 10^seq(10, -2, length = 100)),
            trControl = en_cv)
plot(la_en)
par(mfrow = c(1,2))
plot(la_en$finalModel, xvar = "lambda", ylab="Standardized Coefficients",label = T)
plot(la_en$finalModel, xvar = "dev", ylab="Standardized Coefficients",label = T)

la_en_pred = predict(la_en, la_test)
rmse_en = mean((la_en_pred - la_test$dispframe)^2)
rmse_en #815.7882

par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue")
lines(la_en_pred, type = "l", col = "red")

get_best_result(la_en)

#alpha    lambda     RMSE  Rsquared    MAE    RMSESD RsquaredSD     MAESD
#1     1 0.2848036 28.46236 0.8045469 19.599 0.7961651 0.01003719 0.2394656

####################################################################
##################### Linear Model LA ########################
la_lm = train(dispframe~.,
                 la_train,
                 method = 'lm',
                 trControl = en_cv)
summary(la_lm)
summary(lm_model)
###################### Predictions on xy_train ##############
train$pred = predict(lm_model)

rmse_train <- RMSE(train$pred, train$dispframe)

rmse_train #25.42566 / 25.59475

###################### Predictions on xy_test ##############
test$pred = predict(lm_model, x_test)

rmse_test <- RMSE(test$pred, test$dispframe)
rmse_test  #2727107416

####
############## Compare Models #################################
model_list = list(UA_LinearModel = lm_model, UA_Ridge = ridge_mod, UA_Lasso = lasso_mod)
res = resamples(model_list)
summary(res)
xyplot(res, metric = 'RMSE')




##################### Linear Model ########################
lm = train(dispframe~.,
           la_train,
           method = 'lm',
           trControl = en_cv)
lm$results
lm
summary(lm)
##################### Ridge Regression ####################

ridge = train(dispframe ~ .,
              la_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 0,
                                     lambda = 10^seq(10, -2, length = 100)),
              trControl = en_cv)
plot(ridge)
ridge
plot(ridge$finalModel, xvar = "lambda", label = T)
plot(ridge$finalModel, xvar = "dev", label = T)
plot(varImp(ridge, scale = F))

################### Lasso Regression #######################

lasso = train(dispframe ~ .,
              la_train,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 1,
                                     lambda = 10^seq(10, -2, length = 100)),
              trControl = en_cv)
lasso
plot(lasso$finalModel, xvar = "lambda", label = T)
plot(lasso$finalModel, xvar = "dev", label = T)
varImp(lasso)

model_list = list(LA_LinearModel = lm, LA_Ridge = ridge, LA_Lasso = lasso, LA_EN = la_en)
res = resamples(model_list)
summary(res)
xyplot(res, metric = 'RMSE')

##### plots of predictions against real / actual values
#LA-EN
par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(la_en_pred, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)

#LA-LASSO
par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(la_lasso_pred, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)
#LA-RIDGE
par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(la_ridge_pred, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)

#LA-Linear Model
par(mfrow = c(1,1))
plot(la_test$dispframe, type = "l", lty = 1.8, col = "blue", ylab = "Signal Amplitude", xlab="Time(seconds)")
lines(test$pred, type = "l", col = "red")
legend("topright",legend=c("Actual RTP values", "Predicted RTP values"), col=c("blue", "red"),lty=1.2,cex=0.6)

#################LA using load values for 30, 40, 50, 60, 72
en72 = train(dispframe ~ .,
              t72,
              method = 'glmnet',
              tuneGrid = expand.grid(alpha = 0:1,
                                     lambda = 10^seq(10, -2, length = 100)),
              trControl = en_cv)
get_best_result(en72)

all30 = filter(XY, X19 == 30)


tcp30 = filter(tcp, X19 == 30)
####subset selection
num_vars = ncol(t30) - 1

fit_fwd= regsubsets(dispframe ~ ., data = t30, nvmax = num_vars, method = "backward")
fit_fwd_sum = summary(fit_fwd)
names(fit_fwd_sum)
which.min(fit_fwd_sum$cp)
coef(fit_fwd, which.max(fit_fwd_sum$adjr2))

set.seed(42)
train.control = trainControl(method = "cv", number = 10)
step_all = train(dispframe ~ ., data = train,
                   method = "leapBackward",
                   tuneGrid = data.frame(nvmax = 20),
                   trControl = train.control)
step.model$results
summary(step.model.30)
##################################
t50 = filter(tcp, X19 == 50)

##### EN using LA for the load values 30,40,50,60,72
en72 = train(dispframe ~ .,
             t72,
             method = 'glmnet',
             tuneGrid = expand.grid(alpha = 0:1,
                                    lambda = 10^seq(10, -2, length = 100)),
             trControl = en_cv)

model_LA = list(EN30 = en30, EN40 = en40, EN50 = en50, EN60 = en60, EN72 = en72)
res_LA = resamples(model_LA)
summary(res_LA)
xyplot(res, metric = 'RMSE')
######################################
set.seed(123)
train.control = trainControl(method = "cv", number = 10)
en_la_ss_72 = train(dispframe ~ ., data = t72,
                 method = "leapBackward",
                 tuneGrid = data.frame(nvmax = 30),
                 trControl = train.control)

model_LA_SS = list(EN_30 = en_la_ss_30, EN_40 = en_la_ss_40, EN_50 = en_la_ss_50, EN_60 = en_la_ss_60, EN_72 = en_la_ss_72)
res_LA_SS = resamples(model_LA_SS)
summary(res_LA_SS)
################################-----------------#################################
##################################################
####### Unadjusted learning using EN for dispframes ##############
#custome control
en_cv = trainControl(method = "cv", number = 10)

en2 = train(dispframe ~ .,
            train,
            method = 'glmnet',
            tuneGrid = expand.grid(alpha = 0:1,
                                   lambda = 10^seq(10, -2, length = 100)),
            trControl = en_cv)

en2$bestTune
#  alpha    lambda
#113     1 0.2848036
#####################################################
############# Predict on the EN Model - Dispframes
#RMSE for train
p1 = predict(en2, train)
rmse_train_UA = sqrt(mean((p1-train$dispframe)^2))
rmse_train_UA   #28.30061

#RMSE for test
p2 = predict(en2, test)
rmse_test_UA = sqrt(mean((p2-test$dispframe)^2))
rmse_test_UA  #27.6736
###############################################
plot(train$dispframe, p1)
plot(test$dispframe, p2)
###############################################################
######## Audio Buffer Rate & Net Operations Read ##############
#combine all 
XY_cnsm = cbind(Y_cnsm[,c(2,4)], X0_cnsm)

table(is.na(XY_cnsm))
XY_cnsm = na.omit(XY_cnsm)
############ Data partition
set.seed(103)
id2 = sample(2, nrow(XY_cnsm), replace = T, prob = c(0.6, 0.4))
train_cnsm = XY_cnsm[id2 == 1,]
test_cnsm = XY_cnsm[id2 == 2,]

##############################################################
####### Unadjusted learning using EN for Audio Buffer Rate####
en_audio = train(noAudioPlayed ~ .,
            train_cnsm[,-2],
            method = 'glmnet',
            tuneGrid = expand.grid(alpha = 0:1,
                                   lambda = 10^seq(10, -2, length = 100)),
            trControl = en_cv)

en_audio$bestTune
#alpha     lambda
# 1       0.01747528

#####################################################
############# Predict on the EN Model - Audio
#RMSE for train (Audio Buffer Rate)
pred_audio = predict(en_audio, train_cnsm)
rmse_train_UA_audio = sqrt(mean((pred_audio-train_cnsm$noAudioPlayed)^2))
rmse_train_UA_audio   #19.61647

#RMSE for test
pred_audio_test = predict(en_audio, test_cnsm)
rmse_test_UA_audio = sqrt(mean((pred_audio_test-test_cnsm$noAudioPlayed)^2))
rmse_test_UA_audio #19.40701
#######################################################
get_best_result(en_audio)
#alpha     lambda     RMSE   Rsquared      MAE    RMSESD RsquaredSD
# 1       0.01747528 19.71877 0.08758169 14.31381 0.3798021 0.01550241
#MAESD
#1 0.2163152

##############################################################
####### Unadjusted learning using EN for Net Read Operations####
en_netreadOp = train(NetReadOperations ~ .,
                 train_cnsm[,-1],
                 method = 'glmnet',
                 tuneGrid = expand.grid(alpha = 0:1,
                                        lambda = 10^seq(10, -2, length = 100)),
                 trControl = en_cv)

en_netreadOp$bestTune
#alpha    lambda
#1         0.1629751

#####################################################
############# Predict on the EN Model - Net Read Operations
#RMSE for train (Audio Buffer Rate)
pred_netread = predict(en_netreadOp, train_cnsm)
rmse_train_UA_netread = sqrt(mean((pred_netread-train_cnsm$NetReadOperations)^2))
rmse_train_UA_netread  #154.9323

#RMSE for test
pred_netread_test = predict(en_netreadOp, test_cnsm)
rmse_test_UA_netread = sqrt(mean((pred_netread_test-test_cnsm$NetReadOperations)^2))
rmse_test_UA_netread #149.0355
###############################################
#Compare the EN models for Dispframes, Audio & NetReadOperations
model_list = list(EN_Dispframes = en2, EN_Audio = en_audio, EN_NetReadOp = en_netreadOp)
res = resamples(model_list)
summary(res)
xyplot(res, metric = 'RMSE')

XY[["X19"]]
search()

?gzfile
zz=gzfile('realm-im2015-vod-traces','rt')
dat=read.csv(zz,header=F)

a=fread('realm-im2015-vod-traces')
dt <- fread(input = 'realm-im2015-vod-traces.gz')
??fread             
install.packages("data.table")
