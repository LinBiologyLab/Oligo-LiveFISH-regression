library(ggplot2)
library(reshape)
library(ROCR)
library(pROC)
library(rlang)
library(glmnet)
library(corrplot)





cv_test <- function(ss){
  print(ss)
  set.seed(ss)
allSampleID <- 1:19

testSampleID <- sample(allSampleID,2)
trainSampleID <- setdiff(allSampleID,testSampleID)
print(testSampleID)


Ys <- read.table('data/Groups_SNR_2nd.txt', header=T, row.names = 1)
Xs <- read.table('data/All_Groups_FeaturesSummary_2nd.txt',header=T, row.names=2)[rownames(Ys),2:16]

sel_feature <- c(
  "Guide_Num",
  "GC_Content",
  "Length",
  "Space_2nd",
  "crRNA_2nd",
  "Dis_H3K27acP",
  "Dis_ATAC_P",
  "H3K27ac",
  "ATAC",
  "H3K18ac",
  "H3K56ac",
  "H3K9ac",
  "H3K9me3",
  "H3K4me3", 
  "H3K4me1_siNC"
)

Xtrain <- Xs[trainSampleID,sel_feature]
ytrain <- Ys[trainSampleID,selY]

Xtest <- Xs[testSampleID,sel_feature]
ytest <- Ys[testSampleID,selY]

Xtrain_mean <- apply(Xtrain, 2, mean)
Xtrain_sd <- apply(Xtrain, 2, sd)
Xtrain_scaled <- scale(Xtrain, center = Xtrain_mean, scale = Xtrain_sd)
Xtest_scaled <- scale(Xtest, center = Xtrain_mean, scale = Xtrain_sd)



for(i in 1:1000){
 
an <- 0.1
set.seed(i)
elastic_model <- cv.glmnet(Xtrain_scaled,ytrain,grouped=F,standardize=FALSE,alpha=an)
coef(elastic_model, s = "lambda.min")
ytrain_pred <- predict(elastic_model,Xtrain_scaled,s = "lambda.min")
ytest_pred <- predict(elastic_model,Xtest_scaled,s = "lambda.min")
print(cor(ytrain_pred,ytrain))
coef_lmin = coef(elastic_model, s = "lambda.min")

if(!is.na(cor(ytrain_pred,ytrain))){
  break
}
print(coef(elastic_model, s = "lambda.min"))
coef(elastic_model, s = "lambda.min")
}


allresult <- sum((ytest-ytest_pred)*(ytest-ytest_pred))



Ys <- read.table('data/Groups_SNR_2nd.txt', header=T, row.names = 1)
Xs <- read.table('data/All_Groups_FeaturesSummary_2nd.txt',header=T, row.names=2)[rownames(Ys),2:16]
sel_feature <- c(
  "Guide_Num",
  "GC_Content",
  "Length",
  "Space_2nd",
  "Dis_H3K27acP",
  "ATAC", 
  "H3K4me3", 
  "H3K4me1_siNC"
)

Xtrain <- Xs[trainSampleID,sel_feature]
ytrain <- Ys[trainSampleID,selY]

Xtest <- Xs[testSampleID,sel_feature]
ytest <- Ys[testSampleID,selY]

Xtrain_mean <- apply(Xtrain, 2, mean)
Xtrain_sd <- apply(Xtrain, 2, sd)
#Xtrain_scaled <- apply(Xtrain,2,scale)
Xtrain_scaled <- scale(Xtrain, center = Xtrain_mean, scale = Xtrain_sd)
#Xtrain_scaled <- apply(Xtrain,2,scale), center = Xtrain_mean, scale = Xtrain_sd)
Xtest_scaled <- scale(Xtest, center = Xtrain_mean, scale = Xtrain_sd)

#pdf("giscore_lasso_linear_model_training_lambda_test.pdf")
# m <-lapply(seq(0,1,0.05), function(an){
#   set.seed(900)
# elastic_model <- cv.glmnet(Xtrain_scaled[,1:8],ytrain,grouped=F,standardize=FALSE,alpha=an)
# y_pred <- predict(elastic_model,Xtrain_scaled[,1:8],s = "lambda.min")
# print(i);
# print(cor(y_pred,ytrain))
# print(coef(elastic_model, s = "lambda.min"))
# coef(elastic_model, s = "lambda.min")
# return(cor(y_pred,ytrain))
# })
# dev.off()

for(i in 1:1000){
  
  an <- 0.5
  #cvseed <- 3728
  #cvseed <- sample(1:10000,1)
  set.seed(i)
  elastic_model <- cv.glmnet(Xtrain_scaled,ytrain,grouped=F,standardize=FALSE,alpha=an)
  #plot(elastic_model)
  coef(elastic_model, s = "lambda.min")
  ytrain_pred <- predict(elastic_model,Xtrain_scaled,s = "lambda.min")
  ytest_pred <- predict(elastic_model,Xtest_scaled,s = "lambda.min")
  print(cor(ytrain_pred,ytrain))
  coef_lmin = coef(elastic_model, s = "lambda.min")
  
  if(!is.na(cor(ytrain_pred,ytrain))){
    break
  }
  print(coef(elastic_model, s = "lambda.min"))
  coef(elastic_model, s = "lambda.min")
  #plot(elastic_model)
  #plot(ytrain_pred,ytrain)
}
selresult  <- sum((ytest-ytest_pred)*(ytest-ytest_pred))

return(data.frame(allrs =allresult, sel=selresult))

}



library(parallel)
allrs <- c()
selrs <- c()
for(i in 1:5){
  selY <- c("SNR",
            "SBN",
            "detection_SNR_greater_than_4",
            "detection_SNR_greater_than_5",
            "Detection_Ratios_SBN_greater_than_5")[i]
  rs <- mclapply(1:100,FUN = cv_test,mc.cores = 10)
  df <- do.call(rbind,rs)
  allrs <- c(allrs,mean(df$allrs))
  selrs <- c(selrs, mean(df$sel))
}

write.csv(x = data.frame(allvar=allrs,selvar=selrs),file = 'result/CV_result.csv')


