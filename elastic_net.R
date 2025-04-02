library(ggplot2)
library(glmnet)
library(corrplot)
library(rlang)


selYs <- c("SNR",
          "SBN",
          "detection_SNR_greater_than_4",
          "detection_SNR_greater_than_5",
          "Detection_Ratios_SBN_greater_than_5")

allSampleID <- 1:19
trainSampleID <- c(1:4,6,8:19)
testSampleID <- c(5,7)

Ys <- read.table('data/Groups_SNR_2nd.txt', header=T, row.names = 1)
Xs <- read.table('data/All_Groups_FeaturesSummary_2nd.txt',header=T, row.names=2)[rownames(Ys),2:16]
# Xs$Dis_H3K27acP <- log10(Xs$Dis_H3K27acP)
# Xs$Guide_Num <- log10(Xs$Guide_Num)
corrplot(corr = cor(Xs))
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
Xtrain_scaled <- scale(Xtrain, center = Xtrain_mean, scale = Xtrain_sd)
Xtest_scaled <- scale(Xtest, center = Xtrain_mean, scale = Xtrain_sd)


an <- 0.5
set.seed(cvseed)
elastic_model <- cv.glmnet(Xtrain_scaled,ytrain,grouped=F,standardize=FALSE,alpha=an)
plot(elastic_model)
coef(elastic_model, s = "lambda.min")
ytrain_pred <- predict(elastic_model,Xtrain_scaled,s = "lambda.min")
ytest_pred <- predict(elastic_model,Xtest_scaled,s = "lambda.min")
print(cor(ytrain_pred,ytrain))
print(coef(elastic_model, s = "lambda.min"))
coef(elastic_model, s = "lambda.min")
plot(elastic_model)
plot(ytrain_pred,ytrain)
coef_lmin = coef(elastic_model, s = "lambda.min")




mmm <- as.numeric(ytest_pred)-as.numeric(ytest)
mse_test <- sum(mmm*mmm)

v <- coef(lm(y_true~y_pred,data = data.frame(y_pred=as.numeric(ytrain_pred), 
                                             y_true=ytrain)))


kk <- v[2]
bb <- v[1]
maxx0 <- max(ytrain_pred,ytest_pred)
minx0 <- min(ytrain_pred,ytest_pred)
maxy0 <- max(ytrain,ytest)
miny0 <- min(ytrain,ytest)

maxx <- max(maxx0, (maxy0-bb)/kk)
maxy <- max(maxy0, kk*maxx0+bb)
minx <- min(minx0, (miny0-bb)/kk)
miny <- min(miny0, kk*minx0+bb)

g0<-ggplot(data.frame(y_pred= c(ytrain_pred,ytest_pred), 
                      y_true=c(ytrain,ytest), 
                      type=c(rep('train' , length(ytrain_pred)),rep('test', length(ytest_pred))))) +
  geom_point(aes(x=y_pred,y=y_true,color=type)) + xlab(sprintf('Predicted %s',selY)) + ylab(selY) + 
  geom_abline(intercept = v[1], slope = v[2]) + 
  ggtitle(sprintf("p = %.2e, r = %.3f",
                  cor.test(ytrain,predict(elastic_model, newx = Xtrain_scaled, s = "lambda.min"))[3][[1]], 
                  cor(ytrain,predict(elastic_model, newx = Xtrain_scaled, s = "lambda.min")))) + 
  xlim(minx,maxx)+
  ylim(miny,maxy)+
  theme_bw(base_size=8) + theme(panel.grid = element_blank(), 
                                axis.text = element_text(size = 8),
                                axis.title = element_text(size = 8))


ggsave(sprintf('result/cor_pred_sel_%s.pdf',selY),width = 3, height = 2.3, units = "in")
write.csv(x = data.frame(y_pred= c(ytrain_pred,ytest_pred), 
                         y_true=c(ytrain,ytest), 
                         type=c(rep('train' , length(ytrain_pred)),rep('test', length(ytest_pred)))),
          file = sprintf('result/data_sel_%s.csv',selY)
)



feature_importance_v = c()
mean_sqr_all = mean((predict(elastic_model, newx = Xtrain_scaled, s = "lambda.min") - t(t(ytrain)))^2)
for (selected_feature in 
     names(coef_lmin[c(coef_lmin[,1])!=0,])[2:length(coef_lmin[c(coef_lmin[,1])!=0,])]){
  x_new <- duplicate(Xtrain_scaled, shallow = FALSE)
  x_new[,selected_feature] <- sample(x_new[,selected_feature])
  y_1shuffle = predict(elastic_model, newx = x_new, s = "lambda.min")
  mean_sqr_1shuffle = mean((y_1shuffle- t(t(ytrain)))^2)
  feature_importance = mean_sqr_1shuffle-mean_sqr_all
  feature_importance_v = c(feature_importance_v, feature_importance)
  message(selected_feature,": ",round(feature_importance,4))
}
names(feature_importance_v) <- names(coef_lmin[c(coef_lmin[,1])!=0,])[2:length(coef_lmin[c(coef_lmin[,1])!=0,])]




xx <- names(feature_importance_v)[order(feature_importance_v,decreasing=T)]
g1<-ggplot(data.frame(x=factor(xx,xx),
                      y = feature_importance_v[order(feature_importance_v,decreasing=T)])) + 
  ylab("Feature Importance\n(mean squared error") +
  xlab("Features") +
  geom_bar(stat="identity",aes(x=x,y=y)) + scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size =7)   + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

ggsave(sprintf('result/feature_importance_sel_%s.pdf',selY),width = 2, height = 1.8, units = "in")


print(feature_importance_v)
print(sum((ytest-ytest_pred)*(ytest-ytest_pred)))

