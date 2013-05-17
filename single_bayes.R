library(bnlearn)
library(caret)
data = read.csv('result425.csv',header=F)
data = data[data$V2==1,]


fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,returnResamp = "all")
model1 <- train(V45~V4+V5+V6+V7+V8+V9+V10, data=data[c(-1,-2,-3)],method='nb',trControl = fitControl,
                tuneGrid = data.frame(.fL=1,.usekernel=F))

resampleHist(model1)

pre <- predict(model1)
confusionMatrix(pre,data$Class)





library(bnlearn)
data2 <- discretize(data[c(-1,-2,-3,-45)],method='quantile')
data2$class <- as.factor(data[,45])
bayesnet <- hc(data2)
plot(bayesnet)
fitted <- bn.fit(bayesnet, data2,method='mle')
pre <- predict(fitted,data=data2,node='class')
confusionMatrix(pre,data2$class)