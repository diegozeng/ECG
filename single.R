#load data
div <- div_data(newdata)
train = newdata[div,]
print(table(train$V45))
test = newdata[-div,]
print(table(test$V45))

train = newdata
test = newdata

library(randomForest)
s = randomForest(fs[[1]],data = train)
p = predict(s,test)
	
plot(test$V45,col=1,pch=2)
points(p,col=3,pch=20)


library(ada)
control <- rpart.control(cp = -1, maxdepth = 14,maxcompete = 1,xval = 0)
gen1 <- ada(V45~., data = train[c(-1,-2,-3)], test.x = test[,-45], test.y = test[,45], type = "gentle", control = control, iter = 70)

p <- predict(gen1,test)
confusionMatrix(p,test$V45)

summary(gen1)
varplot(gen1)

library(kernlab)
s = ksvm(fs[[1]],
         data = train,
         type = "C-bsvc", kernel = "rbfdot",
         kpar = list(sigma = 0.1), C = 10,
         prob.model = TRUE)

p = predict(s,newdata=test)
confusionMatrix(p,test$V45)

#lm
s = lm(f1,data = train)
p = predict(s,test)
pos = nrow(health)/nrow(newdata)
p[p>=pos] = 1
p[p<pos] = 0
confusionMatrix(p,test$V45)
plot(test$V45,col=1,pch=2)
points(p,col=3,pch=20)


library(RSofia)

s <- sofia(fs[[1]],data = newdata,loop_type="balanced-stochastic")
p = predict(s,newdata=test)
confusionMatrix(p,test$V45)


