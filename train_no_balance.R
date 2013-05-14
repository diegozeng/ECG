library(tree)
library(e1071)
library(kernlab)
library(randomForest)

data = read.csv('result425.csv',header=F)
newdata = data[data$V2==1,]

dead = newdata[newdata$V45 == 1,]
health = newdata[newdata$V45 == 0,]


f1 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26+V27+V28+V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+V43+V44
f2 <- V45 ~ V4+V5+V6+V7+V8+V9+V10#RR
f3 <- V45 ~ V11+V12+V13+V14+V15+V16+V17#QT
f4 <- V45 ~ V18+V19+V20+V21+V22+V23+V24#QT/RR
f5 <- V45 ~ V25+V26+V27+V28+V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+V43+V44#Reg
f6 <- V45 ~ V25+V26#Reg_lin
f7 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24#RR+QT+QT/RR
f8 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17#RR+QT
f9 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V18+V19+V20+V21+V22+V23+V24#RR+QT/RR
f10 <- V45 ~ V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24#QT+QT/RR
f11 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26#RR+QT+QT/RR+Reg_Lin

fs <- c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11)
rate = matrix(NA, ncol = 10, nrow = 11)
r1 = matrix(NA, ncol = 10, nrow = 11)
r2 = matrix(NA, ncol = 10, nrow = 11)
for(i in 1:10){
	deadlen = nrow(dead)
	halfdead = deadlen/2
	healthlen = nrow(health)
	halfhealth = healthlen /2

	d = sample(1:deadlen,halfdead,replace=T)
	dead1 = dead[d,]
	dead2 = dead[-d,]

	d = sample(1:healthlen,healthlen/2)
	train = rbind(health[d,],dead1)
	#dim(train)
	print(table(train$V45))

	test = rbind(health[-d,],dead2)
	#dim(test)
	print(table(test$V45))

	
	for(j in 1:length(fs)){
		s = ksvm(fs[[j]],
			data = train,
			type = "C-bsvc", kernel = "rbfdot",
			kpar = list(sigma = 0.1), C = 10,
			prob.model = TRUE)

		p = predict(s,test)

		print(summary(p == test$V45))
		t=table(p,test$V45)
		print(t)
		rate[j,i] = length(which(p == test$V45))/nrow(test)
		r1[j,i] = t[1]/(t[1]+t[2])
		r2[j,i] = t[4]/(t[3]+t[4])

		#plot(test$V45,col=1,pch=2)
		#points(p,col=3,pch=20)
	}
}
print(rate)
for(j in 1:11){
	print(mean(rate[j,]))
	print(mean(r1[j,]))
	print(mean(r2[j,]))
}





#single train
s = randomForest(f1,data = train)
p = predict(s,test)
p[p>=0.5] = 1
p[p<0.5] = 0
summary(p == test$V45)
table(p,test$V45)		
plot(test$V45,col=1,pch=2)
points(p,col=3,pch=20)



s = ksvm(f1,data = all,
type = "C-bsvc", kernel = "rbfdot",
kpar = list(sigma = 0.1), C = 10,
cross=5)









