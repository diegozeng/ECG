source('prepare_data.R')

rate = matrix(NA, ncol = 10, nrow = 11)
r1 = matrix(NA, ncol = 10, nrow = 11)
r2 = matrix(NA, ncol = 10, nrow = 11)
for(i in 1:10){
	select = sample(1:len,len/2)
	train = all[-select,]
	dim(train)
	table(train$V45)

	test = all[select,]
	dim(test)
	table(test$V45)

	
	for(j in 1:length(fs)){
		s = ksvm(fs[[j]],
			data = train,
			type = "C-bsvc", kernel = "rbfdot",
			kpar = list(sigma = 0.1), C = 10,
			prob.model = TRUE)

		p = predict(s,test)

		print(summary(p == test$V45))
		t=table(p,test$V45)
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














