library(DMwR)
source('div_data.func.R')
source('formulas.func.R')

fs <- formulas()

data = read.csv('result425.csv',header=F)
newdata = data[data$V2==1,]

newdata$V45 <- as.factor(newdata$V45)
tmp <- SMOTE(fs[[1]], newdata, perc.over = 300,perc.under = 900)
table(tmp$V45)

d <- sample(1:nrow(tmp),nrow(tmp)/2)

train <- tmp[d,]
test <- tmp[-d,]