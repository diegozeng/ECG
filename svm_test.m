data = load('result425.csv');

DATA_SERIAL_NUM = 2;
LABEL = 45;
FEATURES = 4:44;
LABEL_HEALTH = 0;
LABEL_DEAD = 1;


newdata = data(data(:,DATA_SERIAL_NUM) == 1,:);
health = newdata(newdata(:,LABEL == LABEL_HEALTH),:);
dead = newdata(newdata(:,LABEL == LABEL_DEAD),:);

p = cvpartition(newdata(:,45),'Holdout',0.5);
% tabulate(newdata(p.training,45))
% tabulate(newdata(p.test,45))

model = svmtrain(newdata(p.training,LABEL),newdata(p.training,FEATURES),'-s 0 -w1 100 -w0 0.5');

[predict_label, accuracy,d] = svmpredict(newdata(p.test,LABEL), newdata(p.test,FEATURES), model);

accuracy

confusionmat(newdata(p.test,LABEL),predict_label)




