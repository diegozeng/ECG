l = load('list.txt');

parfor i = 1:size(l)
    tmp = num2str(l(i));
    disp(['Get png of data : ',tmp]);
    
    data = load(strcat(tmp,'.mat'));
    data = struct2cell(data);
    data = cell2mat(data);
    set(figure(1),'visible','off');
    
    subplot(3,1,1);plot(data(10000:11000,1));title(strcat(tmp,' No.1 (10000-11000)'));
    subplot(3,1,2);plot(data(10000:11000,2));title(strcat(tmp,' No.2 (10000-11000)'));
    subplot(3,1,3);plot(data(10000:11000,3));title(strcat(tmp,' No.3 (10000-11000)'));
    print(gcf,'-dpng',strcat(tmp,'.png'))   
end