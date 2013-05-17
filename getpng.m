l = load('all_mat_list.txt');

parfor i = 1:425
    tmp = num2str(l(i));
    disp(['Get png of data : ',tmp]);
    
    data = load(strcat(tmp,'.mat'));
    data = struct2cell(data);
    data = cell2mat(data);
	
    set(gcf,'visible','off');
    subplot(3,1,1);plot(data(1:1000,1));title(strcat(tmp,' No.1 (1-1000)'));
    subplot(3,1,2);plot(data(1:1000,2));title(strcat(tmp,' No.2 (1-1000)'));
    subplot(3,1,3);plot(data(1:1000,3));title(strcat(tmp,' No.3 (1-1000)'));
    print(gcf,'-dpng',strcat(tmp,'_1.png'));
	
	set(gcf,'visible','off');  
    subplot(3,1,1);plot(data(10000:11000,1));title(strcat(tmp,' No.1 (10000-11000)'));
    subplot(3,1,2);plot(data(10000:11000,2));title(strcat(tmp,' No.2 (10000-11000)'));
    subplot(3,1,3);plot(data(10000:11000,3));title(strcat(tmp,' No.3 (10000-11000)'));
    print(gcf,'-dpng',strcat(tmp,'_2.png'));
	
	set(gcf,'visible','off');
	subplot(3,1,1);plot(data(20000:21000,1));title(strcat(tmp,' No.1 (20000-21000)'));
    subplot(3,1,2);plot(data(20000:21000,2));title(strcat(tmp,' No.2 (20000-21000)'));
    subplot(3,1,3);plot(data(20000:21000,3));title(strcat(tmp,' No.3 (20000-21000)'));
    print(gcf,'-dpng',strcat(tmp,'_3.png'));
end