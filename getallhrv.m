l = load('ECG_change.txt');

name = l(:,1);
von = l(:,2:5);

parfor i = 83:185
    for j=1:3
        out = HRV(num2str(name(i)),von(i,:),j);
        
        tmp = strcat(num2str(out.name),',',num2str(out.MEAN),',',num2str(out.SDNN),',',num2str(out.SDANN),',',num2str(out.ASDNN ),',',num2str( out.DC ), ',',num2str( out.AC ), ',', num2str( out.TI ));
        
        disp(tmp);
    end
end