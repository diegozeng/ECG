l = load('208.csv');

name = l(:,1);
von = l(:,2:5);
parfor i = 1:208
    tmp = num2str(name(i));
    if(exist(strcat(tmp,'.mat'),'file') == 2)
        j = l(i,6);
        [AC,DC,qtSDNN,radioSDNN,qtMEAN,radioMEAN] = QT(tmp,j,von(i,:),false);
        disp([tmp,',',num2str(j),',',num2str(AC),',',num2str(DC),',',num2str(qtSDNN),...
            ',',num2str(radioSDNN),',',num2str(qtMEAN),',',num2str(radioMEAN)]);
    end
end