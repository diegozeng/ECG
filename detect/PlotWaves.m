function [] = PlotWaves( RR,RRm,RRew,QT,QTm,QTew,s )
    figure
    plot(s,RR(s),'y');
    hold on;
    plot(s,RRm(s),'g');
    hold on;
    plot(s,RRew(s),'r');
    hold on;

    title('RR Waves','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
    xlabel('Time in cell','FontName','Times New Roman','FontSize',14);
    ylabel('Wave','FontName','Times New Roman','FontSize',14,'Rotation',0);

    figure
    plot(s,QT(s),'y');
    hold on;
    plot(s,QTm(s),'g');
    hold on;
    plot(s,QTew(s),'r');
    hold on;

    title('QT Waves','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
    xlabel('Time in cell','FontName','Times New Roman','FontSize',14);
    ylabel('Wave','FontName','Times New Roman','FontSize',14,'Rotation',0);
end

