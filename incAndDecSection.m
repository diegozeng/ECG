%find section where RR inc or dec
function[inc,dec,incnum,decnum] = incAndDecSection(data,Rnum,R,mintime,begin,plot_outcome)
MAX_NUM = 1000;
inc = zeros(MAX_NUM,3);
incnum = 0;
dec = zeros(MAX_NUM,3);
decnum = 0;
pre = begin-1;
% mintime = 7000;
for i = begin:Rnum-begin
    if(data(i-1) < data(i) && data(i) > data(i+1))
        if(R(i,2) - R(pre,2) > mintime)
            decnum = decnum + 1;
            dec(decnum,1) = i - pre;
            dec(decnum,2) = pre;
            dec(decnum,3) = i;
        end
        pre = i;
    elseif(data(i-1) > data(i) && data(i) < data(i+1))
        if(R(i,2) - R(pre,2) > mintime)
            incnum = incnum + 1;
            inc(incnum,1) = i - pre;
            inc(incnum,2) = pre;
            inc(incnum,3) = i;
        end
        pre = i;
    end
end

if plot_outcome
    figure
    plot(data);
    hold on;
    for i = 1:decnum
        s = dec(i,2):dec(i,3);
        plot(s,data(s),'r');
    end
    hold on;
    for i = 1:incnum
        s = inc(i,2):inc(i,3);
        plot(s,data(s),'g');
    end
    hold on;

    [max_dec,didx] = max(dec(1:decnum,1));
     s = dec(didx,2):dec(didx,3);
     plot(s,data(s),'r','LineWidth',4);
     
    [max_inc,iidx] = max(inc(1:incnum,1));
    s = inc(iidx,2):inc(iidx,3);
     plot(s,data(s),'g','LineWidth',4);
    % plot(dec(1:decnum,2),data(dec(1:decnum,2)),'+r');
    % hold on;
    % plot(inc(1:incnum,2),data(inc(1:incnum,2)),'+g');
    title('Inc and Dec','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
    xlabel('Time','FontName','Times New Roman','FontSize',14);
    ylabel('RR','FontName','Times New Roman','FontSize',14,'Rotation',0);
end