function [ R,Rnum ] = RWavesDetect( Sig,cell )
    maxsize= size(Sig,1);
    mindistant = 40;        %对R点的搜索最近距离
    maxdistant = 5000/cell; %对R点的搜索最远距离
    uprange = 3;            %幅值过滤上界
    downrange = 0.5;        %幅值过滤下界
    dis = 0.5;              %斜率过滤指标
    tstep = 3;              %斜率步长
    kran = 0.6;             %初始点阈值

    S = zeros(maxsize,1);%slope

    maxrnum = 150000;
    R = zeros(maxrnum,3);%R wave

    %-----------------------------确定第一个R峰位置------------------------
    S(1+tstep:maxsize-tstep) = ((2*Sig(1+tstep:maxsize-tstep)-Sig(1:maxsize-2*tstep) - Sig(1+2*tstep:maxsize))/2);
    k = max(S(1+tstep:1000));
    for i = 1+tstep:4000
        if S(i) > k*kran,
            break;
        end
    end
    m = -inf;
    if i <= 7
        i = 7;
    end
    k = i - 5;
    for j = i -5:i+5
        if m < Sig(j) && Sig(j) > Sig(j-1) && Sig(j) >= Sig(j+1),
            m = Sig(j);
            k = j;
        end
    end
    start = k;

    t = S(start);%第一个RR波峰的M值
    m = 1;       %用于统计R点和RR周期的关系
    i = start;   %从第一个R波波峰开始扫描
    n = t*tstep; %第一个R波波峰的高度
    while i+maxdistant < maxsize                             %找出R点并算出RR周期存入A
        if S(i)/t < dis||Sig(i) < Sig(i-1)||Sig(i) < Sig(i+1)||S(i)*tstep < n*downrange||S(i)*tstep > n*uprange,%除去伪R波
            [i] = nextR(i, maxdistant, S);
            continue;
        end
        tar = 0;
        for j = i+mindistant:i+maxdistant
            if S(j)*tstep > n*uprange,                      %排除波形过高的噪音波
                tar = 1;
                break;
            elseif Sig(j) > Sig(j-1)&&Sig(j) >= Sig(j+1)&&S(j)/t >= dis&&S(j)*tstep > n*downrange,
                r = j;
                break;
            end
        end
        if j+maxdistant > maxsize,
            break;
        end
        if tar == 1||j == i+maxdistant,
            [i] = nextR(j, maxdistant, S);
            continue;
        end
        R(m, 1) = (r - i);%第m个RR间期的长度
        R(m, 2) = i;%第m个RR间期的起点
        R(m, 3) = r;%第m个RR间期的终点
        m = m + 1;
        t = S(i);
        n = t*tstep;
        i = r;
    end

    Rnum = m - 1;%RR周期的个数
end

 %------------------------------寻找下一个R波----------------------------
function [k] = nextR(i, maxdis, M)
    temp = 0;
    l = i + 1;
    for j = i+1:i+maxdis                                             %寻找下一个R波
        if temp < M(j),
            temp = M(j);
            l = j;
        end
    end
    k = l;
end

