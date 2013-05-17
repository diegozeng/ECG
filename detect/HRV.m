%----------------------HRV分析------------------------------------------
function [hrv] = HRV(name, von, datanum)
%----------------------导入数据----------------------------------------------
s = load(name);                                                   %导出数据
c = struct2cell(s);
N =cell2mat(c);                                              %存储3导的心电信号
h = size(N);                                                    
max = h(1);                                                    %把信号长度存入max当中
for i = 1:3
    if von(i) == 1
        N(:, i) = -N(:, i);                                       %将倒置图像调整过来并存入N中
    end
end
hrv.name = str2num(name);
%----------------------HRV的时域分析---------------------------------------
%----------------------常量------------------------------------------------
Fs = 125;                                                             %采样频率
cell = 1000/Fs;                                                     %时间单位
Rmax = 150000;                                                  %RR周期个数限定
lenghth =30;                                                       %心率段长度(DC/AC)
len = 300;                                                            %Y的最大容量
five = 300000/cell;                                               %五分钟内包含的采样点个数
mindistant = 40;                                                  %对R点的搜索最近距离 
maxdistant = 5000/cell;                                       %对R点的搜索最远距离
uprange = 3;                                                       %幅值过滤上界
downrange = 0.5;                                                %幅值过滤下界
dis = 0.5;                                                             %斜率过滤指标
tstep = 3;                                                            %斜率步长
d = 0.2;                                                               %去除RR间期比的范围
w = 20;                                                                %去除RR间期窗口大小
Tomax = 100;                                                      %最大To个数
kran = 0.6;                                                           %初始点阈值
dn = 13;                                                              %R波到Q波和S波的搜索范围
Trange = 120;                                                      %室性早搏QRS群判断界限
sizeTo = 23;                                                         %To数组的大小
Tsize = 64;                                                          %T波电交替的连续心搏数
TWAstep = 3;                                                      %同步采样半径
Tcell = 20;                                                           %心电压单位
Fc = 40;                                                               %截止频率
FR = 1;                                                                %渐量中值滤波窗口半径
nw = 20;                                                              %噪声测量区域大小
%----------------------数组------------------------------------------------
M = zeros(max, 1);                                   %将第i个点与第i+-tstep个点之间的斜率算出存入M
G = zeros(Rmax, 3);                                        %用来存放RR周期和其对应点的横坐标
B = zeros(Rmax, 1);                                              %在计算DC和AC时，用来存放去除噪声的RR周期
X = zeros(lenghth+1, 1);                                      %用于存储减速周期心率段均值
Z = zeros(lenghth+1, 1);                                      %用于存储加速周期心率段均值
Y = zeros(len, 1);                                                  %用于存储每五分钟的RR周期均值
I = zeros(len, 1);                                                  %用于存放ASDNN
P = zeros(maxdistant, 1);                                     %用于统计RR周期频数
% odd = zeros(Tsize/2, 200);                              %用于存储TWA分析中奇数心搏
% even = zeros(Tsize/2, 200);                             %用于存储TWA分析中偶数心搏
%----------------------要计算的量-------------------------------------------
hrv.MEAN = 0;                                                     %均值
hrv.SDNN = 0;                                                     %总体标准差
hrv.SDANN = 0;                                                   %均值标准差 
hrv.ASDNN = 0;                                                   %标准差均值
hrv.DC = 0;                                                          %心脏减速率
hrv.AC = 0;                                                          %心脏加速率
hrv.TI = 0;                                                            %三角指数
To = zeros(Tomax, sizeTo);                             %用于存放To值以及To之后的前20个RR间期
Ts = zeros(Tomax, 2);                                     %用于存放TS
%TWAM = zeros(Tsize, 2);                                    %用于存放TWA分析中的T波电压平均值
%----------------------计算RR周期------------------------------------------
% N(1:5000000, datanum) = wden(N(1:5000000, datanum), 'heursure', 's', 'one', 5, 'sym8');
% N(5000001:max, datanum) = wden(N(5000001:max, datanum), 'heursure', 's', 'one', 5, 'sym8');
if von(4) == 1,
    zer =  zeros(3000, 3);                                           %对于无用波形的去除
    N(1:3000, :) = zer;
end
[M, start] = position(M, tstep, kran, N, max, datanum); %找到第一个R波波峰，返回波峰和M（双坡导数）
t = M(start);                                             %第一个RR波峰的M值
m = 1;                                                                 %用于统计R点和RR周期的关系
i = start;                                                         %从第一个R波波峰开始扫描
j = i + 1; 
n = t*tstep;                                                         %第一个R波波峰的高度       
while i+maxdistant < max                             %找出R点并算出RR周期存入A
    if M(i)/t < dis||N(i, datanum) < N(i-1, datanum)||N(i, datanum) < N(i+1, datanum)||M(i)*tstep < n*downrange||M(i)*tstep > n*uprange,%除去伪R波
        [i] = nextR(i, maxdistant, M);
        continue;
    end
    tar = 0;
    for j = i+mindistant:i+maxdistant     
        if M(j)*tstep > n*uprange,                      %排除波形过高的噪音波
            tar = 1;
            break;
        elseif N(j, datanum) > N(j-1, datanum)&&N(j, datanum) >= N(j+1, datanum)&&M(j)/t >= dis&&M(j)*tstep > n*downrange,
            r = j;
            break;
        end
    end
    if j+maxdistant > max,
        break;
    end
    if tar == 1||j == i+maxdistant,                         %除去伪R波（包括室性早搏）
        [i] = nextR(j, maxdistant, M);
        continue;
    end
    G(m, 1) = r - i;                                            %第m个RR间期的长度
    G(m, 2) = i;                                                %第m个RR间期的起点
    G(m, 3) = r;                                                %第m个RR间期的终点
    m = m + 1;
    t = M(i); 
    n = t*tstep; 
    i = r;
end
Rnum = m - 1;                                               %RR周期的个数
[G, Rnum, To, Tonum] = myfilter(G, d, Rnum, w, To, sizeTo);%选择何种消噪方式
for i = 1:Rnum
    hrv.MEAN = hrv.MEAN + G(i, 1);
end
hrv.MEAN = hrv.MEAN / Rnum;         %均值的值(未乘cell)
%N_R = Rnum/(m-1);
%----------------------计算DC和AC----------------------------------------------
Rtemp = G(1, 1);                                            %用于存储上一个未被舍去的RR周期值
j = 1;                                                                   %用以统计RR周期的个数
for i = 1:Rnum                                               %将去除噪声的RR周期存入B中
    if G(i, 1)/Rtemp <= 1.05&&G(i, 1)/Rtemp >= 0.95,
        B(j) = G(i, 1);
        Rtemp = G(i, 1);
        j = j + 1;
    end
end
RRnum = j - 1;                                               %计算DC时去噪后的RR周期的个数
Dnum = 0;                                                           %用于存储减速周期的个数
Anum = 0;                                                           %用于存储加速周期的个数
for i = 16:RRnum-1                                        %求减速心率段中每个点对应的总和
    if B(i) > B(i-1),
        for j = 1:31
            X(j) = X(j) + B(i-(16-j));
        end
        Dnum = Dnum + 1;
    end
end
for i = 15:RRnum-1                                        %求加速心率段中每个点对应的总和
    if B(i) < B(i-1),
        for j = 1:31
            Z(j) = Z(j) + B(i-(15-j));
        end
        Anum = Anum + 1;
    end
end
for i = 1:31                                                          %求心率段每个点的平均  
    X(i) = X(i)/Dnum;
end
for i = 1:31                                                          %求心率段每个点的平均  
    Z(i) = Z(i)/Anum;
end
hrv.DC = (X(16) + X(17) - X(15) - X(14))*cell/4;      %DC的值
hrv.AC = (Z(15) + Z(16) - Z(14) - Z(13))*cell/4;      %AC的值
%----------------------计算总体平方差---------------------------------------
for i = 1:Rnum
    hrv.SDNN = hrv.SDNN + (G(i, 1) - hrv.MEAN)^2;
end
hrv.SDNN = sqrt(hrv.SDNN/Rnum)*cell;          %总体平方差的值
%----------------------计算均值---------------------------------------------
hrv.MEAN = hrv.MEAN*cell;
%-----------------------计算均值标准差--------------------------------------
RF = 1;                                                           %用于记录5分钟RR间期总个数
k = 1;                                                                  %RR间期序号，用l记录每个5min开始的RR序号
i = 1;
while i+five-1 < max  
    l = k;
    temp = 0;
    while k < Rnum
        if G(k, 3) > i+five-1
            break;
        end
        temp = temp + G(k, 1);
        k = k + 1;
    end
   if k < Rnum,
       if k == l,
           i = i + five;
           continue;
       end
       Y(RF) = temp / (k - l);
       for j = l:k-1
           I(RF) = I(RF) + (G(j, 1) - Y(RF))^2;
       end
       I(RF) = sqrt(I(RF)/(k - l));
       RF = RF + 1;
   else
       break;
   end
   i = i + five;
end
RF = RF - 1;
temp = 0;
for i = 1:RF
    temp = temp + Y(i);
end
avg = temp/RF;
for i = 1:RF
    hrv.SDANN = hrv.SDANN + (Y(i) - avg)^2;
end
hrv.SDANN = sqrt(hrv.SDANN/RF)*cell;                %均值标准差的值
%------------------------计算标准差均值-------------------------------------
for i = 1:RF
    hrv.ASDNN = hrv.ASDNN + I(i);
end
hrv.ASDNN = hrv.ASDNN*cell / RF; 
%--------------------------三角指数-----------------------------------------
for i = 1:Rnum
    for j = 1:maxdistant
        if G(i, 1) == j,
            P(j) = P(j) + 1;
            break;
        end
    end
end
temp = 0;
for i = 1:maxdistant
    if temp < P(i),
        temp = P(i);
    end
end
hrv.TI = Rnum / temp;
% %--------------------------HRV的频谱分析------------------------------------
% %-----------------------------------待求量-------------------------------------
% hrv.VLF = 0;
% hrv.LF = 0;
% hrv.HF = 0;
% hrv.TP = 0;
% hrv.L_H = 0;
% %-----------------------------频谱分析程序------------------------------------
% E = G(:, 1)*cell;
% C = fft(E, Rnum);
% D = C.*conj(C)/Rnum;
% s = 0:Rnum-1;
% s = s/Rnum; 
% for i = 1:Rnum
%     if (s(i) < 0.0033&&s(i+1)>0.0033)||(s(i) == 0.0033)
%         first = i;
%     elseif (s(i) < 0.04&&s(i+1)>0.04)||(s(i) == 0.04)
%         sr = i;
%     elseif (s(i) < 0.15&&s(i+1)>0.15)||(s(i) == 0.15)
%         third = i;
%     elseif (s(i) < 0.4&&s(i+1)>0.4)||(s(i) == 0.4)
%         four = i;
%         break;
%     end
% end
% for i = first+1:sr
%     hrv.VLF = hrv.VLF + (D(i-1)+D(i))/2/Rnum;
% end
% for i = sr+1:third
%     hrv.LF = hrv.LF + (D(i-1)+D(i))/2/Rnum;
% end
% for i = third+1:four
%     hrv.HF = hrv.HF + (D(i-1)+D(i))/2/Rnum;
% end
% hrv.VLF = sqrt(hrv.VLF);
% hrv.LF = sqrt(hrv.LF);
% hrv.HF = sqrt(hrv.HF);
% hrv.TP = hrv.VLF+hrv.LF+hrv.HF;
% hrv.L_H = hrv.LF/hrv.HF;
% %---------------------------------HRT分析--------------------------------------
% %--------------------------HRT中室性早搏的判定-------------------------------
% %------------------------检测正常波是否T波倒置--------------------------------
% [tp, base] = tposition(N, G(1, 2), G(1, 3), dn, datanum);
% if N(tp, datanum) <= N(tp-5, datanum) && N(tp, datanum) <= N(tp+5, datanum),
%     Tonum = 0;                                               %所有正常T波都倒置的时候不考虑室性早搏
%     To = zeros(Tomax, sizeTo); 
% end
% %----------------------------------开始检测--------------------------------------
% z = zeros(1, sizeTo);
% if Tonum ~= 0,
%     i = 1;
%     while i <= Tonum                        
%         j = To(i, 1);
%         m = N(j, datanum);
%         k = j;
%         while k >= j-dn
%             if N(k, datanum) < m,
%                 m = N(k, datanum);
%                 left = k;
%             end
%             k = k - 1;
%         end
%         m = N(j, datanum);
%         for k = j:j+dn
%             if N(k, datanum) < m,
%                 m = N(k, datanum);
%                 right = k;
%             end
%         end
%         if (right - left)*cell < Trange,                        %将QRS波的时长小于0.12s的早搏排除
%             To(i:Tonum-1, :) = To(i+1:Tonum, :);
%             To(Tonum, :) = z;
%             Tonum = Tonum - 1;
%         else
%             i = i + 1;
%         end
%     end
% end
% %-------------------------------T波倒置检测-------------------------------------
% if Tonum ~= 0,
%     i = 1;
%     while i <= Tonum
%         [tp, base] = tposition(N, To(i, 1), To(i, sizeTo), dn, datanum);
%         if N(tp, datanum) > N(tp-5, datanum) || N(tp, datanum) > N(tp+5, datanum),
%             To(i:Tonum-1, :) = To(i+1:Tonum, :);
%             To(Tonum, :) = z;
%             Tonum = Tonum - 1;    
%         else
%             i = i + 1;
%         end
%     end
% end
% %------------------------------HRT中TS的分析----------------------------------
% if Tonum ~= 0,
%     x = [1; 2; 3; 4; 5];
%     mx = mean(x);
%     sx = std(x);
%     sx = (sx^2)*4;
%     for i = 1:Tonum
%         if To(i, 7) == 0,
%             Ts(i) = -100;
%             continue;
%         end
%         temp = -100;
%         for j = 3:18
%             y = To(i, j:j+4);
%             my = mean(y);
%             if temp < (y*x - 5*mx*my) / sx,
%                 temp = (y*x - 5*mx*my) / sx;
%                 k = j;
%             end
%         end
%         Ts(i, 1) = temp;
%         Ts(i, 2) = k;
%     end
% end
% %------------------------------------TWA分析------------------------------------
% i = 2;
% Tstart = 0;
% while i + Tsize-1 <= Rnum                             %确定TWA分析的开始点
%     Tstart = i;
%     for j = 2:Tsize
%         if G(i, 2) ~= G(i-1, 3),
%             break;
%         end
%     end
%     if j < Tsize,
%         Tstart = Tstart + j - 1;
%         i = Tstart;
%     else
%         Tstart = Tstart - 1;
%         break;
%     end
% end
%----------------------------------TWA频域分析----------------------------------
% if Tstart + Tsize-1 <= Rnum,
%     for i = Tstart:Tstart+Tsize-1
%         [tp, base] = tposition(N, G(i, 2), G(i, 3), dn, datanum);
%         TWAM(i - Tstart + 1, 2) = tp;                
%         TWAM(i - Tstart + 1, 1) = mean(N(tp-TWAstep:tp+TWAstep, datanum)) - base; 
%     end                                                                 %以T波顶点为中心TWAstep为半径的窗口内的平均值存入TWAM
% %     m = mean(TWAM(1:Tsize, 1));
% %     for i = 1:Tsize
% %         TWAM(i, 1) = TWAM(i, 1) - m;
% %     end
% 
% %     E = TWAM(:, 1)*Tcell;
% %     C = fft(E, Tsize);
% %     s = 0:Tsize-1;
% %     s = s / Tsize;
% %     plot(s(1:Tsize), abs(C(1:Tsize)));
%     
% else
%     Tstart = 0;
% end
% %----------------------------------TWA时域分析----------------------------------
% if Tstart + Tsize-1 <= Rnum,
%     for i = Tstart:Tstart+Tsize-1                             %基线漂移校正
%         base = (N(G(i, 2)+dn, datanum) + N(G(i, 3)-dn, datanum)) / 2;
%         s = G(i, 2);
%         e = G(i, 3);
%         N(s:e, datanum) = N(s:e, datanum) - base;
%     end
%     s = G(Tstart, 2);
%     e = G(Tstart + Tsize - 1, 3);
%     FN = N(s:e, datanum);
%     [i, j] = butter(8, Fc*2/Fs, 'low');
%     FN = filter(i, j, FN);                                           %以40Hz为截止频率进行低通滤波
%     m = G(Tstart, 2) - 1;
%     for i = Tstart:Tstart+Tsize-1                             %分成奇偶心搏
%         s = G(i, 2) - m;
%         e = G(i, 3) - m;
%         if mod(i, 2) == 1,
%             odd((i+1)/2, 1:e-s+1) = FN(s:e);
%         else
%             even(i/2, 1:e-s+1) = FN(s:e);
%         end
%     end
%     [odd] = mid_filter(odd, FR);                              %对奇数心搏进行渐量中值修正
%     [even] = mid_filter(even, FR);                           %对偶数心搏进行渐量中值修正
%     odd = sort(odd, 1);                               
%     even = sort(even, 1);
%     midodd = (odd(Tsize/4, :) + odd(Tsize/4+1, :)) / 2;    %算出奇数心搏的中位数波形
%     mideven = (even(Tsize/4, :) + even(Tsize/4+1, :)) / 2;%算出偶数心搏的中位数波形
%     [oddt] = ttposition(midodd, dn, 90);                               %找出奇数中位数心搏的T波峰值
%     [event] = ttposition(mideven, dn, 90);                            %找出偶数中位数心搏的T波峰值
%     t = floor((oddt + event) / 2);                                     %将两者的平均作为TWA测量中心
%     hrv.twa = Tcell*mean(abs(midodd(t-dn:t+dn) - mideven(t-dn:t+dn)));        %计算TWA
%     x = zeros(1, 1:nw-w+1);
%     [K, x] = noise_k(mideven, event, nw, x);         %计算信噪比K
% else
%     Tstart = 0;
% end
%--------------------------主函数用到的子函数--------------------------
%------------------------------中值滤波函数-----------------------------
function [x] = mid_filter(x, n)
[h, l] = size(x);                                                      %求行、列数，适用于对行的数据进行滤波
for i = 1:h
    for j = 1+n:l-n
        y = sort(x(i, j-n:j+n));
        x(i, j) = y(n+1);
    end
end
%------------------------------计算信噪比------------------------------------
function [K, x] = noise_k(N, start, wide, x)
w = 5;                                                                 %求取局部方差时求方差的数据宽度
for i = start:start+wide-w
    x(i-start+1) = var(N(i:i+w-1));
end
x = sort(x);
K = x(wide-w+1) / x(1);
%-------------------------频域分析确定T波顶点位置--------------------------
function [tp, m] = tposition(N, leftR, rightR, dn, datanum)
t = 0;
m = (N(leftR+dn, datanum) + N(rightR-dn, datanum)) / 2;%m代表RR间期中的基线位置
for k = leftR+dn:rightR-dn
    if t < abs(N(k, datanum) - m),
        t = abs(N(k, datanum) - m);
        tp = k;            
    end
end
%-------------------------时域分析确定T波顶点位置--------------------------
function [tp] = ttposition(N, leftR, rightR)
t = 0;
m = (N(leftR) + N(rightR)) / 2;%m代表RR间期中的基线位置
for k = leftR:rightR
    if t < abs(N(k) - m),
        t = abs(N(k) - m);
        tp = k;            
    end
end
%-----------------------------确定第一个R峰位置------------------------
function [M, start] = position(M, tstep, kran, N, max, datanum)
M(1+tstep:max-tstep) = ((2*N(1+tstep:max-tstep, datanum)-N(1:max-2*tstep, datanum)- N(1+2*tstep:max, datanum))/(2*tstep));
k = 0; 
for i = 1+tstep:4000
    if M(i)  > k,
        k = M(i);
    end
end
for i = 1+tstep:4000
    if M(i) > k*kran,
        break;
    end
end
m = -1000;
k = i - 2;
for j = i -2:i+4
    if m < N(j, datanum) && N(j, datanum) > N(j-1, datanum) && N(j, datanum) >= N(j+1, datanum),
        m = N(j, datanum);
        k = j;
    end
end
start = k;
if k == 3,
    M(3) = (N(3, datanum) - N(3+tstep, datanum)) / tstep;
end
%---------------------------------去噪------------------------------------
function [G, Rnum, To, Tonum] = myfilter(G, d, Rnum, w, To, sizeTo)    
Tonum = 0;                                                         %To个数
target = 0;                                                           %判断标记     
n = 0;                                                                  %用于记录早搏后的前20个心动周期
t = 0;                                                                   %记录早搏前后两个异常RR间期的和
i = w+1;
while i+w <= Rnum
    sum = 0;
    for j = i-w:i+w
        sum = G(j, 1) + sum;
    end
    sum = (sum - G(i, 1))/(2*w);                             %计算窗口内的平均RR间期
    if G(i, 1)<sum*(1-d),
        if target == 0,
            target = 1;
            t = G(i, 1);
        else
            target = 0;
            t = 0;
        end
        G(i:Rnum-1, :) = G(i+1:Rnum, :);
        G(Rnum, :) = [0 0 0];
        Rnum = Rnum-1;
    elseif G(i, 1)>sum*(1+d),
        if target == 1 ,
            t = t + G(i, 1);
            if t <= 2*G(i-1, 1),                                    %检查代偿间歇是否小于等于2倍正常RR间期检测
                target = 0;
                ToRRL = G(i-2, 1) + G(i-1, 1);                    %左侧RR和
                ToRRR = G(i+1, 1) + G(i+2, 1);                  %右侧RR和
                Tonum = Tonum + 1;
                To(Tonum, 1) = G(i, 2);
                To(Tonum, sizeTo) = G(i, 3);
                To(Tonum, 2) = (ToRRR - ToRRL) / ToRRL;
                n = 20;                                                    %用于记录早搏后的前20个心动周期 
                t = 0;
            else
                t = 0;
                target = 0;
            end
        end
        G(i:Rnum-1, :) = G(i+1:Rnum, :);
        G(Rnum, :) = [0 0 0];
        Rnum = Rnum-1;
    else
        t = 0;
        i = i + 1;
        target = 0;
        if n > 0,
            To(Tonum, 23 - n) = G(i-1, 1);
            n = n - 1;
        end
    end
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