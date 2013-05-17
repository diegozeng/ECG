%----------------------HRV����------------------------------------------
function [hrv] = HRV(name, von, datanum)
%----------------------��������----------------------------------------------
s = load(name);                                                   %��������
c = struct2cell(s);
N =cell2mat(c);                                              %�洢3�����ĵ��ź�
h = size(N);                                                    
max = h(1);                                                    %���źų��ȴ���max����
for i = 1:3
    if von(i) == 1
        N(:, i) = -N(:, i);                                       %������ͼ���������������N��
    end
end
hrv.name = str2num(name);
%----------------------HRV��ʱ�����---------------------------------------
%----------------------����------------------------------------------------
Fs = 125;                                                             %����Ƶ��
cell = 1000/Fs;                                                     %ʱ�䵥λ
Rmax = 150000;                                                  %RR���ڸ����޶�
lenghth =30;                                                       %���ʶγ���(DC/AC)
len = 300;                                                            %Y���������
five = 300000/cell;                                               %������ڰ����Ĳ��������
mindistant = 40;                                                  %��R�������������� 
maxdistant = 5000/cell;                                       %��R���������Զ����
uprange = 3;                                                       %��ֵ�����Ͻ�
downrange = 0.5;                                                %��ֵ�����½�
dis = 0.5;                                                             %б�ʹ���ָ��
tstep = 3;                                                            %б�ʲ���
d = 0.2;                                                               %ȥ��RR���ڱȵķ�Χ
w = 20;                                                                %ȥ��RR���ڴ��ڴ�С
Tomax = 100;                                                      %���To����
kran = 0.6;                                                           %��ʼ����ֵ
dn = 13;                                                              %R����Q����S����������Χ
Trange = 120;                                                      %�����粫QRSȺ�жϽ���
sizeTo = 23;                                                         %To����Ĵ�С
Tsize = 64;                                                          %T���罻��������Ĳ���
TWAstep = 3;                                                      %ͬ�������뾶
Tcell = 20;                                                           %�ĵ�ѹ��λ
Fc = 40;                                                               %��ֹƵ��
FR = 1;                                                                %������ֵ�˲����ڰ뾶
nw = 20;                                                              %�������������С
%----------------------����------------------------------------------------
M = zeros(max, 1);                                   %����i�������i+-tstep����֮���б���������M
G = zeros(Rmax, 3);                                        %�������RR���ں����Ӧ��ĺ�����
B = zeros(Rmax, 1);                                              %�ڼ���DC��ACʱ���������ȥ��������RR����
X = zeros(lenghth+1, 1);                                      %���ڴ洢�����������ʶξ�ֵ
Z = zeros(lenghth+1, 1);                                      %���ڴ洢�����������ʶξ�ֵ
Y = zeros(len, 1);                                                  %���ڴ洢ÿ����ӵ�RR���ھ�ֵ
I = zeros(len, 1);                                                  %���ڴ��ASDNN
P = zeros(maxdistant, 1);                                     %����ͳ��RR����Ƶ��
% odd = zeros(Tsize/2, 200);                              %���ڴ洢TWA�����������Ĳ�
% even = zeros(Tsize/2, 200);                             %���ڴ洢TWA������ż���Ĳ�
%----------------------Ҫ�������-------------------------------------------
hrv.MEAN = 0;                                                     %��ֵ
hrv.SDNN = 0;                                                     %�����׼��
hrv.SDANN = 0;                                                   %��ֵ��׼�� 
hrv.ASDNN = 0;                                                   %��׼���ֵ
hrv.DC = 0;                                                          %���������
hrv.AC = 0;                                                          %���������
hrv.TI = 0;                                                            %����ָ��
To = zeros(Tomax, sizeTo);                             %���ڴ��Toֵ�Լ�To֮���ǰ20��RR����
Ts = zeros(Tomax, 2);                                     %���ڴ��TS
%TWAM = zeros(Tsize, 2);                                    %���ڴ��TWA�����е�T����ѹƽ��ֵ
%----------------------����RR����------------------------------------------
% N(1:5000000, datanum) = wden(N(1:5000000, datanum), 'heursure', 's', 'one', 5, 'sym8');
% N(5000001:max, datanum) = wden(N(5000001:max, datanum), 'heursure', 's', 'one', 5, 'sym8');
if von(4) == 1,
    zer =  zeros(3000, 3);                                           %�������ò��ε�ȥ��
    N(1:3000, :) = zer;
end
[M, start] = position(M, tstep, kran, N, max, datanum); %�ҵ���һ��R�����壬���ز����M��˫�µ�����
t = M(start);                                             %��һ��RR�����Mֵ
m = 1;                                                                 %����ͳ��R���RR���ڵĹ�ϵ
i = start;                                                         %�ӵ�һ��R�����忪ʼɨ��
j = i + 1; 
n = t*tstep;                                                         %��һ��R������ĸ߶�       
while i+maxdistant < max                             %�ҳ�R�㲢���RR���ڴ���A
    if M(i)/t < dis||N(i, datanum) < N(i-1, datanum)||N(i, datanum) < N(i+1, datanum)||M(i)*tstep < n*downrange||M(i)*tstep > n*uprange,%��ȥαR��
        [i] = nextR(i, maxdistant, M);
        continue;
    end
    tar = 0;
    for j = i+mindistant:i+maxdistant     
        if M(j)*tstep > n*uprange,                      %�ų����ι��ߵ�������
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
    if tar == 1||j == i+maxdistant,                         %��ȥαR�������������粫��
        [i] = nextR(j, maxdistant, M);
        continue;
    end
    G(m, 1) = r - i;                                            %��m��RR���ڵĳ���
    G(m, 2) = i;                                                %��m��RR���ڵ����
    G(m, 3) = r;                                                %��m��RR���ڵ��յ�
    m = m + 1;
    t = M(i); 
    n = t*tstep; 
    i = r;
end
Rnum = m - 1;                                               %RR���ڵĸ���
[G, Rnum, To, Tonum] = myfilter(G, d, Rnum, w, To, sizeTo);%ѡ��������뷽ʽ
for i = 1:Rnum
    hrv.MEAN = hrv.MEAN + G(i, 1);
end
hrv.MEAN = hrv.MEAN / Rnum;         %��ֵ��ֵ(δ��cell)
%N_R = Rnum/(m-1);
%----------------------����DC��AC----------------------------------------------
Rtemp = G(1, 1);                                            %���ڴ洢��һ��δ����ȥ��RR����ֵ
j = 1;                                                                   %����ͳ��RR���ڵĸ���
for i = 1:Rnum                                               %��ȥ��������RR���ڴ���B��
    if G(i, 1)/Rtemp <= 1.05&&G(i, 1)/Rtemp >= 0.95,
        B(j) = G(i, 1);
        Rtemp = G(i, 1);
        j = j + 1;
    end
end
RRnum = j - 1;                                               %����DCʱȥ����RR���ڵĸ���
Dnum = 0;                                                           %���ڴ洢�������ڵĸ���
Anum = 0;                                                           %���ڴ洢�������ڵĸ���
for i = 16:RRnum-1                                        %��������ʶ���ÿ�����Ӧ���ܺ�
    if B(i) > B(i-1),
        for j = 1:31
            X(j) = X(j) + B(i-(16-j));
        end
        Dnum = Dnum + 1;
    end
end
for i = 15:RRnum-1                                        %��������ʶ���ÿ�����Ӧ���ܺ�
    if B(i) < B(i-1),
        for j = 1:31
            Z(j) = Z(j) + B(i-(15-j));
        end
        Anum = Anum + 1;
    end
end
for i = 1:31                                                          %�����ʶ�ÿ�����ƽ��  
    X(i) = X(i)/Dnum;
end
for i = 1:31                                                          %�����ʶ�ÿ�����ƽ��  
    Z(i) = Z(i)/Anum;
end
hrv.DC = (X(16) + X(17) - X(15) - X(14))*cell/4;      %DC��ֵ
hrv.AC = (Z(15) + Z(16) - Z(14) - Z(13))*cell/4;      %AC��ֵ
%----------------------��������ƽ����---------------------------------------
for i = 1:Rnum
    hrv.SDNN = hrv.SDNN + (G(i, 1) - hrv.MEAN)^2;
end
hrv.SDNN = sqrt(hrv.SDNN/Rnum)*cell;          %����ƽ�����ֵ
%----------------------�����ֵ---------------------------------------------
hrv.MEAN = hrv.MEAN*cell;
%-----------------------�����ֵ��׼��--------------------------------------
RF = 1;                                                           %���ڼ�¼5����RR�����ܸ���
k = 1;                                                                  %RR������ţ���l��¼ÿ��5min��ʼ��RR���
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
hrv.SDANN = sqrt(hrv.SDANN/RF)*cell;                %��ֵ��׼���ֵ
%------------------------�����׼���ֵ-------------------------------------
for i = 1:RF
    hrv.ASDNN = hrv.ASDNN + I(i);
end
hrv.ASDNN = hrv.ASDNN*cell / RF; 
%--------------------------����ָ��-----------------------------------------
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
% %--------------------------HRV��Ƶ�׷���------------------------------------
% %-----------------------------------������-------------------------------------
% hrv.VLF = 0;
% hrv.LF = 0;
% hrv.HF = 0;
% hrv.TP = 0;
% hrv.L_H = 0;
% %-----------------------------Ƶ�׷�������------------------------------------
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
% %---------------------------------HRT����--------------------------------------
% %--------------------------HRT�������粫���ж�-------------------------------
% %------------------------����������Ƿ�T������--------------------------------
% [tp, base] = tposition(N, G(1, 2), G(1, 3), dn, datanum);
% if N(tp, datanum) <= N(tp-5, datanum) && N(tp, datanum) <= N(tp+5, datanum),
%     Tonum = 0;                                               %��������T�������õ�ʱ�򲻿��������粫
%     To = zeros(Tomax, sizeTo); 
% end
% %----------------------------------��ʼ���--------------------------------------
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
%         if (right - left)*cell < Trange,                        %��QRS����ʱ��С��0.12s���粫�ų�
%             To(i:Tonum-1, :) = To(i+1:Tonum, :);
%             To(Tonum, :) = z;
%             Tonum = Tonum - 1;
%         else
%             i = i + 1;
%         end
%     end
% end
% %-------------------------------T�����ü��-------------------------------------
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
% %------------------------------HRT��TS�ķ���----------------------------------
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
% %------------------------------------TWA����------------------------------------
% i = 2;
% Tstart = 0;
% while i + Tsize-1 <= Rnum                             %ȷ��TWA�����Ŀ�ʼ��
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
%----------------------------------TWAƵ�����----------------------------------
% if Tstart + Tsize-1 <= Rnum,
%     for i = Tstart:Tstart+Tsize-1
%         [tp, base] = tposition(N, G(i, 2), G(i, 3), dn, datanum);
%         TWAM(i - Tstart + 1, 2) = tp;                
%         TWAM(i - Tstart + 1, 1) = mean(N(tp-TWAstep:tp+TWAstep, datanum)) - base; 
%     end                                                                 %��T������Ϊ����TWAstepΪ�뾶�Ĵ����ڵ�ƽ��ֵ����TWAM
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
% %----------------------------------TWAʱ�����----------------------------------
% if Tstart + Tsize-1 <= Rnum,
%     for i = Tstart:Tstart+Tsize-1                             %����Ư��У��
%         base = (N(G(i, 2)+dn, datanum) + N(G(i, 3)-dn, datanum)) / 2;
%         s = G(i, 2);
%         e = G(i, 3);
%         N(s:e, datanum) = N(s:e, datanum) - base;
%     end
%     s = G(Tstart, 2);
%     e = G(Tstart + Tsize - 1, 3);
%     FN = N(s:e, datanum);
%     [i, j] = butter(8, Fc*2/Fs, 'low');
%     FN = filter(i, j, FN);                                           %��40HzΪ��ֹƵ�ʽ��е�ͨ�˲�
%     m = G(Tstart, 2) - 1;
%     for i = Tstart:Tstart+Tsize-1                             %�ֳ���ż�Ĳ�
%         s = G(i, 2) - m;
%         e = G(i, 3) - m;
%         if mod(i, 2) == 1,
%             odd((i+1)/2, 1:e-s+1) = FN(s:e);
%         else
%             even(i/2, 1:e-s+1) = FN(s:e);
%         end
%     end
%     [odd] = mid_filter(odd, FR);                              %�������Ĳ����н�����ֵ����
%     [even] = mid_filter(even, FR);                           %��ż���Ĳ����н�����ֵ����
%     odd = sort(odd, 1);                               
%     even = sort(even, 1);
%     midodd = (odd(Tsize/4, :) + odd(Tsize/4+1, :)) / 2;    %��������Ĳ�����λ������
%     mideven = (even(Tsize/4, :) + even(Tsize/4+1, :)) / 2;%���ż���Ĳ�����λ������
%     [oddt] = ttposition(midodd, dn, 90);                               %�ҳ�������λ���Ĳ���T����ֵ
%     [event] = ttposition(mideven, dn, 90);                            %�ҳ�ż����λ���Ĳ���T����ֵ
%     t = floor((oddt + event) / 2);                                     %�����ߵ�ƽ����ΪTWA��������
%     hrv.twa = Tcell*mean(abs(midodd(t-dn:t+dn) - mideven(t-dn:t+dn)));        %����TWA
%     x = zeros(1, 1:nw-w+1);
%     [K, x] = noise_k(mideven, event, nw, x);         %���������K
% else
%     Tstart = 0;
% end
%--------------------------�������õ����Ӻ���--------------------------
%------------------------------��ֵ�˲�����-----------------------------
function [x] = mid_filter(x, n)
[h, l] = size(x);                                                      %���С������������ڶ��е����ݽ����˲�
for i = 1:h
    for j = 1+n:l-n
        y = sort(x(i, j-n:j+n));
        x(i, j) = y(n+1);
    end
end
%------------------------------���������------------------------------------
function [K, x] = noise_k(N, start, wide, x)
w = 5;                                                                 %��ȡ�ֲ�����ʱ�󷽲�����ݿ��
for i = start:start+wide-w
    x(i-start+1) = var(N(i:i+w-1));
end
x = sort(x);
K = x(wide-w+1) / x(1);
%-------------------------Ƶ�����ȷ��T������λ��--------------------------
function [tp, m] = tposition(N, leftR, rightR, dn, datanum)
t = 0;
m = (N(leftR+dn, datanum) + N(rightR-dn, datanum)) / 2;%m����RR�����еĻ���λ��
for k = leftR+dn:rightR-dn
    if t < abs(N(k, datanum) - m),
        t = abs(N(k, datanum) - m);
        tp = k;            
    end
end
%-------------------------ʱ�����ȷ��T������λ��--------------------------
function [tp] = ttposition(N, leftR, rightR)
t = 0;
m = (N(leftR) + N(rightR)) / 2;%m����RR�����еĻ���λ��
for k = leftR:rightR
    if t < abs(N(k) - m),
        t = abs(N(k) - m);
        tp = k;            
    end
end
%-----------------------------ȷ����һ��R��λ��------------------------
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
%---------------------------------ȥ��------------------------------------
function [G, Rnum, To, Tonum] = myfilter(G, d, Rnum, w, To, sizeTo)    
Tonum = 0;                                                         %To����
target = 0;                                                           %�жϱ��     
n = 0;                                                                  %���ڼ�¼�粫���ǰ20���Ķ�����
t = 0;                                                                   %��¼�粫ǰ�������쳣RR���ڵĺ�
i = w+1;
while i+w <= Rnum
    sum = 0;
    for j = i-w:i+w
        sum = G(j, 1) + sum;
    end
    sum = (sum - G(i, 1))/(2*w);                             %���㴰���ڵ�ƽ��RR����
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
            if t <= 2*G(i-1, 1),                                    %��������Ъ�Ƿ�С�ڵ���2������RR���ڼ��
                target = 0;
                ToRRL = G(i-2, 1) + G(i-1, 1);                    %���RR��
                ToRRR = G(i+1, 1) + G(i+2, 1);                  %�Ҳ�RR��
                Tonum = Tonum + 1;
                To(Tonum, 1) = G(i, 2);
                To(Tonum, sizeTo) = G(i, 3);
                To(Tonum, 2) = (ToRRR - ToRRL) / ToRRL;
                n = 20;                                                    %���ڼ�¼�粫���ǰ20���Ķ����� 
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
%------------------------------Ѱ����һ��R��----------------------------
function [k] = nextR(i, maxdis, M)
temp = 0;
l = i + 1;
for j = i+1:i+maxdis                                             %Ѱ����һ��R��
     if temp < M(j),
         temp = M(j);
         l = j;
     end
end
k = l;