function [hrv,qtv,r10] = TEST( name,datanum )

data = load(name);
data = struct2cell(data);
data = cell2mat(data);

name = strcat(name,'_',num2str(datanum));
Sig = data(:,datanum);
% if(headstand(datanum) == 1)
% 	Sig = -Sig;
% end

maxsize= size(Sig,1);
%maxsize= 200000;
%disp(['Load ',name,' with size ', num2str(maxsize)]);

Fs = 125;               %采样频率
cell = 1000/Fs;         %时间单位
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
disp(['Rnum: ',num2str(Rnum)]);

if Rnum < 100
    hrv = NaN;
    qtv = NaN;
    r10 = NaN;
    return
end

%----------------------------Detect T wave----------------------------
winsize = 16;
p = 4;
lambada = 6;
Twave = zeros(Rnum,2);

for i = 1:Rnum
    if R(i,1) < 110
        ka = R(i,2) + floor(R(i,1) * 0.15) +18;
        kb = R(i,2) + ceil(R(i,1) * 0.7) - 4;
    else
        ka = R(i,2) + 35;
        kb = R(i,2) + ceil(R(i,1) * 0.2) + 50;
    end
    Twave(i,1) = ka - winsize;
    
    maxk1 = -inf;
    k1 = ka;
    mink2 = inf;
    k2 = ka;
    for k = ka:kb
        sk = 0;
        Ak = 0;
        for j = k - p:k + p
            sk = sk + Sig(j);
        end
        sk = sk / (2 * p + 1);
        for j = k - winsize + 1:k
            Ak = Ak + Sig(j) - sk;
        end
        if Ak > maxk1
            maxk1 = Ak;
            k1 = k;
        end
        if Ak < mink2
            mink2 = Ak;
            k2 = k;
        end
    end
    maxk1 = abs(maxk1);
    mink2 = abs(mink2);
    tmp = maxk1 / mink2;
    if 1 / lambada < tmp && lambada > tmp
        if k1 > k2
            Twave(i,2) = k1;
        else
            Twave(i,2) = k2;
        end
    else
        if maxk1 > mink2
            Twave(i,2) = k1;
        else
            Twave(i,2) = k2;
        end
    end
end
%------------------------ Detect Q wave ----------------------------
Qwave = zeros(Rnum,1);
for i = 1:Rnum
    j =    R(i,2) - 2;
    endj = R(i,2) - 15;
    min = Sig(j+1);
    if(j <= 0 || endj <= 0)
        Qwave(i) = 1;
        continue;
    end
    while j > endj
        if Sig(j) < min
            min = Sig(j);
        else
            break;
        end
        j = j - 1;
    end
    Qwave(i) = j-1;
end

%-------------------    draw    -------------------------
QTwave = zeros(Rnum,1);
QTwave(:) = Twave(:,2) - Qwave;

figure
plot(Sig);
figure
plot(Sig);
hold on;
plot(Twave(:,2),Sig(Twave(:,2)),'+r');
hold on;
plot(Qwave(:),Sig(Qwave(:)),'+g');

figure
plot(R(1:Rnum,1));
figure
plot(QTwave);

R(1:Rnum,1) = meanFilter(R(1:Rnum,1),Rnum,30,0.05);
%R(1:Rnum,1) = medfilt1(R(1:Rnum,1),3,30); 
QTwave(:) = meanFilter(QTwave(:),Rnum,30,0.05);
%QTwave(:) = medfilt1(QTwave(:),3,30); 

Ni = 150;
%EffectiveCoverage = Ni:Rnum - Ni;


RRew = mainFilter(R(:,1),Rnum,Ni);

%QTlow = low_pass_filter(QTwave(:));
QTew = mainFilter(QTwave(:),Rnum,Ni);

radioew = RRew ./ QTew;

radio = R(1:Rnum,1) ./ QTwave;



figure
s = 1000:Rnum-Ni;
plot(s,R(s,1),'y');
hold on;
plot(s,RRew(s),'g');
hold on;
plot(s,QTwave(s),'y');
hold on;
plot(s,QTew(s),'g');
hold on;
plot(s,radio(s),'y');
hold on;
plot(s,radioew(s),'g');

title('Waves compare to before','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
xlabel('Time in cell','FontName','Times New Roman','FontSize',14);
ylabel('Wave','FontName','Times New Roman','FontSize',14,'Rotation',0);

% figure
% plot(RRew(s),QTew(s),'.g','Markersize',2);
% title('RR vs QT','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
% xlabel('RR','FontName','Times New Roman','FontSize',14);
% ylabel('QT','FontName','Times New Roman','FontSize',14,'Rotation',0);


%----------------HRV and QTV-----------------
hrv = RV(R(s,1),Rnum,cell);
disp('HRV:');
disp(print_rv(hrv));

qtv = RV(QTwave(s),Rnum,cell);
disp('QTV:');
disp(print_rv(qtv));

radiov = RV(radio(s),Rnum,cell);
disp('RADIOV:');
disp(print_rv(radiov));

%------------Rgress models-----------------
%RRew = mainFilter(R(:,1),Rnum,Ni);
r10 = regress10(QTwave(s)*cell/1000,RRew(s)*cell/1000,true);

% %------------Inc and  Dec intervals-----------
% RRlow = low_pass_filter(R(1:Rnum,1));
% [inc,dec,incnum,decnum] = incAndDecSection(RRlow,Rnum,R,10000,1000,true);
% incnum
% decnum

% [max_inc,iidx] = max(inc(1:incnum,1));
% disp(['Max inc ',num2str((R(inc(iidx,3),2) - R(inc(iidx,2),2))*cell/1000),' seconds']);
% figure
% is = inc(iidx,2):inc(iidx,3);
% plot(RRew(is),QTwave(is),'.r');

% [max_dec,didx] = max(dec(1:decnum,1));
% disp(['Max dec ',num2str((R(dec(didx,3),2) - R(dec(didx,2),2))*cell/1000),' seconds']);
% hold on;
% ds = dec(didx,2):dec(didx,3);
% plot(RRew(ds),QTwave(ds),'.g');

% title('Max inc and dec interval','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
% xlabel('RR','FontName','Times New Roman','FontSize',14);
% ylabel('QT','FontName','Times New Roman','FontSize',14,'Rotation',0);

% % figure
% % for i = 1:incnum  
% %     s = inc(i,2):inc(i,3);hold on;
% %     plot(RRew(s),QTew(s),'.r','Markersize',2);
% % end 
% % for i = 1:decnum  
% %     s = dec(i,2):dec(i,3);hold on;
% %     plot(RRew(s),QTew(s),'.g','Markersize',2);
% % end











function[data] = low_pass_filter(wave)
Fs = 1;  % Sampling Frequency
N  = 2;     % Order
Fc = 0.001;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');
data=filter(Hd,wave);

function[data] = band_pass_filter(wave)
Fs = 1;  % Sampling Frequency

N   = 2;     % Order
Fc1 = 0.03;  % First Cutoff Frequency
Fc2 = 0.05;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');
data=filter(Hd,wave);


function[WAVEew] = mainFilter(wave,Rnum,Ni)
%Ni = 150;
u = 2 / (1 + Ni);
K = 1 / (1 - (1-u)^Ni);
WAVEew = zeros(Rnum,1);

param = K * u * ((1-u) .^ (0:Ni-1));

for i = Ni:Rnum
    WAVEew(i) = sum(wave(i-Ni+1:i) .* param(:) );
%     for j = -Ni + 1:0
%         WAVEew(i) = WAVEew(i) + K * u * (1-u)^-j * wave(i+j);
%     end
end

% @wave => wave to be filtered
% @Rnum => RR number
% @w    => window size
function[wave] = meanFilter(wave,Rnum,w,level)
%w = 30;
%level = 0.2
avg = zeros(Rnum - w,1);
for i = 1:w
    avg = avg + wave(i:Rnum + i - w - 1);
end
avg = avg / w;
up = 1 + level;
down = 1 - level;
for i = w:Rnum - w
    a = avg(i - w + 1);
    if(wave(i) > a * up || wave(i) < a * down)
        wave(i) = a;
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