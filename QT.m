function [AC,DC,qtSDNN,radioSDNN,qtMEAN,radioMEAN] = QT( name,datanum,headstand,showGraph )

data = load(name);
data = struct2cell(data);
data = cell2mat(data);

Sig = data(:,datanum);
if(headstand(datanum) == 1)
	Sig = -Sig;
end

%maxsize= size(Sig);
maxsize= 200000;
%disp(['Load ',name,' with size ', num2str(maxsize)]);

Fs = 125;%采样频率
cell = 1000/Fs;%时间单位
mindistant = 40;%对R点的搜索最近距离
maxdistant = 5000/cell;%对R点的搜索最远距离
uprange = 3;%幅值过滤上界
downrange = 0.5;%幅值过滤下界
dis = 0.5;%斜率过滤指标
tstep = 3;%斜率步长
kran = 0.6;%初始点阈值

S = zeros(maxsize,1);%slope
R = zeros(maxsize,3);%R wave

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

%------------If R number is too small, fail------------------
if(Rnum < 10)
	AC = NaN;
	DC = NaN;
	qtSDNN = NaN;
	radioSDNN = NaN;
    qtMEAN = NaN;
    radioMEAN = NaN;
	return;
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

%---------------------Draw maker of R wave-----------------------------
if showGraph
    plot(Sig(1:4000),'b');
    hold on;
    % tmp = R(:,2);
    % plot(tmp,data(tmp,datanum),'md','MarkerSize',5);
    for i = 1:40
        %plot(R(i,2),data(R(i,2),datanum),'md','MarkerSize',8);
        plot(Twave(i,1),Sig(Twave(i,1)),'md','MarkerSize',6);
        plot(Twave(i,2),Sig(Twave(i,2)),'md','MarkerSize',8);
       
        plot(Qwave(i),Sig(Qwave(i)),'md','MarkerSize',6);
        %line([R(i,2)-3 R(i,2)-3],[-20 20],'Color','r');
        %line([R(i,2)-10 R(i,2)-10],[-20 20],'Color','r');
    end
end

%-------------------------calc AC and DC----------------------------

QTwave = zeros(Rnum,2);
QTwave(:,1) = Twave(:,2) - Qwave;

pre = QTwave(1,1);
for i = 2:Rnum
    if QTwave(i,1) > pre
        QTwave(i,2) = 1;
    elseif QTwave(i,1) < pre
        QTwave(i,2) = -1;
    else
        QTwave(i,2) = 0;
    end
    pre = QTwave(i,1);
end
% plot(tmp);
% plot(QTwave(:,1));
% hold on;
% for i = 2:Rnum
%     if QTwave(i,2) == 1
%         plot(i,QTwave(i,1),'marker','square','color','r','markersize',4);
%     elseif QTwave(i,2) == -1
%         plot(i,QTwave(i,1),'marker','square','color','g','markersize',4);
%     end
% end


X = zeros(31,1);
Z = zeros(31,1);
ACnum = 0;
DCnum = 0;
for i = 16:Rnum-16
    if QTwave(i,2) == -1
        X = X + QTwave(i - 15:i+15,1);
        ACnum = ACnum + 1;
    elseif QTwave(i,2) == 1
        Z = Z + QTwave(i - 15:i+15,1);
        DCnum = DCnum + 1;
    end
end
X = X / ACnum;
Z = Z / ACnum;
DC = (X(16) + X(17) - X(15) - X(14))*cell/4;      %DC的值
AC = (Z(16) + Z(17) - Z(15) - Z(14))*cell/4;      %AC的值

%-------------------------calc SDNN---------------------------
radio = zeros(Rnum,1);
for i = 1:Rnum
    radio(i) = QTwave(i,1) / R(i,1);
end

qtMEAN = mean(QTwave(:,1));
radioMEAN = mean(radio);
qtSDNN= std(QTwave(:,1),1);
radioSDNN= std(radio,1);

% subplot(3,1,1);plot(radio);
% subplot(3,1,2);plot(R(1:Rnum,1));
% subplot(3,1,3);plot(QTwave(:,1));

% off = 0;
% for i = 1:Rnum
%     if(R(i,1) == 110)
%         off = R(i,1) - QTwave(i,1);
%         break;
%     end
% end
% plot(R(1:Rnum,1));hold on;
% plot(QTwave(:,1) + off,'Color','r');

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






