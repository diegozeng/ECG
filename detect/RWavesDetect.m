function [ R,Rnum ] = RWavesDetect( Sig,cell )
    maxsize= size(Sig,1);
    mindistant = 40;        %��R��������������
    maxdistant = 5000/cell; %��R���������Զ����
    uprange = 3;            %��ֵ�����Ͻ�
    downrange = 0.5;        %��ֵ�����½�
    dis = 0.5;              %б�ʹ���ָ��
    tstep = 3;              %б�ʲ���
    kran = 0.6;             %��ʼ����ֵ

    S = zeros(maxsize,1);%slope

    maxrnum = 150000;
    R = zeros(maxrnum,3);%R wave

    %-----------------------------ȷ����һ��R��λ��------------------------
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

    t = S(start);%��һ��RR�����Mֵ
    m = 1;       %����ͳ��R���RR���ڵĹ�ϵ
    i = start;   %�ӵ�һ��R�����忪ʼɨ��
    n = t*tstep; %��һ��R������ĸ߶�
    while i+maxdistant < maxsize                             %�ҳ�R�㲢���RR���ڴ���A
        if S(i)/t < dis||Sig(i) < Sig(i-1)||Sig(i) < Sig(i+1)||S(i)*tstep < n*downrange||S(i)*tstep > n*uprange,%��ȥαR��
            [i] = nextR(i, maxdistant, S);
            continue;
        end
        tar = 0;
        for j = i+mindistant:i+maxdistant
            if S(j)*tstep > n*uprange,                      %�ų����ι��ߵ�������
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
        R(m, 1) = (r - i);%��m��RR���ڵĳ���
        R(m, 2) = i;%��m��RR���ڵ����
        R(m, 3) = r;%��m��RR���ڵ��յ�
        m = m + 1;
        t = S(i);
        n = t*tstep;
        i = r;
    end

    Rnum = m - 1;%RR���ڵĸ���
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
end

