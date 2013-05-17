function [ Twave ] = TWavesDetect( R,Sig,Rnum )
%Detect T wave
%   params: 
%       R    => origin R detect matrix
%       Sig  => origin signal
%       Rnum
%   return: Twave, 2xRnum matrix
%       first is begin of T wave, but useless
%       second is end of Twave
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
end

