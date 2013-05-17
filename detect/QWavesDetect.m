function [ Qwave ] = QWavesDetect( R,Sig,Rnum )
%Detect Q wave
%   params: 
%       R    => origin R detect matrix
%       Sig  => origin signal
%       Rnum
%   return: Qwave
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
end

