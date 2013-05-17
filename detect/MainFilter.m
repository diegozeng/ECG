function[WAVEew] = MainFilter(wave,Rnum,Ni)
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
end