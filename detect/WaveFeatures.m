function[rv] = WaveFeatures(wave,Rnum,cell)
rv.MEAN = mean(wave)*cell;
rv.SDNN = std(wave)*cell;
fiveMinutes = 5000 / cell;
len = numel(wave);
F = zeros(len,2);
Fnum = 0;
i = 1;
while i < len
    s = 0;
    prei = i;
    while s < fiveMinutes
        if i >= len
            s = 0;
            break;
        end
        s = s + wave(i);
        i = i + 1;      
    end
    if s >= fiveMinutes
        Fnum = Fnum + 1;
        F(Fnum,1) = mean(wave(prei:i));
        F(Fnum,2) = std(wave(prei:i));
    end
end
rv.SDANN = std(F(1:Fnum,1))*cell;
rv.ASDNN = mean(F(1:Fnum,2))*cell;

%TI
t = tabulate(wave);
m = max(t(:,2));
rv.TI = Rnum / m;

%AC and DC
waveState = zeros(len,1);
pre = wave(1);
for i = 2:len
    if wave(i) > pre
        waveState(i) = 1;
    elseif wave(i) < pre
        waveState(i) = -1;
    else
        waveState(i) = 0;
    end
    pre = wave(i);
end

X = zeros(31,1);
Z = zeros(31,1);
ACnum = 0;
DCnum = 0;
for i = 16:len-16
    if waveState(i) == -1
        X = X + wave(i - 15:i+15);
        ACnum = ACnum + 1;
    elseif waveState(i) == 1
        Z = Z + wave(i - 15:i+15);
        DCnum = DCnum + 1;
    end
end
X = X / ACnum;
Z = Z / ACnum;
rv.DC = (X(16) + X(17) - X(15) - X(14))*cell/4;
rv.AC = (Z(16) + Z(17) - Z(15) - Z(14))*cell/4;
