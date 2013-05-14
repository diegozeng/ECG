Fs = 1;  % Sampling Frequency

N  = 2;     % Order
Fc = 0.03;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');

data=filter(Hd,R(1:Rnum,1));

plot(R(1:Rnum,2)*0.008,R(1:Rnum,1));
hold on;
plot(R(1:Rnum,2)*0.008,data,'color','r','LineWidth',3);