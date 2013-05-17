function [hrv,qtv,radiov,heart_radio_avg,r10] = TEST( Sig )
Fs = 125;               %采样频率
cell = 1000/Fs;         %时间单位

tic;
[ R,Rnum ] = RWavesDetect( Sig,cell );
disp(['Rnum: ',num2str(Rnum)]);
disp(['R detect time: ',num2str(toc)]);

if Rnum < 2000
    hrv = NaN;
    qtv = NaN;
    radiov = NaN;
    heart_radio_avg = NaN;
    r10 = NaN;
    return
end

Twave = TWavesDetect( R,Sig,Rnum );
Qwave = QWavesDetect( R,Sig,Rnum );

QTwave = zeros(Rnum,1);
QTwave(:) = Twave(:,2) - Qwave;

% PlotDetectOutcome(Sig,Twave(:,2),Qwave);
RRm = MeanFilter(R(1:Rnum,1),Rnum,20,0.2);
%R(1:Rnum,1) = medfilt1(R(1:Rnum,1),3,30); 
QTm = MeanFilter(QTwave(:),Rnum,20,0.2);
%QTwave(:) = medfilt1(QTwave(:),3,30); 

Ni = 150;


RRew = MainFilter(R(1:Rnum,1),Rnum,Ni);

%QTlow = low_pass_filter(QTwave(:));
QTew = MainFilter(QTwave(:),Rnum,Ni);
% radioew = RRew ./ QTew;

radio = RRm ./ QTm;

s = 1000:Rnum-Ni;

PlotWaves( R(:,1),RRm,RRew,QTwave,QTm,QTew,s );

% figure
% plot(RRew(s),QTew(s),'.g','Markersize',2);
% title('RR vs QT','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
% xlabel('RR','FontName','Times New Roman','FontSize',14);
% ylabel('QT','FontName','Times New Roman','FontSize',14,'Rotation',0);


%----------------WaveFeatures-----------------
hrv = WaveFeatures(RRm(s),Rnum,cell);
disp('HRV:');
disp(PrintFeatures(hrv));

qtv = WaveFeatures(QTm(s),Rnum,cell);
disp('QTV:');
disp(PrintFeatures(qtv));

radiov = WaveFeatures(radio(s),Rnum,cell);
disp('RADIOV:');
disp(PrintFeatures(radiov));

heart_radio_avg = mean(1 ./ (R(s,1))) / cell * 1000 * 60;
disp(['Heart avg radio: ',num2str(heart_radio_avg) ,' or ', num2str(1 / hrv.MEAN * 1000 * 60)]);


r10 = Regress10(QTm(s)*cell/1000,RRew(s)*cell/1000,true);


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

end









function[data] = low_pass_filter(wave)
    Fs = 1;  % Sampling Frequency
    N  = 2;     % Order
    Fc = 0.001;  % Cutoff Frequency

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
    Hd = design(h, 'butter');
    data=filter(Hd,wave);
end

function[data] = band_pass_filter(wave)
    Fs = 1;  % Sampling Frequency

    N   = 2;     % Order
    Fc1 = 0.03;  % First Cutoff Frequency
    Fc2 = 0.05;  % Second Cutoff Frequency

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
    Hd = design(h, 'butter');
    data=filter(Hd,wave);
end