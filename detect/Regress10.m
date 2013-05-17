%lin      qt = b + a * rr
%hyp      qt = b + a / rr     => qt = b + a * rr .^ -1
%par      qt = b * rr ^ a     => ln(qt) = ln(b) + a * ln(rr) 
%log      qt = b + a * ln(rr)
%shlog    qt = ln(b + a * rr) => e^qt = b + a * rr
%exp      qt = b + a * e ^ -rr
%atan     qt = b + a * arctag(rr)
%htan     qt = b + a * tgh(rr)
%ahs      qt = b + a * arcsinh(rr)
%ahc      qt = b + a * arccosh(rr+1)
function[r] = Regress10(qt,rr,plot_outcome)
lnrr = log(rr);
o = ones(size(rr));
r.lin   = regress(qt,      [o rr]);         %qt = b + a * rr
r.hyp   = regress(qt,      [o 1 ./ rr]);    %qt = b + a * rr .^ -1
r.par   = regress(log(qt), [o lnrr]);       %ln(qt) = ln(b) + a * ln(rr) 
r.par(1) = exp(r.par(1));

r.log   = regress(qt,      [o lnrr]);       %qt = b + a * ln(rr)
r.shlog = regress(exp(qt), [o rr]);         %e^qt = b + a * rr
r.exp   = regress(qt,      [o exp(-rr)]);   %qt = b + a * e ^ -rr
r.atan  = regress(qt,      [o atan(rr)]);   %qt = b + a * arctag(rr)
r.htan  = regress(qt,      [o tanh(rr)]);   %qt = b + a * tgh(rr)
r.ahs   = regress(qt,      [o asinh(rr)]);  %qt = b + a * arcsinh(rr)
r.ahc   = regress(qt,      [o acosh(rr+1)]);%qt = b + a * arccosh(rr+1)

disp(['Lin:   ',num2str(r.lin(1)),   ',',num2str(r.lin(2))]);
disp(['Hyp:   ',num2str(r.hyp(1)),   ',',num2str(r.hyp(2))]);
disp(['Par:   ',num2str(r.par(1)),   ',',num2str(r.par(2))]);
disp(['Log:   ',num2str(r.log(1)),   ',',num2str(r.log(2))]);
disp(['shlog: ',num2str(r.shlog(1)), ',',num2str(r.shlog(2))]);
disp(['exp:   ',num2str(r.exp(1)),   ',',num2str(r.exp(2))]);
disp(['atan:  ',num2str(r.atan(1)),  ',',num2str(r.atan(2))]);
disp(['htan:  ',num2str(r.htan(1)),  ',',num2str(r.htan(2))]);
disp(['ahs:   ',num2str(r.ahs(1)),   ',',num2str(r.ahs(2))]);
disp(['ahc:   ',num2str(r.ahc(1)),   ',',num2str(r.ahc(2))]);

if plot_outcome   
    figure
    plot(rr,qt,'.','Markersize',5);
    hold on;plot(rr, r.lin(1) + r.lin(2) .* rr,         'g','LineWidth',2);%qt = b + a * rr
    hold on;plot(rr, r.hyp(1) + r.hyp(2) ./ rr,         'g','LineWidth',2);%qt = b + a * rr .^ -1
    hold on;plot(rr, r.par(1) * rr .^ r.par(2),         'g','LineWidth',2);%ln(qt) = ln(b) + a * ln(rr) 
    hold on;plot(rr, r.log(1) + r.log(2) .* lnrr,       'g','LineWidth',2);%qt = b + a * ln(rr)
    hold on;plot(rr, log(r.shlog(1) + r.shlog(2) .* rr),'g','LineWidth',2);%e^qt = b + a * rr
    hold on;plot(rr, r.exp(1) + r.exp(2) .* exp(-rr),   'g','LineWidth',2);%qt = b + a * e ^ -rr
    hold on;plot(rr, r.atan(1) + r.atan(2) .* atan(rr), 'g','LineWidth',2);%qt = b + a * arctag(rr)
    hold on;plot(rr, r.htan(1) + r.htan(2) .* tanh(rr), 'g','LineWidth',2);%qt = b + a * tgh(rr)
    hold on;plot(rr, r.ahs(1) + r.ahs(2) .* asinh(rr),  'g','LineWidth',2);%qt = b + a * arcsinh(rr)
    hold on;plot(rr, r.ahc(1) + r.ahc(2) .* acosh(rr+1),'g','LineWidth',2);%qt = b + a * arccosh(rr+1)

    [mr,mi] = min(rr);
    rrend = rr(mi);
    text(rrend, r.lin(1) + r.lin(2) .* rrend,         'lin');
    text(rrend, r.hyp(1) + r.hyp(2) ./ rrend,         'hyp');
    text(rrend, r.par(1) * rrend .^ r.par(2),         'par');
    text(rrend, r.log(1) + r.log(2) .* lnrr(mi),      'log');
    text(rrend, log(r.shlog(1) + r.shlog(2) .* rrend),'shlog');
    text(rrend, r.exp(1) + r.exp(2) .* exp(-rrend),   'exp');
    text(rrend, r.atan(1) + r.atan(2) .* atan(rrend), 'atan');
    text(rrend, r.htan(1) + r.htan(2) .* tanh(rrend), 'htan');
    text(rrend, r.ahs(1) + r.ahs(2) .* asinh(rrend),  'ahs');
    text(rrend, r.ahc(1) + r.ahc(2) .* acosh(rrend+1),'ahc');

    title('Regress models','FontName','Times New Roman','FontWeight','Bold','FontSize',16);
    xlabel('RR(s)','FontName','Times New Roman','FontSize',14);
    ylabel('QT(s)','FontName','Times New Roman','FontSize',14,'Rotation',0);
end