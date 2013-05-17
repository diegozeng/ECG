function [] = PlotDetectOutcome( Sig,T,Q,s )
    if nargin == 3
        s = 1:size(Q,1);
    end
    figure
    plot(Sig);
    hold on;
%     if nargin >= 4
%         R = R(s);
%         plot(R,Sig(R),'dm');
%         hold on;
%     end
    T = T(s);
    plot(T,Sig(T),'+r');
    hold on;
    Q = Q(s);
    plot(Q,Sig(Q),'+g');
end

