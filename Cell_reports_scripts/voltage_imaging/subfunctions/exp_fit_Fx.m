function [ fitbaseline, coeff]=exp_fit_Fx(v,FS)
    %% fits an expoential to estimate the photobleaching
    v1=v;
    v1([FS-10:FS*2+20])= mean([ v([FS-40:FS-20 2*FS+20:2*FS+40 ])]);
    v1(1)=v1(2);
    F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3));%+ x(3)*exp(- xdata./x(4))  ;
    x0 = [mean(v1) 40 1.5   ] ;
    OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    t = (1:length(v))./FS;
    tsel=1:length(t);
    [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);
    fitbaseline=xunc(1)+xunc(2)*exp(-t./xunc(3));
    coeff=xunc;
% figure,plot(t,fit_baseline)
%       hold on,plot(t,v)
