clear all
ISI=[]; ISI2=[];
wind2=2500;n=0;
clear allA allP
for iter=1:500
    iter
    
    clear fall fall2 gg  sig_inp
    tspan=5500;
    gg=zeros(7,tspan*100);
    II=fastsmooth(randn(1,tspan*100),10).*10;%.+-0 + sin([0.001*0.01:0.001*0.01:tspan./1000].*2*pi.*69).*1;
  
    ST=tspan.*100;%tspan.*100;
    R1=700.*100;
    R2=300.*100;
    mMOD=4;minval=-0.2;
    tstart=2000.*100;
    II2= zeros(1,ST)+minval;
    IX1=linspace(minval,mMOD,R1);
    IX2=linspace(mMOD,minval,R2);
    IX=fastsmooth([IX1,IX2],5,1,1);
    II2(tstart:+tstart+length(IX)-1)=IX;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    II2= fastsmooth(II2,15000,1,1);
    ST=tspan.*100;%tspan.*100;
    beg=0.1; endv=1;
    modA=0.5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ph=-pi;  frs=linspace(0.2,50,ST).*.01;
    FREQ=8+randn(1).*.3; %%%
    frs=linspace(FREQ,FREQ,ST).*.01;
    %frs=linspace(1000,15,ST);%.*.01;
    for id=2:ST
        ph(id)=  ph(id-1) + 0.001*(2.*pi).*frs(id);
    end
    I3=sin(ph).*0.3;

    
    %%%%%%%%%%%%%%%%
    s1=randn(1,length(I3)+10000);
    FS=1000*100;
    Fn = FS/2;FB=[ 60 80];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    I4=((filtfilt(B,A,  s1')));
    I4=zscore( I4(5000:end-4999))./2;
    s1=randn(1,length(I3));
    FS=1000*100;
    Fn = FS/2;FB=[ 5 9];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    %%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%
    gg(:,1)=[-65 0 0 0  -65  0 0];
    for ind=1:tspan*100 %1000      
        [f,varargout] = dXdT_HH_3(1,gg(:,ind),5,II2(ind)+I3(ind)+II(ind)+I4(ind),0.4);
        sig_inp(ind)=II2(ind)+I3(ind)+II(ind)./10;
        gg(:,ind+1)=gg(:,ind)+0.01*f;
        fall(ind)=gg(1,ind);
        fall2(ind)=gg(5,ind);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%5
    fall2=fastsmooth(fall,100,1,1);
    fall2=fall2(1:100:end);
    resultS=spike_detect_SNR_sim3_1K(fall2',5,5,7);
    
    clear spikTyp
    wind=10;
    for tria = 1:length(resultS.spike_idx)
        s=resultS.spike_idx{tria}  ;
        for in = 1:length(s)
            if length(s) >2
                if in==1
                    d=(s(in+1)-  s(in));
                    if d <wind;   spikTyp{tria}(in)= 1 ;
                    else; spikTyp{tria}(in)= 0 ;end
                    
                elseif in == length(s)
                    d=( s(in)-s(in-1) );
                    if d <wind;   spikTyp{tria}(in)= 1 ;
                    else; spikTyp{tria}(in)= 0 ;end
                else
                    d1=(s(in+1)-  s(in));         d=( s(in)-s(in-1) );
                    if d <wind |d1 <wind ;   spikTyp{tria}(in)= 1 ;
                    else; spikTyp{tria}(in)= 0 ;end
                    
                end
            else
                spikTyp{tria}=[];
            end
        end
    end
    
    clear spikTyp5
    spikTyp5=spikTyp;
    for id=1:length(spikTyp)
        z= find(spikTyp{id}==1);
        s_order=zeros(1,length(z));
        sptim2= resultS.spike_idx{id}( spikTyp{id}==1 );
        s_order(1)=1;bstart=1;bx=0;
        for sind=2:length(s_order)
            other_s=find(sptim2(sind)-sptim2(bstart+bx)<14  & sptim2(sind)-sptim2(bstart+bx)>0  );
            if ~isempty(other_s)
                s_order(sind)= sind-bstart+1;bx=bx+1;
            else
                bstart=sind;
                s_order(bstart) =1;bx=0;
            end
        end
        spikTyp5{id}(z)=s_order;
    end
    
    
    spikTyp6=spikTyp;
    FS=1000;
    for id=1:length(spikTyp)
        z= find(spikTyp{id}>-7);
        sptim2= resultS.spike_idx{id}( spikTyp{id}==1 );
        
        bursts_onsets = find(diff(sptim2)>(30.*(FS./1000)))+1;
        ff=find(spikTyp{id}==1)-1;;
        for x=1:length(bursts_onsets)
            if any(bursts_onsets(x)==ff)
                found_sp=find((sptim2>sptim2(bursts_onsets(x))+1)&  sptim2<= sptim2(bursts_onsets(x))+(40.*(FS./1000)));
                if length(found_sp)>0
                    spikTyp6{id}(z(bursts_onsets(x)))=2;
                    for j=length(found_sp)
                        spikTyp6{id}(z(found_sp(j)))=3;%+j;
                    end;end
            end
        end
    end
    
    
    sA=angle(hilbert(zscore(gg(7,2:100:end))));
    sA1=angle(hilbert(zscore(I3(1:100:end))));
    
    
    savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_model\'
    pheight=150;
    if iter==1
        figure('COlor','w','Position', [ 300 400 300 pheight],'Renderer', 'painters')
        plot(fall2,'k');hold on, plot(I3(1:100:end).*12+10,'g');plot(zscore(gg(6,1:100:end)).*4-10,'m')
        plot(zscore(gg(7,1:100:end)).*4-30,'b')
        %plot(frs(1:100:end).*100,'Color',[0.5 0.5 0.5])
        hold on,plot(resultS.spike_idx{1}(spikTyp{1}==1),ones(1,length(find(spikTyp{1}==1)))-10,'.r','Markersize',10)
        hold on,plot(resultS.spike_idx{1}(spikTyp{1}==0),ones(1,length(find(spikTyp{1}==0)))-11,'.b','Markersize',10)
        axis tight
        plot(II2(1:100:end).*12+0,'b')
        %plot(circ_dist(sA1,sA),'k')
        plot(sA1,'k')
        ylim([-65 35])
    end
    roast1= zeros(1,length(resultS.roaster));
    roast1(resultS.spike_idx{1}(spikTyp{1}==1))=1;
    roast2= zeros(1,length(resultS.roaster));
    roast2(resultS.spike_idx{1}(spikTyp{1}==0))=1;
    roast3= zeros(1,length(resultS.roaster));
    roast3(resultS.spike_idx{1}(spikTyp{1}==1))=1;
    %%%%%%%%%%%%%%%%
    
    sA2=angle(hilbert(zscore(II2(1:100:end))));
    Va=sA2;
    deltpeak  =zeros(1,length(Va));
    timer=50;
    for x1=1:length(Va)
        timer=timer+1;
        if Va(x1) >=0. &Va(x1) < 0.01 & timer>50
            deltpeak(x1)= 1; timer=0; end
    end
    tims=find(deltpeak);
    tims=tstart./100+800;
    
    
    
    sAX=circ_dist(sA1,pi);
    ph1_r= sAX.*roast1;
    ph2_r= sAX.*roast2;
    fg=s;
    
    for ind=1:length(tims)
        if tims(ind)>wind2  & tims(ind) +wind2 <length(roast1)
            n=n+1;
            allA(:,1,n)=roast2(tims(ind)-wind2:tims(ind)+wind2);
            allA(:,2,n)=roast1(tims(ind)-wind2:tims(ind)+wind2);
            allP(:,1,n)=  ph2_r(tims(ind)-wind2:tims(ind)+wind2);
            allP(:,2,n)=  ph1_r(tims(ind)-wind2:tims(ind)+wind2);
            
        end
    end
    
    ISI=[ISI;diff(find(roast3))'];
    ISI2=[ISI2;diff(find(roast2))'];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\model_PF\'
pheight=150;

wind=wind2;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
%plot([-wind:wind], nanmean(allA(:,1,:),3)./2,'b','Color',[0.8 0.8 0.99])
%hold on,plot([-wind:wind],nanmean(allA(:,2,:),3)./2-0.05,'r','Color',[0.99 0.8 0.8])
plot([-wind:wind],zscore(fastsmooth(nanmean(allA(:,1,:),3),200,1,1)),'b','Linewidth',2);hold on,
plot([-wind:wind],zscore(fastsmooth(nanmean(allA(:,2,:),3),200,1,1)-0.0),'r','Linewidth',2)
axis tight
xlim([ -1400 800])
   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'firing_PF.pdf'])

winsel=1500:3000;clear CM
batch_tr=10;
iter=floor(size(allP,3)./batch_tr)
clear allI
for iter1=1:iter
xt=[1+(batch_tr* (iter1-1)):(batch_tr* (iter1))]
SS=fastsmooth(nanmean(allA(winsel,1,xt),3),200,1,1).*1000;
BS=fastsmooth(nanmean(allA(winsel,2,xt),3),200,1,1).*1000;
AS=fastsmooth(nanmean(nanmean(allA(winsel,:,xt),2),3),200,1,1).*1000;
%[n1 n2]=max(fastsmooth(nanmean(nanmean(allA,2),3),200,1,1));
CM(iter1,1)=sum(winsel'.* (SS./sum(SS)))-sum(winsel'.* (AS./sum(AS)));
CM(iter1,2)=sum(winsel'.* (BS./sum(BS)))-sum(winsel'.* (AS./sum(AS)))
end
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(CM(:,2)- CM(:,1),[1.3 ],'ViolinColor', [ 0.4 0.1 0.6; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[]);ylim([-200 100])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CM_PF.pdf'])
[h,p,ci,stats] =ttest(CM(:,2),CM(:,1))

%%%%
batch_tr=10;
iter=floor(size(allP,3)./batch_tr)
clear allI
for iter1=1:iter
xt=[1+(batch_tr* (iter1-1)):(batch_tr* (iter1))]
SS=fastsmooth(nanmean(allA(:,1,xt),3),150,1,1);
BS=fastsmooth(nanmean(allA(:,2,xt),3),150,1,1);
SS=abs(SS(1:10:end-100).*1000);
BS=abs(BS(1:10:end-100).*1000);
p_prob=ones(1,length(SS))./length(SS);
 meanf= sum(p_prob.*SS')  ;
 meanfB= sum(p_prob.*BS')  ;
     allI(iter1,1)= nansum ( p_prob.*((SS')./meanf).* log2((SS')./meanf));
allI(iter1,2)= nansum ( p_prob.*((BS')./meanfB).* log2((BS')./meanfB))
end

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(allI(:,2)- allI(:,1),[1.3 ],'ViolinColor', [ 0.9 0.4 0.6; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[])
axis tight;%ylim([-0.55 0.55])
[h,p,ci,stats] =ttest(allI(:,1),allI(:,2))
 %1.5088e-05
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'INFO_PF.pdf'])


allF(1)= max(SS)./mean(SS(1:150));
allF(2)= max(BS)./mean(BS(1:150));

% figure,plot([-wind2:wind2], nanmean(allA(:,1,:),3)./2,'b','Color',[0.8 0.8 0.8])
% hold on,plot([-wind2:wind2],nanmean(allA(:,2,:),3)./2-0.1,'r','Color',[0.8 0.8 0.8])
% plot([-wind2:wind2],fastsmooth(nanmean(allA(:,1,:),3),100,1,1),'b','Linewidth',1)
% plot([-wind2:wind2],fastsmooth(nanmean(allA(:,2,:),3),100,1,1)-0.1,'r','Linewidth',1)
% axis tight





% figure
% for id=1:size(allP,3)
% plot(allP(:,1,id),'b.')
% hold on,
% plot(allP(:,1,id)+2*pi,'b.')
% end
% ylim([-pi pi+2*pi]);xlim([0 2*wind2])
% allP(allP==0)=NaN;
% for id=1:size(allP,3)
% plot(allP(:,2,id),'r.')
% hold on,plot(allP(:,2,id)+2*pi,'r.')
% end
% ylim([-pi pi+2*pi])
% xlim([0 2*wind2])


allP(allP==0)=NaN;
figure('COlor','w','Position', [ 300 400 200 150],'Renderer', 'painters')
subplot(1,1,1)
for id=1:size(allP,3)
plot(allP(:,1,id),'b.')
hold on,
plot(allP(:,1,id)+2*pi,'b.')
end
%plot(fastsmooth(nanmean(allA(:,1,:),3),100,1,1).*500+3*pi,'b','Linewidth',1)
ylim([-pi pi+1.5*pi]);%xlim([0 1*wind])
xlim([1500 3000])
fx= tims(ind)-wind2:tims(ind)+wind2;InputPF=zscore(II2(1:100:end));
,plot(InputPF(fx),'k')
plot(zscore(fastsmooth(nanmean(allA(:,1,:),3),200,1,1)),'b');xlim([1500 3000])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PH_RAST_SS_PF.pdf'])
allP(allP==0)=NaN;

figure('COlor','w','Position', [ 300 400 200 150],'Renderer', 'painters')
subplot(1,1,1)
for id=1:size(allP,3)
plot(allP(:,2,id),'r.')
hold on,
plot(allP(:,2,id)+2*pi,'r.')
end
ylim([-pi pi+1.5*pi]);%xlim([0 1*wind])
,plot(InputPF(fx),'k')
plot(zscore(fastsmooth(nanmean(allA(:,2,:),3),200,1,1)),'r');xlim([1500 3000])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PH_RAST_BS_PF.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear SA
batch_tr=10;
iter=floor(size(allP,3)./batch_tr)
for iter1=1:iter
      POS=[]; positions=1500:3000;
ANG=[];
for tr=[1+(batch_tr* (iter1-1)):(batch_tr* (iter1))]
vv=allP(positions,1,tr);
ANG=[ANG;vv(~isnan(vv))];
POS=[POS; positions(~isnan(vv))'];
end
SA(iter1)=circ_corrcl(ANG,POS);
end
clear BA
for iter1=1:iter
      POS=[]; positions=1500:3000;
ANG=[];
for tr=[1+(batch_tr* (iter1-1)):(batch_tr* (iter1))]
vv=allP(positions,2,tr);
ANG=[ANG;vv(~isnan(vv))];
POS=[POS; positions(~isnan(vv))'];
end
BA(iter1)=circ_corrcl(ANG,POS);
end;

  figure('COlor','w','Position', [ 300 400 100 150],'Renderer', 'painters') 
violinplot2((BA(1,:)- SA(1,:))',[1.3 ],'ViolinColor', [ 0.4 0.3 0.4; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[])
axis tight;%ylim([-0.55 0.55])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PHASE_PREC_CODE.pdf'])
[h,p,ci,stats] = ttest(BA(1,:),SA(1,:))
% df=49,   9.9009e-10
%%%%%%%%%%%  4.2386e-11%%%%%%%%%%%%%

% 
% 
%  [ n1 n2]=hist(ISI,[2:5:2000]); xlim([1 300])
%   [ n3 n4]=hist(ISI2,[2:5:2000]); xlim([1 300])
%   
%   figure,plot(n2,n1)
%   hold on,plot(n4,n3)
%   xlim([1 300])
%   