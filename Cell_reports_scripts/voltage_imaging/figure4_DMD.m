cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\Hip_DMD\')
clear all
FS=500;
ses=dir('*.mat');
% remove 608442_fov3_4 neurons-all individual_allresults_indi.mat
%remove wide
% 14 low

trlength= [ 5000,  10000 , 5000,10000,  5000,10000, 10000,10000,10000,10000,10000,10000,10000]-9;
COH1=[];COH2=[];
rate=[];
for sesd=[1:13]
    
    load(ses(sesd).name);
    sc=0.1625*2;
    clear irm
    for id2= 1:length(allresults.roi)
        [ x y]=find(allresults.roi{id2});
        irm(:,id2)= round(mean([x , y]));end
    
    %% distance between neurons
    clear DIST
    for id2= 1:length(allresults.roi)
        for id3= 1:length(allresults.roi)
            
            DIST(id2,id3)=sqrt(sum((irm(:,id2).*sc -irm(:,id3).*sc).^2));
        end
    end
    
    j=0;
      vraw=allresults.orig_trace(:,j+1:end);
    v=allresults.trace_ws(:,j+1:end);
    clear allV  vsig
    mm=0;
    for id=1:size(v,1)
        mm=mm+1;
        Fn = 500/2;FB=[ 68 72];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        LFPg= ((filtfilt(B,A,    v(id,:))));
        vsig(mm,:)= v(id,:)-LFPg;
         LFPg= ((filtfilt(B,A,    vraw(id,:))));
        vr= vraw(id,:)-LFPg;
          Fn = 500/2;FB=[ 4 11];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        LFPg= ((filtfilt(B,A,    vr)));
        vraw(id,:)=LFPg;
        %  allV(:,id)=  angle(hilbert(LFPg));
    end
    
    lfp=[];
    
    lfp.trial{1}= vsig;
     lfp.trial2{1}=  vraw;
    lfp.time{1}= (1:size(vsig,2))./500;
    for id=1:size(vsig,1)
        lfp.label{id}= num2str(id);
    end
    

         [cwt_out,frs]=runcwtDMD(lfp, [3 180],FS,trlength(sesd));
     freq2.freq=frs;
     
          wavA = abs(squeeze(cwt_out));
      wavD = angle(squeeze(cwt_out));
  %  wavD = angle(squeeze(freq2.fourierspctrm));
    %     figure, imagesc(squeeze(wavA(3,:,:)))
    
    
    rate=[rate;sum(allresults.roaster,2)./((size(allresults.roaster,2))./500)]
    %clear   dataSFC                                   %v=   v-LFPg;
    mj=0;
    for burst=[ 0 1];
        mj=mj+1;
        
        
        SP=allresults.roaster2;  %% spike trains
        % artefact removal
        if sesd==3
            SP(:,22000:25000)=0;
        elseif sesd==2
            SP(:,11000:18000)=[];
        elseif sesd==9
            SP(:,1000:4000)=[];
        elseif sesd==13
            SP(:,1600:10000)=0;
        end
        %%
        
        clear allSNR
        
        for NN=1:size(allresults.spike_snr,1)
            Z=[];
            for x=1:size(allresults.spike_snr,2)
                Z=[Z;allresults.spike_snr{NN,x}];
            end
            allSNR(NN)=mean(Z);
        end
        
       
        %%%
         allSTA=[];
        clear allM allCM  allAM allM2
        mm=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for NN=1:size(allresults.spike_snr,1)  %Main loop
            NN
            
            G=[];  Allang=[];
            for x=1:size(allresults.spike_snr,2)
                G=[G;allresults.spike_snr{NN,x}];
            end
            
            
            j=0;
            sptrig=SP(NN,1:end-j);
            
            if length(find(sptrig))./(length(sptrig)./FS)>0.%length(find(sptrig))>4
            
           %%%% finding SS and CS
                wind=7; % 500Hz, means minimum 10ms (5) interspike interval
                tria=1;clear spikTyp
                s= find(sptrig)  ;
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
           
                
                sptrig=SP(NN,1:end-j);
                
          
                
                end
            
            %%%
           
            spindex=find( SP(NN,1:end-j));
            spikTyp5=spikTyp;
for id=1:length(spikTyp)
   z= find(spikTyp{id}==1);
sptim2=  spindex( spikTyp{id}==1 );

bursts_offets = find(diff(sptim2)>30);
bursts_onsets = find(diff(sptim2)>30)+1
spikTyp5{id}(z(bursts_offets))=2;
spikTyp5{id}(z(bursts_onsets))=3;
end

spikTyp5=spikTyp;
for id=1:length(spikTyp)
       z= find(spikTyp{id}==1);
sptim2=  spindex( spikTyp{id}==1 );
   s_order=zeros(1,length(z));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikTyp6=spikTyp;
for id=1:length(spikTyp)
   z= find(spikTyp{id}>-7);
sptim2=spindex( spikTyp{id}>-7);

bursts_onsets = find(diff(sptim2)>(30.*(FS./1000)))+1;
 ff=find(spikTyp{id}==1)-1;;
for x=1:length(bursts_onsets)
if any(bursts_onsets(x)==ff)
            found_sp=find((sptim2>sptim2(bursts_onsets(x)+1))&  sptim2<= sptim2(bursts_onsets(x)+1)+(40.*(FS./1000)) & (spikTyp{id}==1));
if length(found_sp)>0
spikTyp6{id}(z(bursts_onsets(x)+1))=2;
for jj=1:length(found_sp)
spikTyp6{id}(z(found_sp(jj)))=2+jj;
end;end
end


end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear spikTyp7
spikTyp7=spikTyp;
for tr=1:length(spikTyp)
         clear allW allW2; mm=0;  
sptim=spindex(spikTyp5{tr}==0);z=find(spikTyp5{tr}==0) ;
if ~isempty(z)
volt=zscore(lfp.trial2{tr}(1,:)-fastsmooth(lfp.trial2{tr}(1,:),1600,1,1));
volt2=zscore(lfp.trial{tr}(1,:)-fastsmooth(lfp.trial{tr}(1,:),1600,1,1));
wind=round(100.*(FS./1000));wind2=round(50.*(FS./1000));
for ind=1:length(sptim)
    if sptim(ind)-wind2 >0 & sptim(ind)+wind <length(volt)
       mm=mm+1;
      allW(:,ind)= volt(sptim(ind)-wind2:  sptim(ind)+wind);
       allW2(:,ind)= volt2(sptim(ind)-wind2:  sptim(ind)+wind); 
  
    else
          allW(:,ind)= zeros(1,length(-wind:wind2)).*NaN;
       allW2(:,ind)=  zeros(1,length(-wind:wind2)).*NaN;
   
    end
end
spADP= nanmean(allW,2);ax=(-wind2:wind).*(1000/FS);
windSA=round(25.*(FS./1000));windSB=round(5.*(FS./1000));
windS2A=round(5.*(FS./1000));windS2B=round(25.*(FS./1000));
base=mean(mean(allW2(wind2-windSA:wind2-windSB,:),1));
 A= nanmean(allW2(wind2+windS2A:wind2+windS2B,:),1)-nanmean(allW2(wind2-windSA:wind2-windSB,:),1);
 sp_height= nanmedian(nanmedian(allW(wind2+1,:),1));
 ADP=A./sp_height;avADPSS=nanmean(ADP);%ADP_s=[ADP_s,ADP];
end
     clear allW allW2; mm=0;  
%sptim=dataM.spikes2{tr}(spikTyp5{tr}==1);z=find(spikTyp5{tr}==1) ; 
sptim=spindex(spikTyp5{tr}==1);z=find(spikTyp5{tr}==1) ;

if ~isempty(z)
volt=zscore(lfp.trial2{tr}(1,:)-fastsmooth(lfp.trial2{tr}(1,:),1600,1,1));
volt2=zscore(lfp.trial{tr}(1,:)-fastsmooth(lfp.trial{tr}(1,:),1600,1,1));
wind=round(100.*(FS./1000));wind2=round(50.*(FS./1000));
for ind=1:length(sptim)
    if sptim(ind)-wind2 >0 & sptim(ind)+wind <length(volt)
       mm=mm+1;
      allW(:,ind)= volt(sptim(ind)-wind2:  sptim(ind)+wind);
       allW2(:,ind)= volt2(sptim(ind)-wind2:  sptim(ind)+wind); 
  
    else
          allW(:,ind)= zeros(1,length(-wind:wind2)).*NaN;
       allW2(:,ind)=  zeros(1,length(-wind:wind2)).*NaN;
   
    end
end
spADP= nanmean(allW,2);ax=(-wind2:wind).*(1000/FS);
windSA=round(25.*(FS./1000));windSB=round(5.*(FS./1000));
windS2A=round(5.*(FS./1000));windS2B=round(25.*(FS./1000));
base=mean(mean(allW2(wind2-windSA:wind2-windSB,:),1));
 A= nanmean(allW2(wind2+windS2A:wind2+windS2B,:),1)-nanmean(allW2(wind2-windSA:wind2-windSB,:),1);
 sp_height= nanmedian(nanmedian(allW(wind2+1,:),1));
 ADP=A./sp_height;avADP=nanmean(ADP);%ADP_b=[ADP_b,ADP];
 %%%%%%%
 if ~isnan(avADP)
z2= find(spikTyp5{tr}>0);
sorder=spikTyp5{tr}(z2);sordersel=zeros(1,length(sorder));
bonsets=find(sorder==1);
for j=1:length(bonsets)
    if j<length(bonsets)
        if ADP(j)>.15
           sordersel(bonsets(j):bonsets(j+1)-1)=1;
        end
    else
          if ADP(j)>.15
           sordersel(bonsets(j):end)=1;
        end
    end
end
spikTyp7{tr}(z2)=sordersel;
 end;end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%% SS or CS
      j=0;
       sptrig=SP(NN,1:end-j);       sptrig2=SP(NN,1:end-j);
         s= find(sptrig)  ;     s2= find(sptrig2)  ;
         sptrigCRIT=sptrig;
         sptrigCRIT(s(spikTyp{1}==1))=0;
         sptrigCRIT2=sptrig;
         sptrigCRIT2(s(spikTyp7{1}~=1 ))=0;
                if burst==1
                    sptrig(s(spikTyp7{1}~=1  ))=0;
                else
                    sptrig(s(spikTyp{1}==1))=0;
                end
  sptrig2(s2(spikTyp{1}==0  ))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lwin=100;clear STA
rwin=100; sdata{1}=sptrig;
if  sum(sdata{1})>5
[ STA]= STA_computeDMD(lfp.trial2,sdata,lwin,rwin);
else
STA=zeros(size(lfp.trial2{1},1),length(-lwin:rwin) ).*NaN;
end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if nanmean(G)>0  % G=SNR
                mm=mm+1;
    allSTA=[allSTA ;nanmean(STA(DIST(:,NN)>30,:),1)];           
                thetasel=find(freq2.freq>=4  & freq2.freq<=9);
                
                allP=[]; allC=[]; allP2=[]; 
                
                for id= 1:size(wavD,1)
                    sptrig2=SP(id,1:end-j); % spike train, non-selected
                    clear spikTyp2 % ?
                    s2= find(sptrig2)  ;
                    
                    if id~=NN & DIST(id,NN) >=40 & length(find(sptrigCRIT(1:end)))>1 & length(find(sptrigCRIT2(1:end)))>1 %& (length(find(sptrig2))./(length(sptrig2)./FS)>0.) %  & (sum(sptrig))>4 %)  %length(s2)>10
                        
                        
                        M=squeeze(exp(1i.*wavD(id,:,logical(sptrig(1:end)))));
                        Z=abs(nanmean(M,2));
                        T=   Z.^2;NT= size(M,2) ;
                        Z= (((1/(NT-1))*((T.*NT)-1)));
                       spiketim= logical(sptrig(1:end)) ;
                         M=squeeze(exp(1i.*wavD(id,:,  spiketim(end:-1:1)  )));
                        Z2=abs(nanmean(M,2));
                        T=   Z2.^2;NT2= size(M,2) ;
                        Z2= (((1/(NT2-1))*((T.*NT2)-1)));
                      
                        if NT>=10
                            allP(:,id)=Z;  allP2(:,id)=Z2;
                        else
                            allP(1:length(Z),id)=NaN; allP2(1:length(Z2),id)=NaN;
                        end
                        
                        
                        %% preferred phase
%                         A1=angle(M         );     A1=  circ_mean(A1(thetasel,:))        ;

    
                        %%%%%
                        %% ?
                        sptrig2=SP(id,1:end-j);
                        clear spikTyp2
                        s2= find(sptrig2)  ;
                        for in = 1:length(s2)
                            if length(s2) >2
                                if in==1
                                    d=(s2(in+1)-  s2(in));
                                    if d <wind;   spikTyp2{tria}(in)= 1 ;
                                    else; spikTyp2{tria}(in)= 0 ;end
                                    
                                elseif in == length(s2)
                                    d=( s2(in)-s2(in-1) );
                                    if d <wind;   spikTyp2{tria}(in)= 1 ;
                                    else; spikTyp2{tria}(in)= 0 ;end
                                else
                                    d1=(s2(in+1)-  s2(in));         d=( s2(in)-s2(in-1) );
                                    if d <wind |d1 <wind ;   spikTyp2{tria}(in)= 1 ;
                                    else; spikTyp2{tria}(in)= 0 ;end
                                    
                                end
                            else
                                spikTyp2{tria}=[];
                            end
                        end
                        
%                         if burst==1
%                             sptrig2(s2(spikTyp2{1}==1))=0;
%                         else
%                             sptrig2(s2(spikTyp2{1}==0))=0;
%                         end
                        %% cross_correlation
                        if sum(sptrig)>5  & sum(sptrig2)>5
                            [cw,lags]=xcorr(zscore(nanfastsmooth(sptrig,255,1,1)),zscore(nanfastsmooth(sptrig2,255,1,1)),100,'Coeff');
                        allC(:,id)=(cw);
                        else
                           allC(:,id)=ones(1,201).*NaN;  
                        end
                        
                    else
                        %  allP(:,id)=NaN;
                    end
                    
                end  % trials 
                
                %figure,plot(allP)
                if ~isempty(allP)
                    allP(allP==0)=NaN;allC(allC==0)=NaN;  % Remove 0 due to empty parts of the matrix
                      allP2(allP2==0)=NaN;
                    allM(1:size(allP,1),mm)=nanmean(allP,2);  % Phase locking value (PLV) adjusted for N
                    allM2(1:size(allP2,1),mm)=nanmean(allP2,2); 
                    allCM(:,mm)=nanmean(allC,2);  % cross corr
                    allN(mm)=size(allP,2);
                    
                    if burst==0
                        COH1=[COH1,   allP-allP2];
else
    COH2=[COH2,   allP-allP2];
end
                    
                end
            %    [nn1 nn2]=hist(Allang,[-pi:pi/16:pi]);
               % allAM(:,mm)=fastsmooth(nn1,2,1,1);
            end
            end   
        end  % MAIN NEURON LOOP
        
        try
            allM(allM==0)=NaN;   allM2(allM2==0)=NaN;
            dataSFC{mj,sesd}=allM-allM2; %% MAIN PLV data matrix
            allCM(allCM==0)=NaN;
            dataCORR{mj,sesd}=allCM; %% MAIN Xcorr data matrix
            
            allAM(allAM==0)=NaN;
            dataA{mj,sesd}=allAM; %% Main preferred angle matrix
            
            dataSTA{mj,sesd}=allSTA;
        end  % 
    end  % SS and CS LOOP    
end % SESSION LOOP

A=[];B=[];C=[];D=[]; E=[];F=[];
for id=1:size(dataSFC,2)
   if  ~isempty(dataSFC{1,id})&~isempty(dataSFC{2,id})
       try
               B=[B,  dataSFC{2,id}];
    A=[A,  dataSFC{1,id}]; %PLV

    
    C=[C,  dataCORR{1,id}]; %COrr
    D=[D,  dataCORR{2,id}];
       
    E=[E;  dataSTA{1,id}]; % STA (spike triggered average)
    F=[F;  dataSTA{2,id}];end
   end
end

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_figure2\'
pheight=150;
A=COH1; B=COH2;
      figure('COlor','w','Position',[300 300 200 pheight],'Renderer', 'painters'),
   fill_error_area2(freq2.freq,nanmean(B,2),nanstd(B,[],2)./sqrt(sum(~isnan(B(20,:)))), [ 0.5 0.5 0.5]);
    plot(freq2.freq,nanmean(A(:,:),2), 'b');axis tight
        fill_error_area2(freq2.freq,nanmean(A,2),nanstd(A,[],2)./sqrt(sum(~isnan(A(20,:)))), [ 0.5 0.5 0.5]);
         plot(freq2.freq,nanmean(B,2), 'r','Linewidth',1);axis tight; hold on
   plot(freq2.freq,nanmean(A,2), 'b','Linewidth',1);axis tight
 
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath ['DMD_PLV_VM' '.pdf']])
savefig(gcf, [ savepath ['DMD_PLV_VM' '.fig']])


% 
% 
% lags2=lags.*2;
%       figure('COlor','w','Position',[300 300 250 200],'Renderer', 'painters'),
%    plot(lags2,nanmean(D,2), 'r','Linewidth',2);axis tight; hold on
%    plot(lags2,nanmean(C,2), 'b','Linewidth',2);axis tight
%    fill_error_area2(lags2,nanmean(D,2),nanstd(D,[],2)./sqrt(sum(isnan(D(100,:)))), [ 0.5 0.5 0.5]);
%     plot(lags2,nanmean(C(:,:),2), 'b','Linewidth',2);axis tight
%         fill_error_area2(lags2,nanmean(C,2),nanstd(C,[],2)./sqrt(sum(isnan(C(100,:)))), [ 0.5 0.5 0.5]);
% %

 fr=freq2.freq>4 & freq2.freq<12
V1=nanmean(B(fr,:));
 V2=nanmean(A(fr,:));
 V2(isnan(V1))=[];;V1(isnan(V1))=[];;
  V1(isnan(V2))=[];;V2(isnan(V2))=[];;
  fr=freq2.freq>30 & freq2.freq<=90
V3=nanmean(B(fr,:));
 V4=nanmean(A(fr,:));
  V4(isnan(V3))=[];;V3(isnan(V3))=[];;
  V3(isnan(V4))=[];;V4(isnan(V4))=[];;

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2((V1-V2)',[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.25 0.4])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_THETA__BARplotMATCH.pdf'])
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2((V3-V4)',[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
 ylim([-0.08 0.08])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_GAM__BARplotMATCH.pdf'])



%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DMD_SFC_VM_BARplot.pdf'])
%savefig(gcf, [ savepath 'DMD_SFC_VM_BARplot.fig'])

[h,p,ci,stats] =ttest(V1-V2)
% df=1025, 

[h,p,ci,stats] =ttest(V3-V4)
% 2.2927e-08
%%%%%%%%%%%%%%%%%%%%%%%
ax=(-lwin:rwin)./FS;
     figure('COlor','w','Position',[300 300 200 pheight],'Renderer', 'painters'),
   plot(ax,nanmean(F,1), 'r');axis tight; hold on
   plot(ax,nanmean(E,1), 'b');axis tight
    fill_error_area2(ax,nanmean(F,1),nanstd(F,[],1)./sqrt(sum(~isnan(F(:,12)))), [ 0.5 0.5 0.5]);
        fill_error_area2(ax,nanmean(E,1),nanstd(E,[],1)./sqrt(sum(~isnan(E(:,12)))), [ 0.5 0.5 0.5]);
axis tight
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DMD_STA_filtt.pdf'])
savefig(gcf, [ savepath 'DMD_STA_filttt.fig'])


%        
%  figure('COlor','w','Position', [ 300 400 440 200])
%  fr=freq2.freq>3 & freq2.freq<=10; V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
%  fr=freq2.freq>30 & freq2.freq<=80; V3=nanmean(B(fr,:))'; V4=nanmean(A(fr,:))';
% fr=freq2.freq>90 & freq2.freq<=180; V5=nanmean(B(fr,:))'; V6=nanmean(A(fr,:))'
% V1mF=nanmean(V1);V2mF=nanmean(V2);  V3m=nanmean(V3);V4m=nanmean(V4);  V5m=nanmean(V5);V6m=nanmean(V6); 
%        M=[V1'  , V2'  , V3'  , V4', V5'  , V6'];
% boxplot( M   ,[ ones(length(V1),1);ones(length(V2),1).*2;ones(length(V3),1).*3;ones(length(V4),1).*4,;ones(length(V5),1).*5;ones(length(V6),1).*6],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
% hold on,            [h,p,ci,stats] = ttest2(V1,V2);
%              V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
% %plot([ 1 2],[V1m V2m], '.--k','Markersize',15);hold on, errorbar([1 2], [ V1m V2m], [V1s V2s],'.k')
%  fr=freq2.freq>30 & freq2.freq<=90;
%        V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
%              V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
% %plot([ 3 4],[V1m V2m], '.--k','Markersize',15);hold on, errorbar([3 4], [ V1m V2m], [V1s V2s],'.k')
%  fr=freq2.freq>90 & freq2.freq<=180;
%        V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
%              V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
% %plot([ 5 6],[V1m V2m], '.--k','Markersize',15);hold on, errorbar([5 6], [ V1m V2m], [V1s V2s],'.k')
% set(gca,'Xtick', [ 1 2 3 4 5 6],'Xticklabel', {'CS' ; 'SS';'CS' ; 'SS';'CS' ; 'SS'})
% xlim([0.5 6.5])

%plot([ 1 3 5],[V1mF V3m V5m], '.r','Markersize',15);
%plot([ 1 3 5]+1,[V2mF V4m V6m], '.b','Markersize',15);

% nanmean(allM)
% std(allM)./sqrt(length(allM))







