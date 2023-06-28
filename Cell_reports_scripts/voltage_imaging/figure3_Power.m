

%addpath(genpath('Z:\EricLowet\'))
addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))

clear all
mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
cd(mainpath);

foldnames = {'DC'; 'DC_multi_level_trial';'place'; 'RES_L'; 'DC_L'; 'other'; 'RES_C'; 'focal_wide';'ramp';'DC_DMD'; 'hip_osc'; 'hip_osc40'} ;
clear allISI allPHASE allCOH   allBURST  lfp_q allname  allangs
allSNR=[];  allRATE=[]; LFPpre=[]; allB=[];allSPEC=[];
loc_D =[];pN=[];pB=[];
mr=0; mmm=0;   % 2  13 21  25 30  35 42 49 54
for foldn= [1]%:length(foldnames);   %% folders
    cd(mainpath)
    cd(foldnames{foldn})
    ses=dir('*.mat');
    %problem with LFP 36 37 39
   % 76 out
   %74?
    for ih=[1:length(ses)]  %% sessions    , place 2  3  4  7 14
        
        Cpath= ses(ih).name
        load(Cpath)     
        
        if isfield(result,'resultS')
        
        %% is LFP present
        if isfield(result,'lfp')
            if isfield(result.lfp,'sampling_rate')
                for neuron=1:size(result.traces,2) %length(result.REL_SP)
                    LFPpre=[LFPpre , 1];   end; else
                for neuron=1:length(result.REL_SP)
                    LFPpre=[LFPpre , 0];
                end;
            end;else
            for neuron=1:length(result.REL_SP)
                LFPpre=[LFPpre , 0];    end; end;
               
        
        %%%%%%%%%%%
        %% count neurons
        for neuron=1:size(result.traces,2) %length(result.REL_SP)
            loc_D = [loc_D ; [ foldn ih]];
        end
        
        %% check sampling rate
        if isfield(result,'lfp')
            if isfield(result.lfp,'sampling_rate')
                try;  FS=result.lfp.sampling_rate; end
                try FS=result.lfp.sampling_rate{1};  end
            else     FS=828;
            end; else    FS=828;   end;
        
        
        if       length(result.resultS)  == length(unique(result.trial_vec    ))    % check trial correspondence
            try                
                if   1%   length(result.resultS)  == length((result.lfp.trial    ))      %% optional additional check             
                    
                    %%
                    for neuron=1%:size(result.traces,2)   %neuron within a session
                        
                        SNR=[];RATE=[];  %% SNR and RATE aggregated over trials
                        for id=1:length(result.resultS)
                            if ~isempty(result.resultS{id})
                                clear A B
                                % for neuron=1:length(result.resultS{id}.spike_snr)
                                A=   nanmean(result.resultS{id}.spike_snr{neuron}) ;
                                B=sum(result.resultS{id}.roaster(neuron,:))./(size(result.resultS{id}.roaster,2)./FS);
                                %  end;
                                SNR= [ SNR,  A'   ];
                                RATE= [ RATE, B'   ];
                            end
                        end
                        
                        mmm=mmm+1;  %% overall neuron indicator !!
                        
                        try
                            
                                                                if isfield(result,'mis')
                                if ~isempty(result.mis)
                                    if length(result.mis)==1
                                        result.lfp.trial_order=[ 1: result.mis-1 result.mis+1: length(result.lfp.trial)];
                                    else
                                        %%%%%%%%%%%%%%%%%
                                    end
                                end
                                                                end
                            
                            mm=0; clear nall allF allopto
                            spang=[];dataM=[];
                            for tt=unique(result.trial_vec)   %% trial iteration!!!!
                                if ~isempty(result.resultS{tt})
                                    try
                                        SpikeTimes=result.resultS{tt}.spike_idx{neuron};%./FS;
                                     %   Vm_trace=result.traces(result.trial_vec==tt,neuron) ;
                                        Vm_trace=result.resultS{tt}.trace_ws(neuron,:);
                                          Vm_trace= Vm_trace-fastsmooth(Vm_trace,2000,1,1);
                                            Fn = FS/2;FB=[ 86 89];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,    Vm_trace)));
                                           Vm_trace=   Vm_trace-LFPg;
%                                                  Fn = FS/2;FB=[ 18 20];
%                                          [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
%                                          LFPg= ((filtfilt(B,A,    Vm_trace)));
%                                            Vm_trace=   Vm_trace-LFPg;
% %                                        
                                          
                                        Vm_trace2= result.resultS{tt}.orig_trace(neuron,:);
                                        SpikeSNR= result.resultS{tt}.spike_snr{neuron};
                                        %%%%%%%%%%%%%%%%%
                                        [n nax]=hist(diff(SpikeTimes./FS),[0.00:0.002:0.5]);   %% spike intervals distribution!!
                                        mm=mm+1; %% trial indicator
                                        nall(:,mm)= n/sum(n);  %% normalized spike intervals distribution
                                        
                                        %% get lfp trial order if present
                      
                                        tr_order=[];
                                        if isfield(result.lfp,'trial_order')
                                            if ~isempty(result.lfp.trial_order)
                                                tr_order=result.lfp.trial_order;
                                            end
                                        end
                                        %%%%%%%%%%%%%%%%
                                        %% get LFP and opto
                                        if ~isempty(tr_order)
                                            if size( result.lfp.trial{tr_order(1)},1) ==1
                                                LFPtr=result.lfp.trial{tr_order(tt)}(1,:);
                                                opto=result.lfp.opto{tr_order(tt)};
                                            else
                                                LFPtr=result.lfp.trial{tr_order(tt)}(3,:);
                                                opto=result.lfp.opto{tr_order(tt)};
                                            end
                                        else
                                            if size( result.lfp.trial{1},1) ==1
                                                LFPtr=result.lfp.trial{tt}(1,:);
                                                opto=result.lfp.opto{(tt)};
                                            else
                                                LFPtr=result.lfp.trial{tt}(3,:);
                                                opto=result.lfp.opto{(tt)};
                                            end
                                        end
                                        
                                    %     LFPtr=zscore(result.lfp.opto{(tt)});
                                          
                                        if ih==34 | ih==29  & foldn ==1
                                            LFPtr=-LFPtr;end
                                        %%%%%%%%%%%%%
                                        %% denoise
                                         Fn = FS/2;FB=[ 58 62];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                       %   LFPtr= LFPtr-LFPg;
                                          
                                               Fn = FS/2;FB=[ 117 123];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                     %   LFPtr= LFPtr-LFPg;
                                                     Fn = FS/2;FB=[ 177 183];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                       %% LFPtr= LFPtr-LFPg;
                                        
                                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            Fn = FS/2;FB=[ 30 55];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                        
                                     
                                        LFPphase= angle(hilbert(filtfilt(B,A,  LFPtr)));
                                        Vmphase= angle(hilbert(filtfilt(B,A,  Vm_trace)));
%                                         
                                        %LFPphase(700:end)=NaN;
                                        %LFPphase(end-700:end)=NaN;
                                        % LFPphase((opto> prctile(opto,50)))=NaN;
                                        opto= opto-mean(opto(1:400));opto2=opto;
                                        opto(1:400)=0; opto(2700:end)=0;
                                       optostart= find(fastsmooth(opto,1,1,1)>mean(opto));
                                       if ~isempty(optostart)
                                      optostart=optostart(1);
                                      
                                       else
                                           if FS <600   % quick solution for some missin opto trials
                                         optostart=500; 
                                           else
                                            optostart=950;    
                                           end
                                       end
                                      optostart2=   optostart;
                                       optostart=0;
                                       
                                       
%                                         if max(SpikeTimes) <= length(LFPphase)  % spike phases!!
%                                             spang=[ spang , LFPphase(SpikeTimes  )];
%                                         else
%                                             spang=[ spang , LFPphase(SpikeTimes( (SpikeTimes) <= length(LFPphase) )  )  ];
%                                         end
                                        %%%%%%%%%%%%%%%%
                                        %%%%%%%%%%
                                 
                                        
                                        
                                        
                                        if length(Vm_trace) <= length(LFPtr) 
                                      dataM.trial{mm}(1,:) = zscore(double(Vm_trace));
                                       dataM.trial{mm}(2,:) =  zscore(double(LFPtr(1:length(Vm_trace))));
                                      %       dataM.trial{mm}(1,:) = zscore(double(Vm_trace));
                                     %   dataM.trial{mm}(2,:) =  zscore(double(LFPtr(length(LFPtr)-length(Vm_trace)+1:end)))
                                        
                                        dataM.trial2{mm}(1,:) = (double(Vm_trace2));
                                        else
                                            dataM.trial{mm}(1,:) =zscore(double(Vm_trace(1:length(LFPtr))));
                                        dataM.trial{mm}(2,:) =  zscore(double(LFPtr)); 
                                          dataM.trial2{mm}(1,:) =(double(Vm_trace2(1:length(LFPtr))));
                                        end
                                       % dataM.trial{mm}(1,:)=
                                        dataM.opto{mm}= opto2; 
                                   
                                          dataM.optostart{mm}=optostart; 
                                        if length(LFPphase)>=length(Vmphase)
                                        dataM.phaseLFP{mm}= LFPphase(1:length(Vmphase)); 
                                        dataM.phaseVm{mm}= Vmphase; else
                                        dataM.phaseLFP{mm}= LFPphase; 
                                        dataM.phaseVm{mm}= Vmphase(1:length(LFPphase));
                                        end
                                          if max(SpikeTimes) <= length(LFPphase) 
                                         dataM.phasespike{mm}= LFPphase(SpikeTimes  );else
                                          dataM.phasespike{mm}=  LFPphase(SpikeTimes( (SpikeTimes) <= length(LFPphase) )  ) ;
                                          end
                                          dataM.phasespikeVm{mm}= Vmphase(SpikeTimes  );
                                        dataM.spikes{mm} = (SpikeTimes-optostart);
                                       dataM.spikeSNR{mm} =   SpikeSNR;
                                         dataM.spikes2{mm} = SpikeTimes;
                                        dataM.time{mm}= ((1:length(Vm_trace))-optostart)./FS;
                                        dataM.roaster{mm}= result.resultS{tt}.roaster(neuron,:);
                                       allopto(1:length(opto),mm)= opto2;
                                            dataM.optostartM=   optostart2;
                                    end
                                end
                            end  %% trials
                         
                              dataM.label= {'A'; 'B'};
          
      %%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     


clear spikTyp
wind=12;
for tria = 1:length(dataM.spikes)
   s= dataM.spikes{tria}  ;
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










if   1%nanmean(SNR,2) >0 & allBURSTF > 0 &  (length(A))>50 & length(B)>5  %avADP>0.15 &
if 1 


    
    

           warning off
  cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
     cfg.channel= 'all'%; %chans=cfg.channel;
    cfg.foi= [3:0.5:63];
     cfg.toi=dataM.time{2}(1:1:end) ;
     cfg.width =6;
    cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
freq2 = ft_freqanalysis(cfg, dataM);
    


%   figure,plot(freq2.freq,squeeze( nanmean(nanmean(abs(freq2.fourierspctrm(:,2,:,:)).^2,1),4)))

    end

   

 if 1
   mr=mr+1
      detG(mr)=1;
%%
                            %% putting everything in final structures %%%%%%%%
               %    allCorr(:,mr)= CC;         
                            if isfield(result.lfp,'preproc_v')
                                lfp_q(mr)= 1;
                            else
                                lfp_q(mr)= 0;
                            end 
                           
                               allFS(mr)=FS;
                               allIND{mr}=ih;
                      
                              allaM(:,1,mr)= squeeze( nanmean(nanmean(abs(freq2.fourierspctrm(:,1,:,:)).^2,1),4));%nanmean(abs(aM),3);
                              allaM(:,2,mr)= squeeze( nanmean(nanmean(abs(freq2.fourierspctrm(:,2,:,:)).^2,1),4));
                        end                                            
                       
                    end
                end
            end;end;end;end;end
    
    end
end
%% RESULT ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%
% allaM3(allaM3==0)=NaN;allaM3B(allaM3B==0)=NaN;
% allaM3R(allaM3R==0)=NaN;allaMB3R(allaMB3R==0)=NaN;
wind=150;
ax=(-wind:wind)./830;


savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_figure1\'
pheight=160;

fsel1=freq2.freq > 5 &freq2.freq<10;
fsel2=freq2.freq > 2 &freq2.freq<5;
%selC= find(squeeze(allaM(end-2,2,:)) <3);
selC= find((squeeze(nanmean(allaM(fsel1,2,:),1))./squeeze(nanmean(allaM(fsel2,2,:),1)) >1)  );
    figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
%,plot(freq2.freq,squeeze(allaM(:,1,selC)),'k');hold on,
plot(freq2.freq,nanmean(squeeze(allaM(:,2,selC)),2),'Color', [ 0.2 0.7 0.5],'Linewidth', 1);
  fill_error_area2(freq2.freq,nanmean(squeeze(allaM(:,2,selC)),2), nanstd(squeeze(allaM(:,2,selC)),[],2)./sqrt(length(selC)), [ 0.5 0.5 0.5]);
selC= find((squeeze(nanmean(allaM(fsel1,2,:),1))./squeeze(nanmean(allaM(fsel2,2,:),1)) <1)  );
plot(freq2.freq,nanmean(squeeze(allaM(:,2,selC)),2),'Color', [ 0.7 0.5 0.2],'Linewidth', 1);
  fill_error_area2(freq2.freq,nanmean(squeeze(allaM(:,2,selC)),2), nanstd(squeeze(allaM(:,2,selC)),[],2)./sqrt(length(selC)), [ 0.5 0.5 0.5]);

  axis tight;xlim([3 25])
%set(gca,'Xscale','log')


if 0
    figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
    ,subplot(1,1,1),imagesc(ax,freq2.freq, (nanmean(allaM,3).*repmat((freq2.freq.^0.5)',1, size(data,3)))    )
   axis xy;colormap(jet)
   %set(gca,'Clim',[ 10 65])
      figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
  subplot(1,1,1), imagesc(ax,freq2.freq, (nanmean(allaMB,3).*repmat((freq2.freq.^0.5)',1, size(data,3)))    )
   axis xy;colormap(jet)
 % set(gca,'Clim',[ 10 65])
end
  
  allaMB2=allaMB;%bsxfun(@minus,allaMB,nanmean(nanmean(allaMB(:,40:140,:),1),2));
    allaM2=allaM;%bsxfun(@minus,allaM,nanmean(nanmean(allaMB(:,40:140,:),1),2));
   SC=(freq2.freq.^0.5)';tsel=40:140;
 figure('COlor','w','Position',[ 300 400 250 pheight],'Renderer', 'painters'),
  plot(freq2.freq,nanmean(nanmean(allaM2(:,tsel,:),3),2).*SC,'r'); hold on,
  fill_error_area2(freq2.freq,mean(nanmean(allaM2(:,tsel,:),3),2).*SC,(nanstd(nanmean(allaM2(:,tsel,:),2).*SC,[],3))./sqrt(size(allaM,3)), [ 0.5 0.5 0.5]);
  plot(freq2.freq,nanmean(nanmean(allaMB2(:,tsel,:),3),2).*SC,'b')
    fill_error_area2(freq2.freq,mean(nanmean(allaMB2(:,tsel,:),3),2).*SC,(nanmean(nanstd(allaMB2(:,tsel,:).*SC,[],3),2))./sqrt(size(allaMB,3)), [ 0.5 0.5 0.5]);
  axis tight
  set(gca,'Xscale','log')
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Power_CHangeplot.pdf'])
savefig(gcf, [ savepath 'Power_CHangeplot.fig'])

  
  
    B=squeeze(nanmean(allaM2(:,40:140,:),2)).*SC ;A=squeeze(nanmean(allaMB2(:,40:140,:),2)).*SC;
 figure('COlor','w','Position', [ 300 400 250 pheight])
 fr=freq2.freq>3 & freq2.freq<=10; V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
 fr=freq2.freq>30 & freq2.freq<=80; V3=nanmean(B(fr,:))'; V4=nanmean(A(fr,:))';
fr=freq2.freq>90 & freq2.freq<=180; V5=nanmean(B(fr,:))'; V6=nanmean(A(fr,:))'
V1mF=nanmean(V1);V2mF=nanmean(V2);  V3m=nanmean(V3);V4m=nanmean(V4);  V5m=nanmean(V5);V6m=nanmean(V6); 
       M=[V2'  , V1' ];
boxplot( M   ,[ ones(length(V1),1);ones(length(V2),1).*2;],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
hold on,          %  [h,p,ci,stats] = ttest2(V1,V2);
set(gca,'Xtick', [ 1 2 ],'Xticklabel', {'SS';'CS'})
xlim([0.5 2.5])
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Power_CHangeBARquant.pdf'])
savefig(gcf, [ savepath 'Power_CHangeBARquant.fig'])

[h,p,ci,stats] = ttest(V1,V2)
  
  
  B=squeeze(nanmean(allaM2(:,40:140,:),2)) ;A=squeeze(nanmean(allaMB2(:,40:140,:),2));
 figure('COlor','w','Position', [ 300 400 440 200])
 fr=freq2.freq>3 & freq2.freq<=10; V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
 fr=freq2.freq>30 & freq2.freq<=80; V3=nanmean(B(fr,:))'; V4=nanmean(A(fr,:))';
fr=freq2.freq>90 & freq2.freq<=180; V5=nanmean(B(fr,:))'; V6=nanmean(A(fr,:))'
V1mF=nanmean(V1);V2mF=nanmean(V2);  V3m=nanmean(V3);V4m=nanmean(V4);  V5m=nanmean(V5);V6m=nanmean(V6); 
       M=[V1'  , V2'  , V3'  , V4', V5'  , V6'];
boxplot( M   ,[ ones(length(V1),1);ones(length(V2),1).*2;ones(length(V3),1).*3;ones(length(V4),1).*4,;ones(length(V5),1).*5;ones(length(V6),1).*6],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
hold on,            [h,p,ci,stats] = ttest2(V1,V2);
             V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
%plot([ 1 2],[V1m V2m], '.k','Markersize',15);
hold on, errorbar([1 2], [ V1m V2m], [V1s V2s],'.k')
 fr=freq2.freq>30 & freq2.freq<=90;
       V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
             V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
%plot([ 3 4],[V1m V2m], '.k','Markersize',15);
hold on, errorbar([3 4], [ V1m V2m], [V1s V2s],'.k')
 fr=freq2.freq>100 & freq2.freq<=180;
       V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
             V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
%plot([ 5 6],[V1m V2m], '.k','Markersize',15);
hold on, errorbar([5 6], [ V1m V2m], [V1s V2s],'.k')
set(gca,'Xtick', [ 1 2 3 4 5 6],'Xticklabel', {'CS' ; 'SS';'CS' ; 'SS';'CS' ; 'SS'})
xlim([0.5 6.5])
plot([ 1 3 5],[V1mF V3m V5m], '.r','Markersize',15);
plot([ 1 3 5]+1,[V2mF V4m V6m], '.b','Markersize',15);

  
   fr=freq2.freq>3 & freq2.freq<=10; V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
 fr=freq2.freq>30 & freq2.freq<=80; V3=nanmean(B(fr,:))'; V4=nanmean(A(fr,:))';
fr=freq2.freq>90 & freq2.freq<=200; V5=nanmean(B(fr,:))'; V6=nanmean(A(fr,:))'

[h,p,ci,stats] =ttest(V1,V2)
 
[h,p,ci,stats] =ttest(V3,V4)
[h,p,ci,stats] =ttest(V5,V6)


%    SC=(freq2.freq.^0.0)';tsel=20:120;
%  figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
%   plot(log(freq2.freq(30:end)),log(nanmean(nanmean(allaMB(30:end,tsel,:),3),2)),'r'); hold on,
%   axis tight

%   
  
       figure,subplot(1,1,1),imagesc(ax,freq2.freq, (nanmean((allaM- allaMB)./(allaM+allaMB),3))    )
   axis xy;colormap(jet)
 
%     
%     figure('COlor','w'),
%     subplot(2,1,1),plot(freq2.freq,nanmean(nanmean(allaM2(:,10:end,1:1:end),3),2),'Linewidth',2);axis tight
%     title('LFP-Vm coherence')
% subplot(2,1,2),
%    plot(freq2.freq,nanmean(allaM3B,2),'Linewidth',2);axis tight; hold on
%     plot(freq2.freq,nanmean(allaMB3R,2),'Linewidth',2);axis tight
%        title('SP-LFP coherence')
   
       figure('COlor','w','Renderer', 'painters'),
   plot(freq2.freq,nanmean(allaM3(:,:),2), 'r','Linewidth',2);axis tight; hold on
   plot(freq2.freq,nanmean(allaM3B(:,:),2), 'b','Linewidth',2);axis tight
   fill_error_area2(freq2.freq,nanmean(allaM3(:,:),2),nanstd(allaM3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
    plot(freq2.freq,nanmean(allaM3B(:,:),2), 'b','Linewidth',2);axis tight
        fill_error_area2(freq2.freq,nanmean(allaM3B,2),nanstd(allaM3B,[],2)./sqrt(size(allaM3B,2)), [ 0.5 0.5 0.5]);
set(gca,'Xscale','log')




%     j=23
%        figure('COlor','w','Renderer', 'painters'),
%    plot(freq2.freq,nanmean(allaM3(:,j),2), 'r','Linewidth',2);axis tight; hold on
%    plot(freq2.freq,nanmean(allaM3B(:,j),2), 'b','Linewidth',2);axis tight
%    fill_error_area2(freq2.freq,nanmean(allaM3(:,j),2),nanstd(allaM3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
%     plot(freq2.freq,nanmean(allaM3B(:,j),2), 'b','Linewidth',2);axis tight
%         fill_error_area2(freq2.freq,nanmean(allaM3B(:,j),2),nanstd(allaM3B,[],2)./sqrt(size(allaM3B,2)), [ 0.5 0.5 0.5]);
% %set(gca

      figure('COlor','w','Renderer', 'painters'),subplot(2,1,1)
   plot(freq2.freq,fastsmooth(nanmean(allaM3,2),1,1,1), 'r','Linewidth',2);axis tight; hold on
    fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3,2),1,1,1),nanstd(allaM3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
   ,subplot(2,1,2), plot(freq2.freq,fastsmooth(nanmean(allaM3B,2),1,1,1), 'b','Linewidth',2);axis tight
        fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3B,2),1,1,1),nanstd(allaM3B,[],2)./sqrt(size(allaM3B,2)), [ 0.5 0.5 0.5]);
%set(gca,'Xscale','log')
       %title('SP-Vm coherence')
       
       
       
%        
%              figure('COlor','w','Renderer', 'painters'),subplot(2,1,1)
%    plot(freq2.freq,fastsmooth(nanmean(allaM3_1,2),3,1,1), 'r','Linewidth',2,'COlor',[ 0.3 0 0 ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3_1,2),3,1,1),nanstd(allaM3_1,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
%  plot(freq2.freq,fastsmooth(nanmean(allaM3_2,2),3,1,1), 'r','Linewidth',2,'COlor',[ 0.6 0 0 ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3_2,2),3,1,1),nanstd(allaM3_2,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
%  plot(freq2.freq,fastsmooth(nanmean(allaM3_3,2),3,1,1), 'r','Linewidth',2,'COlor',[ 0.9 0 0 ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3_3,2),3,1,1),nanstd(allaM3_3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
% subplot(2,1,2)
%    plot(freq2.freq,fastsmooth(nanmean(allaM3B_1,2),3,1,1), 'r','Linewidth',2,'COlor',[0 0 0.3  ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3B_1,2),3,1,1),nanstd(allaM3B_1,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
%  plot(freq2.freq,fastsmooth(nanmean(allaM3B_2,2),3,1,1), 'r','Linewidth',2,'COlor',[0  0 0.6  ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3B_2,2),3,1,1),nanstd(allaM3B_2,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
%  plot(freq2.freq,fastsmooth(nanmean(allaM3B_3,2),3,1,1), 'r','Linewidth',2,'COlor',[ 0 0 0.9  ]);axis tight; hold on
%     fill_error_area2(freq2.freq,fastsmooth(nanmean(allaM3B_3,2),3,1,1),nanstd(allaM3B_3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
% 
% 
% %        figure('COlor','w','Position', [ 300 500 240 200])
% %        fr=freq2.freq>3 & freq2.freq<=9;
% %        V1=mean(allaM3(fr,:))'; V2=mean(allaM3B(fr,:))';
% %              [h,p,ci,stats] = ttest2(V1,V2)
% % boxplot( [V1;V2]  ,[ ones(length(V1),1);ones(length(V2),1).*2],  'notch','on','colors',[ 0.5 0 0], 'symbol','.k')
% % xlim([ 0.5 2.5]);title([ 'p= ' num2str(p)])
% %       figure('COlor','w','Position', [ 300 500 240 200])
% %        fr=freq2.freq>34 & freq2.freq<=49;
% %        V1=mean(allaM3(fr,:))'; V2=mean(allaM3B(fr,:))'; [h,p,ci,stats] = ttest2(V1,V2)
% % boxplot( [V1;V2]  ,[ ones(length(V1),1);ones(length(V2),1).*2],  'notch','on','colors',[ 0.5 0.5 0], 'symbol','.k')
% % xlim([ 0.5 2.5]);title([ 'p= ' num2str(p)])
% %   figure('COlor','w','Position', [ 300 500 240 200])
% %        fr=freq2.freq>130 & freq2.freq<=300;
% %        V1=mean(allaM3(fr,:))'; V2=mean(allaM3B(fr,:))'; [h,p,ci,stats] = ttest2(V1,V2)
% % boxplot( [V1;V2]  ,[ ones(length(V1),1);ones(length(V2),1).*2],  'notch','on','colors',[ 0 0.5 0.5], 'symbol','.k')
% % xlim([ 0.5 2.5]);title([ 'p= ' num2str(p)])
% 
       
 figure('COlor','w','Position', [ 300 500 240 200])
       fr=freq2.freq>3 & freq2.freq<=9;
       V1=nanmean(allaM3(fr,:))'; V2=nanmean(allaM3B(fr,:))';
             [h,p,ci,stats] = ttest2(V1,V2);
             V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
bar([ 1 2],[V1m V2m], 'b');hold on, errorbar([1 2], [ V1m V2m], [V1s V2s],'.k')
title([ 'p= ' num2str(p)]);set(gca,'Xtick', [ 1 2],'Xticklabel', {'Burst' ; 'Non burst'})
%  figure('COlor','w','Position', [ 300 500 240 200])
%     fr=freq2.freq>33 & freq2.freq<=60;
%        V1=nanmean(allaM3(fr,:))'; V2=nanmean(allaM3B(fr,:))';
%              [h,p,ci,stats] = ttest2(V1,V2);
%              V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
% bar([ 1 2],[V1m V2m], 'r');hold on, errorbar([1 2], [ V1m V2m], [V1s V2s],'.k')
% title([ 'p= ' num2str(p)]);set(gca,'Xtick', [ 1 2],'Xticklabel', {'Burst' ; 'Non burst'})
% figure('COlor','w','Position', [ 300 500 240 200])
%     fr=freq2.freq>160 & freq2.freq<=349;
%        V1=nanmean(allaM3(fr,:))'; V2=nanmean(allaM3B(fr,:))';
%              [h,p,ci,stats] = ttest2(V1,V2);
%              V1m=nanmean(V1);   V2m=nanmean(V2);  V1s=nanstd(V1)./sqrt(length(V1));   V2s=nanstd(V2)./sqrt(length(V2));
% bar([ 1 2],[V1m V2m], 'm');hold on, errorbar([1 2], [ V1m V2m], [V1s V2s],'.k')
% title([ 'p= ' num2str(p)]);set(gca,'Xtick', [ 1 2],'Xticklabel', {'Burst' ; 'Non burst'})
% 
% 
% 
% %            figure('COlor','w'),
% %    plot(freq2.freq,circ_mean(allaM4,[],2),'Linewidth',2);axis tight; hold on
% %     plot(freq2.freq,circ_mean(allaM4B,[],2),'Linewidth',2);axis tight
% %        title('SP-LFP coherence')
% %        
% % %        
% % fr=3;
% % figure,subplot(1,2,1),rose(allaM4(fr,allaM4(fr,:)~=0))
% % subplot(1,2,2),,rose(allaM4B(fr,allaM4B(fr,:)~=0))
% 
% %             figure('COlor','w'),
% %    plot(freq2.freq,nanmean(allaM3-allaM3R,2),'Linewidth',2);axis tight; hold on
% %     plot(freq2.freq,nanmean(allaM3B-allaMB3R,2),'Linewidth',2);axis tight
% %        title('SP-LFP coherence')
% % %        
% % %        
% %          figure('COlor','w'),
% %            plot(freq2.freq,nanmean(allaM3(:,24),2),'Linewidth',2);axis tight
% 
% %       figure,
% %      
% %      for ind=1:size(allaM,3)
% %          subplot(8,8,ind),
% %          imagesc(ax,freq2.freq, smoothn(nanmean(allaM2(:,:,ind),3).*repmat((freq2.freq.^0.)',1, size(data,3)),10)    )
% %    axis xy;colormap(jet);title(num2str(allIND{ind}));end
