
% 28
addpath(genpath('Z:\EricLowet\'))
addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))

clear all
mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
cd(mainpath);

foldnames = {'DC'; 'DC_multi_level_trial';'place'; 'RES_L'; 'DC_L'; 'other'; 'RES_C'; 'focal_wide';'ramp';'DC_DMD';'DC_level'} ;
clear allISI allPHASE allCOH   allBURST  lfp_q allname  allangs
allSNR=[];  allRATE=[]; LFPpre=[]; allB=[];allSPEC=[];allB2=[];
loc_D =[];
BTP1=[];STP1=[];BTF1=[];STF1=[];
IFmod_E=[];PLV_E=[];
IFmod_E2=[];STIFD=[];PLA_E=[];
mr=0; mmm=0;mr2=0;
for foldn= [ 11  ]  %:length(foldnames);   %% folders
    cd(mainpath)
    cd(foldnames{foldn})
    ses=dir('*.mat');

    for ih= [1:length(ses)]  %% sessions  % *1* 5 6 8 33 39 *10*  1  2 two conditions  % 32 multiple cond?  & 16?
        
        Cpath= ses(ih).name;
        load(Cpath)     
        
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
        
        
        if     1%  length(result.resultS)  == length(unique(result.trial_vec    ))    % check trial correspondence
            try                
                if   1%   length(result.resultS)  == length((result.lfp.trial    ))      %% optional additional check             
                    
                    %%
                    for neuron=1:size(result.traces,2)   %neuron within a session
                        
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
                         %   result.lfp.trial_order=[ 1:31 33:42]
                            
                                               %% get lfp trial order if present
                                        tr_order=[];
                                        if isfield(result.lfp,'trial_order')
                                            if ~isempty(result.lfp.trial_order)
                                                tr_order=result.lfp.trial_order;
                                            end
                                        end
                                             if isfield(result,'mis')
                                if ~isempty(result.mis)
                                    if length(result.mis)==1
                                        result.lfp.trial_order=[ 1: result.mis-1 result.mis+1: length(result.lfp.trial)];
                                    else
                                        %%%%%%%%%%%%%%%%%
                                    end
                                end;end
                            
                            

                            
                            mm=0; clear nall allF allopto
                            spang=[];dataM=[];
                            for tt=unique(result.trial_vec); %unique(result.trial_vec)   %% trial iteration!!!!
                                if ~isempty(result.resultS{tt})
                                    try
                                        SpikeTimes=result.resultS{tt}.spike_idx{neuron};%./FS;
                                        Vm_trace= result.resultS{tt}.trace_ws(neuron,:);
                                        
                                        
                                                 Fn = FS/2;FB=[ 86 89];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,    Vm_trace)));
                                           Vm_trace=   Vm_trace-LFPg;
                                        
                                         Vm_trace2= result.resultS{tt}.orig_trace(neuron,:);
                                        
                                        %%%%%%%%%%%%%%%%%
                                        
                                        %% get lfp trial order if present
                       
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
                                            elseif size( result.lfp.trial{1},1) ==2
                                                LFPtr=result.lfp.trial{tt}(2,:);
                                                opto=result.lfp.opto{(tt)};
                                            
                                            else
                                                LFPtr=result.lfp.trial{tt}(3,:);
                                                opto=result.lfp.opto{(tt)};
                                            end
                                        end
                                           opto3=opto;
                                          opto=result.lfp.opto{(3)};
                                        %%%%%%%%%%%%%
                                        %% extract phases
                                         Fn = FS/2;FB=[ 58 62];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                          LFPtr= LFPtr-LFPg;
                                           
                                         Fn = FS/2;FB=[ 117 123];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                        LFPtr= LFPtr-LFPg;
                                                     Fn = FS/2;FB=[ 177 183];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                        LFPtr= LFPtr-LFPg;
                                                            Fn = FS/2;FB=[ 238 243];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                        LFPtr= LFPtr-LFPg;
                                                   Fn = FS/2;FB=[ 219 223];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                        LFPtr= LFPtr-LFPg;
                                     
                                         Fn = FS/2;FB=[ 4 9];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                        LFPphase= angle(hilbert(filtfilt(B,A,  Vm_trace)));
%                                         LFPphase2= angle(hilbert(filtfilt(B,A,  Vm_trace)));
%                                         
                                        %LFPphase(700:end)=NaN;
                                        %LFPphase(end-700:end)=NaN;
                                        % LFPphase((opto> prctile(opto,50)))=NaN;
                                        opto= opto-mean(opto(1:400));opto2=opto;
                                        opto(1:400)=0;
                                       optostart= find(medfilt1(opto,3,1)>mean(opto));
                                       if ~isempty(optostart)
                                      optostart=optostart(1);
                                      
                                       else
                                           if FS <600   % quick solution for some missin opto trials
                                         optostart=500; 
                                           else
                                            optostart=950;    
                                           end
                                       end
                                       
                                       if foldn==10
                                           optostart=828;
                                       end
                                        %optostart=408
                                        
                                        
                                        %                                         if max(SpikeTimes) <= length(LFPphase)  % spike phases!!
                                        %                                             spang=[ spang , LFPphase(SpikeTimes  )];
                                        %                                         else
                                        %                                             spang=[ spang , LFPphase(SpikeTimes( (SpikeTimes) <= length(LFPphase) )  )  ];
                                        %                                         end
                                        
                                        [n nax]=hist(diff(SpikeTimes(SpikeTimes>optostart & SpikeTimes<optostart+FS.*1.4  )./FS),[0.00:0.002:0.5]);   %% spike intervals distribution!!
                                        mm=mm+1; %% trial indicator
                                        nall(:,mm)= n/sum(n);  %% normalized spike intervals distribution
                                        
                                         [n nax]=hist(diff(SpikeTimes(SpikeTimes>0 & SpikeTimes<optostart  )./FS),[0.00:0.002:0.5]);   %% spike intervals distribution!!
                                      %  mm=mm+1; %% trial indicator
                                        nall2(:,mm)= n/sum(n);  %% normalized spike intervals distribution
                                        
                                        %%%%%%%%%%%%%%%%
                                        %%%%%%%%%%
                                        if length(Vm_trace) <= length(LFPtr)
                                            dataM.trial{mm}(1,:) = (double(Vm_trace));
                                            dataM.trial{mm}(2,:) =  zscore(double(LFPtr(1:length(Vm_trace))));
                                            dataM.trial2{mm}(1,:) = zscore(double(Vm_trace2));
                                            
                                        else
                                            dataM.trial{mm}(1,:) = (double(Vm_trace(1:length(LFPtr))));
                                        dataM.trial{mm}(2,:) =  zscore(double(LFPtr));  
                                               dataM.trial2{mm}(1,:) = zscore(double(Vm_trace2(1:length(LFPtr))));
                                 
                                        end
                                       % dataM.trial{mm}(1,:)=
                                        dataM.opto{mm}=  opto3; 
                                        dataM.optostart=optostart;
                                        dataM.phaseVm{mm}= LFPphase; 
                                        dataM.spikes{mm} = SpikeTimes-optostart;
                                          dataM.spikes2{mm} = SpikeTimes;
                                        dataM.time{mm}= ((1:length(Vm_trace))-optostart)./FS;
                                       allopto(1:length(opto3),mm)=  opto3;
                                        dataM.roaster{mm}= result.resultS{tt}.roaster(neuron,:);
                                    end
                                end
                            end  %% trials
                       
                              dataM.label= {'A'; 'B'};
      
      %%    

     lW=3.5;
hW=9;
 lW2=7;
  hW2=17;
Fa=[];
for id=1:length(dataM.trial)
FF=ft_preproc_bandpassfilter(zscore(dataM.trial{id}(1,:)),FS,[  lW2 hW2]);FF=abs(hilbert(FF));
Fa=[Fa, FF];
end
mF=mean(Fa); mS=std(Fa);
%  
      
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
% 
clear spikTyp2
for tria = 1:1:length(dataM.spikes)
   s= dataM.spikes2{tria}  ;
   FF=ft_preproc_bandpassfilter(zscore(dataM.trial{tria}(1,:)),FS,[ lW2 hW2]);FF=abs(hilbert(FF));
Fref=(FF-mF)./mS;
  for in = 1:length(s)
      if Fref(s(in)) >0
   spikTyp2{tria}(in)= 1  ;
      else
     spikTyp2{tria}(in)= 0   ;   
      end
               
  end
end

clear spikTyp5
spikTyp5=spikTyp;
for id=1:length(spikTyp)
   z= find(spikTyp{id}==1);
   s_order=zeros(1,length(z));
sptim2=dataM.spikes{id}( spikTyp{id}==1 );
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

for id=1:length(spikTyp)
   z= find(spikTyp{id}>-7);
sptim2=dataM.spikes{id}( spikTyp{id}>-7);

bursts_onsets = find(diff(sptim2)>(30.*(FS./1000)))+1;
 ff=find(spikTyp{id}==1)-1;;
for x=1:length(bursts_onsets)
if any(bursts_onsets(x)==ff)
            found_sp=find((sptim2>sptim2(bursts_onsets(x)+1))&  sptim2<= sptim2(bursts_onsets(x)+1)+(40.*(FS./1000)) & (spikTyp{id}==1)');
if length(found_sp)>0
spikTyp6{id}(z(bursts_onsets(x)+1))=2;
for j=1:length(found_sp)
spikTyp6{id}(z(found_sp(j)))=2+j;
end;end
end


end
end


%%%%%%%%%%%%%
clear spikTyp7
spikTyp7=spikTyp;
for tr=1:length(dataM.spikes2)
     clear allW allW2; mm=0;  
sptim=dataM.spikes2{tr}(spikTyp5{tr}==1);z=find(spikTyp5{tr}==1) ; 
if ~isempty(z)
volt=zscore(dataM.trial2{tr}(1,:)-fastsmooth(dataM.trial2{tr}(1,:),1600,1,1));
volt2=zscore(dataM.trial{tr}(1,:)-fastsmooth(dataM.trial{tr}(1,:),1600,1,1));
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
 A= mean(allW2(wind2+windS2A:wind2+windS2B,:),1)-mean(allW2(wind2-windSA:wind2-windSB,:),1);
 sp_height= median(median(allW(wind2+1,:),1));
 ADP=A./sp_height;avADP=nanmean(ADP);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ik=1:size(dataM.roaster,2)
    ROAST1=zeros(1,length(dataM.roaster{ik}));try
ROAST1(dataM.spikes2{ik}( logical((spikTyp7{ik}==1   ) ) ))=1;
 % ROAST1(dataM.spikes2{ik}( logical((spikTyp{ik}==1  ) ) ))=1;  
%ROAST1(dataM.spikes2{ik}( logical((spikTyp{ik}==1   &  (spikTyp7{ik}>0.15  &   spikTyp7{ik}<300 ) )) ))=1;
dataM.roasterT1{ik}=ROAST1;end
end

for ik=1:size(dataM.roaster,2)
    ROAST1=zeros(1,length(dataM.roaster{ik}));try
        ROAST1(dataM.spikes2{ik}( logical((spikTyp{ik}==0  ) ) ))=1; 
% ROAST1(dataM.spikes2{ik}( logical(~(spikTyp{ik}==1   &  spikTyp2{ik}>0) ) ))=1;

%ROAST1(dataM.spikes2{ik}(logical(~(spikTyp{ik}==1   &  (spikTyp7{ik}>0.15  &   spikTyp7{ik}<300 )) ) ))=1;
dataM.roasterT2{ik}=ROAST1;end
end


  
                  ISI=nanmean(nall,2);
                            FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                 
                               allBURST5= sum(ISI(FI))./sum(ISI(FL));

 
baseN=[]; stimN=[];
for id2= 1:length( dataM.trial)
  
 baseN=  [  baseN ;find(dataM.spikes{id2}./FS>-0.8  & (dataM.spikes{id2}./FS)< -0.1)];
 stimN= [stimN; find((dataM.spikes{id2}./FS)>0.1  & (dataM.spikes{id2}./FS)< 1)];

end

 mr2=mr2+1

 A=[]; B=[];
for id2=1
nn=0;clear allF
for ind=1:1:length(dataM.spikes)
nn=nn+1;try
t= ((dataM.spikes{ind}./1000)>-3);
A=[A,dataM.spikes{ind}(spikTyp7{ind}==1 &  t')'];
B=[B,dataM.spikes{ind}(spikTyp7{ind}==0 &  t')'];
end
end
end
%%%%%%%%%%%%%%%
 
 
 
 
if  nanmean(SNR,2) >4 & FS>700  & length(A)>=10 &  length(B)>=10 %  &     allBURST5>3.5%
  
   NT=length(dataM.trial);
   
   clear M1 M2 
   for i=1:NT
   M1(1:length(dataM.trial{i}),i)=dataM.trial{i}(1,:);
    M2(1:length(dataM.opto{i}),i)=dataM.opto{i}(1,:);
   end
   
   
Basesel=dataM.time{3} <0 & dataM.time{3}>-0.7;
Stimsel=dataM.time{3} >0.1 & dataM.time{3}<1.4;

    Stopto=nanmean(allopto(Stimsel,:),1);
     Bopto =  nanmean(allopto(Basesel,:),1);
   
 optolevels=zscore(Stopto-Bopto);
 
%  [sP, C]=kmeans( optolevels',3);
%  [C1, C2]=sort(C);
%  nf=0;clear cond_tr
%  for cond=C2'
%      nf=nf+1;
%   cond_tr{nf}=find(sP==cond);
%  end
 
 if 1
cond_tr{1}=[1:3:NT];
 cond_tr{2}=[2:3:NT];cond_tr{3}=[3:3:NT];
end

  if ih==6
cond_tr{1}=[1:3:40 41:3:NT];
 cond_tr{2}=[1:3:40 42:3:NT];cond_tr{3}=[1:3:40 43:3:NT];
end


               warning off
  cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
     cfg.channel= 'all'%; %chans=cfg.channel;
    cfg.foi= [3:3:180];
     cfg.toi=dataM.time{2}(1:1:end) ;
     cfg.width =5;
    cfg.t_ftimwin =[ones(1,length(cfg.foi))*0.2];
freq2 = ft_freqanalysis(cfg, dataM);

%    [cwt_out,frs]=runcwt(dataM, [4 180],FS);
%      freq2.freq=frs;
% 
%      wavA = abs(squeeze(cwt_out));
%       wavD = angle(squeeze(cwt_out));
 wavD = angle(squeeze(freq2.fourierspctrm));
 wavA = abs(squeeze(freq2.fourierspctrm));


cond_tr{1}=[1:length(dataM.roasterT1)];
cond_tr{2}=[1:length(dataM.roasterT1)];
cond_tr{3}=[1:length(dataM.roasterT1)];
  CH=1;
     %%%%%%%%%%%%
     [PLV_b1]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{1},1,FS,dataM.optostart) ;
       [PLV_b2]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{2},1,FS,dataM.optostart) ;
       [PLV_b3]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{3},1,FS,dataM.optostart) ;

         [PLV_s1]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{1},1,FS,dataM.optostart) ;
       [PLV_s2]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{2},1,FS,dataM.optostart) ;
       [PLV_s3]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{3},1,FS,dataM.optostart) ;
    
       [PLVB_b1]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{1},0,FS,dataM.optostart) ;
       [PLVB_b2]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{2},0,FS,dataM.optostart) ;
       [PLVB_b3]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{3},0,FS,dataM.optostart) ;

              [PLVB_s1]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{1},0,FS,dataM.optostart) ;
             [PLVB_s2]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{2},0,FS,dataM.optostart) ;
       [PLVB_s3]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{3},0,FS,dataM.optostart) ;

       
       
  CH=2;
     %%%%%%%%%%%%
     [PLV_b1L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{1},1,FS,dataM.optostart) ;
       [PLV_b2L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{2},1,FS,dataM.optostart) ;
       [PLV_b3L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{3},1,FS,dataM.optostart) ;

         [PLV_s1L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{1},1,FS,dataM.optostart) ;
       [PLV_s2L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{2},1,FS,dataM.optostart) ;
       [PLV_s3L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{3},1,FS,dataM.optostart) ;

       
           [PLVB_b1L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{1},0,FS,dataM.optostart) ;
       [PLVB_b2L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{2},0,FS,dataM.optostart) ;
       [PLVB_b3L]= spike_field_ppc_condopto(wavD,dataM.roasterT1,CH,cond_tr{3},0,FS,dataM.optostart) ;

              [PLVB_s1L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{1},0,FS,dataM.optostart) ;
             [PLVB_s2L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{2},0,FS,dataM.optostart) ;
       [PLVB_s3L]= spike_field_ppc_condopto(wavD,dataM.roasterT2,CH,cond_tr{3},0,FS,dataM.optostart) ;

if  0
  figure,plot(PLV_b1); hold on,  
    plot(PLV_b2); plot(PLV_b3)
    
     figure,plot(PLV_s1); hold on,  
    plot(PLV_s2); plot(PLV_s3)
    
end

Vm1= nanmean(M1(:,cond_tr{1}),2);
Vm2= nanmean(M1(:,cond_tr{2}),2);
Vm3= nanmean(M1(:,cond_tr{3}),2);


StimO=dataM.time{3} >0 & dataM.time{3}<0.2;
if ih<=5 | ih==20 | ih==27  | ih==23  | ih==26  %verify 23
Vml=[mean(Vm1(StimO)) mean(Vm2(StimO)) mean(Vm3(StimO))];
[Z1 Z2]=sort(Vml);
clear cond_tr2
for u=1:length(Z2)
cond_tr2{u}=cond_tr{Z2(u)};
end
cond_tr=cond_tr2;
end

PowM1V= squeeze(nanmean(wavA(:,1,:,:),1));%;squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{1},1,:,:)),1));
PowM2V= squeeze(nanmean(wavA(:,1,:,:),1));%squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{2},1,:,:)),1));
PowM3V=squeeze(nanmean(wavA(:,1,:,:),1));% squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{3},1,:,:)),1));

PowM1L= squeeze(nanmean(wavA(:,2,:,:),1));%squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{1},2,:,:)),1));
PowM2L= squeeze(nanmean(wavA(:,2,:,:),1));%squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{2},2,:,:)),1));
PowM3L= squeeze(nanmean(wavA(:,2,:,:),1));%squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{3},2,:,:)),1));

%PowM3= squeeze(nanmean(abs(freq2.fourierspctrm( cond_tr{3},2,:,:)),1));




%figure,imagesc(PowM3-PowM1);colormap(jet);axis xy

% %%
Vm1= nanmean(M1(:,cond_tr{1}),2);
Vm2= nanmean(M1(:,cond_tr{2}),2);
Vm3= nanmean(M1(:,cond_tr{3}),2);

timed= -0.8:0.005:1.7;%  cfg.toi;%dataM.time{3}(1:1:end) ;%-0.8:0.001:1.78;
clear spM
for hj=1:length(dataM.roaster)
    for op=1:length(timed)-1
hh=dataM.roaster{hj}(find(dataM.time{hj} >=timed(op)   & dataM.time{hj} <timed(op+1) ));
spM(op,hj)=nanmean(hh);
    end
end



clear spA
for hj=1:length(dataM.roasterT1)
    for op=1:length(timed)-1
        if ~isempty(dataM.roasterT1{hj})
hh=dataM.roasterT1{hj}(find(dataM.time{hj} >=timed(op)   & dataM.time{hj} <timed(op+1) ));
spA(op,hj)=nanmean(hh);end
    end
end

clear spB
for hj=1:length(dataM.roasterT2)
   for op=1:length(timed)-1
         if ~isempty(dataM.roasterT2{hj})
hh=dataM.roasterT2{hj}(find(dataM.time{hj} >=timed(op)   & dataM.time{hj} <timed(op+1) ));
spB(op,hj)=nanmean(hh);end
    end
end






%%
if 1

    mr=mr+1

                            %% putting everything in final structures %%%%%%%%
                   detG(mr)=1;                
                            if isfield(result.lfp,'preproc_v')
                                lfp_q(mr)= 1;
                            else
                                lfp_q(mr)= 0;
                            end 
                            allN(mr)= length(spang);
                            ISI=nanmean(nall,2);
                            FI=(nax> 0.002 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                            allISI(:,mr)= ISI(2:end-1);
                             ISI2=nanmean(nall2,2);
                           % FI=(nax> 0.002 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                           % allISI(:,mr)= ISI(2:end-1);
                               allBURST(:,mr)= sum(ISI(FI))./sum(ISI(FL));
                            [nn nn2]= hist(spang,30);
                          %  allSPEC(:,mr)= fastsmooth(nn./sum(nn),3,1,1);
                               allFS(mr)=FS;
                      %    allFreqs{mr} = freq1;
                            allname{mr}=  Cpath;
                            %allangs{mr}= spang;
                            alloptoM{mr}= allopto;
                            allTFR{mr}=   squeeze(nanmean(wavA(:,1,:,:),1));
                             allTFR2{mr}=   squeeze(nanmean(wavA(:,2,:,:),1));
                            %% SFC
allSFC_VM_B(:,1,mr)= PLV_b1;
 allSFC_VM_B(:,2,mr)= PLV_b2;
 allSFC_VM_B(:,3,mr)= PLV_b3;
allSFC_VM_S(:,1,mr)= PLV_s1;
 allSFC_VM_S(:,2,mr)= PLV_s2;
 allSFC_VM_S(:,3,mr)= PLV_s3;
 
 allSFC_LFP_B(:,1,mr)= PLV_b1L;
 allSFC_LFP_B(:,2,mr)= PLV_b2L;
 allSFC_LFP_B(:,3,mr)= PLV_b3L;
allSFC_LFP_S(:,1,mr)= PLV_s1L;
 allSFC_LFP_S(:,2,mr)= PLV_s2L;
 allSFC_LFP_S(:,3,mr)= PLV_s3L;
 
 %%%%% BASELINE SFC
 allSFCB_VM_B(:,1,mr)= PLVB_b1;
 allSFCB_VM_B(:,2,mr)= PLVB_b2;
 allSFCB_VM_B(:,3,mr)= PLVB_b3;
allSFCB_VM_S(:,1,mr)= PLVB_s1;
 allSFCB_VM_S(:,2,mr)= PLVB_s2;
 allSFCB_VM_S(:,3,mr)= PLVB_s3;
 
 allSFCB_LFP_B(:,1,mr)= PLVB_b1L;
 allSFCB_LFP_B(:,2,mr)= PLVB_b2L;
 allSFCB_LFP_B(:,3,mr)= PLVB_b3L;
allSFCB_LFP_S(:,1,mr)= PLVB_s1L;
 allSFCB_LFP_S(:,2,mr)= PLVB_s2L;
 allSFCB_LFP_S(:,3,mr)= PLVB_s3L;
 
 %
%5%%
spM(isnan(spM))=0;spA(isnan(spA))=0;spB(isnan(spB))=0;
                        sm=10;
                            v1=nanfastsmooth(nanmean(spM(:, cond_tr{1}),2),sm,3,1).*FS;
                            allSPm(1:length(v1),mr,1)= v1;
                                v1=nanfastsmooth(nanmean(spM(:, cond_tr{2}),2),sm,3,1).*FS;
                            allSPm(1:length(v1),mr,2)= v1;
                                v1=nanfastsmooth(nanmean(spM(:, cond_tr{3}),2),sm,3,1).*FS;
                            allSPm(1:length(v1),mr,3)= v1;
                            
                            %%
                               v1=nanfastsmooth(nanmean(spA(:, cond_tr{1}),2),sm,3,1).*FS;
                            allSPmA(1:length(v1(5:end-5)),mr,1)= v1(5:end-5);
                                v1=nanfastsmooth(nanmean(spA(:, cond_tr{2}),2),sm,3,1).*FS;
                            allSPmA(1:length(v1(5:end-5)),mr,2)= v1(5:end-5);
                                v1=nanfastsmooth(nanmean(spA(:, cond_tr{3}),2),sm,3,1).*FS;
                            allSPmA(1:length(v1(5:end-5)),mr,3)= v1(5:end-5);
                            
                                        v1=nanfastsmooth(nanmean(spB(:, cond_tr{1}),2),sm,3,1).*FS;
                            allSPmB(1:length(v1(5:end-5)),mr,1)= v1(5:end-5);
                                v1=nanfastsmooth(nanmean(spB(:, cond_tr{2}),2),sm,3,1).*FS;
                            allSPmB(1:length(v1(5:end-5)),mr,2)= v1(5:end-5);
                                v1=nanfastsmooth(nanmean(spB(:, cond_tr{3}),2),sm,3,1).*FS;
                            allSPmB(1:length(v1(5:end-5)),mr,3)= v1(5:end-5);
                            
                            
                            allVm(1:length(Vm1),1,mr)=Vm1;
                              allVm(1:length(Vm1),2,mr)=Vm2;
                                allVm(1:length(Vm1),3,mr)=Vm3;
                                
 
                              firing_ratio(mr)= nanmean(stimN)./nanmean(baseN);
  allIND(mr)=ih;allSNR(mr)=nanmean(SNR,2);
                            
                                optostarts(mr)=optostart;                                    
                        allB= [allB ;sum(ISI(FI))./sum(ISI(FL))];
                           allB2= [allB2 ;sum(ISI2(FI))./sum(ISI2(FL))];
                     %   allSNR= [ allSNR; nanmean(SNR,2)];
                        allRATE= [ allRATE; nanmean(RATE,2)];
                        %%%%%%%%%%%%
                        allPowM1V(:,1:size(PowM1V,2),mr)= PowM1V;
                         allPowM2V(:,1:size(PowM1V,2),mr)= PowM2V;
                         allPowM3V(:,1:size(PowM1V,2),mr)= PowM3V;
                           allPowM1L(:,1:size(PowM1V,2),mr)= PowM1L;
                         allPowM2L(:,1:size(PowM1V,2),mr)= PowM2L;
                         allPowM3L(:,1:size(PowM1V,2),mr)= PowM3L;
                        
%%%%%%%%%%%%%%%%%%%%%
mr
mr

end
                    end
                end
            end;end;end;end;end
    
end
%% RESULT ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_fig_3\'
pheight=150;

V1=nanmean(nanmean(allSFC_VM_B(:,:,:),2),3);
V1s=nanstd(nanmean(allSFC_VM_B(:,:,:),2),[],3)./sqrt(size(allSFC_VM_B,3));;
V2=nanmean(nanmean(allSFCB_VM_B,2),3);
V2s=nanstd(nanmean(allSFCB_VM_B,2),[],3)./sqrt(size(allSFC_VM_B,3));;
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
,plot(freq2.freq,V1,'r','Linewidth',1.5); 
hold on,plot(freq2.freq,V2,'Color', [ 0.5 0.5 0.5]);
fill_error_area2(freq2.freq,V2,V2s,[ 0.5 0.5 0.5])
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
axis tight
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_SFC_VM_CS.pdf'])
savefig(gcf, [ savepath 'DCopto_SFC_VM_CS.fig'])
 %%
V1=nanmean(nanmean(allSFC_VM_S(:,:,:),2),3);
V1s=nanstd(nanmean(allSFC_VM_S(:,:,:),2),[],3)./sqrt(size(allSFC_VM_S,3));;
V2=nanmean(nanmean(allSFCB_VM_S,2),3);
V2s=nanstd(nanmean(allSFCB_VM_S,2),[],3)./sqrt(size(allSFC_VM_S,3));;
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
,plot(freq2.freq,V1,'b','Linewidth',1.5); 
hold on,plot(freq2.freq,V2,'Color', [ 0.5 0.5 0.5]);
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
fill_error_area2(freq2.freq,V2,V2s,[ 0.5 0.5 0.5])
axis tight
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_SFC_VM_SS.pdf'])
savefig(gcf, [ savepath 'DCopto_SFC_VM_SS.fig'])
%%%%
 





%%%%
V1=nanmean(nanmean(allSFC_LFP_B(:,3,:),2),3);
V1s=nanstd(nanmean(allSFC_LFP_B(:,3,:),2),[],3)./sqrt(size(allSFC_VM_S,3));;
V2=nanmean(nanmean(allSFCB_LFP_B(:,3,:),2),3);
V2s=nanstd(nanmean(allSFCB_LFP_B,2),[],3)./sqrt(size(allSFC_VM_S,3));;
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
fill_error_area2(freq2.freq,V2,V2s,[ 0.5 0.5 0.5])
,plot(freq2.freq,V1,'r','Linewidth',1); 
hold on,plot(freq2.freq,V2,'Linewidth',1,'Color', [ 0.5 0.5 0.5]);
axis tight;ylim([-0.02 0.15])
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_CS.pdf'])


V1=nanmean(nanmean(allSFC_LFP_S(:,1,:),2),3);
V1s=nanstd(nanmean(allSFC_LFP_S(:,1,:),2),[],3)./sqrt(size(allSFC_VM_S,3));;
V2=nanmean(nanmean(allSFCB_LFP_S,2),3);
V2s=nanstd(nanmean(allSFCB_LFP_S,2),[],3)./sqrt(size(allSFC_VM_S,3));;
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
fill_error_area2(freq2.freq,V2,V2s,[ 0.5 0.5 0.5])
,plot(freq2.freq,V1,'b','Linewidth',1); 
hold on,plot(freq2.freq,V2,'Linewidth',1,'Color', [ 0.5 0.5 0.5]);
axis tight;ylim([-0.02 0.15])
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_SS.pdf'])
ylim([-0.01 0.04])
xlim([ 20 70])

 figure('COlor','w','Position',[300 300 90 90],'Renderer', 'painters')
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
fill_error_area2(freq2.freq,V2,V2s,[ 0.5 0.5 0.5])
,plot(freq2.freq,V1,'b','Linewidth',1); 
hold on,plot(freq2.freq,V2,'Linewidth',1,'Color', [ 0.5 0.5 0.5]);
axis tight;
set(gca,'Xscale','log')
ylim([-0.008 0.025])
xlim([ 20 70])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_SSzoom.pdf'])

 fr=freq2.freq>30 & freq2.freq<90
V1=squeeze(nanmean(nanmean(allSFC_LFP_B(fr,1,:),2),1))';
 V2=squeeze(nanmean(nanmean(allSFCB_LFP_B(fr,1,:),2),1))';
V3=squeeze(nanmean(nanmean(allSFC_LFP_S(fr,1,:),2),1))';
 V4=squeeze(nanmean(nanmean(allSFCB_LFP_S(fr,1,:),2),1))';
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
violinplot2((V1-V2)',[1.25 ],'ViolinColor', [ 0.9 0.3 0.1; 0 0 0.9])
hold on,
violinplot2((V3-V4)',[1.85 ],'ViolinColor', [0.2 0.3 0.9; 0 0 0.9])
line([ 0.8 2.3], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5);axis tight
set(gca,'Xtick',[],'Xticklabel',[])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_GAM__violinSS.pdf'])
[h,p,ci,stats] = ttest(V1-V2)
% df=17,  0.9472
[h,p,ci,stats] = ttest(V3-V4)
% 0.0475

 fr=freq2.freq>3 & freq2.freq<12
V1=squeeze(nanmean(nanmean(allSFC_LFP_B(fr,1,:),2),1))';
 V2=squeeze(nanmean(nanmean(allSFCB_LFP_B(fr,1,:),2),1))';
 %fr=freq2.freq>30 & freq2.freq<90
V3=squeeze(nanmean(nanmean(allSFC_LFP_S(fr,1,:),2),1))';
 V4=squeeze(nanmean(nanmean(allSFCB_LFP_S(fr,1,:),2),1))';
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
violinplot2((V1-V2)',[1.25 ],'ViolinColor', [ 0.9 0.3 0.1; 0 0 0.9])
violinplot2((V3-V4)',[1.85 ],'ViolinColor', [ 0.2 0.3 0.9; 0 0 0.9])
line([ 0.8 2.3], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5);axis tight
set(gca,'Xtick',[],'Xticklabel',[])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_GAM__violinCS.pdf'])
[h,p,ci,stats] = ttest(V1-V2)
% df=17,  0.0343
[h,p,ci,stats] = ttest(V3-V4)
% df=25;  0.0159
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

  allPowM1V(  allPowM1V==0)=NaN;
 allPowM2V(  allPowM2V==0)=NaN;
 allPowM3V(  allPowM3V==0)=NaN;
  allPowM1L(  allPowM1L==0)=NaN;
 allPowM2L(  allPowM2L==0)=NaN;
 allPowM3L(  allPowM3L==0)=NaN;
 
 allSPmA(allSPmA==0)=NaN;
  allSPmB(allSPmB==0)=NaN;
 %figure,plot(squeeze(nanmean(allSPm,2)))
  ax=timed(5:end-6)+0.005;
  bsel=find(timed<=0);
AA= bsxfun(@minus, allSPmA, nanmean(nanmean(allSPmA(bsel,:,:),1),3));
AC= bsxfun(@plus, allSPmA, nanmean(nanmean(allSPmA(bsel,:,:),1),3));
%AA=AA./AC;
BB= bsxfun(@minus, allSPmB, nanmean(nanmean(allSPmB(bsel,:,:),1),3));
BC= bsxfun(@plus, allSPmB, nanmean(nanmean(allSPmB(bsel,:,:),1),3));
%BB=BB./BC;
 
 NT=nansum(nansum(nanmean(abs(allSPmA(:,:,:)),3),1)>0);
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
 fill_error_area2(ax,nanmean(AA(:,:,1),2) ,nanstd(AA(:,:,1),[],2)./sqrt(NT),[0.6 0.6 0.6])  
 fill_error_area2(ax,nanmean(AA(:,:,2),2) ,nanstd(AA(:,:,2),[],2)./sqrt(NT),[0.6 0.6 0.6])  
  fill_error_area2(ax,nanmean(AA(:,:,3),2) ,nanstd(AA(:,:,3),[],2)./sqrt(NT),[0.6 0.6 0.6])  
   plot(ax,squeeze(nanmean(AA(:,:,1),2)),'COlor',[ 0.3 0 0],'Linewidth',2);hold on,
  plot(ax,squeeze(nanmean(AA(:,:,2),2)),'COlor',[ 0.6 0 0],'Linewidth',2);hold on,
   plot(ax,squeeze(nanmean(AA(:,:,3),2)),'COlor',[ 0.9 0 0],'Linewidth',2);hold on,
 axis tight;ylim([-1 8])
   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_TCS_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_CS_PLOT.fig'] )

  NT=nansum(nansum(nanmean(abs(allSPmB(:,:,:)),3),1)>0);
 figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters')
 fill_error_area2(ax,nanmean(BB(:,:,1),2) ,nanstd(BB(:,:,1),[],2)./sqrt(NT),[0.6 0.6 0.6])  
 fill_error_area2(ax,nanmean(BB(:,:,2),2) ,nanstd(BB(:,:,2),[],2)./sqrt(NT),[0.6 0.6 0.6])  
  fill_error_area2(ax,nanmean(BB(:,:,3),2) ,nanstd(BB(:,:,3),[],2)./sqrt(NT),[0.6 0.6 0.6]) ;hold on, 
   plot(ax,squeeze(nanmean(BB(:,:,1),2)),'COlor',[ 0. 0 0.3],'Linewidth',2)
  plot(ax,squeeze(nanmean(BB(:,:,2),2)),'COlor',[ 0. 0 0.6],'Linewidth',2);hold on,
   plot(ax,squeeze(nanmean(BB(:,:,3),2)),'COlor',[ 0. 0 0.9],'Linewidth',2);hold on,
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_SS_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_SS_PLOT.fig'] )
 axis tight;ylim([-1 8])
for h=1:2
    if h==1
tsel=find(ax>0 & ax<0.2);else
tsel=find(ax>0.2 & ax<1.2);end
 Am= squeeze(nanmean(nanmean(AA(tsel,:,:),1),2));
 As= squeeze(nanstd(nanmean(AA(tsel,:,:),1),[],2))./sqrt(size(AA,2));
 Bm= squeeze(nanmean(nanmean(BB(tsel,:,:),1),2));
 Bs= squeeze(nanstd(nanmean(BB(tsel,:,:),1),[],2))./sqrt(size(BB,2));
 Aval=0.4
  figure('COlor','w','Position',[300 300 150 pheight],'Renderer', 'painters')
  b1=bar(1,Am(1),'Facecolor',[ 0.3 0 0]); hold on,set(b1,'FaceAlpha',Aval)
  b1=bar(2,Am(2),'Facecolor',[ 0.6 0 0]); set(b1,'FaceAlpha',Aval)
  b1=bar(3,Am(3),'Facecolor',[ 0.9 0 0]); set(b1,'FaceAlpha',Aval)
  errorbar([1 2 3], Am, As,'.k')
    b1=bar(5,Bm(1),'Facecolor',[ 0 0 0.3]); hold on,set(b1,'FaceAlpha',Aval)
  b1=bar(6,Bm(2),'Facecolor',[ 0 0 0.6]); set(b1,'FaceAlpha',Aval)
  b1=bar(7,Bm(3),'Facecolor',[ 0 0 0.9]); set(b1,'FaceAlpha',Aval)
  errorbar([5:7], Bm, Bs,'.k')
  if h==1
     print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_TRANS_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_TRANS_PLOT.fig'] )
  else
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_SUST_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_SUST_PLOT.fig']);end
end
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   ax=timed(5:end-6)+0.005;
  bsel=find(timed<=0);
AA= allSPmA;%bsxfun(@minus, allSPmA, nanmean(nanmean(allSPmA(bsel,:,:),1),3));
%AA=AA./AC;
BB=allSPmB;% bsxfun(@minus, allSPmB, nanmean(nanmean(allSPmB(bsel,:,:),1),3));
%B 
for h=1:2
    if h==1
tsel=find(ax>0 & ax<0.2);else
tsel=find(ax>0.2 & ax<1.2);end
tselB=find(ax>-0.7 & ax<0);
 Am= squeeze(nanmean(nanmean((AA(tsel,:,:)),1),3));
 AmB= squeeze(nanmean(nanmean((AA(tselB,:,:)),1),3));
 Bm= squeeze(nanmean(nanmean((BB(tsel,:,:)),1),3));
BmB= squeeze(nanmean(nanmean((BB(tselB,:,:)),1),3));
 Aval=0.4
  figure('COlor','w','Position',[300 300 200 pheight],'Renderer', 'painters')
violinplot2((Am-AmB)',[1.25 ],'ViolinColor', [ 0.9 0.3 0.1; 0 0 0.9])
violinplot2((Bm-BmB)',[1.85 ],'ViolinColor', [ 0.2 0.3 0.9; 0 0 0.9])
line([ 0.8 2.3], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5);axis tight
set(gca,'Xtick',[],'Xticklabel',[])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_GAM__violinCS.pdf'])

  if h==1
     print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_TRANS_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_TRANS_PLOT.fig'] )
  else
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'DCopto_FIRING_SUST_PLOT.pdf'])
savefig(gcf, [ savepath 'DCopto_FIRING_SUST_PLOT.fig']);end
end


%% STAT for FIGURE 3 e and f  (firing rate)
  tsel=find(ax>0 & ax<0.2);
tsel=find(ax>0.2 & ax<1.2);
tselB=find(ax>-0.5 & ax<0);
V1=squeeze(nanmean(nanmean(nanmean(AA(tsel,:,:),3),1),1));
%V2=squeeze(nanmean(nanmean(nanmean(AA(tselB,:,:),3),1),1));
V2=squeeze(nanmean(nanmean(nanmean(BB(tsel,:,:),3),1),1));
[h,p,ci,stats] = ttest(V1,V2)

% df=26,  1.5966e-05


%%%%%TFR
clear allM allM2
for id=1:length(allTFR)
   
       A=  allTFR{id}(:,:);
     B=nanmean(A(:,100:800),2);
 allM(:,1:size(A,2),id)=    smooth2a(bsxfun(@minus,A,B)./bsxfun(@plus,A,B),1,35) ;
 allM2(:,1:size(A,2),id)=    smooth2a(A,1,35) ;
end
allM(allM==0)=NaN;allM2(allM2==0)=NaN;
%        figure('COlor','w','Position', [ 300 400 250 120],'Renderer', 'painters')
% imagesc((1:size(allM,2))./828-1.2,freq2.freq,nanmean(allM,3) );%ylim([2 20])
% axis xy;colormap(jet);set(gca,'Clim',[-0.1 0.1])
% xlim([ -0.7 1.7])

fsel=freq2.freq>30 & freq2.freq<90;
%fsel=freq2.freq>5 & freq2.freq<14;
V1=squeeze(nanmean(nanmean(allM(fsel,900:1600,:),2),1));
tsel=find(ax>0.1 & ax<1.5);tselB=find(ax>-1 & ax<-0.1);
S1=squeeze(nanmean(nanmean(nanmean(AA(tsel,:,:),3),1),1));S2=squeeze(nanmean(nanmean(nanmean(BB(tsel,:,:),3),1),1));
S1b=squeeze(nanmean(nanmean(nanmean(AA(tselB,:,:),3),1),1));S2b=squeeze(nanmean(nanmean(nanmean(BB(tselB,:,:),3),1),1));
h1=(S2-S2b)./(S2+S2b);h2=(S1-S1b)./(S1+S1b);
h1(isnan(h2))=NaN;h2(isnan(h1))=NaN;
xdata1=h1;ydata1=V1;ydata1(isnan(xdata1))=[];xdata1(isnan(xdata1))=[];
    figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(xdata1,ydata1,'.k','Markersize',15);fitResults1 = polyfit(xdata1,ydata1',1);
yplot1 = polyval(fitResults1,xdata1);hold on,plot(xdata1,yplot1,'COlor', [ 0.4 0.4 0.8],'Linewidth',1);
xlim([-1 1])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Scatter_GAM_SS.pdf'])
savefig(gcf, [ savepath 'Scatter_GAMSS.fig'])
xdata1=h2;
ydata1=(V1);ydata1(isnan(xdata1))=[];xdata1(isnan(xdata1))=[];
    figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(xdata1,ydata1,'.k','Markersize',15);fitResults1 = polyfit(xdata1,ydata1',1);
yplot1 = polyval(fitResults1,xdata1);hold on,plot(xdata1,yplot1,'COlor', [ 0.8 0.4 0.4],'Linewidth',1);
xlim([-1 1])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Scatter_GAM_CS.pdf'])
savefig(gcf, [ savepath 'Scatter_GAMCS.fig'])

[B,BINT,R,RINT,stats] = regress((V1),[ones(length(S2-S2b),1),(((S2-S2b)./(S2+S2b))')])
[B,BINT,R,RINT,stats] = regress((V1),[ones(length(S1-S1b),1),(((S1-S1b)./(S1+S1b))')])



%fsel=freq2.freq>30 & freq2.freq<90;
fsel=freq2.freq>7 & freq2.freq<15;
V1=squeeze(nanmean(nanmean(allM(fsel,900:1600,:),2),1));
tsel=find(ax>0.1 & ax<1.5);tselB=find(ax>-1 & ax<-0.1);
S1=squeeze(nanmean(nanmean(nanmean(AA(tsel,:,:),3),1),1));
S2=squeeze(nanmean(nanmean(nanmean(BB(tsel,:,:),3),1),1));
S1b=squeeze(nanmean(nanmean(nanmean(AA(tselB,:,:),3),1),1));
S2b=squeeze(nanmean(nanmean(nanmean(BB(tselB,:,:),3),1),1));
h1=(S2-S2b)./(S2+S2b); %SS
h2=(S1-S1b)./(S1+S1b); %CS
h1(isnan(h2))=NaN;h2(isnan(h1))=NaN;
xdata1=h1; % SS
ydata1=V1;ydata1(isnan(xdata1))=[];xdata1(isnan(xdata1))=[];
    figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(xdata1,ydata1,'.k','Markersize',15);fitResults1 = polyfit(xdata1,ydata1',1);
yplot1 = polyval(fitResults1,xdata1);hold on,plot(xdata1,yplot1,'COlor', [ 0.4 0.4 0.8],'Linewidth',1);
xlim([-1 1])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Scatter_Theta_SS.pdf'])
savefig(gcf, [ savepath 'Scatter_ThetaSS.fig'])
xdata1=h2;%CS
ydata1=(V1);ydata1(isnan(xdata1))=[];xdata1(isnan(xdata1))=[];
    figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(xdata1,ydata1,'.k','Markersize',15);fitResults1 = polyfit(xdata1,ydata1',1);
yplot1 = polyval(fitResults1,xdata1);hold on,plot(xdata1,yplot1,'COlor', [ 0.8 0.4 0.4],'Linewidth',1);
xlim([-1 1])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Scatter_Theta_CS.pdf'])
savefig(gcf, [ savepath 'Scatter_ThetaCS.fig']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,BINT,R,RINT,stats] = regress((V1),[ones(length(S2-S2b),1),(((S2-S2b)./(S2+S2b))')])
[B,BINT,R,RINT,stats] = regress((V1),[ones(length(S1-S1b),1),(((S1-S1b)./(S1+S1b))')])

% 
tsel=find(ax>0.1 & ax<1.5);tselB=find(ax>-0.9 & ax<-0.1);
S1=squeeze(nanmean(nanmean(nanmean(AA(tsel,:,:),3),1),1));
S2=squeeze(nanmean(nanmean(nanmean(BB(tsel,:,:),3),1),1));
S1b=squeeze(nanmean(nanmean(nanmean(AA(tselB,:,:),3),1),1));
S2b=squeeze(nanmean(nanmean(nanmean(BB(tselB,:,:),3),1),1));
h1=(S2-S2b)./(S2+S2b);
h2=(S1-S1b)./(S1+S1b);
nx=0;wid=5;
clear allFIT
for f=1:size(allM,1)
    nx=nx+1;
%fsel=freq2.freq>f & freq2.freq<f+wid;

V1=squeeze(nanmean(nanmean(allM(f,900:1600,:),2),1));
[B,BINT,R,RINT,stats] = regress(zscore(V1),[ones(length(h1),1),((h1)')]);
allFIT(nx,1)=stats(1);
[B,BINT,R,RINT,stats] = regress(zscore(V1),[ones(length(h1),1),((h2)')]);
allFIT(nx,2)=stats(1);
end
    figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
plot(freq2.freq,allFIT(:,1),'Linewidth',1.5,'COlor',[ 0.1 0.1 0.8]); hold on,
plot(freq2.freq,allFIT(:,2),'Linewidth',1.5,'COlor',[ 0.8 0.1 0.1]);axis tight
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Rsquare_spec.pdf'])
savefig(gcf, [ savepath 'Rsquare_spec.fig']);

% M1=squeeze(nanmean(allSFC_LFP_S(:,:,:),2));
% M2=squeeze(nanmean(allSFCB_LFP_S(:,:,:),2));
% figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
% plot(freq2.freq,nanmean(M1(:,(S1./S2)< prctile((S1./S2),50)),2),'b','Linewidth',1);
% hold on,
% plot(freq2.freq,nanmean(M1(:,(S1./S2)> prctile((S1./S2),50)),2),'r','Linewidth',1);
% set(gca,'Xscale','log')

% %%%%%%%%%%%%%%%%%%%
% PLOT histo of CS and SS opto-induced distribution
%  figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
% b=violin([h2(~isnan(h2))'   h1(~isnan(h1))'    ],'Facecolor',[ 0.5 0.1 0])
% 
% set(b(1),'Facecolor',[0.6 0 0])
% set(b(2),'Facecolor',[0 0 0.6])
%%%%%%%%%

 figure('COlor','w','Position', [ 300 400 125 pheight],'Renderer', 'painters')
b=violinplot([h2(~isnan(h2))'  ],[ 1 ],'ViolinColor',[0.6 0 0]   )
ylim([ -1 1]); xlim([ 0.5 1.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Violingplot_CS.pdf'])


 figure('COlor','w','Position', [ 300 400 125 pheight],'Renderer', 'painters')
b=violinplot([h1(~isnan(h1))'  ],[ 1 ],'ViolinColor',[0. 0 0.6]   )
ylim([ -1 1]);xlim([ 0.5 1.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Violingplot_SS.pdf'])
% 


V1=nanmean(nanmean(allM(:,900:2000,:),2),3);
V1s=nanstd(nanmean(allM(:,900:2000,:),2),[],3)./sqrt(size(allM,3));
    figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
,plot(freq2.freq,V1,'k')
fill_error_area2(freq2.freq,V1,V1s,[ 0.5 0.5 0.5])
,plot(freq2.freq,V1,'k')
line([ freq2.freq(1) freq2.freq(end)],[ 0 0],'COlor',[ 0 0 0]) 
axis tight,set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Opow_spec.pdf'])
savefig(gcf, [ savepath 'Opow_spec.fig']);



%     figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
% 
% for ses=1:36;
%    subplot(6,6,ses)
% V2=(squeeze(nanmean(nanmean((allM2(:,100:800,ses)),2),3)))';
% V1=(squeeze(nanmean(nanmean((allM2(:,900:1400,ses)),2),3)))';
% plot(freq2.freq,V1.*freq2.freq.^0.5,'Linewidth',1.5,'COlor',[ 0.1 0.7 0.1]); hold on,
% plot(freq2.freq,V2.*freq2.freq.^0.5,'Linewidth',1.5,'COlor',[ 0.1 0.1 0.1]); hold on,
% axis tight
% set(gca,'Xscale','log')
% end
% 
% 
% 
