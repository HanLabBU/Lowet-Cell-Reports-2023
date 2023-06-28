

addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))
addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\hip_scripts\'))
%
clear all
mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
cd(mainpath);

foldnames = {'DATA_FIG4'} ;
clear allISI allPHASE allCOH   allBURST  lfp_q allname  allangs
allSNR=[];  allRATE=[]; LFPpre=[]; allB=[];allSPEC=[];
loc_D =[];pN=[];pB=[];allIF=[];
allS_N=[] ;allB_N=[];POWsel=[];POWunsel=[];  allS_BA=[];

mr=0; mmm=0;   %
foldn=1;
cd(mainpath)
cd(foldnames{foldn})
ses=dir('*.mat');
for ih=[1:length(ses)]  %%
    Cpath= ses(ih).name
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
    
    
    
    ADP_s=[];  ADP_b=[];
    
    
    %ompute overall SNR and RATE
    for neuron=1:size(result.traces,2)   %neuron within a session
        
        
        SNR=[];RATE=[];  %% SNR and RATE aggregated over trials
        for id=1:length(result.resultS)
            if ~isempty(result.resultS{id})
                clear A B
                % for neuron=1:length(result.resultS{id}.spike_snr)
                A=   nanmean(result.resultS{id}.spike_snr{neuron}) ;
                B=sum(result.resultS{id}.roaster(neuron,10:FS-10))./(size(result.resultS{id}.roaster(10:FS-10),2)./FS);
                %  end;
                SNR= [ SNR,  A'   ];
                RATE= [ RATE, B'   ];
            end
        end
        
        
        mmm=mmm+1;  %% overall neuron indicator !!
        
        
        if isfield(result,'mis')  % Deal with a missing trial
            if ~isempty(result.mis)
                if length(result.mis)==1
                    result.lfp.trial_order=[ 1: result.mis-1 result.mis+1: length(result.lfp.trial)];
                else
                end
            end
        end
        
        mm=0; clear nall allF allopto
        spang=[];dataM=[];
        gP=[]; tP=[];
        for tt=unique(result.trial_vec)   %% trial iteration!!!!
            
            if ~isempty(result.resultS{tt}) & tt <= length(result.lfp.trial)
                if nanmean(result.resultS{tt}.spike_snr{neuron}) >4
                    
                    
                    SpikeTimes=result.resultS{tt}.spike_idx{neuron};%./FS;
                    Vm_trace=result.resultS{tt}.trace_ws(neuron,:);
                    Vm_trace(1:10)= mean(Vm_trace(11:15));
                    if length(Vm_trace)<3500
                        try
                            [ fitbaseline, coeff]=exp_fit_Fx(Vm_trace',round(FS)); %remove photobleaching
                            Vm_trace=(Vm_trace-fitbaseline);end
                        
                    else
                        try
                            [ fitbaseline, coeff]=exp_fit_Fx(Vm_trace',round(FS)); %remove photobleaching
                            Vm_trace=(Vm_trace-fitbaseline); end
                    end
                    
                    Fn = FS/2;FB=[ 86 89]; % Remove camera noise
                    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                    LFPg= ((filtfilt(B,A,    Vm_trace)));
                    Vm_trace=   Vm_trace-LFPg;
                    FB=[ 3 12]; % theta frequency range
                    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                    LFPg= ((filtfilt(B,A,    Vm_trace)));
                    Vm_trace3=   LFPg;
                    %
                    
                    Vm_trace2= result.resultS{tt}.orig_trace(neuron,:);
                    
                    if length(Vm_trace2)<3500
                        try
                            [ fitbaseline, coeff]=exp_fit_Fx(Vm_trace2',round(FS)); %remove photobleaching
                            Vm_trace2=(Vm_trace2-fitbaseline);end
                    end
                    
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
                    
                    %%%%%%%%%%%%%
                    %% Line noise denoise
                    Fn = FS/2;FB=[ 58 62]; % 
                    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                    LFPg= ((filtfilt(B,A,  LFPtr)));
                    LFPtr= LFPtr-LFPg;
                    
                    %% USED as quality index
                    Fn = FS/2;FB=[ 1 3]; % low delta
                    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                    LFPg1= ((filtfilt(B,A,  LFPtr)));
                    Fn = FS/2;FB=[ 4 10]; % theta
                    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                    LFPTP= ((filtfilt(B,A,  LFPtr)));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    tP=[tP; mean(abs(hilbert(LFPTP))) ];
                    gP=[gP; mean(abs(hilbert(LFPg1))) ];
                  
                    
                 
                    opto= opto-mean(opto(1:400));opto2=opto;
                    opto(1:400)=0; opto(2700:end)=0;
                    optostart= find(fastsmooth(opto,1,1,1)>mean(opto));
                    if ~isempty(optostart)
                        optostart=optostart(1); 
                    else
                        if FS <600   % solution for some missing opto trials
                            optostart=500;
                        else
                            optostart=950;
                        end
                    end
                    optostart2=   optostart; optostart=0;
                    %%%%%%%%%%%%%%%%
                    %%%%%%%%%%
                    
                    if length(Vm_trace) <= length(LFPtr)
                        dataM.trial{mm}(1,:) = (double(Vm_trace));
                        dataM.trial{mm}(2,:) =  zscore(double(LFPtr(1:length(Vm_trace))));
                        
                        dataM.trial2{mm}(1,:) = (double(Vm_trace2)); %orig_trace
                        dataM.trial3{mm}(1,:) = zscore(double(Vm_trace3)); % Vm theta phase
                    else
                        dataM.trial{mm}(1,:) =(double(Vm_trace(1:length(LFPtr))));
                        dataM.trial{mm}(2,:) =  zscore(double(LFPtr));
                        dataM.trial2{mm}(1,:) =(double(Vm_trace2(1:length(LFPtr))));
                        dataM.trial3{mm}(1,:) = zscore(double(Vm_trace3(1:length(LFPtr))));
                    end
                    dataM.opto{mm}= opto2;
                    dataM.optostart{mm}=optostart;
                    dataM.spikes{mm} = (SpikeTimes-optostart);
                    dataM.spikeSNR{mm} =   SpikeSNR;
                    dataM.spikeamplitude{mm}=result.resultS{tt}.spike_amplitude{neuron};
                    dataM.spikes2{mm} = SpikeTimes;
                    dataM.time{mm}= ((1:length(Vm_trace))-optostart)./FS;
                    dataM.roaster{mm}= result.resultS{tt}.roaster(neuron,:);
                    allopto(1:length(opto),mm)= opto2;
                    dataM.optostartM=   optostart2; 
                    
                end  %% trial
            end
        end
    
        dataM.label= {'A'; 'B'};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear spikTypTET
        wind=9; % 828Hz, approx 1.2ms time stamps
        for tria = 1:length(dataM.spikes)
            s= dataM.spikes{tria}  ;
            for in = 1:length(s)
                if length(s) >2
                    if in==1
                        d=(s(in+1)-  s(in));
                        if d <wind;   spikTypTET{tria}(in)= 1 ;
                        else; spikTypTET{tria}(in)= 0 ;end
                        
                    elseif in == length(s)
                        d=( s(in)-s(in-1) );
                        if d <wind;   spikTyp{tria}(in)= 1 ;
                        else; spikTypTET{tria}(in)= 0 ;end
                    else
                        d1=(s(in+1)-  s(in));         d=( s(in)-s(in-1) );
                        if d <wind |d1 <wind ;   spikTypTET{tria}(in)= 1 ;
                        else; spikTypTET{tria}(in)= 0 ;end
                        
                    end
                else
                    spikTypTET{tria}=[];
                end
            end
        end
        
        clear spikTypS
        wind=12;
        for tria = 1:length(dataM.spikes)
            s= dataM.spikes{tria}  ;
            for in = 1:length(s)
                if length(s) >2
                    if in==1
                        d=(s(in+1)-  s(in));
                        if d <wind;   spikTypS{tria}(in)= 1 ;
                        else; spikTypS{tria}(in)= 0 ;end
                        
                    elseif in == length(s)
                        d=( s(in)-s(in-1) );
                        if d <wind;   spikTypS{tria}(in)= 1 ;
                        else; spikTypS{tria}(in)= 0 ;end
                    else
                        d1=(s(in+1)-  s(in));         d=( s(in)-s(in-1) );
                        if d <wind |d1 <wind ;   spikTypS{tria}(in)= 1 ;
                        else; spikTypS{tria}(in)= 0 ;end
                        
                    end
                else
                    spikTypS{tria}=[];
                end
            end
        end
        
        clear spikTyp
        wind=12; % 828Hz, approx 1.2ms time stamps
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear spikTyp7
        spikTyp7=spikTyp;
        for tr=1:length(dataM.spikes2)
            clear allW allW2; mm=0;
            sptim=dataM.spikes2{tr}(spikTyp5{tr}==0);z=find(spikTyp5{tr}==0) ;
            if ~isempty(z)
                volt=zscore(dataM.trial2{tr}(1,:)-fastsmooth(dataM.trial2{tr}(1,:),1400,1,1));
                volt2=zscore(dataM.trial{tr}(1,:)-fastsmooth(dataM.trial{tr}(1,:),1400,1,1));
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
                ADP=A./sp_height;avADPSS=nanmean(ADP);ADP_s=[ADP_s,ADP];
            end
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
                A= nanmean(allW2(wind2+windS2A:wind2+windS2B,:),1)-nanmean(allW2(wind2-windSA:wind2-windSB,:),1);
                sp_height= nanmedian(nanmedian(allW(wind2+1,:),1));
                ADP=A./sp_height;avADP=nanmean(ADP);ADP_b=[ADP_b,ADP];
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
        
        [spikTypB2]=spike_burst_n(dataM,spikTyp7,FS,1);
        [spikTypB3]=spike_burst_n(dataM,spikTyp7,FS,2);
        [spikTypB4]=spike_burst_n(dataM,spikTyp7,FS,3);
        [spikTypB5]=spike_burst_n(dataM,spikTyp7,FS,4);
        
        
        dataM.roasterB=burst_spike_select(dataM,spikTyp5,spikTyp7,-1,1,2) ;
        dataM.roasterBTET=burst_spike_select(dataM,spikTypTET,spikTyp7,1,-1,2) ;
        dataM.roasterBON=burst_spike_select(dataM,spikTyp5,spikTyp7,-1,1,2) ;
        dataM.roasterBONS=burst_spike_select(dataM,spikTyp5,spikTyp7,1,1,1) ;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SHUFFLING trial relation %%%
        nj=0;
        for ik=[2:size(dataM.roaster,2) 1] %%
            ROAST1=zeros(1,length(dataM.roaster{ik}));try
               timsel=((dataM.spikes{ik}>10))';
                ROAST1(dataM.spikes{ik}(spikTyp7{ik}==1  & timsel))=1;
                nj=nj+1;
                dataM.roasterBON1{nj}=ROAST1(end:-1:1);end
        end
        nj=0;
        for ik=[2:size(dataM.roaster,2) 1]
            ROAST1=zeros(1,length(dataM.roaster{ik}));try
                     timsel=((dataM.spikes{ik}>10))';
                ROAST1(dataM.spikes{ik}(spikTyp{ik}==0    & timsel))=1;
                nj=nj+1;
                dataM.roasterBON2{nj}=ROAST1(end:-1:1);end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SELECTION OF SB AS A FUCNTION OF SB NUMBER
        dataM.roasterB2=burst_spike_select(dataM,spikTypB2,spikTyp7,1,1,2) ;
        dataM.roasterB3=burst_spike_select(dataM,spikTypB3,spikTyp7,1,1,2) ;
        dataM.roasterB4=burst_spike_select(dataM,spikTypB4,spikTyp7,1,1,2) ;
        dataM.roasterB5=burst_spike_select(dataM,spikTypB5,spikTyp7,1,1,2) ;
        
        
        %% SS
        dataM.roasterS=burst_spike_select(dataM,spikTypS,spikTyp7,0,0,1) ;
        %%%%
        %burst 1st spike
         dataM.roasterB2f=burst_spike_select(dataM,spikTypB2,spikTyp7,2,1,1) ;
         dataM.roasterB3f=burst_spike_select(dataM,spikTypB3,spikTyp7,2,1,1) ;
         dataM.roasterB4f=burst_spike_select(dataM,spikTypB4,spikTyp7,2,1,1) ;
        for ik=1:size(dataM.roaster,2)
            ROAST1=zeros(1,length(dataM.roaster{ik}));try
                ROAST1(dataM.spikes{ik}(spikTypB5{ik}==2 | spikTypB4{ik}==2 ))=1;
                dataM.roasterB45f{ik}=ROAST1;end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A=[]; B=[];C=[]; %% NUMBER OF SB AND SS for the given neuron
        for id2=1
            nn=0;clear allF;ttim=0;
            for ind=1:length(dataM.spikes)
                nn=nn+1;try
                    t= ((dataM.spikes{ind}>10) );
                    A=[A,dataM.spikes{ind}(spikTyp7{ind}==1 & spikTyp5{ind}==1  &  t')']
                    B=[B,dataM.spikes{ind}(spikTyp{ind}==0 & t')'];
                    C=[C,dataM.spikes{ind}(spikTyp7{ind}==1   &  t')']
                    
                    ttim=ttim+length(t);
                end
            end
        end
        %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ISI=nanmean(nall,2);,
        FI=(nax> 0.005 & nax<0.015);FL=(nax> 0.035 & nax<0.055);
        allBURSTF= sum(ISI(FI))./sum(ISI(FL));
        
        
        if    length(A)>10 &  length(B)>10  & (nanmean(tP)./nanmean(gP))>1.1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [cwt_out,frs]=runcwt(dataM, [3 180],FS);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freq2.freq=frs;
            wavA = abs(squeeze(cwt_out));
            wavD = angle(squeeze(cwt_out));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            mr=mr+1
            detG(mr)=1;
            %%
            %% putting everything in final structures %%%%%%%%
            if isfield(result.lfp,'preproc_v')
                lfp_q(mr)= 1;
            else
                lfp_q(mr)= 0;
            end
            
            ISI=nanmean(nall,2);
            FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
            allISI(:,mr)= ISI(2:end-1);
            allBURST(:,mr)= sum(ISI(FI))./sum(ISI(FL));
            
            allFS(mr)=FS;
            allIND{mr}=ih;
            %     allspADP(:,mr)= spADP;
            allavADP(mr)= nanmean(ADP_b)  ;
            allavADPSS(mr)=  nanmean(ADP_s);
            allname{mr}=  Cpath;
            alloptoM{mr}= allopto;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%
            %%%%%%%%%%%%%%
            CH=1;sp_shift=7;
            %%%%%%%%%%%%
            [PLV_all]= spike_field_ppc_adj(wavD,dataM.roaster,CH,sp_shift) ;
            [PLV_b2]= spike_field_ppc_adj(wavD,dataM.roasterB2,CH,sp_shift) ;
            [PLV_b3]= spike_field_ppc_adj(wavD,dataM.roasterB3,CH,sp_shift) ;
            [PLV_b4]= spike_field_ppc_adj(wavD,dataM.roasterB4,CH,sp_shift) ;
            [PLV_b5]= spike_field_ppc_adj(wavD,dataM.roasterB5,CH,sp_shift) ;
            [PLV_s]= spike_field_ppc_adj(wavD,dataM.roasterS,CH,sp_shift) ;
            [PLV_b]= spike_field_ppc_adj(wavD,dataM.roasterBON,CH,sp_shift) ;
            [PLV_b_1]= spike_field_ppc_adj(wavD,dataM.roasterBON1,CH,sp_shift) ;
            [PLV_b_2]= spike_field_ppc_adj(wavD,dataM.roasterBON2,CH,sp_shift) ;;
            [PLV_bs]= spike_field_ppc_adj(wavD,dataM.roasterBONS,CH,sp_shift) ;
        
         
            [ PLV_sC]= spike_field_ppc_adjCOH(wavD,dataM.roasterS,CH,sp_shift,wavA);
            [ PLV_bC]= spike_field_ppc_adjCOH(wavD,dataM.roasterB,CH,sp_shift,wavA);
            [ PLV_sCSH]= spike_field_ppc_adjCOH(wavD,dataM.roasterBON2,CH,sp_shift,wavA);
            [ PLV_bCSH]= spike_field_ppc_adjCOH(wavD,dataM.roasterBON1,CH,sp_shift,wavA);
            
            allPLVs(:,1,mr)= PLV_all;
            allPLVs(:,2,mr)= PLV_s;
            allPLVs(:,3,mr)= PLV_b;
            allPLVs(:,4,mr)= PLV_b2;
            allPLVs(:,5,mr)= PLV_b3;
            allPLVs(:,6,mr)= PLV_b4;
            allPLVs(:,7,mr)= PLV_b5;
            allPLVs(:,8,mr)= PLV_b_1;
            allPLVs(:,9,mr)= PLV_b_2;
            allPLVs(:,12,mr)= (PLV_bs);%+PLV_bs2)./2;;
            try
                allPLVsLFP(:,13,mr)= PLV_bCSH;
                allPLVsLFP(:,14,mr)= PLV_sCSH;
                allPLVsLFP(:,10,mr)= PLV_bC;%PLV_Match(:,1);
                allPLVsLFP(:,11,mr)= PLV_sC;%PLV_Match(:,2);
                
            end
            CH=2;sp_shift=0;
            %%%%%%%%%%%%,sp_shift
            [PLV_all]= spike_field_ppc_adj(wavD,dataM.roaster,CH,sp_shift) ;
            [PLV_b2]= spike_field_ppc_adj(wavD,dataM.roasterB2,CH,sp_shift) ;
            [PLV_b3]= spike_field_ppc_adj(wavD,dataM.roasterB3,CH,sp_shift) ;
            [PLV_b4]= spike_field_ppc_adj(wavD,dataM.roasterB4,CH,sp_shift) ;
            [PLV_b5]= spike_field_ppc_adj(wavD,dataM.roasterB5,CH,sp_shift) ;
            [PLV_s]= spike_field_ppc_adj(wavD,dataM.roasterS,CH,sp_shift) ;
            [ PLV_sC]= spike_field_ppc_adjCOH(wavD,dataM.roasterS,CH,sp_shift,wavA);
            [PLV_b]= spike_field_ppc_adj(wavD,dataM.roasterBON,CH,sp_shift) ;
            [ PLV_bC]= spike_field_ppc_adjCOH(wavD,dataM.roasterB,CH,sp_shift,wavA);
            [PLV_b_1]= spike_field_ppc_adj(wavD,dataM.roasterBON1,CH,sp_shift) ;
            [PLV_b_2]= spike_field_ppc_adj(wavD,dataM.roasterBON2,CH,sp_shift) ;
            [ PLV_sCSH]= spike_field_ppc_adjCOH(wavD,dataM.roasterBON2,CH,sp_shift,wavA);
            [ PLV_bCSH]= spike_field_ppc_adjCOH(wavD,dataM.roasterBON1,CH,sp_shift,wavA);
            [PLV_Match]= spike_field_ppc_adjMatch(wavD,dataM.roasterBONS,dataM.roasterS,CH,sp_shift) ;
            [PLV_bs]= spike_field_ppc_adj(wavD,dataM.roasterBONS,CH,sp_shift) ;
       
            
            [PLV_bTET]= spike_field_ppc_adj(wavD,dataM.roasterBTET,CH,sp_shift);
            [PLV_bTETSH]= spike_field_ppc_adj(wavD(end:-1:1,:,:,end:-1:1),dataM.roasterBTET,CH,sp_shift)   ;
            allPLVsLFP(:,1,mr)= PLV_all;
            allPLVsLFP(:,2,mr)= PLV_s;
            allPLVsLFP(:,3,mr)= PLV_b;
            allPLVsLFP(:,4,mr)= PLV_b2;
            allPLVsLFP(:,5,mr)= PLV_b3;
            allPLVsLFP(:,6,mr)= PLV_b4;
            allPLVsLFP(:,7,mr)= PLV_b5;
            allPLVsLFP(:,8,mr)= PLV_b_1;
            allPLVsLFP(:,9,mr)= PLV_b_2;
            allPLVsLFP(:,12,mr)= (PLV_bs);
            try
                allPLVsLFP(:,13,mr)= PLV_bCSH;
                allPLVsLFP(:,14,mr)= PLV_sCSH;
                allPLVsLFP(:,10,mr)= PLV_bC;%PLV_Match(:,1);
                allPLVsLFP(:,11,mr)= PLV_sC;%PLV_Match(:,2);
                
                
            end
            sp_shift=0;
            [PLV_LFPVM_B]= spike_field_ppc_adjVM(wavD,dataM.roasterBONS,1,sp_shift) ;
            [PLV_LFPVM_S]= spike_field_ppc_adjVM(wavD,dataM.roasterS,1,sp_shift) ;
            [PLV_LFPVM_B1]= spike_field_ppc_adjVMSH(wavD,dataM.roasterBONS,1,sp_shift) ;
            [PLV_LFPVM_S1]= spike_field_ppc_adjVMSH(wavD,dataM.roasterS,1,sp_shift) ;
            [PLV_LFPVM_BB]= spike_field_ppc_adjVM(wavD,dataM.roasterB45f,1,sp_shift) ;
            
            M= abs(squeeze(nanmean(nanmean(exp(1i.*circ_dist(wavD(:,1,:,:),wavD(:,2,:,:))),1),4)));
            
            allVmLFP(:,mr)= M;
            allPLV_VMLFP(:,1,mr)= PLV_LFPVM_B;
            allPLV_VMLFP(:,2,mr)= PLV_LFPVM_S;
            allPLV_VMLFP(:,3,mr)= PLV_LFPVM_B1;
            allPLV_VMLFP(:,4,mr)= PLV_LFPVM_S1;
            allPLV_VMLFP(:,5,mr)= PLV_LFPVM_BB;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CH=1;sp_shift=7;
            %%%%%%%%%%%%
            [PLV_b2]= spike_field_ppc_adj(wavD,dataM.roasterB2f,CH,sp_shift) ;
            [PLV_b3]= spike_field_ppc_adj(wavD,dataM.roasterB3f,CH,sp_shift) ;
            [PLV_b4]= spike_field_ppc_adj(wavD,dataM.roasterB4f,CH,sp_shift) ; 
            [PLV_s]= spike_field_ppc_adj(wavD,dataM.roasterS,CH,sp_shift) ;
            allPLVsf(:,1,mr)= PLV_s;
            allPLVsf(:,2,mr)= PLV_b2;
            allPLVsf(:,3,mr)= PLV_b3;
            allPLVsf(:,4,mr)= PLV_b4;
            % allPLVsf(:,5,mr)= PLV_b5;
            
            CH=2;sp_shift=0;
            %%%%%%%%%%%%
            [PLV_b2]= spike_field_ppc_adj1stF(wavD,dataM.roasterB2f,CH,sp_shift) ;
            [PLV_b3]= spike_field_ppc_adj1stF(wavD,dataM.roasterB3f,CH,sp_shift) ;
            [PLV_b4]= spike_field_ppc_adj1stF(wavD,dataM.roasterB4f,CH,sp_shift) ;
            % [PLV_b5]= spike_field_ppc_adj1stF(wavD,dataM.roasterB5f,CH,sp_shift) ;
            [PLV_s]= spike_field_ppc_adj1stF(wavD,dataM.roasterS,CH,sp_shift) ;
            allPLVsflfp(:,1,mr)= PLV_s;
            allPLVsflfp(:,2,mr)= PLV_b2;
            allPLVsflfp(:,3,mr)= PLV_b3;
            allPLVsflfp(:,4,mr)= PLV_b4;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            allB_N=[allB_N, length(A)];
            allS_N=[allS_N, length(B)];
            allS_BA=[allS_BA, length(C)];
         
            mr
            
        end
        allB= [allB ;sum(ISI(FI))./sum(ISI(FL))];
        allSNR= [ allSNR; nanmean(SNR,2)];
        allRATE= [ allRATE; nanmean(RATE,2)];
    end
end



%% RESULT ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%
allPLVs(allPLVs==0)=NaN;
allPLVsLFP(allPLVsLFP==0)=NaN;


%%%%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\fig4_volt\'
pheight=150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BB=allavADP;%
sm=1;frs=freq2.freq;
selN=    allB_N>20;% 
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
fill_error_area2(frs, fastsmooth(nanmean(allPLVsLFP(:,11,selN)-allPLVsLFP(:,[ 14 ],selN),3),sm,1,1),nanstd(allPLVsLFP(:,11,:),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs, fastsmooth(nanmean(nanmean(allPLVsLFP(:,[10   ],selN)-allPLVsLFP(:,[ 13 ],selN),2),3),sm,1,1),nanstd(nanmean(allPLVsLFP(:,[10   ],selN),2),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
plot(frs,  fastsmooth((nanmean(allPLVsLFP(:,11,selN)-allPLVsLFP(:,[ 14 ],selN),3)),sm,1,1),'b','Linewidth',1)   ; hold on,
plot(frs,  fastsmooth(nanmean(allPLVsLFP(:,[ 10 ],selN),3)-nanmean(allPLVsLFP(:,[ 13 ],selN),3),sm,1,1),'r','Linewidth',1)
axis tight;set(gca,'Xscale','log')

figure('COlor','w','Position', [ 300 400 90 90],'Renderer', 'painters')
fill_error_area2(frs, fastsmooth(nanmean(allPLVsLFP(:,11,selN)-allPLVsLFP(:,[ 14 ],selN),3),sm,1,1),nanstd(allPLVsLFP(:,11,:),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs, fastsmooth(nanmean(nanmean(allPLVsLFP(:,[10   ],selN)-allPLVsLFP(:,[ 13 ],selN),2),3),sm,1,1),nanstd(nanmean(allPLVsLFP(:,[10   ],selN),2),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
plot(frs,  fastsmooth((nanmean(allPLVsLFP(:,11,selN)-allPLVsLFP(:,[ 14 ],selN),3)),sm,1,1),'b')   ; hold on,
plot(frs,  fastsmooth(nanmean(allPLVsLFP(:,[ 10 ],selN),3)-nanmean(allPLVsLFP(:,[ 13 ],selN),3),sm,1,1),'r')
axis tight;xlim([ 20 90])
set(gca,'Xscale','log')

%% SB 1st
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
fill_error_area2(frs, fastsmooth(nanmean(allPLVsLFP(:,2,selN),3),sm,1,1),nanstd(allPLVsLFP(:,2,:),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs, fastsmooth(nanmean(allPLVsLFP(:,[12   ],selN),3),sm,1,1),nanstd(nanmean(allPLVsLFP(:,[12   ],selN),2),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
plot(frs,  fastsmooth(nanmean(allPLVsLFP(:,2,selN),3),sm,1,1),'b')   ; hold on,
plot(frs,  fastsmooth(nanmean(allPLVsLFP(:,[ 12 ],selN),3),sm,1,1),'r')
axis tight;set(gca,'Xscale','log')


fr=freq2.freq>3 & freq2.freq<12
V1=   squeeze(nanmean(allPLVsLFP(fr,12,:),1));
V2=   squeeze(nanmean(allPLVsLFP(fr,2,:),1));
fr=freq2.freq>30 & freq2.freq<=90
V3=   squeeze(nanmean(allPLVsLFP(fr,12,:),1));
V4=   squeeze(nanmean(allPLVsLFP(fr,2,:),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.1 0.1])

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.04 0.04])



%% SB 1st
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
fill_error_area2(frs, fastsmooth(nanmean(allPLVs(:,2,selN),3),sm,1,1),nanstd(allPLVs(:,2,:),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs, fastsmooth(nanmean(allPLVs(:,[12   ],selN),3),sm,1,1),nanstd(nanmean(allPLVs(:,[12   ],selN),2),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
plot(frs,  fastsmooth(nanmean(allPLVs(:,2,selN),3),sm,1,1),'b')   ; hold on,
plot(frs,  fastsmooth(nanmean(allPLVs(:,[ 12 ],selN),3),sm,1,1),'r')
axis tight
set(gca,'Xscale','log')

fr=freq2.freq>3 & freq2.freq<12
V1=   squeeze(nanmean(allPLVs(fr,12,:),1));
V2=   squeeze(nanmean(allPLVs(fr,2,:),1));
fr=freq2.freq>30 & freq2.freq<=90
V3=   squeeze(nanmean(allPLVs(fr,12,:),1));
V4=   squeeze(nanmean(allPLVs(fr,2,:),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.45 0.45])

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.2 0.2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fr=freq2.freq>3 & freq2.freq<=9
V1=   squeeze(nanmean(allPLVsLFP(fr,10,:),1));
V2=   squeeze(nanmean(allPLVsLFP(fr,11,:),1));
fr=freq2.freq>30 & freq2.freq<=90
V3=   squeeze(nanmean(allPLVsLFP(fr,10,:),1));
V4=   squeeze(nanmean(allPLVsLFP(fr,11,:),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.4 0.4])

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.15 0.15])

%%%%%%%%%%%%%%%% VM  %%%%%%%%%%%%%%%%%%%

frs=freq2.freq;
%selN= BB>prctile(BB,0);% allavADP>=0.0;
frs=freq2.freq;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
fill_error_area2(frs,nanmean(allPLVs(:,2,selN),3),nanstd(allPLVs(:,2,:),[],3)./sqrt(size(allPLVs,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs,nanmean(allPLVs(:,[12 ],selN),3),nanstd(nanmean(allPLVs(:,[3 ],:),2),[],3)./sqrt(size(allPLVs,3)),[ 0.5 0.5 0.5])
plot(frs,  (nanmean(allPLVs(:,2,selN),3)),'b')   ; hold on,
plot(frs,  (nanmean(nanmean(allPLVs(:,[12  ],selN),2),3)),'r')
axis tight
set(gca,'Xscale','log')

allPLV_VMLFP(allPLV_VMLFP==0)=NaN;
BB=allavADP-allavADPSS;
gg=0;%selN=   BB>prctile(BB,0);%allBURST>5;
frs=freq2.freq;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
fill_error_area2(frs,(nanmean(allPLV_VMLFP(:,2,  selN)-allPLV_VMLFP(:,[4  ],  selN),3)),nanstd(allPLV_VMLFP(:,2, selN)-allPLV_VMLFP(:,[4  ], selN),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
fill_error_area2(frs,  (nanmean(allPLV_VMLFP(:,1,  selN)-allPLV_VMLFP(:,[3  ],  selN),3)),nanstd(allPLV_VMLFP(:,1, selN)-allPLV_VMLFP(:,[3  ], selN),[],3)./sqrt(size(allPLVsLFP,3)),[ 0.5 0.5 0.5])
plot(frs,  (nanmean(allPLV_VMLFP(:,2,  selN)-allPLV_VMLFP(:,4, selN),3)),'b','Linewidth',1)   ; hold on,
plot(frs,   (nanmean(allPLV_VMLFP(:,1,  selN)-allPLV_VMLFP(:,3,  selN),3)),'r','Linewidth',1)
axis tight;ylim([ -0.01 0.1])
set(gca,'Xscale','log')

fr=freq2.freq>3 & freq2.freq<12
V1=   squeeze(nanmean(allPLV_VMLFP(fr,1,selN)-allPLV_VMLFP(fr,3,selN),1));
V2=   squeeze(nanmean(allPLV_VMLFP(fr,2,selN)-allPLV_VMLFP(fr,4,selN),1));
fr=freq2.freq>30 & freq2.freq<=90
V3=   squeeze(nanmean(allPLV_VMLFP(fr,1,selN)-allPLV_VMLFP(fr,3,selN),1));
V4=   squeeze(nanmean(allPLV_VMLFP(fr,2,selN)-allPLV_VMLFP(fr,4,selN),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.15 0.15])
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.03 0.03])



fr=freq2.freq>3 & freq2.freq<12
V4=   squeeze(nanmean(allPLVsf(fr,4,:),1));
V3=   squeeze(nanmean(allPLVsf(fr,3,:),1));
V1=   squeeze(nanmean(allPLVsf(fr,2,:),1));
V2=   squeeze(nanmean(allPLVsf(fr,1,:),1));
figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.3 0 0.5; 0 0 0.9])
violinplot2(V3- V2,[1.8 ],'ViolinColor', [ 0.5 0 0.3; 0 0 0.9])
violinplot2(V4- V2,[2.3 ],'ViolinColor', [ 0.99 0 0.; 0 0 0.9])
%violinplot2(V5- V2,[2.8 ],'ViolinColor', [ 0.999 0 0; 0 0 0.9]);axis tight
hold on,hold on, line([ 0.7 2.8], [0 0],'Color', [ 0.2 0.2 0.2])
M=[V1-V2, V3-V2, V4-V2];M1=[ ones(length(V1),1)*0.8,ones(length(V1),1)*1.8,ones(length(V1),1)*2.7];
[coef, pval] =corr(M1(~isnan(M)),M(~isnan(M)))
Xvals=[ones(length(M(~isnan(M))),1),M1(~isnan(M))];
[B,BINT,R,RINT,stats] = regress((M(~isnan(M))),Xvals)

yplot1 = polyval(B(end:-1:1),Xvals(:,2));
hold on,plot(Xvals(:,2),yplot1,'k','Linewidth',1,'Color',[ 0.4 0.4 0.4])

fr=freq2.freq>3 & freq2.freq<12
nsel= squeeze(nanmean(allPLVsLFP(fr,1,:),1))>0.005;allB_N>10;%;
%V5=   squeeze(nanmean(allPLVsflfp(fr,5, nsel),1));
V4=   squeeze(nanmean(allPLVsflfp(fr,4, nsel),1));
V3=   squeeze(nanmean(allPLVsflfp(fr,3, nsel),1));
V1=   squeeze(nanmean(allPLVsflfp(fr,2, nsel),1));
V2=   squeeze(nanmean(allPLVsflfp(fr,1, nsel),1));
figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.3 0 0.5; 0 0 0.9])
violinplot2(V3- V2,[1.8 ],'ViolinColor', [ 0.5 0 0.3; 0 0 0.9])
violinplot2(V4- V2,[2.3 ],'ViolinColor', [ 0.99 0 0.; 0 0 0.9])
%violinplot2(V5- V2,[2.8 ],'ViolinColor', [ 0.999 0 0; 0 0 0.9]);axis tight
hold on,hold on, line([ 0.7 2.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.05 0.12])
M=[V1-V2, V3-V2, V4-V2];M1=[ ones(length(V1),1)*0.8,ones(length(V1),1)*1.8,ones(length(V1),1)*2.7];
[coef, pval] =corr(M1(~isnan(M)),M(~isnan(M)))
Xvals=[ones(length(M(~isnan(M))),1),M1(~isnan(M))];
[B,BINT,R,RINT,stats] = regress((M(~isnan(M))),Xvals)
%  0.0593    4.4730    0.0379    0.9540
yplot1 = polyval(B(end:-1:1),Xvals(:,2));
hold on,plot(Xvals(:,2),yplot1,'k','Linewidth',1,'Color',[ 0.4 0.4 0.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selN=    allB_N>20;
fr=freq2.freq>3 & freq2.freq<12
V1=   squeeze(nanmean(allPLVsLFP(fr,10,selN)-allPLVsLFP(fr,13,selN),1));
V2=   squeeze(nanmean(allPLVsLFP(fr,11,selN)-allPLVsLFP(fr,14,selN),1));
fr=freq2.freq>30 & freq2.freq<90
V3=   squeeze(nanmean(allPLVsLFP(fr,10,selN)-allPLVsLFP(fr,13,selN),1));
V4=   squeeze(nanmean(allPLVsLFP(fr,11,selN)-allPLVsLFP(fr,14,selN),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.15 0.15])

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.03 0.03])


fr=freq2.freq>3 & freq2.freq<=12
V1=   squeeze(nanmean(allPLVs(fr,3,selN),1));
V2=   squeeze(nanmean(allPLVs(fr,2,selN),1));
fr=freq2.freq>30 & freq2.freq<=90
V3=   squeeze(nanmean(allPLVs(fr,3,selN),1));
V4=   squeeze(nanmean(allPLVs(fr,2,selN),1));
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
hold on,hold on, line([ 0.7 1.8], [0 0],'Color', [ 0.2 0.2 0.2])
ylim([-0.6 0.6])

figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V3- V4,[1.3 ],'ViolinColor', [ 0.3 0.4 0.5; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
ylim([-0.15 0.15])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


