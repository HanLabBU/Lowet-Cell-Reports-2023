

%addpath(genpath('Z:\EricLowet\'))
%addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))

clear all
mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
cd(mainpath);

foldnames = {'DC'; 'DC_multi_level_trial';'place'; 'RES_L'; 'DC_L'; 'other'; 'RES_C'; 'focal_wide';'ramp';'DC_DMD'; 'hip_osc'; 'hip_osc40'} ;
clear allISI allPHASE allCOH   allBURST  lfp_q allname  allangs
allSNR=[];  allRATE=[]; LFPpre=[]; allB=[];allSPEC=[];
loc_D =[];pN=[];pB=[];
   allSTAopt=[];   allSTAB=[];   allSTAs=[];
    allBSR=[];allBBR=[];allSSR=[];allSBR=[];
%%% FOR PLOTTING 40Hz-
mr=0; mmm=0;   % 2  13 21  25 30  35 42 49 54
for foldn= [ 11]% 3 8]%:length(foldnames);   %% folders
    cd(mainpath)
    cd(foldnames{foldn})
    ses=dir('*.mat');
    %problem with LFP 36 37 39
    for ih=[1:length(ses)]  %% sessions    , place 2  3  4  7 14
        
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
            else     FS=820;
            end; else    FS=820;   end;
        
        
        if       length(result.resultS)  == length(unique(result.trial_vec    ))    % check trial correspondence
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
                                        
                                         LFPtr=zscore(fastsmooth(result.lfp.opto{(tt)},3,3,1));
                                          LFPtr=LFPtr+ randn(1,length(LFPtr)).*0;
                                        if ih==34 | ih==29  & foldn ==1
                                            LFPtr=-LFPtr;end
                                        %%%%%%%%%%%%%
                                        %% denoise
                                         Fn = FS/2;FB=[ 58 62];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                          LFPtr= LFPtr-LFPg;
                                          
                                               Fn = FS/2;FB=[ 117 123];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  LFPtr)));
                                     %  LFPtr= LFPtr-LFPg;
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
                                        
                                        dataM.trial2{mm}(1,:) = zscore(double(Vm_trace2));
                                        else
                                            dataM.trial{mm}(1,:) =zscore(double(Vm_trace(1:length(LFPtr))));
                                        dataM.trial{mm}(2,:) =  zscore(double(LFPtr)); 
                                          dataM.trial2{mm}(1,:) = zscore(double(Vm_trace2(1:length(LFPtr))));
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
                                    end
                                end
                            end  %% trials
                         
                              dataM.label= {'A'; 'B'};
          
      %%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%
lW=3.5
hW=9
 lW2=26
  hW2=46
Fa=[];
for id=1:length(dataM.trial2)
FF=ft_preproc_bandpassfilter(zscore(dataM.trial{id}(1,:)),FS,[  lW2 hW2]);FF=abs(hilbert(FF));
Fa=[Fa, FF];end
mF=mean(Fa); mS=std(Fa);

Fl=[];
for id=1:length(dataM.trial2)
FF=ft_preproc_bandpassfilter((dataM.trial{id}(2,:)),FS,[ lW hW]);FF=abs(hilbert(FF));
Fl=[Fl, FF];end
lF=mean(Fa); lS=std(Fa);

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
ROAST1(dataM.spikes{ik}( spikTyp7{ik}==1  ))=1;
dataM.roasterT1{ik}=ROAST1;end
end

for ik=1:size(dataM.roaster,2)
    ROAST1=zeros(1,length(dataM.roaster{ik}));try
ROAST1(dataM.spikes{ik}( spikTyp{ik}==0    ))=1;
dataM.roasterT2{ik}=ROAST1;end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[]; B=[];C=[];D=[];
for id2=1
nn=0;clear allF
for ind=1:1:length(dataM.spikes)
nn=nn+1;try
%plot( dataM.spikes{ind}(spikTyp{ind}==1), dataM.phasespike{ind}(spikTyp{ind}==1),'.r'); hold on,
t= ((dataM.spikes{ind}./1000)>-3);

%t= ((dataM.spikes{ind}./1000)<0    |  (dataM.spikes{ind}./1000)>3);
%t= ((dataM.spikes{ind}./1000)>0.2  & (dataM.spikes{ind}./1000)<1.6);
A=[A,dataM.phasespike{ind}(spikTyp{ind}==0 &  t')]
B=[B,dataM.phasespike{ind}(spikTyp{ind}==1 & t')];
C=[C,  dataM.spikeSNR{ind}(spikTyp2{ind}==0 & t',1)'];
D=[D,  dataM.spikeSNR{ind}(spikTy2{ind}==1 & t',1)'];
end
end
end
 %A=A(C>0);
 %C=C(C>0);
 %B=B(D>3);
 %D=D(D>3);
p= exp(1i.*A) ;P=abs(mean(p)); NT=length(A);T=   P.^2;P= (((1/(NT-1))*((T.*NT)-1)));P(P<0)=0;
p= exp(1i.*B) ;P2=abs(mean(p)); NT2=length(B);T=   P2.^2;P2= (((1/(NT2-1))*((T.*NT2)-1)));;P2(P2<0)=0;
          ISI=nanmean(nall,2);
                            FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                               allBURSTF= sum(ISI(FI))./sum(ISI(FL));


if   1%nanmean(SNR,2) >0 & allBURSTF >0 % & (length(A) +length(B))>10
if 1 

   lW=4.3
hW=8.9
clear angleSig
for id=1:length(dataM.trial2)
FF=ft_preproc_bandpassfilter(zscore(dataM.trial{id}(2,:)),FS,[  lW hW]);
angleSig{id}=angle(hilbert(FF));
end
   lW2=25
hW2=40
Fl=[];clear gamSig
for id=1:length(dataM.trial2)
FF=ft_preproc_bandpassfilter((dataM.trial{id}(2,:)),FS,[ lW2 hW2]);FF=abs(hilbert(FF));
gamSig{id}= (((FF)));end
Fx=[];Fx2=[];
    for xx=1:length(gamSig)
   Fx=[Fx;    gamSig{xx}'];%circ_dist(gamSig{xx},angleSig{xx}) ];
     Fx2=[Fx2;    angleSig{xx}'];
    end
    

    
    
    [cwt_out,frs]=runcwt(dataM, [2 140],FS);
     freq2.freq=frs;
     %figure,imagesc([],f,abs(wt))
     %figure,imagesc(squeeze(nanmean(wavA(2,2,:,:),1)))
     
  %   figure, plot(f, zscore(nanmean(abs(wt).^2,2))); hold on,plot(fr,zscore(nanmean(squeeze(nanmean(wavA(3,1,:,:),1)),2)))
     
     wavA = abs(squeeze(cwt_out));
      wavD = angle(squeeze(cwt_out));
    


  lwin=50;rwin=70;
  CH=1;
  [ STA_opto]= STA_compute_opto(dataM.opto,dataM.opto,CH,lwin,rwin);
  
    [ STA_B]= STA_compute_opto(dataM.roasterT1,dataM.opto,CH,lwin,rwin);

    [ STA_S]= STA_compute_opto(dataM.roasterT2,dataM.opto,CH,lwin,rwin);

     allSTAopt=[allSTAopt , STA_opto];
   allSTAB=[allSTAB , STA_B];
    allSTAs=[allSTAs , STA_S];
%   
%     figure,subplot(1,2,1),imagesc([],freq2.freq, (squeeze(nanmean(abs(data(2,:,:,:)),4))).*repmat((freq2.freq.^0.5)',1, size(data,3))    )
%    axis xy;colormap(jet)
%   subplot(1,2,2), imagesc([],freq2.freq, (squeeze(nanmean(abs(data(1,:,:,:)),4))).*repmat((freq2.freq.^0.5)',1, size(data,3))    )
%    axis xy;colormap(jet)
% %        %   
SBR=[]; BBR=[];
for ij=1:length(dataM.roasterT1)
SBR=[SBR,sum(dataM.roasterT1{ij}(round([FS*1.55:FS*2.5]-1)))];
BBR=[BBR,sum(dataM.roasterT1{ij}(round([FS*0.05:FS*1]-1)))];
end

allSBR=[allSBR,nanmean(SBR)];
allBBR=[allBBR,nanmean(BBR)];

SSR=[]; BSR=[];
for ij=1:length(dataM.roasterT2)
SSR=[SSR,sum(dataM.roasterT2{ij}(round([FS*1.55:FS*2.5]-1)))];
BSR=[BSR,sum(dataM.roasterT2{ij}(round([FS*0.05:FS*1]-1)))];
end

allSSR=[allSSR,nanmean(SSR)];
allBSR=[allBSR,nanmean(BSR)];


    end

%size(freq1.fourierspctrm);
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
                            allN(mr)= length(spang);
                            ISI=nanmean(nall,2);
                            FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                            allISI(:,mr)= ISI(2:end-1);
                               allBURST(:,mr)= sum(ISI(FI))./sum(ISI(FL));
                            [nn nn2]= hist(spang,30);
                          %  allSPEC(:,mr)= fastsmooth(nn./sum(nn),3,1,1);
                               allFS(mr)=FS;
                               allIND{mr}=ih;
                               
                       %   allFreqs{mr} = freq1;
                            allname{mr}=  Cpath;
                            %allangs{mr}= spang;
                                     allTFR{mr}=   squeeze(nanmean(wavA(:,1,:,:),1));     
                            alloptoM{mr}= allopto;
                              allaM(:,:,mr)= (squeeze(nanmean(nanmean(abs(wavA(:,2,:,:)),4),1)));%nanmean(abs(aM),3);
                              AA=exp(squeeze(1i.*wavD(:,1,:,:))); BB=exp(1i.*squeeze(wavD(:,2,:,:)));
                              AA=angle(AA./BB);
                               MM= abs(nanmean(nanmean(exp(1i.*AA),3),1));
                              allaM2(:,:,mr)=MM;
                              
                        
                        CH=1;sp_shift=0;
     %%%%%%%%%%%%
          [PLV_all]= spike_field_ppc_adj(wavD,dataM.roaster,CH,sp_shift) ;
     [PLV_all]= spike_field_ppc_adj(wavD,dataM.roaster,CH,sp_shift) ;
      [PLV_s]= spike_field_ppc_adj(wavD(:,:,:,:),dataM.roasterT2,CH,sp_shift) ;
     [PLV_b]= spike_field_ppc_adj(wavD(:,:,:,:),dataM.roasterT1,CH,sp_shift) ;
       [PLV_s_SH]= spike_field_ppc_adjSHUFFLE_v2(wavD(:,:,:,:),dataM.roasterT2,CH,sp_shift) ;
       [PLV_b_SH]= spike_field_ppc_adjSHUFFLE_v2(wavD(:,:,:,:),dataM.roasterT1,CH,sp_shift) ;
                     
                        allaM3B(:,mr)=PLV_s-PLV_s_SH;%PLV_s;
                        allaM3(:,mr)=PLV_b-PLV_b_SH;%PLV_b;
                                 
                                 
                              allSP(:,mr)=fastsmooth(nanmean(data3,2),25,1,1);
                               allSP(:,mr)= allSP(:,mr)-mean( allSP(:,mr));
                                     allSP2(:,mr)=fastsmooth(nanmean(data4,2),25,1,1);
                               allSP2(:,mr)= allSP2(:,mr)-mean( allSP2(:,mr));
                        %      firing_ratio(mr)= mean(stimN)./mean(baseN);
                              mr
                              
                      %         allSFCM{mr}= nanmean(allSFC,3);
 %allSFCM2{mr}= nanmean(allSFC2,3);
                        end                                            
                        allB= [allB ;sum(ISI(FI))./sum(ISI(FL))];
                        allSNR= [ allSNR; nanmean(SNR,2)];
                        allRATE= [ allRATE; nanmean(RATE,2)];
                    end
                    end
                end
            end;end;end;end;
    
end
%% RESULT ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%
allaM3(allaM3==0)=NaN;allaM3B(allaM3B==0)=NaN;
%llaM3R(allaM3R==0)=NaN;allaMB3R(allaMB3R==0)=NaN;
wind=120;
ax=(-wind:wind)./830;

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_fig_3\'
pheight=150;


   figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
,plot(freq2.freq,nanmean(nanmean(allaM2,3),2),'k','Linewidth',2)
axis tight


  tr=7;VX=dataM.trial2{tr};OPTO1=dataM.trial{tr}(2,:);
      figure('COlor','w','Position',[300 300 500 200],'Renderer', 'painters')
     ,plot((1:length(dataM.trial2{tr}))./FS-1,VX-fastsmooth(VX,2000,1,1),'k');hold on,
      plot(find(dataM.roasterT1{tr})./FS-1, ones(1,length(find(dataM.roasterT1{tr})))+3.5, 'r.','Markersize',12)
       plot(find(dataM.roasterT2{tr})./FS-1, ones(1,length(find(dataM.roasterT2{tr})))+3.5, 'b.','Markersize',12)
       OPTO1=  single(diff(OPTO1)>0.4); %OPTO1(OPTO1<0.5)=NaN;
       ,plot((2:length(dataM.trial2{tr}))./FS-1,(OPTO1)+5,'b-','COlor',[ 0.5 0.5 1]);hold on,
       axis tight
     %xlim([-0.5 1.6])
     
     
ax=(-lwin:rwin)./FS
      figure('COlor','w','Position',[300 300 200 pheight],'Renderer', 'painters')
m=0;for ind=1:size(allSTAs,2)
    if length(find(allSTAs(:,ind)))>0
        m=m+1;
plot(ax(find(allSTAs(:,ind))), ones(1,length(find(allSTAs(:,ind)))) +m,'b.')
    hold on,end
end;axis tight
optmod=nanmedian(allSTAopt,2);optmod=optmod./max(optmod);
hold on,plot(ax,optmod.*(m-2),'k','Linewidth',2.5,'COlor',[0.3 0.8 1])
%xlim([-0.015 0.02])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_opto_THET_opt_trigSS.pdf'])
savefig(gcf, [ savepath 'OSC_opto_THET_opt_trigSS.fig'])

      figure('COlor','w','Position',[300 300 200 pheight],'Renderer', 'painters')
m=0;for ind=1:size(allSTAB,2)
    if length(find(allSTAB(:,ind)))>0
        m=m+1;
plot(ax(find(allSTAB(:,ind))), ones(1,length(find(allSTAB(:,ind)))) +m,'r.')
    hold on,end
end;axis tight
optmod=nanmedian(allSTAopt,2);optmod=optmod./max(optmod);
hold on,plot(ax,optmod.*(m-2),'k','Linewidth',2.5,'COlor',[0.3 0.8 1])
%xlim([-0.015 0.02])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_opto_THET_opt_trigBS.pdf'])
savefig(gcf, [ savepath 'OSC_opto_THET_opt_trigBS.fig'])
% 
% 
%     figure,subplot(3,1,1),imagesc(ax,freq2.freq, (nanmean(allaM,3).*repmat((freq2.freq.^0.65)',1, size(data,3)))    )
%    axis xy;colormap(jet)
%   subplot(3,1,2), imagesc(ax,freq2.freq,smoothn(nanmean(allaM2,3),15)   )
%    axis xy;colormap(jet)
%     subplot(3,1,3),plot(ax,nanmean(allSP,2));axis tight
%    hold on,plot(ax,nanmean(allSP2,2));axis tight

%     
%     figure('COlor','w'),
%     subplot(2,1,1),plot(freq2.freq,nanmean(nanmean(allaM2(:,10:end,1:1:end),3),2),'Linewidth',2);axis tight
%     title('LFP-Vm coherence')
% subplot(2,1,2),
%    plot(freq2.freq,nanmean(allaM3B,2),'Linewidth',2);axis tight; hold on
%     plot(freq2.freq,nanmean(allaMB3R,2),'Linewidth',2);axis tight
%        title('SP-LFP coherence')
       
       figure('COlor','w','Position',[300 300 250 200],'Renderer', 'painters'),
  
    fill_error_area2(freq2.freq,nanmean(allaM3,2),nanstd(allaM3,[],2)./sqrt(size(allaM3,2)), [ 0.5 0.5 0.5]);
            fill_error_area2(freq2.freq,nanmean(allaM3B,2),nanstd(allaM3B,[],2)./sqrt(size(allaM3B,2)), [ 0.5 0.5 0.5]);
    plot(freq2.freq,nanmean(allaM3,2), 'r','Linewidth',2);axis tight; hold on
    plot(freq2.freq,nanmean(allaM3B,2), 'b','Linewidth',2);axis tight
set(gca,'Xscale','log')


fsel= find(freq2.freq>6  &freq2.freq<10);
V1=(nanmean(allaM3(fsel,:),1))'
V2=(nanmean(allaM3B(fsel,:),1))'
  fsel= find(freq2.freq>38  &freq2.freq<42);
V3=(nanmean(allaM3(fsel,:),1))'
V4=(nanmean(allaM3B(fsel,:),1))'
figure('COlor','w','Position', [ 300 400 190 pheight],'Renderer', 'painters')
line([ 0.8 2.3], [ -0.01 -0.01],'COlor', [ 1 1 1 ],'Linewidth',0.5)
violinplot2(V1-V2,[1.2 ],'ViolinColor', [ 0.3 0.5 0.9; 0 0 0.9])
violinplot2(V3-V4,[1.9 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9]);axis tight
line([ 0.9 1.5], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
line([ 1.6 2.2], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
ylim([-0.3 .3]);set(gca,'Xticklabel',[],'Xtick',[])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_8Hzopto_THET_PLV_BAR.pdf'])
savefig(gcf, [ savepath 'OSC_opto8HZ_THET_PLV_BAR.fig'])
   [h,p,ci,stats] = ttest(V3,V4)
  %  0.3818
   [h,p,ci,stats] = ttest(V1,V2)
   %     1.0482e-04,  df=8
   
   
%    fsel= find(freq2.freq>35  &freq2.freq<45);
% V1=(nanmean(allaM3(fsel,:),1))'
% V2=(nanmean(allaM3B(fsel,:),1))'
% figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
% violinplot2(V1-V2,[1.3 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9])
% %violinplot2(V2,[1.8 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
% line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
% ylim([-0.3 .3])
%  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_opto_GAM_PLV_BAR.pdf'])
% %savefig(gcf, [ savepath 'OSC_8HZopto_GAM_PLV_BAR.fig'])
%    [h,p,ci,stats] = ttest(V1,V2)
%    
   
   
    V4=   allSSR';
 V3=   allBSR';
 V1=   allBBR';
 V2=   allSBR';;
 
   
   
    V4=   allSSR';
 V3=   allBSR';
 V1=   allBBR';
 V2=   allSBR';;
 
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
  line([ 0.8 2.3], [ -0.01 -0.01],'COlor', [ 1 1 1 ],'Linewidth',0.5)
violinplot2(V2-V1,[1.2 ],'ViolinColor', [ 0.9 0 0; 0 0 0.9])
violinplot2(V4-V3,[1.9 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
line([ 0.9 1.5], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
line([ 1.6 2.2], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
axis tight;set(gca,'Xticklabel',[],'Xtick',[])
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_opto_THET_FIRING_BAR.pdf'])
%savefig(gcf, [ savepath 'OSC_opto_40Hz_FIRING_BAR.fig'])

 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OSC_opto_THET_FIRING_BAR.pdf'])
savefig(gcf, [ savepath 'OSC_opto_THET_FIRING_BAR.fig'])


[h,p,ci,stats] = ttest(V1,V2)
%SB-SSSB-SSGammaTheta

%%TFR
%%%%%TFR
clear allM
for id=1:length(allTFR)
   
       A=  allTFR{id}(:,:);
     B=nanmean(A(:,100:800),2);
 allM(:,1:size(A,2),id)=    smooth2a(bsxfun(@minus,A,B)./bsxfun(@plus,A,B),1,55) ;
end
allM(allM==0)=NaN;
       figure('COlor','w','Position', [ 300 400 250 120],'Renderer', 'painters')
imagesc(((1:size(allM,2))./828)-1,freq2.freq,nanmean(allM,3) );ylim([3 100])
axis xy;colormap(jet);set(gca,'Clim',[-0.2 0.2])
   hold on, ,plot((1:length(OPTO1))./FS-1,(OPTO1./0.04)+15.5,'b-','COlor',[ 0.3 0.8 1],'Linewidth',1);hold on,
xlim([ -0.5 2.85])


