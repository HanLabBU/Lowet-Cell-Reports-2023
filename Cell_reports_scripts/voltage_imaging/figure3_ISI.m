

%addpath(genpath('Z:\EricLowet\'))
addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))
addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\hip_scripts\'))

clear all
mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
cd(mainpath);

foldnames = {'DC'; 'DC_multi_level_trial';'place'; 'RES_L'; 'DC_L'; 'other'; 'RES_C'; 'focal_wide';'ramp';'DC_DMD'; 'hip_osc'; 'hip_osc40'} ;
clear allISI allPHASE allCOH   allBURST  lfp_q allname  allangs
allSNR=[];  allRATE=[]; LFPpre=[]; allB=[];allSPEC=[];
loc_D =[];pN=[];pB=[];
mr=0; mmm=0;   % 2  13 21  25 30  35 42 49 54
for foldn= [1 3 8]% 3  8 10 11 12]%:length(foldnames);   %% folders
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
                                B=sum(result.resultS{id}.roaster(neuron,50:800 ))./(size(result.resultS{id}.roaster(neuron,50:800),2)./FS);
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
                                            Fn = FS/2;FB=[ 5 35];
                                         [B, A] = butter(3, [min(FB)/Fn max(FB)/Fn]);
                        
                                     
                                        LFPphase= angle(hilbert(filtfilt(B,A,  LFPtr)));
                                        Vmphase= angle(hilbert(filtfilt(B,A,  Vm_trace)));
                                      Vmpow= abs(hilbert(filtfilt(B,A,  Vm_trace)));
  
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
                                        dataM.phaseVm{mm}= Vmpow; else
                                        dataM.phaseLFP{mm}= LFPphase; 
                                        dataM.phaseVm{mm}= Vmpow(1:length(LFPphase));
                                        end
                                          if max(SpikeTimes) <= length(LFPphase) 
                                         dataM.phasespike{mm}= LFPphase(SpikeTimes  );else
                                          dataM.phasespike{mm}=  LFPphase(SpikeTimes( (SpikeTimes) <= length(LFPphase) )  ) ;
                                          end
                                          dataM.phasespikeVm{mm}= Vmpow(SpikeTimes  );
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
windSA=round(20.*(FS./1000));windSB=round(5.*(FS./1000));
windS2A=round(5.*(FS./1000));windS2B=round(20.*(FS./1000));
base=nanmean(nanmean(allW2(wind2-windSA:wind2-windSB,:),1));
 A= nanmean(allW2(wind2+windS2A:wind2+windS2B,:),1)-nanmean(allW2(wind2-windSA:wind2-windSB,:),1);
 sp_height= nanmedian(nanmedian(allW(wind2+1,:),1));
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









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for ik=1:size(dataM.roaster,2)
%     ROAST1=zeros(1,length(dataM.roaster{ik}));try
% ROAST1(dataM.spikes{ik}(spikTyp6{ik}>1& spikTyp6{ik}<5  & spikTyp7{ik}>0.6 & spikTyp7{ik}<41.9  ))=1;
% dataM.roasterT1{ik}=ROAST1;end
% end
% 
% for ik=1:size(dataM.roaster,2)
%     ROAST1=zeros(1,length(dataM.roaster{ik}));try
% ROAST1(dataM.spikes{ik}( spikTyp6{ik}>1& spikTyp6{ik}<5 & spikTyp7{ik}<0.1))=1;
% dataM.roasterT2{ik}=ROAST1;end
% end

for ik=1:size(dataM.roaster,2)
    ROAST1=zeros(1,length(dataM.roaster{ik}));try
ROAST1(dataM.spikes{ik}(spikTyp7{ik}==1 ))=1;
dataM.roasterT1{ik}=ROAST1;end
end

for ik=1:size(dataM.roaster,2)
    ROAST1=zeros(1,length(dataM.roaster{ik}));try
ROAST1(dataM.spikes{ik}( ~(spikTyp7{ik}==1)  ))=1;
dataM.roasterT2{ik}=ROAST1;end
end

% for ik=1:size(dataM.roaster,2)
%     ROAST1=zeros(1,length(dataM.roaster{ik}));try
% ROAST1(dataM.spikes{ik}(spikTyp6{ik}==0 & spikTyp7{ik}>0.3 & spikTyp7{ik}<41.9  ))=1;
% dataM.roasterT1{ik}=ROAST1;end
% end
% 
% for ik=1:size(dataM.roaster,2)
%     ROAST1=zeros(1,length(dataM.roaster{ik}));try
% ROAST1(dataM.spikes{ik}( spikTyp6{ik}==0 & spikTyp7{ik}<0))=1;
% dataM.roasterT2{ik}=ROAST1;end
% end
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
% C=[C,  dataM.spikeSNR{ind}(spikTyp2{ind}==0 & t',1)'];
% D=[D,  dataM.spikeSNR{ind}(spikTy2{ind}==1 & t',1)'];
end
end
end
 %A=A(C>0);
 %C=C(C>0);
 %B=B(D>3);
 %D=D(D>3);
% p= exp(1i.*A) ;P=abs(mean(p)); NT=length(A);T=   P.^2;P= (((1/(NT-1))*((T.*NT)-1)));P(P<0)=0;
% p= exp(1i.*B) ;P2=abs(mean(p)); NT2=length(B);T=   P2.^2;P2= (((1/(NT2-1))*((T.*NT2)-1)));;P2(P2<0)=0;
          ISI=nanmean(nall,2);
                            FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
                               allBURSTF= sum(ISI(FI))./sum(ISI(FL));


if  nanmean(SNR,2) >4   & FS>400% &length(A)>3 &  length(B)>3
if 1 
    
    
    
    
Vx1=[];Vx2=[];
for id=1:length(result.resultS)
    try
Vx=result.resultS{id}.spike_idx{1}( spikTyp{id}>-1);
Vx1= [Vx1;(-Vx(1:end-2)+Vx(2:end-1)).*(1000/FS)];
Vx2= [Vx2;(-Vx(2:end-1)+Vx(3:end)).*(1000/FS)];end
end

VxD=[];;
for id=1:length(result.resultS)
    try
Vx=result.resultS{id}.spike_idx{1}( spikTyp{id}>-1);
VxD= [VxD;(-Vx(1:end-1)+Vx(2:end)).*(1000/FS)];
end
end

VxP=[];VxP2=[];
for id=1:length(result.resultS)
   % try
Vx=result.resultS{id}.spike_idx{1}( spikTyp{id}>-1);
tg=dataM.phaseVm{id};
tg2=nanfastsmooth(dataM.trial{id}(1,:),30,1,1);

VxP= [VxP, tg(Vx(1:end-1))];
VxP2=[VxP2,tg2(Vx(1:end-1))- tg2(abs(Vx(1:end-1)-40)+1)];
%;end
end


Vx1B=[];Vx2B=[];
for id=1:length(result.resultS)
    try
Vx=result.resultS{id}.spike_idx{1}( ~(spikTyp{id}==1));
Vx1B= [Vx1B;(-Vx(1:end-2)+Vx(2:end-1)).*(1000/FS)];
Vx2B= [Vx2B;(-Vx(2:end-1)+Vx(3:end)).*(1000/FS)];end
end

[n12 n2]=hist(Vx2B,[2:2:400 ]);
[n1 n2]=hist(Vx2,[2:2:400 ]); %n12=n1+n12;
nt=sum(n1+n12);
% 

   figure('COlor','w','Position', [ 300 200 250 200],'Renderer', 'painters') ;
   subplot(1,2,1)
scatter1 = scatter(Vx1+randn(1,length(Vx1)).*0.3,Vx2+randn(1,length(Vx2)).*0.3,8,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
scatter1.MarkerFaceAlpha = .7;
scatter1.MarkerEdgeAlpha = .1;hold on,
scatter1 = scatter(Vx1B'+randn(1,length(Vx1B)).*0.3,Vx2B'+randn(1,length(Vx2B)).*0.3,8,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
set(gca,'Xscale','log','Yscale','log')
axis tight;axis([ 2 1500 2 1500])
line([ 2 1500],[2 1500],'COlor',[ 0.8 0 0],'Linewidth',1)
box on
   subplot(1,2,2)
hist(Vx1,[1:2:1200],'k')
xlim([ 2 1000])
set(gca,'Xscale','log')
% intb1=[];
% for id=1:length(result.resultS)
%     try
%     for onsets=    find(spikTyp5{id}==1 & spikTyp7{id}==1)
%         onsets2= find(spikTyp5{id}>1);
%          onsets2=   onsets2(find(onsets2==onsets+1));
%          if ~isempty(onsets2)
%  intb1=[intb1,  diff(result.resultS{id}.spike_idx{1}([onsets onsets2])).*(1000/FS) ]; end   
%     end;end;end
% 
% intb2=[];
% for id=1:length(result.resultS)
%     try
%     for onsets=    find(spikTyp6{id}==3 & spikTyp2{id}==1)
%         onsets2= find(spikTyp6{id}==4);
%          onsets2=   onsets2(find(onsets2==onsets+1));
%          if ~isempty(onsets2)
%  intb2=[intb2,  diff(result.resultS{id}.spike_idx{1}([onsets onsets2])).*(1000/FS) ]; end   
%     end;end;end    
% 
% intb3=[];
% for id=1:length(result.resultS)
%     try
%     for onsets=    find(spikTyp6{id}==4 & spikTyp2{id}==1)
%         onsets2= find(spikTyp6{id}==5);
%          onsets2=   onsets2(find(onsets2==onsets+1));
%          if ~isempty(onsets2)
%  intb3=[intb3,  diff(result.resultS{id}.spike_idx{1}([onsets onsets2])).*(1000/FS) ]; end   
%     end;end;end   
% 
% intb4=[];
% for id=1:length(result.resultS)
%     try
%     for onsets=    find(spikTyp6{id}==5 & spikTyp2{id}==1)
%         onsets2= find(spikTyp6{id}==6);
%          onsets2=   onsets2(find(onsets2==onsets+1));
%          if ~isempty(onsets2)
%  intb4=[intb4,  diff(result.resultS{id}.spike_idx{1}([onsets onsets2])).*(1000/FS) ]; end   
%     end;end;end  

k1=[];k2=[];k3=[];k4=[];k5=[];;k6=[];
for id=1:length(result.resultS)
    try
        t= ((dataM.spikes{id}./1000)<1 | (dataM.spikes{id}./1000)>2.5 )';
 k1=[k1,find(spikTyp{id}>-1& t)]; 
 k2=[k2,find(spikTyp7{id}==1 & t)];
%  k3=[k3,find(spikTyp6{id}==2& spikTyp2{id}==1)];
%   k4=[k4,find(spikTyp6{id}==3& spikTyp2{id}==1)];
%     k5=[k5,find(spikTyp6{id}==4& spikTyp2{id}==1)];
%       k6=[k6,find(spikTyp6{id}>4& spikTyp2{id}==1)];
    end
end
%      abs( mean( (Fx./mean(Fx)).*exp(1i.*Fx2)))
%    figure,%hist(Fx)
%     for id=1:size(angleSig,2)
%         plot(angleSig{id}(1:10:end),gamSig{id}(1:10:end),'.k');hold on,end
%     
 


%   
%     figure,subplot(1,2,1),imagesc([],freq2.freq, (squeeze(nanmean(abs(data(1,:,:,:)),4))).*repmat((freq2.freq.^0.5)',1, size(data,3))    )
%    axis xy;colormap(jet)
%   subplot(1,2,2), imagesc([],freq2.freq, (squeeze(nanmean(abs(dataB(1,:,:,:)),4))).*repmat((freq2.freq.^0.5)',1, size(data,3))    )
%    axis xy;colormap(jet)
% %        
 
%     figure,subplot(1,2,1),imagesc([],freq2.freq, abs(squeeze(nanmean((data(2,:,:,:)),4))).*repmat((freq2.freq.^0)',1, size(data,3))    )
%    axis xy;colormap(jet)
%   subplot(1,2,2), imagesc([],freq2.freq, abs(squeeze(nanmean((data(1,:,:,:)),4))).*repmat((freq2.freq.^0)',1, size(data,3))    )
%    axis xy;colormap(jet)
   
    end

   
    
 

  
%      figure,
%     subplot(1,1,1)
%     imagesc([],freq2.freq,abs(nanmean((aM(:,:,1:1:end)),3))); axis xy

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
                            allADP(mr)=nanmean(ADP);
%                             ISI=nanmean(nall,2);
%                             FI=(nax> 0.005 & nax<0.025);FL=(nax> 0.035 & nax<0.055);
%                             allISI(:,mr)= ISI(2:end-1);
%                                allBURST(:,mr)= sum(ISI(FI))./sum(ISI(FL));
%                             [nn nn2]= hist(spang,30);
                          %  allSPEC(:,mr)= fastsmooth(nn./sum(nn),3,1,1);
                               allFS(mr)=FS;
                               allIND{mr}=ih;
                          %     allspADP(:,mr)= spADP;
                         %       allavADP(mr)= avADP;
                       %   allFreqs{mr} = freq1;
                            allname{mr}=  Cpath;
                            %allangs{mr}= spang;
                           % alloptoM{mr}= allopto;
                          
                            allISI_1(:,mr)=n1;
                             allISI_2(:,mr)=n12;
                           allNT(mr)=nt;
                         
%                            allINTB(1,mr)=mean(intb1);
%                             allINTB(2,mr)=mean(intb2);
%                              allINTB(3,mr)=mean(intb3);
%                               allINTB(4,mr)=mean(intb4);
%                                  
                              allProbs(1,mr)= length(k1);
                         allProbs(2,mr)= length(k2);
%                          allProbs(3,mr)= length(k3);
%                          allProbs(4,mr)= length(k4);
%                          allProbs(5,mr)= length(k5);
%                          allProbs(6,mr)= length(k6);
                              
                         %     allSP(:,mr)=fastsmooth(nanmean(data3,2),25,1,1);
                           %    allSP(:,mr)= allSP(:,mr)-mean( allSP(:,mr));
                            %         allSP2(:,mr)=fastsmooth(nanmean(data4,2),25,1,1);
                          %     allSP2(:,mr)= allSP2(:,mr)-mean( allSP2(:,mr));
                        %      firing_ratio(mr)= mean(stimN)./mean(baseN);
                              mr
                              
                      %         allSFCM{mr}= nanmean(allSFC,3);
 %allSFCM2{mr}= nanmean(allSFC2,3);
                        end                                            
                    %    allB= [allB ;sum(ISI(FI))./sum(ISI(FL))];
                     %   allSNR= [ allSNR; nanmean(SNR,2)];
                  allRATE= [ allRATE; nanmean(RATE,2)];
                    end
                    end
                end
            end;end;end;end;
    
    end
end
%% RESULT ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%
% allaM3(allaM3==0)=NaN;allaM3B(allaM3B==0)=NaN;
% allaM3R(allaM3R==0)=NaN;allaMB3R(allaMB3R==0)=NaN;
wind=150;
ax=(-wind:wind)./830;
% 5.25,Hz 0.63

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\nat_figure1\'
pheight=160;




   overall_prob_burst=  mean(allProbs(2,:)./allProbs(1,:),2);
   overall_prob_burst_STD=  std(allProbs(2,:)./allProbs(1,:),[],2)./sqrt(size(allProbs,2));

Neuron_b=nansum(((allProbs(2,:))./nanmean(allProbs(1,:),1))>0.40)
Neuron_t=size(allProbs,2);
Neuron_b./Neuron_t
%figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
p1=   overall_prob_burst;%length(find(allProbs(2,:)./mean(allProbs(1,:),1)>0.15));
p2=1-   overall_prob_burst;%length(find(allProbs(2,:)./mean(allProbs(1,:),1)<=0.15));
 figure('COlor','w','Position',[ 300 400 250 pheight],'Renderer', 'painters'),
h=pie(  [ 28.19 100-28.19],[ 0 1] )
% h=pie() output is a vector of alternating patch and text handles. 
% Isolate the patch handles
patchHand = findobj(h, 'Type', 'Patch'); 
% Set the color of all patches using the nx3 newColors matrix
set(patchHand, {'FaceColor'}, mat2cell([0.6 0  0;0  0 0.6], ones(2,1), 3))
% Or set the color of a single wedge
patchHand(2).FaceAlpha = 0.3;
patchHand(1).FaceAlpha = 0.3;
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Pie_prob_neuron.pdf'])
savefig(gcf, [ savepath 'Pie_prob_neuron.fig'])

if 0
x = [1 3 0.5 2.5 2];
explode = [0 1 0 0 0];
figure
pie(x,explode)
colormap jet
end

%    figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
% ,plot(freq2.freq,nanmean(nanmean(allaM2,3),2),'k','Linewidth',2)
% axis tight


    figure('COlor','w','Position',[ 300 400 250 pheight],'Renderer', 'painters'),
    bar(n2,sum(allISI_2+allISI_1,2)./sum(allNT),'Facecolor',[0.5 0.5 0.5],'FaceAlpha',0.5); hold on,
    bar(n2,sum(allISI_1,2)./sum(allNT),'Facecolor',[1 0 0],'FaceAlpha',0.5);axis tight
    xlim([ 3 300])
set(gca,'Xscale','log') 
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'ISI_figure.pdf'])
savefig(gcf, [ savepath 'ISI_figure.fig'])




if  0
    ,subplot(1,1,1),imagesc(ax,freq2.freq, (nanmean(allaM,3).*repmat((freq2.freq.^0.5)',1, size(data,3)))    )
   axis xy;colormap(jet)
   set(gca,'Clim',[ 10 65])
      figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
  subplot(1,1,1), imagesc(ax,freq2.freq, (nanmean(allaMB,3).*repmat((freq2.freq.^0.5)',1, size(data,3)))    )
   axis xy;colormap(jet)
  set(gca,'Clim',[ 10 65])
   
  
  allaMB2=allaMB;%bsxfun(@minus,allaMB,nanmean(nanmean(allaMB(:,40:140,:),1),2));
    allaM2=allaM;%bsxfun(@minus,allaM,nanmean(nanmean(allaMB(:,40:140,:),1),2));

   SC=(freq2.freq.^0.55)';tsel=40:140;
 figure('COlor','w','Position',[ 300 400 250 200],'Renderer', 'painters'),
  plot(freq2.freq,nanmean(nanmean(allaM2(:,tsel,:),3),2).*SC,'r'); hold on,
  fill_error_area2(freq2.freq,mean(nanmean(allaM2(:,tsel,:),3),2).*SC,(nanstd(nanmean(allaM2(:,tsel,:),2).*SC,[],3))./sqrt(size(allaM,3)), [ 0.5 0.5 0.5]);
  plot(freq2.freq,nanmean(nanmean(allaMB2(:,tsel,:),3),2).*SC,'b')
    fill_error_area2(freq2.freq,mean(nanmean(allaMB2(:,tsel,:),3),2).*SC,(nanmean(nanstd(allaMB2(:,tsel,:).*SC,[],3),2))./sqrt(size(allaMB,3)), [ 0.5 0.5 0.5]);
  axis tight
  set(gca,'Xscale','log')
  
  
  
    B=squeeze(nanmean(allaM2(:,40:140,:),2)).*SC ;A=squeeze(nanmean(allaMB2(:,40:140,:),2)).*SC;
 figure('COlor','w','Position', [ 300 400 240 200])
 fr=freq2.freq>3 & freq2.freq<=10; V1=nanmean(B(fr,:))'; V2=nanmean(A(fr,:))';
 fr=freq2.freq>30 & freq2.freq<=80; V3=nanmean(B(fr,:))'; V4=nanmean(A(fr,:))';
fr=freq2.freq>90 & freq2.freq<=180; V5=nanmean(B(fr,:))'; V6=nanmean(A(fr,:))'
V1mF=nanmean(V1);V2mF=nanmean(V2);  V3m=nanmean(V3);V4m=nanmean(V4);  V5m=nanmean(V5);V6m=nanmean(V6); 
       M=[V2'  , V1' ];
boxplot( M   ,[ ones(length(V1),1);ones(length(V2),1).*2;],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
hold on,          %  [h,p,ci,stats] = ttest2(V1,V2);
set(gca,'Xtick', [ 1 2 ],'Xticklabel', {'SS';'CS'})
xlim([0.5 2.5])

  
  
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
       
end