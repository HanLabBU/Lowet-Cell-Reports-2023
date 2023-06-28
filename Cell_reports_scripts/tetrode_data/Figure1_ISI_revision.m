cd('C:\Users\hanlab\Dropbox (BOSTON UNIVERSITY)\rat_hip_data\')

clear all
figure('COlor','w','Position', [ 300 200 250 200],'Renderer', 'painters') ;

ses=dir('*.mat')
mr=0;mr2=0;
clear allZ allZ2
for ij=[  8 9  1:6 11:length(ses)]

  % clearvars Tet* S S2 spikes FF placeF placeF2  xydata   spikeTable 
   clearvars -except ses ij allplace mj allCC allZ allZ2 mr  mr2  allISI_1  allISI_2 allNT     allfiring
    load(ses(ij).name)

    
clear S Vx

if  ij <5
  Vx=Tet18LFP{1}(1:end,2);%-Tet18LFP{2}(1:end,2);

elseif ij==5
     Vx=Tet12LFP{2}(1:end,2);% 
elseif ij==6
      Vx=Tet12LFP{2}(1:end,2);% 
elseif ij ==7
    Vx=Tet3LFP{2}(1:end,2);% 
    elseif ij ==8
    Vx=Tet18LFP{2}(1:end,2);% 
        elseif ij ==9
    Vx=Tet5LFP{2}(1:end,2);% 
        elseif ij ==11
    Vx=Tet5LFP{2}(1:end,2);% 
      elseif ij ==12
    Vx=Tet22LFP{2}(1:end,2);% 
end

vsig=Vx;;
clearvars Tet*


lfp=[];

lfp.trial{1}= vsig';
lfp.time{1}= (1:size(vsig,1))./1000;
for id=1:size(vsig,2)
lfp.label{id}= num2str(id);
end


selT=round(xydata(:,1).*1000);
%%%%%%%%%%%%

speedsel=selT((velocity(10:end-100)<20));
speedsel(speedsel==0)=[];
% for ix=1:length(speedsel)
% wavD2(:,speedsel(ix)-30:speedsel(ix)+30)=NaN;
% end
for ind=1:length(spikes)
wind=10;
tria=1;clear spikTyp
   s= spikes{ind}(:).*1000  ;
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
    

bsp=round(spikes{ind}(spikTyp{1}==1).*1000);bsp=bsp(1:end-10);
 ssp=round(spikes{ind}(spikTyp{1}==0).*1000);  ssp=ssp(1:end-10); 
  sp=round(spikes{ind}(spikTyp{1}>-8).*1000);  sp=sp(1:end-10); 
[n12 n2]=hist(diff(bsp),[2:2:800 ]);
[n1 n2]=hist(diff(sp),[2:2:800 ]); %n12=n1+n12;
nt=sum(n1);
if (length(sp)./max(sp)).*1000  <10
mr=mr+1;
     allISI_1(:,mr)=n12;
         allISI_2(:,mr)=n1;
            allNT(mr)=nt;
           allfiring(mr)=(length(sp)./max(sp)).*1000; 
           
           
           
           Vx=(bsp);
Vx1=(-Vx(1:end-2)+Vx(2:end-1))';%Vx1(Vx1<2)=[];
Vx2= (-Vx(2:end-1)+Vx(3:end))';;%Vx2(Vx2<2)=[];
  Vx=(ssp);
Vx1S=(-Vx(1:end-2)+Vx(2:end-1))';%Vx1S(Vx1S<2)=[];
Vx2S= (-Vx(2:end-1)+Vx(3:end))';;%Vx2S(Vx2S<2)=[];
gapN=80;
Vx1=Vx1(1:gapN:end);Vx2=Vx2(1:gapN:end);
Vx1S=Vx1S(1:gapN:end);Vx2S=Vx2S(1:gapN:end);
if length(Vx1)>100
    Vx1=Vx1(1:100);Vx2=Vx2(1:100);end
if length(Vx1S)>100
    Vx1S=Vx1S(1:100);Vx2S=Vx2S(1:100);end
scatter1 = scatter(Vx1+randn(1,length(Vx1)).*0.3,Vx2+randn(1,length(Vx2)).*0.3,8,'MarkerFaceColor','r','MarkerEdgeColor','r'); 
scatter1.MarkerFaceAlpha = .7;
scatter1.MarkerEdgeAlpha = .1;hold on,
set(gca,'Xscale','log','Yscale','log');axis tight
scatter1 = scatter(Vx1S+randn(1,length(Vx1S)).*0.3,Vx2S+randn(1,length(Vx2S)).*0.3,8,'MarkerFaceColor',[0 0 0.5],'MarkerEdgeColor',[0 0 0.5]); 
scatter1.MarkerFaceAlpha = .4;
scatter1.MarkerEdgeAlpha = .1;hold on,
set(gca,'Xscale','log','Yscale','log');axis tight
axis([ 2 700 2 700])
pause(0.1)
end
end
end


% 1.32Hz, 0.15
intsel=find(n2< 11);
    figure('COlor','w','Position',[ 300 400 200 160],'Renderer', 'painters'),
    bar(n2,nanmean(allISI_2./allNT,2),'Facecolor',[0.5 0.5 0.5],'FaceAlpha',0.5); hold on,
    bar(n2(intsel),nanmean(allISI_2(intsel,:)./allNT,2),'Facecolor',[1 0 0],'FaceAlpha',0.5);axis tight
    xlim([ 2 300])
set(gca,'Xscale','log') 



