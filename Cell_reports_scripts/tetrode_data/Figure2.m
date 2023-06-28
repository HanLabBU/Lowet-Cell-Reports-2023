

cd('C:\Users\hanlab\Dropbox (BOSTON UNIVERSITY)\rat_hip_data\processed\')

load('run1_placeF_TFR_outboundnew','allTFR','allTFRA','allSp','allI','allID','firing_rate1','firing_rate3','firing_rate2','frs','allplace','allCM','allMR')
 allTFRo=allTFR; allTFRAo=allTFRA;
 allplaceo=allplace;
 allSpo=allSp;
 allIo= abs(allI);
 allIDo= allID;
 firing_rate1o= firing_rate1;
 firing_rate2o= firing_rate2;
 firing_rate3o= firing_rate3;
 allCMo= allCM./1;
 allMRo=allMR;
 load('run1_placeF_TFR_inboundnew','allTFR','allTFRA','allSp','allI','allID','firing_rate1','firing_rate3','firing_rate2','frs','allplace','allCM','allMR')
%%%%%%%%%%%%%%%%%%%%%%%%%

for id=1:length(firing_rate3)
   if  allMR(id,3) > allMRo(id,3)  
   else
     firing_rate3(id)  =firing_rate3o(id) ;
      firing_rate1(id)  =firing_rate1o(id) ;
       firing_rate2(id)  =firing_rate2o(id) ;
       allTFR(:,:,id,:)= allTFRo(:,1:185,id,:);
          allTFRA(:,:,id,:)= allTFRAo(:,1:185,id,:);
       allplace(:,id,:)= allplaceo(1:185,id,:);
       allSp(id,:)= allSpo(id,:);
      allI(id,:)= abs(allIo(id,:));
       allCM(id,:)=allCMo(id,:)./1;
   end
    
end
 allCM= allCM./4;
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\natfig4_pf\'

pheight=160;
if 0
clear errorF fitval errorF1 fitval1 errorF2 fitval2
for ix=1:size(allplace,2)
    ix
x=1:size(allplace,1);
y= allplace(:,ix,3);y=y-min(y);y=y./max(y);
y1= allplace(:,ix,1);y1=y1-min(y1);y1=y1./max(y1);
y2= allplace(:,ix,2);y2=y2-min(y2);y2=y2./max(y2);
     pinit=[0,1,50 ,20];
        LB = [0,0.1 ,1, 1];
     UB = [0.5,1 ,200, 50];
  [pbest,delta_p,ybest,delta_y,chi0,exitflag,output]=easyfit(x,y,pinit,@gausfun, LB,UB);%,'-inf','+inf','np');
errorF(ix,1)=(sum((y-ybest).^2 ));%);pause(0.001)
rr=corrcoef(y,ybest);
errorF(ix,2)=rr(1,2).^2;
fitval(:,ix)=pbest;
  [pbest1,delta_p,ybest1,delta_y,chi0,exitflag,output]=easyfit(x,y1,pinit,@gausfun, LB,UB);%,'-inf','+inf','np');
errorF1(ix,1)=(sum((y1-ybest1).^2 ));%);pause(0.001)
fitval1(:,ix)=pbest1;
rr=corrcoef(y1,ybest1);
errorF1(ix,2)=rr(1,2).^2;
  [pbest2,delta_p,ybest2,delta_y,chi0,exitflag,output]=easyfit(x,y2,pinit,@gausfun, LB,UB);%,'-inf','+inf','np');
errorF2(ix)=(sum((y2-ybest2).^2 ));%);pause(0.001)
fitval2(:,ix)=pbest2;
rr=corrcoef(y2,ybest2);
errorF2(ix,2)=rr(1,2).^2;

end
end
allI=abs(allI);
M=allplace(:,:,3);[ p1 p2]=max(M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modIndex=allSp(:,1)'<1 & allI(:,3)'>0.2 &  firing_rate3<20 %& allMR(:,1)'>0.0 & allMR(:,2)'>0.0% & p2>3 & p2<170% & (nanmean(nanmean(allplace(:,:,1),3),1)>5)'   %& (fitval(3,:)>40   & fitval(3,:)<90  )';
M=allplace(:,modIndex,3);
M1=allplace(:,modIndex,1);
% for ij=1:size(M,2)
%    M1(:,ij)=fastsmooth(M1(:,ij),3,3,1); 
% end
% for ij=1:size(M,2)
%    M2(:,ij)=fastsmooth(M2(:,ij),3,3,1); 
% end
M2=allplace(:,modIndex,2);
clear p5
for ij=1:size(M,2)
 p5(ij)= mean(find( M1(:,ij)>prctile(M1(:,ij),98) ));
end
[ p1 p2]=max(M1);[ p3 p6]=max(M1)
[ p3 p4]=max(M2)
[ b1 b2]=sort(p2)

figure('COlor','w','Position',[ 300 300 200 pheight],'Renderer', 'painters')
MM=bsxfun(@minus, M1(:,b2)  , prctile(M1(:,b2),30));
imagesc(MM');%bsxfun(@rdivide  ,MM, mean(MM))')
imagesc(zscore(M1(:,b2))')
set(gca,'CLim',[ -3 3])
%line([0 length(b2)], [0 length(b2)],'Linewidth',2,'Linestyle','--')
colormap(bone);
axis xy
 %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PF_M1.pdf'])
%savefig(gcf, [ savepath 'PF_M1.fig'])

figure('COlor','w','Position',[ 300 300 200 pheight],'Renderer', 'painters')
imagesc(zscore(M2(:,b2))')
set(gca,'CLim',[ -3  3])
%line([0 length(b2)], [0 length(b2)],'Linewidth',2,'Linestyle','--')
colormap(bone);
axis xy
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PF_M2.pdf'])
%savefig(gcf, [ savepath 'PF_M2.fig'])V1=allI(modIndex,1);V2=allI(modIndex,2);
V1=allI(modIndex,1);V2=allI(modIndex,2);
F1=firing_rate1(modIndex);F2=firing_rate2(modIndex);
%V1=V1((F1-F2)>=0);V2=V2((F1-F2)>=0);
   M5=[V1  , V2 ]; 
%INFORMATION
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.9 0.4 0.6; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'allI.pdf'])
%savefig(gcf, [ savepath 'allI.fig'])
 [h,p,ci,stats] =ttest(V1,V2)
%MATCH,  tstat: 8.2913, df: 103, 4.4853e-13


%MATCH,  tstat: 8.2913, df: 103, 4.4853e-1
V1=allSp(modIndex,1);V2=allSp(modIndex,2);
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.4 0.7 0.6; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'allSp.pdf'])
% 
 [h,p,ci,stats] = ttest(V1,V2)
%MATCH, tstat: -10.3310,  df: 134, 0
V1=allCM(modIndex,1)-allCM(modIndex,3);V2=allCM(modIndex,2)-allCM(modIndex,3) ;
cutoff= 10;
V1(V1>cutoff)=NaN;V1(V1<-cutoff)=-NaN;V2(V2>cutoff)=NaN;V2(V2<-cutoff)=-NaN;
   M5=[V1  , V2 ];
figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters')
violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.4 0.1 0.6; 0 0 0.9])
line([ 0.7 1.8], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',0.5)
axis tight;set(gca,'Xtick',[],'Xticklabel',[])
axis tight;%ylim([-15 15])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'allCM.pdf'])
%savefig(gcf, [ savepath 'allCM.fig']);
 [h,p,ci,stats] =ttest(V1,V2)
 %90,tstat: -3.6808,p=      3.9588e-04
 %MATCH, df: 134, tstat: -4.5004,  1.4574e-05
% [h,p,ci,stats] =ttest((allCM(modIndex,1)./allCM(modIndex,2)) -1)
 %% -22.36cm , +- SE 8.7
 
% % df=146, 0.011
% wd=12;clear FPF FPF2 PP
% mr=0;
% for id=1:length(b1)   
%     if b1(id)-wd>0 & b1(id)+wd<size(M,1)
%         mr=mr+1;
% FPF(:,mr)= (M1( b1(id)-wd:b1(id)+wd,b2(id)));    
% FPF2(:,mr)=(M2( b1(id)-wd:b1(id)+wd,b2(id)));    
% PP(:,mr)= zscore(M(b1(id)-wd:b1(id)+wd,b2(id)));
%     end
% end
% ax=-wd:wd;
% FPF2=bsxfun(@minus, FPF2,min(FPF2,1))
% FPF2=bsxfun(@rdivide, FPF2,max(FPF2,[],1))
% FPF=bsxfun(@minus, FPF,min(FPF,1))
% FPF=bsxfun(@rdivide, FPF,max(FPF,[],1))
% 
% %legend all SS CS
%  % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_PF.pdf'])
% %savefig(gcf, [ savepath 'FR_PF.fig'])
% 
% 
% 
%  %  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Burstprob.pdf'])
% %savefig(gcf, [ savepath 'Burstprob.fig'])
% 
% 
% %%%%%%%%%%%%%%%%%
% 
% 
% aTF1=allTFR(:,:,modIndex,:);
% 
% clear PF PF2
% mr=0;
% for id=1:length(b1)   
%     if b1(id)>wd & b1(id)+wd<size(aTF1,2)
%         mr=mr+1;
% PF(:,:,mr)= ( aTF1(:, b1(id)-wd:b1(id)+wd,b2(id),1)  );    
% PF2(:,:,mr)=( aTF1(:, b1(id)-wd:b1(id)+wd,b2(id),2)  );   
% 
%     end
% end
% 
% %%%%%%%%%%%%%
% clear allCL
% for a1=1:size(PF,3)
%     for a2=1:size(PF,1);
%   allCL(a1,a2)= (corr(PF(a2,:,a1)',-(ax)'));
% 
%     end
% end
% 
% figure,plot(nanmean(allCL,1))
% 
% 
% 
% %%%%%%%%%%%%%%%%%%
% refFR2=nanmean((FPF2),2);
% refFR2=refFR2-min(refFR2);
% refFR2=fastsmooth((refFR2)./max(refFR2),1,1,1)
% refFR=nanmean((FPF),2);refFR=refFR-min(refFR);
% refFR=fastsmooth((refFR-min(refFR))./max(refFR),1,1,1)
% 
% figure
% plot(ax,refFR,'r','Linewidth',1.5);hold on,
% plot(ax,refFR2,'b-','Linewidth',1.5)
% 
% freq2.freq=frs;
% th= squeeze(nanmean(nanmean(PF(end-5:end,:,:),1),2)) <0.91 ;
% ax=-wd:wd;fsel1= freq2.freq>5 & freq2.freq <12; fsel2= freq2.freq>60 & freq2.freq <90; 
% figure('COlor','w','Position',[ 300 300 400 pheight])
% subplot(1,2,1),imagesc(ax,freq2.freq,smoothn(nanmean(PF(:,:,th),3),.5));axis xy;colormap(jet)
% set(gca,'Clim',[ -0. 0.18]);ylim([1 40])
% hold on,xlabel('place field center (cm) ');ylabel('Frequency Hz')
% %plot(ax,zscore(nanmean(FPF,2))*50+50,'m','Linewidth',1.5);plot(ax,zscore(nanmean(FPF2,2))*50+50,'w--')%plot(ax,zscore(nanmean(FPF,2))*50+49.2,'w-');
% title('complex spikes')
% subplot(1,2,2),imagesc(ax,freq2.freq,smoothn(nanmean(PF2(:,:,th),3),.5));axis xy;colormap(jet)
% hold on;%plot(ax,zscore(nanmean(FPF2,2))*50+50,'m','Linewidth',1.5);plot(ax,zscore(nanmean(FPF,2))*50+50,'w--')%plot(ax,zscore(nanmean(FPF2,2))*50+49.2,'w-')
% set(gca,'Clim',[ -0. 0.18]);ylim([1 40])
% xlabel('place field center (cm)');ylabel('Frequency Hz');title('single spikes')
% 
% 
