function [ PLV_output,PLV_phase]= spike_field_ppc_adjMatchTET(wavD,spikes,spikes2,CH,sp_shift)
if size(wavD,4)>1
    [tr chs frs t]= size(wavD);
    
elseif size(wavD,3)>1
    [tr  frs t]= size(wavD);
else
     [ frs t]= size(wavD); tr=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat=[]; mat2=[];
for ind=1
    s=spikes;s=s-sp_shift; s(s<=10)=[]; s(s> length(wavD))=[]; % spike times, shifted to avoid spike contamination issue
     s2=spikes2;s2=s2-sp_shift; s2(s2<=10)=[]; s2(s2> length(wavD))=[];
    mat=[mat,    wavD(:,s)];   mat2=[mat2,   wavD(:,s2)];
end

NT= sum(~isnan(mat(1,:))); % number of observations
NT2= sum(~isnan(mat2(1,:))); % number of observations


if NT >10
    
    Z=abs(nanmean(mat,2));  % MVL
    
   PLV_phase=angle(nanmean(mat,2));
    for ix=1:length(Z);
        ZZ1=squeeze(mat(ix,~isnan(mat(ix,:))));
          ZZ2=squeeze(mat2(ix,~isnan(mat2(ix,:))));
        NT= sum(~isnan(mat(ix,:)));
           NT2= sum(~isnan(mat2(ix,:)));
        [ x1 x2]=min([ NT NT2])
if x2==1
   h1= randperm(NT2); ZZ2=ZZ2(h1(1:NT));
else
  h1= randperm(NT); ZZ1=ZZ1(h1(1:NT2)); 
end
NT=x1;
        T=   abs(nanmean(ZZ1)).^2;
        PLV_output(ix,1)=  abs(nanmean(ZZ1));%(((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
         T=   abs(nanmean(ZZ2)).^2;
        PLV_output(ix,2)=  abs(nanmean(ZZ2));%(((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
        
        if NT<10
            PLV_output(ix)=NaN; %ones(size(wavD,3),1).*
        end
    end
    
else
    PLV_output=ones(size(wavD,3),1).*NaN;
     PLV_phase=ones(size(wavD,3),1).*NaN;
end