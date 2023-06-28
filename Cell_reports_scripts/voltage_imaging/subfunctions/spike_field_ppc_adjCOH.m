function [ PLV_output,PLV_phase]= spike_field_ppc_adjCOH(wavD,spikes,CH,sp_shift,wavA)
if size(wavD,4)>1
    [tr chs frs t]= size(wavD);
    
elseif size(wavD,3)>1
    [tr  frs t]= size(wavD);
else
     [ frs t]= size(wavD); tr=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat=[]; matA=[];
for ind=1:length(spikes)
    if size(wavD,4)>1
        M=exp(1i.*((squeeze(wavD(ind,CH,:,:)))));
        M1=(squeeze(wavA(ind,CH,:,:)));
elseif size(wavD,3)>1
        M=exp(1i.*((squeeze(wavD(ind,:,:)))));
         M1=(squeeze(wavA(ind,:,:)));
    else
        M=exp(1i.*((squeeze(wavD(:,:)))));
                 M1=(squeeze(wavA(:,:)));
    end
    
    s=find(spikes{ind});s=s-sp_shift; s(s<=10)=[]; % spike times, shifted to avoid spike contamination issue
    
    mat=[mat,    M(:,s)];   matA=[matA,    M1(:,s)];
end

NT= sum(~isnan(mat(1,:))); % number of observations
if NT >10
    
   % Zold=abs(nanmean(mat,2));  % MVL
  matA=  bsxfun(@rdivide, matA, nansum(matA,2));
        Z=abs(nansum(matA.*mat,2)); 
   PLV_phase=angle(nanmean(mat,2));
    for ix=1:length(Z);
        NT= sum(~isnan(mat(ix,:)));
        
        T=   Z(ix).^2;
        PLV_output(ix)=  (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
        
        if NT<10
            PLV_output(ix)=NaN; %ones(size(wavD,3),1).*
        end
    end
    
else
    PLV_output=ones(size(wavD,3),1).*NaN;
     PLV_phase=ones(size(wavD,3),1).*NaN;
end