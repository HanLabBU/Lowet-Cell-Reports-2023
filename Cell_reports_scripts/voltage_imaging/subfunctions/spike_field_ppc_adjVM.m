function [ PLV_output,PLV_phase]= spike_field_ppc_adj(wavD,spikes,CH,sp_shift)
if size(wavD,4)>1
    [tr chs frs t]= size(wavD);
    
elseif size(wavD,3)>1
    [tr  frs t]= size(wavD);
else
     [ frs t]= size(wavD); tr=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat=[];
for ind=1:length(spikes)%tr
    if size(wavD,4)>1
        M=(((squeeze(wavD(ind,1,:,:)))));
        M2=(((squeeze(wavD(ind,2,:,:)))));
elseif size(wavD,3)>1
        M=exp(1i.*((squeeze(wavD(ind,:,:)))));
     %   M2=exp(1i.*((squeeze(wavD(ind,3,:,:)))));

    else
        M=exp(1i.*((squeeze(wavD(:,:)))));
    end
    
    s=find(spikes{ind});s=s-sp_shift; s(s<=20)=[]; % spike times, shifted to avoid spike contamination issue
     s(s>=length(M)-20)=[];
    mat=[mat,    exp(1i.*(circ_dist(M2(:,s),M(:,s))))];
end

NT= sum(~isnan(mat(1,:))); % number of observations
if NT >10
    
    Z=abs(nanmean(mat,2));  % MVL
    
   PLV_phase=angle(nanmean(mat,2));
    for ix=1:length(Z);
        NT= sum(~isnan(mat(ix,:)));
        
        T=   Z(ix).^2;
        PLV_output(ix)= (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
        
        if NT<10
            PLV_output(ix)=NaN; %ones(size(wavD,3),1).*
        end
    end
    
else
    PLV_output=ones(size(wavD,3),1).*NaN;
     PLV_phase=ones(size(wavD,3),1).*NaN;
end