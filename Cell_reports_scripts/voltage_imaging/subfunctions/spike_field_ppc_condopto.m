function [ PLV_output]= spike_field_ppc_condopto(wavD,spikes,CH,trs,styp,FS,optostart)
if size(wavD,4)>1
[tr chs frs t]= size(wavD);
else
[tr  frs t]= size(wavD);
end
mat=[];
for ind=trs
    if size(wavD,4)>1
         M=exp(1i.*((squeeze(wavD(ind,CH,:,:))))); 
    else
     M=exp(1i.*((squeeze(wavD(ind,:,:))))); end
 if styp==1
      s=find(spikes{ind});s=s-7; 
      s(s<(optostart))=[];s(s>optostart+1*FS)=[]; s(s<10)=[]; 
 else
     s=find(spikes{ind});s=s-7; s(s>(optostart-50) & s<(optostart+1.6*FS))=[];%s(s>1.5*FS)=[]; 
     s(s<10)=[];
 end
      
  mat=[mat,    M(:,s)];
end
 NT1= sum(~isnan(mat(10,:)));
if NT1>10
  Z=abs(nanmean(mat,2));  % MVL mean vector length
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
    
    
end
