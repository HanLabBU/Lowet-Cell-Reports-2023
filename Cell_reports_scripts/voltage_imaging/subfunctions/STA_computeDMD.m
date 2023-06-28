function [ STA]= STA_computeDMD(SIG,spikes,lwin,rwin)
STA=[];
tr=length(SIG);n=0;
for ind=1:tr
    trace=SIG{ind}(:,:);
   trigs=find( spikes{ind});
   NT=length(trigs);
   for id=1:NT
       if trigs(id)-lwin>0 & trigs(id)+rwin< length(trace)
           n=n+1;
   STA(:,:,n)=      (trace(:,trigs(id)-lwin: trigs(id)+rwin));
       
       end
   
   end
    
    
end

STA=nanmean(STA,3);
STA=bsxfun(@minus,STA, nanmean(STA(:,1:40),2));
if n<5
   STA=STA.*NaN; 
end
