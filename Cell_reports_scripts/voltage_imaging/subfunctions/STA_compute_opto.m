function [ STA]= STA_compute_opto(SIG,spikes,CH,lwin,rwin)

tr=length(SIG);n=0;
for ind=1:tr
    trace=SIG{ind}(CH,:);
   trigs=find(diff(spikes{ind})>0.1);%find( spikes{ind});
   NT=length(trigs);
   for id=1:NT
       if trigs(id)-lwin>0 & trigs(id)+rwin< length(trace)
           n=n+1;
   STA(:,n)=      (trace(trigs(id)-lwin: trigs(id)+rwin));
       
       end
   
   end
    
    
end

%STA=nanmean(STA,2);
%STA=STA-nanmean(STA(1:50));
if n<10
   STA=STA.*NaN; 
end
