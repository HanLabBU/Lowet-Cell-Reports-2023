function  [dataM]=burst_spike_select(dataM,spikTypA,spikTypB,nA,nB) 


for ik=1:size(dataM.roaster,2)
            ROAST1=zeros(1,length(dataM.roaster{ik}));try
                ROAST1(dataM.spikes{ik}(spikTypA{ik}==nA & ( spikTypB{ik}==nB  )   ))=1;
                dataM.roaster{ik}=ROAST1;end
end
        
