function [spikTyp6]=spike_burst_n(dataM,spikTyp,FS,num)

spikTyp6=spikTyp;
for id=1:length(spikTyp)
    z= find(spikTyp{id}>-7);
    sptim2=dataM.spikes{id}( z);
    bursts_onsets = find(diff(sptim2)>(30.*(FS./1000)))+1;
    ff=find(spikTyp{id}==1);
    for x=1:length(bursts_onsets)
        if any(bursts_onsets(x)==ff)
            found_sp=find((sptim2>sptim2(bursts_onsets(x)))&  sptim2<= sptim2(bursts_onsets(x)+1)+(30.*(FS./1000)) & (spikTyp{id}==1)');
            if num <3
                if length(found_sp)==num
                    spikTyp6{id}(z(bursts_onsets(x)))=2;
                    for j=1:length(found_sp)
                        spikTyp6{id}(z(found_sp(j)))=3;%+j;
                    end;end
            else
                if length(found_sp)>=num
                    spikTyp6{id}(z(bursts_onsets(x)))=2;
                    for j=1:length(found_sp)
                        spikTyp6{id}(z(found_sp(j)))=3;%+j;
                    end;end
            end
            
            
        end
    end
end