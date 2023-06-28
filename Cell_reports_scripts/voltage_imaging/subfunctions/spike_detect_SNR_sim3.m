function  result=spike_detect_SNR_v2(traces,up_threshold, down_threshold, sm)

addpath(genpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\Scripts\'))
% result=spike_detect_SNR(traces)
% 
 
FS=1000; % sampling frequency
%addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\Hua-an Tseng\Code\other\eegfilt\'))
plot_yes=0;
  if nargin<3 || isempty(down_threshold) 
        down_threshold = 4;
    end
    
    if nargin<2 || isempty(up_threshold) 
        up_threshold = 4.5;
    end

  event_parameter.moving_window = 501;
    event_parameter.pre_peak_data_point = 1;
    event_parameter.post_peak_data_point = 1;
    event_parameter.event_moving_window = 501; % data points
       
    event_parameter.noise_threshold = 7;
    event_parameter.noise_pre_extension = 0; % data points
    event_parameter.noise_post_extension = 3; % data points
    event_parameter.noise_extension = 3; % data points
    event_parameter.noise_moving_window = 11; % data points
    
    event_parameter.snr_threshold = 0;
%     event_parameter.snr_window = 100;
    event_parameter.refine_threshold = 4; % standard deviation
    event_parameter.down_threshold = down_threshold;
    event_parameter.up_threshold = up_threshold;
event_parameter.moving_window=501;
result=[];
for neuron=1:size(traces,2)

 current_traceOrig = (traces(:,neuron));
 
 %current_traceOrig(1:10)= mean(current_traceOrig(11:end));
 
 current_trace= current_traceOrig-fastsmooth(current_traceOrig,1000,1,1);
              event.idx=[];
            event.amplitude=[];
            event.snr=[];
            
         %   f_trace = eegfilt(current_trace',FS,4,0)';
            f_trace= current_trace-fastsmooth(current_trace,sm,1,1);
 current_trace= f_trace;
  u_f_trace = get_upper_trace(f_trace,event_parameter.moving_window);
            
            d_u_f_trace = diff(u_f_trace);
            d_u_f_trace = [0;d_u_f_trace];

%%%
 l_f_trace = get_lower_trace(f_trace,event_parameter.moving_window);
            d_l_f_trace = diff(l_f_trace);
            d_l_f_trace = [0;d_l_f_trace];
            
 event_parameter.up_threshold_value = event_parameter.up_threshold*nanstd(d_l_f_trace);
            event_parameter.down_threshold_value = event_parameter.down_threshold*nanstd(d_l_f_trace);
            
 trace_val2=f_trace;
   up_idx_list = find( trace_val2>(event_parameter.up_threshold_value));
trace_wp=    current_traceOrig;
rem_sel=[up_idx_list; up_idx_list-1 ;up_idx_list+1 ;up_idx_list+2 ; up_idx_list-2  ;up_idx_list+3 ; up_idx_list-3  ];
   rem_sel(rem_sel<0)=[];
   rem_sel(rem_sel>length(trace_wp))=[];
   try
   trace_wp( rem_sel) =NaN;
   end
   trace_wp=trace_wp -nanfastsmooth(trace_wp,sm,1,1);
   trace_wp(isnan(trace_wp))=NaN;
          l_f_trace = get_lower_trace( trace_wp,event_parameter.moving_window);
   
            l_f_trace2=  l_f_trace;
            l_f_trace2(isnan(l_f_trace2))=[];
            
          
            l_f_trace2=[l_f_trace2;l_f_trace2.*-1];
            current_trace_noise = 1*nanstd(l_f_trace2); 

%%%%%%%%%%%

            event.event_parameter.up_threshold_value = event_parameter.up_threshold_value;           
            event.event_parameter.down_threshold_value = event_parameter.down_threshold_value;
%%%%%%%%%%%%%%%%%%
  noise_idx_list = find_noise_idx(f_trace,event_parameter.moving_window,event_parameter.noise_moving_window,event_parameter.noise_threshold,event_parameter.noise_extension);
%            result.roi(roi_idx).noise_idx = noise_idx_list;
            
            denoise_trace = current_trace;
            for idx=1:numel(noise_idx_list)
                current_noise_idx = noise_idx_list(idx);
                if current_noise_idx-event_parameter.noise_pre_extension>0
                    denoise_trace(current_noise_idx-event_parameter.noise_pre_extension:min(current_noise_idx+event_parameter.noise_post_extension,length(current_trace))) = nan;
                end
            end
%%%%%%%%%%
     pre_d_trace = d_u_f_trace;
            if event_parameter.pre_peak_data_point>0
                for idx=1:event_parameter.pre_peak_data_point
                    shifted_d_trace = [zeros(idx,1);d_u_f_trace(1:end-idx)];
                    shifted_d_trace(shifted_d_trace<0) = 0;
                    pre_d_trace = pre_d_trace+shifted_d_trace;
                end
            end
            
            post_d_trace = d_u_f_trace;
            if event_parameter.post_peak_data_point>0
                for idx=1:event_parameter.post_peak_data_point
                    shifted_d_trace = [d_u_f_trace(idx+1:end);zeros(idx,1)];
                    shifted_d_trace(shifted_d_trace>0) = 0;
                    post_d_trace = post_d_trace+shifted_d_trace;
                end
            end

%%%%%%%%%%%%%
 trace_val=pre_d_trace;
   up_idx_list = find(pre_d_trace>(nanmean(d_u_f_trace)+event_parameter.up_threshold_value));
   
   %%%%%%%%%%%%%%%%%%%%%%
 trace_wp= f_trace;
 rem_sel=[up_idx_list; up_idx_list-1 ;up_idx_list+1 ;up_idx_list+2 ; up_idx_list-2  ;up_idx_list+3 ; up_idx_list-3  ];
   rem_sel(rem_sel<0)=[];
   rem_sel(rem_sel>length(trace_wp))=[];
   try
   trace_wp( rem_sel) =[];
   end
  current_trace_noise = std(trace_wp);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% up_idx_list= find((current_trace)>(4.5*current_trace_noise) )  ;
 for up_idx=up_idx_list'
             %   if up_idx>2 & (up_idx+1)<=numel(d_u_f_trace) & d_u_f_trace(up_idx)>0 %& d_u_f_trace(up_idx+1)<0 %& post_d_trace(up_idx+1)<(nanmean(d_u_f_trace)-event_parameter.down_threshold_value)% & ~isnan(denoise_trace(up_idx))
               if up_idx>2 & (up_idx+1)<=numel(d_u_f_trace) &  d_u_f_trace(up_idx)>0% &  ~isnan(denoise_trace(up_idx)) %&  post_d_trace(up_idx+4)<(nanmean(d_u_f_trace)-event_parameter.down_threshold_value)
           
              peak_intensity = current_trace(up_idx);
                    pre_peak_intensity_1 = current_trace(up_idx-1);
                    pre_peak_intensity_2 = current_trace(up_idx-2);
                     post_peak_intensity_1 = current_trace(up_idx+1);
                  peak_V=trace_val(up_idx);
                valnear= find(abs( up_idx- up_idx_list)<=2 &   abs( up_idx- up_idx_list)>0);
  if isempty(valnear); vthres=0; else;  vthres=max(trace_val(up_idx_list(valnear)));end
           current_signal_intensity = max(peak_intensity-pre_peak_intensity_1,peak_intensity-pre_peak_intensity_2);
                    current_snr = current_signal_intensity/current_trace_noise;
            if peak_V > vthres;% >  post_peak_intensity_1  & peak_intensity > pre_peak_intensity_1 
             if current_snr>=event_parameter.snr_threshold
                        event.idx = cat(1,event.idx,up_idx);
                        event.amplitude = cat(1,event.amplitude,current_signal_intensity);
                        event.snr = cat(1,event.snr,current_snr);
                    end;end
                    
                end
            end





            event.roaster = zeros(size(current_trace));
            event.roaster(event.idx) = 1;   
            event.roaster2 = zeros(size(current_trace));
            event.roaster2(up_idx_list) = 1;  
            event.trace_noise = current_trace_noise;
            event.snr_threshold = event_parameter.snr_threshold;

            
    tracews=current_traceOrig;current_trace=current_traceOrig;
   for sind=1:length(event.idx)
       if event.idx(sind) > 3  & event.idx(sind)< length(current_trace)-3
     tracews( event.idx(sind)-2:event.idx(sind)+2)= mean(current_trace([event.idx(sind)-3 event.idx(sind)+3]));
     %  tracews( event.idx(sind)-1)= mean(tracews([event.idx(sind)-2 event.idx(sind)-2]));
      %  tracews(event.idx(sind)+1)= mean(tracews([event.idx(sind)+2 event.idx(sind)+2]));
       %  tracews( event.idx(sind))= mean(tracews([event.idx(sind)-1  event.idx(sind) event.idx(sind)+1]));
       
       end
   end           
            
            
result.orig_trace(neuron,:) =current_traceOrig;
result.denoise_trace(neuron,:) =denoise_trace;
result.hp_trace(neuron,:)=f_trace;
result.roaster(neuron,:) = event.roaster;
result.roaster2(neuron,:) = event.roaster2;
result.trace_ws(neuron,:) =    tracews;

result.spike_snr{neuron,1} =  event.snr ;

result.spike_amplitude{neuron,1} = event.amplitude' ;
result.spike_idx{neuron,1} =  event.idx  ;
result.trace_noise(neuron,1)= event.trace_noise;
end


if plot_yes==1
rast=result.roaster;
rast(rast==0)=NaN;
% figure,
% for ind=1:size(rast,1)
% plot(rast(ind,:) + ind ,'.k'); hold on,
% end

figure('COlor','w')
subplot(1,3,1:2)
for ind=1:size(result.orig_trace,1)
plot(((result.orig_trace(ind,:)./result.trace_noise(ind)))./15+ ind,'k'); hold on,
hold on,plot((rast(ind,:)+result.trace_noise(ind))./39 + ind ,'.r','Markersize',10); hold on,
end;axis tight;xlabel('Time'); ylabel('neuron')
subplot(1,3,3)
clear SNR_val
for ind=1:length(result.spike_snr)
SNR_val(ind)= mean(result.spike_snr{ind})
end
bar(SNR_val');xlabel('neuron'); ylabel('SNR')




end


function lower_trace = get_lower_trace(current_trace,trace_moving_window)

    m_trace = movmean(current_trace,trace_moving_window);
    lower_trace = current_trace;
    % replace the part below moving average with moving average
    idx = find(lower_trace>m_trace);
    lower_trace(idx)=NaN;%m_trace(idx);
end


function upper_trace = get_upper_trace(current_trace,trace_moving_window)

    m_trace = movmean(current_trace,trace_moving_window);
    upper_trace = current_trace;
    % replace the part below moving average with moving average
    idx = find(upper_trace<m_trace);
    upper_trace(idx)=m_trace(idx);
end

function noise_idx_list = find_noise_idx(current_trace,trace_moving_window,noise_moving_window,noise_threshold,noise_extension)

    m_trace = movmean(current_trace,trace_moving_window);
    lower_current_trace = current_trace;
    % replace the part above moving average with moving average
    idx = find(lower_current_trace>m_trace);
    lower_current_trace(idx)=m_trace(idx);

    movstd_lower_current_trace = movstd(lower_current_trace,noise_moving_window);
    noise_idx_list = find(isoutlier(movstd_lower_current_trace,'gesd')==1);
    %noise_idx_list = find(movstd_lower_current_trace>(mean(movstd_lower_current_trace)+noise_threshold*std(movstd_lower_current_trace)));

    % connect noise index
    noise_idx_list = sort(noise_idx_list);
    d_noise_idx_list = diff(noise_idx_list);
    noise_extension_idx = find(d_noise_idx_list>1 & d_noise_idx_list<noise_extension);
    if ~isempty(noise_extension_idx)
        for idx=1:numel(noise_extension_idx)
            current_idx = noise_extension_idx(idx);
            noise_idx_list = cat(1,noise_idx_list,[noise_idx_list(current_idx):noise_idx_list(current_idx+1)]');
        end
    end

    noise_idx_list = unique(noise_idx_list);

end

end