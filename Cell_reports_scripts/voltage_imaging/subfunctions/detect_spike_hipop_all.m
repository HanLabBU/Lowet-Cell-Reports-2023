
 addpath(genpath('Z:\EricLowet\'))
 % addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))
addpath(genpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\main_analysis\'))
  path1='C:\hippo_opto\'
 
%   rec=31;
% ses=['rec' num2str(rec)]; 
% 
% mouse='611135';%613199'; 
% mouse='613198';%613199'; 
% mouse='611112';%613199'; 
% %mouse='606716';%613199'; 
% %mouse='611185';%613199'; 
% cell=2
% prepath='\\engnas.bu.edu\research\eng_research_handata\';%EricLowet\hippo_opto_main\RESULTS\place
% %%%%%%%%%%%%%%%%%%%
% %Cpath= [prepath 'EricLowet\hippo_opto_main\RESULTS\place' '\'  mouse '_rec'  num2str(rec) '_' 'fov' num2str(cell) '.mat'   ];
% Cpath= [prepath 'EricLowet\hippo_opto_main\RESULTS\DC' '\'  mouse '_rec'  num2str(rec) '_' 'fov' num2str(cell) '.mat'   ];

 %addpath(genpath('Z:\EricLowet\'))
 
 mainpath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\hippo_opto_main\RESULTS\';
 cd(mainpath);
 
 foldnames = {'DC'; 'DC_multi_level_trial';'place'; 'RES_L'; 'DC_L'; 'other'; 'RES_C'; 'focal_wide',;'ramp'; 'DC_DMD';'hip_osc'} ;
 
 
 
 for foldn= [7]%1:length(foldnames);
 cd(mainpath)
 cd(foldnames{foldn})
 ses=dir('*.mat');

 
 
for ih= 1: length(ses)
    ih
Cpath= ses(ih).name;
load(Cpath)

if 1%~isfield(result, 'resultS')
clear resultS
mk=0;
for ind=unique(result.trial_vec)
 % mk=mk+1;
 if    ~isempty(result.traces(result.trial_vec==ind,:))
%resultS{ind}=spike_detect_SNR_hipop(result.traces(result.trial_vec==ind,:));
resultS{ind}=spike_detect_SNR_sim3(result.traces(result.trial_vec==ind,:),4,4,7);

 else
  resultS{ind}= NaN;   
 end
end
%%%%%%%%%%

  result.resultS=resultS;
     save(Cpath,'result')


end

 end
 end

 
 
%  
% ses=dir('*.mat')
% 
% for ih= 1: length(ses)
%     ih
% Cpath= ses(ih).name;
% load(Cpath)
% 
% %%%%%%%%%%%%
% clear resultS
% for ind=unique(result.trial_vec)
%     ind
%  if    ~isempty(result.traces(result.trial_vec==ind,:))
% resultS{ind}=spike_detect_SNR_hipop(result.traces(result.trial_vec==ind,:));
%  else
%   resultS{ind}= NaN;   
%  end
% end
% %%%%%%%%%%
% 
%   result.resultS=resultS;
%      save(Cpath,'result')
% 
% end