function [resp_t, resp_dur, ppr, ppr2] = ppanalysis(prestim, stim_time, light_start, light_dur, Hz, rast, all_light)
% prestim = amount of time (in seconds) before visual stimulus onset
% stim_time = length of visual stimulus (in sec)
% light_start = time after visual stimulus when light starts (in sec)
% light_dur = duration of light period to analyze (in sec)
% Hz = frequency condition to analyze
% rast = raster for unit (num_trials x trial length (in ms))

% 10/29/19 - changed from means to sums

binsize = .005; % 5 ms bins
edges = 0:binsize:stim_time+prestim;
times = -prestim:binsize:stim_time;
[~,psth] = make_psth_v2(.005,edges,1:length(all_light),rast,all_light);     % use # of spikes rather than firing rate
light_conds = unique(all_light);
which_cond = find(light_conds==Hz);
light_bins = (prestim+light_start)/binsize+1:(prestim+light_start+light_dur)/binsize;   % bins from which to extract threshold (light period)

% thresh = mean(psth(which_cond,:))+2*std(psth(which_cond,:));
thresh = mean(psth(light_conds==0,light_bins))+3*std(psth(light_conds==0,light_bins));       % set threshold as 3stds above mean, during light period in NO LIGHT condition
hz_times = light_start:1/Hz:light_start+light_dur-(1/Hz);     % [.5:.1:1.4]
for t = 1:length(hz_times)  % count spikes in 50ms after each pulse
    resps{t} = find(psth(which_cond,times>=hz_times(t) & times<hz_times(t)+.05)>=thresh)+find(times>=hz_times(t),1,'first')-1;  
end

% if sum(cellfun(@isempty,resps))>=5 && isempty(resps{1})  % if there was no significant response after at least 50% of pulses (including first one), dont count!
if sum(cellfun(@isempty,resps(1:2)))==2 && sum(cellfun(@isempty,resps))>=Hz/2 % if neither of first two pulses yielded a significant response AND half or less of pulses didn't elicit a significant response
    resp_t = nan(1,length(hz_times));
    ppr = nan(1,length(hz_times)-1);
    ppr2 = nan(1,length(hz_times)-1);
    resp_dur = nan(1,length(hz_times));
else
    if isempty(resps{1})
    %     resps{1} = find(times==.510);      % hardcoded for light started .5s after visstim and pulse duration=10ms
%         resp_vals(1) = thresh;  % set as threshold (otherwise pprs will get really big, but threshold is the max possible value (otherwise resps{1} would not be empty)
        resp_vals(1) = sum(psth(which_cond,times>=hz_times(1) & times<hz_times(1)+.05));
%         resp_vals(1) = max(psth(which_cond,times>=hz_times(1) & times<hz_times(1)+.05));    % use maximum because, if anything, this should overestimate ppr
        if resp_vals(1)==0; resp_vals(1) = thresh; end      % if there were NO spikes in window after first pulse, set to threshold so that you don't end up dividing by 0
        resp_t(1) = nan;
        resp_dur(1) = 0;    % no bins passed significance threshold after first pulse
    else
        resp_t(1) = times(resps{1}(1))-hz_times(1);
%         resp_vals(1) = mean(psth(which_cond,resps{1}(1):resps{1}(end)));        % take mean peak response (could also be sum or max...)
        resp_vals(1) = sum(psth(which_cond,times>=hz_times(1) & times<hz_times(1)+.05));
        resp_dur(1) = length(resps{1});
    end
    resps(cellfun(@isempty,resps)) = {nan};
    % resps = find(psth(3,find(times>=.5,1,'first'):end)>=thresh)+find(times>=.5,1,'first')-1;
    for t=2:Hz
        if ~isnan(resps{t})
%              resp_vals(t) = mean(psth(which_cond,resps{t}(1):resps{t}(end)));    % take mean peak response (could also be sum or max...)
             resp_vals(t) = sum(psth(which_cond,times>=hz_times(t) & times<hz_times(t)+.05));    % take mean peak response (could also be sum or max...)
    %          if t>1
    %              ppr(t-1) = resp_vals(t)/resp_vals(1);
    %          end
            resp_dur(t) = length(resps{t}); % how many bins passed the significance threshold
            resp_t(t) = times(resps{t}(1))-hz_times(t);
        else
    %         ppr(t-1) = nan;
%             resp_vals(t) = thresh;      % once again, highest possible response
            resp_vals(t) = sum(psth(which_cond,times>=hz_times(t) & times<hz_times(t)+.05)); % sum response during 50ms time window (even if it didn't cross threshold)
            resp_dur(t) = 0;    % no bins passed significance threshold after pulse t
            resp_t(t) = nan;
        end
        if t>1
             ppr(t-1) = resp_vals(t)/resp_vals(1);
             ppr2(t-1) = resp_vals(t)/resp_vals(2); % ratio compared to SECOND pulse
        end
    end
end
% if ~isempty(resps{1}) || ~isempty(resps{2})
% %     reps = find(diff(resps)<10);
% %     resp_vals = zeros(1,length(resps)-length(reps));
% %     goods = ~ismember(1:length(resps),reps);
% %     ct = 1;
% %     for i=1:length(resp_vals)
% %         if ismember(ct,reps)
% %             resp_vals(i) = max(psth(3,resps(ct):resps(find(goods(ct+1:end),1,'first')+ct)));
% %             ct = ct+(find(goods(ct+1:end),1,'first'))+1;
% %         else
% %             resp_vals(i) = psth(3,resps(ct));
% %             ct = ct+1;
% %         end
% %     end
%     for i = 1:2
%         if ~isempty(resps{i})
%             resp_vals(i) = max(psth(3,resps{i}(1):resps{i}(end)));
%         else
%             resp_vals(i) = thresh;
%         end
%     end
%     all_resps = cell2mat(resps);
%     resp_t = times(all_resps(1));
% %     if sum(diff(resps))>=160 && sum(diff(resps))<190      % TEMP - at least two 20ms intervals between peaks (b/c 10hz) - find better way?
%     if length(cell2mat(resps))>=5       % if it followed at least half of the pulses
%         ppr = resp_vals(2)/resp_vals(1);     % "paired-pulse" ratio
%     else
%         ppr=nan;
%     end
% else
%     ppr = nan;
%     resp_t = nan;
% end
% % figure; plot(zscore(psth(3,:)'));

return