function [resp_t, resp_dur, ppr] = ppanalysis_v2(psth, binsize, light_dur, Hz, all_light)
% difference from v1 is that you input the psth from which you want to calculate paired-pulse ratios, instead of creating it from inside this function
% psth = 1 x #bins vector of average spike counts from each bin. calculates
%   responses from WHOLE psth, so psth should only include relevant time
%   points
% binsize = binsize of psth (in sec)
% light_dur = duration of light period to analyze (in sec)
% Hz = frequency condition to analyze
% all_light = 1 x #trials vector where the entry indicates light condition

thresh = mean(psth)+3*std(psth);       % but now the threshold isn't from no light condition...REVISIT (and revert to v1 for now)
hz_times = 0:1/Hz:light_dur-(1/Hz);     % [.5:.1:1.4]
times = 0:binsize:light_dur;
for t = 1:length(hz_times)  % count spikes in 50ms after each pulse
    resps{t} = find(psth(times>=hz_times(t) & times<hz_times(t)+.05)>=thresh)+find(times>=hz_times(t),1,'first')-1;   % hardcoded for 10hz cond - assumes it's 3rd lightcond!
end

if sum(cellfun(@isempty,resps(1:2)))==2 && sum(cellfun(@isempty,resps))>=Hz/2 % if neither of first two pulses yielded a significant response AND half or less of pulses didn't elicit a significant response
    resp_t = nan(1,length(hz_times));
    ppr = nan(1,length(hz_times)-1);
    resp_dur = nan(1,length(hz_times));
else
    if isempty(resps{1})
           resp_vals(1) = max(psth(times>=hz_times(1) & times<hz_times(1)+.05));    % use maximum because, if anything, this should overestimate ppr
        if resp_vals(1)==0; resp_vals(1) = thresh; end      % if there were NO spikes in window after first pulse, set to threshold so that you don't end up dividing by 0
        resp_t(1) = nan;
        resp_dur(1) = 0;    % no bins passed significance threshold after first pulse
    else
        resp_t(1) = times(resps{1}(1))-hz_times(1);
        resp_vals(1) = mean(psth(resps{1}(1):resps{1}(end)));        % take mean peak response (could also be sum or max...)
        resp_dur(1) = length(resps{1});
    end
    resps(cellfun(@isempty,resps)) = {nan};
    for t=2:Hz
        if ~isnan(resps{t})
             resp_vals(t) = mean(psth(resps{t}(1):resps{t}(end)));    % take mean peak response (could also be sum or max...)
            resp_dur(t) = length(resps{t}); % how many bins passed the significance threshold
            resp_t(t) = times(resps{t}(1))-hz_times(t);
        else
            resp_vals(t) = mean(psth(times>=hz_times(t) & times<hz_times(t)+.05)); % mean response during 50ms time window (even if it didn't cross threshold)
            resp_dur(t) = 0;    % no bins passed significance threshold after pulse t
            resp_t(t) = nan;
        end
        if t>1
             ppr(t-1) = resp_vals(t)/resp_vals(1);
         end
    end
end



return