function [burstrate,burstnum,burstresprate,avtonicrate,avburstrate] = calc_bursting(raster,timeinds,plotting)
% MAK 3/2/20 - 
% raster = num_trials x time raster of spiketimes (each bin=1ms)
% timeinds = indices for timeframe in which to analyze bursting (e.g., 1001:2000ms)
% plotting = 1 if you want ISI plots, else 0

visRast = raster(:,timeinds(1)-100:timeinds(end)+4);  % -100 AND +4 b/c otherwise first and last burst spikes in trials would be excluded    
[rows, cols] = find(visRast);
[trialn,ord] = sort(rows);   %sort by trials
sptimes = cols(ord);    % spike times (relative to timeinds(1) by trial)
isis = [sptimes(1); diff(sptimes)]; % interspike intervals by trial (use ms index of first spike for it's isi, since we don't actually know when prior spike occurred
% an isi for every spike - e.g., isis(5) indicates how many ms PRECEDED the 5th spike
trial_trans = logical([0; diff(trialn)]);  % indexes of spikes that start a new trial (previously isis<=0, but misses some trial transitions!
isis(trial_trans) = sptimes(trial_trans); % again, use ms index of first spike in a trial for it's isi (isis<=0 mark transitions between trials)
burst_mid = find(isis(1:end-1)<=4 & sptimes(1:end-1)<=size(visRast,2)-4); % spikes in middle and end of bursts (preceeded by <=4ms ISIs, and NOT in 4ms pad at end of trial)
burst_st = find(isis(1:end-1)>=100 & isis(2:end)<=4 & sptimes(1:end-1)<=size(visRast,2)-4 & trialn(1:end-1)==trialn(2:end));   % spikes at start of bursts (followed by <=4ms ISI but preceded by >=100ms) - NOT in 4ms pad at end of trial and can't be last spike in trial
bursts = union(burst_mid,burst_st); % all burst indices
while sum(~ismember(burst_mid-1,bursts))    % fixed 3/2/20 - previously only did this once, which still left in very small number incorrect burst spikes
    premid = burst_mid(~ismember(burst_mid-1,bursts));  % burst_mid spikes whose preceding spikes were NOT burst spikes
    bursts(ismember(bursts,premid)) = [];
    burst_mid(ismember(burst_mid,premid)) = [];
end
burstrate = length(bursts)/(sum(sum(visRast(:,101:end-4)))); % how many of all spikes in lightstim period were part of bursts (across all trials)
    
% calculate number of spikes in bursts
burst_st(burst_st>length(isis)-4) = [];  % if start of burst is too close to end, won't be able to count # spikes in that burst
burstmat = cumsum([burst_st ones(length(burst_st),25)],2);
burstmat(burstmat>length(isis)) = length(isis);
burstlength = nan(1,length(burst_st));
for b=1:size(burstmat,1)    % for each start of a burst
    if ~isempty(find(isis(burstmat(b,2:end))>4,1,'first'))
        burstlength(b) = find(isis(burstmat(b,2:end))>4,1,'first');   % how many spikes before >4ms ISI
    end
end
burstnum = nanmean(burstlength);    % save unit's mean burstlength
burstresprate = length(burst_st)/(sum(sum(visRast(:,101:end-4)))-length(burst_mid));    % how many burst RESPONSES (as opposed to indiv burst spikes) out of all responses ("tonic responses" = indiv spikes)

% calculate trial-averaged "firing rates" for burst vs. tonic spikes
tonics = trialn(~ismember(1:length(isis),bursts));
for n=1:size(raster,1)
    tonic_tri(n) = sum(tonics==n);
    burst_tri(n) = sum(trialn(bursts)==n);
end
avtonicrate = mean(tonic_tri); % average rate of tonic spikes/s across trials
avburstrate = mean(burst_tri); % average rate of burst spikes/s across trials

%make ISI plot: (only with spikes during analysis window, excluding 1000ms and 4ms pads
if plotting
    preisis = isis(1:end-1);
    postisis = isis(2:end);
    postisis(trial_trans(2:end)&preisis>=100&postisis<=4)=nan;  % exclude post-isis for last spikes in each trial that would otherwise look like first spikes in bursts (since postisi is actually coming from first spike in next trial)
    incl_spks = sptimes>100&sptimes<=size(visRast,2)-4;
    good_preisis = preisis(incl_spks(1:end-1)); % only plot isis in timeinds 
    good_postisis = postisis(incl_spks(1:end-1));
    good_bursts = find(ismember(find(incl_spks),bursts));
    figure;plot(log10(good_preisis),log10(good_postisis),'o')
    xlim([0 3])
    ylim([0 3])
    xax = get(gca,'xlim');
    yax = get(gca,'ylim');
    xticks([log10(1) log10(10) log10(100) log10(1000)])
    set(gca,'xticklabel',[1 10 100 1000])
    yticks([log10(1) log10(10) log10(100) log10(1000)])
    set(gca,'yticklabel',[1 10 100 1000])
    xlabel('Pre-ISI(ms)','Fontsize',18)
    ylabel('Post-ISI(ms)','Fontsize',18)
    line([log10(100) xax(end)],[log10(4) log10(4)],yax,'Color','k','LineStyle','--','linewidth',2)
    line([log10(100) log10(100)],[0 log10(4)],yax,'Color','k','LineStyle','--','linewidth',2)
    yax = get(gca,'ylim');
    line([log10(4) log10(4)],yax,'Color','k','LineStyle','--','linewidth',2)
    hold on;
    set(gca,'fontsize',18,'linewidth',2)
    plot(log10(good_preisis(good_bursts)),log10(good_postisis(good_bursts)),'k.','markersize',16)
end

return