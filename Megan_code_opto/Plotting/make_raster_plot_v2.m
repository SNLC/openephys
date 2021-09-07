function make_raster_plot_v2(spike_raster,prestim,totaltime,light_trialtypes,light_start,pulse_dur)

% spike_raster = num_trials x time (totaltime*1000, so in 1000Hz samp rate)
    % matrix of 1s and 0s indicating when spikes occurred
% prestim = time before visual stimulus onset (in sec)
% totaltime = total time of trial(in sec)
% light_trialtypes = 1xnumtrials vector defining each trial's light
    % condition
% light_start = time when the light started (in sec) for each light cond
% pulse_dur = duration of light pulse (in sec)
    % pulse_dur in trains experiments)

num_trials = size(spike_raster,1);
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2; 0 .8 1; 0 0 1]; % for graphing purposes (first is black, last is green)
% color_mat = [0 0 0; .9 0 .3; 0.50, 0.0780, 0.10]; % for halo (red)
% color_mat = [0 0 0; 0 0.5 .4; 0 0.5 0.05];
% color_mat = [0 0 0; 0.05 .4 1];
% color_mat = [0 0 0; .166 .674 .188];   %L6
% color_mat = [0 0 0; 0 .2 .9];   %L5
% color_mat = [0 0 0;.85 .325 .098];   % LGN=orange
% color_mat = [0 0 0;.494 .184 .556];   % LP = purple

time_vec = linspace(-prestim,totaltime-prestim,size(spike_raster,2));      % for x-axis 
cond_colors = ones(1,num_trials);      % just in case you have multiple conditions


% if light experiment, add blue shading on raster
if find(light_trialtypes)         % it was an opto experiment if there's more than one "light" condition (one condition would be nolight)
    [cond,idx] = sort(light_trialtypes);    % sort trials by light conditions
    diff_conds = unique(cond);
    new_spike_rast = spike_raster(idx,:);             % reordered according to light condition
    new_spike_rast(~new_spike_rast) = nan;        % so that you don't plot absence of spikes as y=0
    for c = 1:length(diff_conds)
        cond_colors(cond==diff_conds(c))=c;
    end
    start_patch = find(cond>0,1,'first');     % first trial with light
    end_patch = num_trials;         % assumes all trials where cond > 0 had light
    x1 = light_start-prestim;       % when the light starts

    for c = 1:length(diff_conds)
        start_patch = find(cond==diff_conds(c),1,'first')-1;
        end_patch = find(cond==diff_conds(c),1,'last')-1;
        if sum(diff_conds>10)         % if it was a trains experiment
            for p = 1:(light_trialtypes(idx(end_patch)))
                space = 1/(light_trialtypes(idx(end_patch)));
%                 if (p-1)*space < (totaltime-prestim-x1(c-1))
                    patch([x1(c-1)+((p-1)*space) x1(c-1)+((p-1)*space) x1(c-1)+pulse_dur(c)+((p-1)*space) x1(c-1)+pulse_dur(c)+((p-1)*space) x1(c-1)+((p-1)*space)],[start_patch end_patch end_patch start_patch start_patch], [0.8 0.8 0.9], 'LineStyle', 'none', 'FaceAlpha',.75)
%                 end
            end
        elseif diff_conds(c)         % any other light condition
            patch([x1(c-sum(diff_conds==0)) x1(c-sum(diff_conds==0)) x1(c-sum(diff_conds==0))+pulse_dur(c) x1(c-sum(diff_conds==0))+pulse_dur(c) x1(c-sum(diff_conds==0))],[start_patch end_patch end_patch start_patch start_patch], [0.8 0.8 0.9], 'LineStyle', 'none', 'FaceAlpha',.75)
%             patch([x1(c-sum(diff_conds==0)) x1(c-sum(diff_conds==0)) x1(c-sum(diff_conds==0))+pulse_dur(c) x1(c-sum(diff_conds==0))+pulse_dur(c) x1(c-sum(diff_conds==0))],[start_patch end_patch end_patch start_patch start_patch], [.9 0 .3], 'LineStyle', 'none', 'FaceAlpha',.10) %for halo
        end
        if length(diff_conds)>1
            line([-prestim totaltime-prestim], [end_patch end_patch]','Color','r','LineStyle','--')
        end
    end
    hold on
else
    new_spike_rast = spike_raster;      %non-opto experiment
    new_spike_rast(~new_spike_rast) = nan; 
end

if exist('xrange','var')        
    rasterrange = [find(time_vec==.4):find(time_vec==.8)];
else
    rasterrange = 1:size(new_spike_rast,2);
end

% draw raster
for t = 1:num_trials     % for each trial
    new_spike_rast(t,find(new_spike_rast(t,:))) = new_spike_rast(t,find(new_spike_rast(t,:)))+(t-1);        % plot 'ones' in vector so that y = trial number
    plot(time_vec(rasterrange),new_spike_rast(t,rasterrange),'.','Color',color_mat(cond_colors(t),:),'MarkerSize',6)    
    hold on
end

ylim([0 num_trials])
xlim([-prestim totaltime-prestim])
xlabel('Time from visual stimulus onset (sec)','fontsize',14)
ylabel('Trial (by condition)','fontsize',14)

return