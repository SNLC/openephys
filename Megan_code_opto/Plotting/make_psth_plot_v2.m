function make_psth_plot_v2(psth,binsize,prestim,stimtime,totaltime,all_light_trials,light_start,light_dur)

% psth = numconds x length(edges)-1 matrix of spike counts (psth) or firing rate(psth_FR, in spikes/s)
    % can be made by make_psth_v2(binsize,edges,which_trials,spike_raster,all_light_trials)
% binsize = size of bins in which to count spikes (in sec)
% prestim = time before visual stimulus onset (in sec)
% stimtime = duration of visual stimulus (in sec)
% totaltime = total time of trial(in sec)
% all_light_trials = 1xnumtrials vector defining each trial's light
    % condition
% light_start = time when the light started (in sec) (if not                                                                                                                                                                                                 an opto
    % experiment, should be []) - for each light condition
% light_dur = duration of light pulse (in sec) (if not an opto
    % experiment, should be [])
    
color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2;0 .8 1; 0 0 1]; % for graphing purposes (first is black, last is green)
% color_mat = [0 0 0; .9 0 .3; 0.50, 0.0780, 0.10]; % for halo (red)
% color_mat = [0 0 0; 0 0.5 .4; 0 0.5 0.05];
% color_mat = [0 0 0; .05 .4 1];
% color_mat = [0 0 0; .166 .674 .188];   %L6
% color_mat = [0 0 0; 0 .2 .9];   %L5
% color_mat = [0 0 0; .85 .325 .098];   % LGN=orange
% color_mat = [0 0 0; .494 .184 .556];   % LP = purple

edges_stim = [-prestim:binsize:(totaltime-prestim-binsize)]'; % x signifies the timepoint of the START of the bin
for c = 1:size(psth,1)
    plot(edges_stim,psth(c,:),'color',color_mat(c,:),'linewidth',2)
    hold on
end

lightconds = unique(all_light_trials);      % different light conditions

xlim([-prestim totaltime-prestim-binsize])  % because points mark the START of the bin
set(gca,'XMinorTick','on')
yax = get(gca,'YLim');
line([0 0], [0 yax(2)]','Color','k','LineStyle','--','LineWidth',2)
% line([stimtime stimtime], [0 yax(2)]','Color','r','LineWidth',2)
xlabel('Time (sec)','fontsize',14)
ylabel('spikes/sec','Fontsize',14)

% draw when light turned on and off
if length(lightconds)>1         % if multiple light conditions
    if length(unique(light_start)) >1  || length(unique(round(light_dur.*1000))) > 2 % if multiple different light start times or duration times, dotted lines to indicate diff start times
        for ii = 1:length(light_start)
            line([light_start(ii)-prestim light_start(ii)-prestim], [0 yax(2)]','Color',color_mat(ii+1,:),'LineStyle','--','LineWidth',2)
            line([light_start(ii)-prestim+light_dur(ii+sum(lightconds==0)) light_start(ii)-prestim+light_dur(ii+sum(lightconds==0))], [0 yax(2)]','Color',color_mat(ii+1,:),'LineStyle','--','LineWidth',2)
        end
    else        % if only one start time, draw patch
        for c = 1:length(light_start)
            x1 = light_start - prestim;     % when the light starts
            patch_start(c) = edges_stim(find(x1(c)-edges_stim>0,1,'last'));     % in case light doesn't evenly start at the beginning of a bin - start the light patch at the earliest bin with any light in it
%             xx = [patch_start(c) patch_start(c) patch_start(c)+light_dur(c+sum(lightconds==0)) patch_start(c)+light_dur(c+sum(lightconds==0)) patch_start(c)];
            xx = [patch_start(c) patch_start(c) patch_start(c)+light_dur(c+sum(lightconds==0)) patch_start(c)+light_dur(c+sum(lightconds==0)) patch_start(c)];
            yy = [0 yax(2) yax(2) 0 0];
            patch(xx, yy, -1 * ones(size(xx)), [0.8 0.8 0.9], 'LineStyle', 'none', 'FaceAlpha',.75)
%              patch(xx, yy, -1 * ones(size(xx)), [.9 0 .3], 'LineStyle', 'none', 'FaceAlpha',.10) %for halo

        end
    end
end
ylim(yax)
    