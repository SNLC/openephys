function [unitinfo,FRs,tuning,waveforms] = unit_analysis_noopto(unit,field_trials,unit_times,params,rez,normalizing,makeplots,fig_dir)

% modified from intan_unit_analysis.m by MAK on 8/8/2017. 
% modified to be more efficient and to postpone most stats until later
% to make things move faster.

% inputs:
% unit = number corresponding to the cluster's name given by Phy (e.g. 27)
% field_trials = num_trials x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% unit_times = samples corresponding to spike times of 'unit' cluster.
    % Still in ORIGINAL sampling rate 
    % extracted as unit_times = spike_times(clusters==unit) 
    % (spike_times and clusters are .npy files directly from Phy)
% params = structure contraining experiment parameters (from analyzer file,
    % get_exp_params.m and get_lightstim_v2.m)
% rez = .mat file output from kilosort (for extracting waveforms)
    % e.g., load(sprintf('%s\\rez.mat',exp_path))
% normalizing = what to use as baseline (e.g. 'blanks' or 'prestim' -
    % anything else is taken as no baseline)
% makeplots = 1 if you want to print plots, 0 if not
% fig_dir = directory where to save figures (leave empty if makeplots=0)
    % ex. fig_dir = sprintf('H:\\LPproject\\LPresults\\%s\\Figures',exp_name)
    
%% Which plots do you want?
plot_rast = 1;      % raster plot
plot_psthV = 1;     % psth of visual trials
plot_psthB = 1;     % psth of blank trials
plot_rasterZ =0;     % plot zoomed in visual trial psth
plot_psthPref = 1;  % psth of preferred orientation trials
plot_psthPrefSFTF = 0;  % psth of trials at preferred SF and TF
plot_psthstdSFTF = 0;  % psth of trials at standard SF (.04) and TF (2Hz)
plot_run = 1;       % bar graph of running v. stationary
plot_FRs = 1;       % bar graph of firing rates at different trial points
plot_ori = 1;       % orientation tuning curve
plot_SF = 0;        % bar plot of FR at diff spatial frequencies
plot_TF = 0;        % bar plot of FR at diff temporal frequencies
plot_waveforms = 1; % plot the unit's waveform


%% rename some inputs
% trial info from params
exp_name = params.exp_name;         % e.g. 'T22_ramp'
exp_type = params.exp_type;
amp_sr = params.amp_sr;
nchs = params.nchs;
trial_type = params.trial_type;     % nxm matrix, n=number of trials, m=number of IVs
IVs = params.IVs;                   % independent variables
prestim_ms = params.prestim*1000;
stimtime_ms = params.stimtime*1000;
poststim_ms = params.poststim*1000;
total_time_ms = prestim_ms+stimtime_ms+poststim_ms;
onset_ms = params.onset*1000;               % amount of time from start of visual stimulus considered "onset"

%% set up fake "light" params to make code work for non-opto experiments (temp...)
all_light = zeros(1,size(field_trials,1));
lighttime = 1000;   % TEMP - manually set to 1s analysis window (light isn't actually on...)
params.lighttime = 1; % TEMP - manually set to 1s analysis window (light isn't actually on...)
params.av_light_start = [];
params.pulse_dur = [];
params.light_dur = [];

%% extract spike times 
unit_times_ds = floor(unit_times./(amp_sr/1000));   % change sample #s to account for downsampling
unit_times_ds = unit_times_ds + 1; % has to be +1 because spike_times starts at 0, but the min possible field_trials value could be 1
spike_raster = make_raster(unit_times_ds,field_trials,1000,total_time_ms/1000);      % trials x time matrix of spike times (0s and 1s)

%% count spikes during periods of interest
disp('Counting spikes...')
num_trials = size(field_trials,1);
% window = round([2*prestim_ms+1 2*prestim_ms+1000]); % TEMP: set analysis time window to [1001:2000]
window = [prestim_ms+251 prestim_ms+stimtime_ms-250];   % 1/14/21 - set to match inact experiments
ev_lighttime = (diff(window)+1)/1000;

spikes_prestim = sum(spike_raster(:,1:prestim_ms),2);
spikes_ev = sum(spike_raster(:,window(1):window(2)),2);
spikes_ev_half = sum(spike_raster(:,window(2)-round(lighttime/2)+1:window(2)),2);      % changed to SECOND HALF 4/6/19
spikes_ev_lighton = sum(spike_raster(:,window(1):window(1)-1+onset_ms),2);      % currently using first av_light_start time as default
spikes_ev_early = sum(spike_raster(:,window(1)+onset_ms:window(1)-1+round(lighttime/2)),2);    % from 100ms after light onset to halfway through light duration (currently using first av_light_start time as default)
spikes_onset = sum(spike_raster(:,prestim_ms+1:prestim_ms+onset_ms),2);

%% Calculate important firing rates
disp('Calculating firing rates...')
% identify variables of interest
orivar = find(strcmp(IVs,'ori'));
oriconds = unique(trial_type(:,orivar));     % get different orientations, including blank
oriinds = oriconds<=360;          % get different orientations, EXCLUDING blanks (which are indicated by 999)
oris = oriconds(oriinds);
runvar = find(strcmp(IVs,'running'));
runconds = unique(trial_type(:,runvar));        % and running vs stationary (should always be two levels)
lightvar = find(strcmp(IVs,'light_bit'));
lightconds = unique(trial_type(:,lightvar));     % get different levels of light variable

% calculate FRs for VISUAL and BLANK trials separately at each light
% condition
% 10/31/17: change to stationary trials only
for i = 1:length(lightconds)                  % get firing rate, across all other variables
    [spikerate_prestim(i), spikerateSE_prestim(i)] = calc_firing_rates(spikes_prestim,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,runvar)==0)),params.prestim); 
    [spikerate_onset(i), spikerateSE_onset(i)] = calc_firing_rates(spikes_onset,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,runvar)==0)),params.onset);
    if find(strcmp(IVs,'s_freq'))               % TEMPORARY - if SF/TF experiment, "spikerate_visual" vals are for standard SF=.04, TF=2Hz ONLY
        SFvar = find(strcmp(IVs,'s_freq'));
        TFvar = find(strcmp(IVs,'t_period'));
        [spikerate_visual(:,i), spikerateSE_visual(:,i)] = calc_firing_rates(spikes_ev,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)&(trial_type(:,SFvar)==.04)&(trial_type(:,TFvar)==30)),ev_lighttime);
        [spikerate_visual_half(i), spikerateSE_visual_half(i)] = calc_firing_rates(spikes_ev_half,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)&(trial_type(:,SFvar)==.04)&(trial_type(:,TFvar)==30)),params.lighttime/2);
        [spikerate_visual_lighton(i), spikerateSE_visual_lighton(i)] = calc_firing_rates(spikes_ev_lighton,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)&(trial_type(:,SFvar)==.04)&(trial_type(:,TFvar)==30)),params.onset);
        [spikerate_visual_early(i), spikerateSE_visual_early(i)] = calc_firing_rates(spikes_ev_early,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)&(trial_type(:,SFvar)==.04)&(trial_type(:,TFvar)==30)),params.lighttime/2-params.onset);
    else
        [spikerate_visual(:,i), spikerateSE_visual(:,i)] = calc_firing_rates(spikes_ev,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)),ev_lighttime);
        [spikerate_visual_half(i), spikerateSE_visual_half(i)] = calc_firing_rates(spikes_ev_half,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)),params.lighttime/2);
        [spikerate_visual_lighton(i), spikerateSE_visual_lighton(i)] = calc_firing_rates(spikes_ev_lighton,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)),params.onset);
        [spikerate_visual_early(i), spikerateSE_visual_early(i)] = calc_firing_rates(spikes_ev_early,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)<=360)&(trial_type(:,runvar)==0)),params.lighttime/2-params.onset);
    end
    [spikerate_blank(:,i), spikerateSE_blank(:,i)] = calc_firing_rates(spikes_ev,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)>360)&(trial_type(:,runvar)==0)),ev_lighttime);
    [spikerate_blank_half(i), spikerateSE_blank_half(i)] = calc_firing_rates(spikes_ev_half,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)>360)&(trial_type(:,runvar)==0)),params.lighttime/2);
    [spikerate_blank_lighton(i), spikerateSE_blank_lighton(i)] = calc_firing_rates(spikes_ev_lighton,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)>360)&(trial_type(:,runvar)==0)),params.onset);
    [spikerate_blank_early(i), spikerateSE_blank_early(i)] = calc_firing_rates(spikes_ev_early,find((trial_type(:,lightvar) == lightconds(i))&(trial_type(:,orivar)>360)&(trial_type(:,runvar)==0)),params.lighttime/2-params.onset);

end

% if SF and TF conditions
if find(strcmp(IVs,'s_freq'))
    SFs = unique(trial_type(:,SFvar));
    SFs = SFs(SFs<999);
    for i = 1:length(SFs)
        for ii = 1:length(lightconds)
            [SF_FR(i,ii) SF_FR_SE(i,ii)] = calc_firing_rates(spikes_ev(:,1),find((trial_type(:,lightvar)==lightconds(ii))&(trial_type(:,SFvar)==SFs(i))&(trial_type(:,runvar)==0)),ev_lighttime);     % in case of multiple start times, using evoked firing rate from FIRST start time condition; only in STATIONARY trials
        end
    end
end
if find(strcmp(IVs,'t_period'))
    TFs = unique(trial_type(:,TFvar));
    TFs = TFs(TFs<999);
    for i = 1:length(TFs)
        for ii = 1:length(lightconds)
            [TF_FR(i,ii) TF_FR_SE(i,ii)] = calc_firing_rates(spikes_ev(:,1),find((trial_type(:,lightvar)==lightconds(ii))&(trial_type(:,TFvar)==TFs(i))&(trial_type(:,runvar)==0)),ev_lighttime);         % in case of multiple start times, using evoked firing rate from FIRST start time condition; only in STATIONARY trials
        end
    end
end

% next, make 3D matrix of evoked FRs according to ori, lightcond, and run
% variables

for lc = 1:length(lightconds)    
    lc_trials{lc} = find(trial_type(:,lightvar) == lightconds(lc));     % make cell in case number of trials per lightcond were unequal
    for o = 1:length(oriconds)
        oricond_trials = find(trial_type(:,orivar) == oriconds(o));    % create matrix of which trials were at given orientation
        trials_per_ori(lc,o) = length(intersect(oricond_trials,lc_trials{lc}));   
        for r = 1:length(runconds)
            run_trials = find(trial_type(:,runvar) == runconds(r));
%             [spikerate_bycond(o,lc,r), spikerateSE_bycond(o,lc,r)] = calc_firing_rates(spikes_ev_half',intersect(intersect(lc_trials{lc},oricond_trials),run_trials),params.lighttime/2);   % currently using FIRST HALF of light period
              [spikerate_bycond(o,lc,r), spikerateSE_bycond(o,lc,r)] = calc_firing_rates(spikes_ev(:,1),intersect(intersect(lc_trials{lc},oricond_trials),run_trials),ev_lighttime);   % currently using WHOLE light period
      end
    end
end

% finally, normalize firing rates (by 'prestim','blanks',or 'none')
if strcmp(normalizing,'prestim')
    baseline = calc_firing_rates(spikes_prestim,find(trial_type(:,runvar)==0),params.prestim);      % across all STATIONARY trials
elseif strcmp(normalizing,'blanks')
    baseline = spikerate_bycond(oriconds>360,lightconds==0,runconds==0);   % baseline FR defined as evoked FR from BLANK,STATIONARY,NOLIGHT trials
else
    baseline = 0;
end
spikerate_bycond_norm = spikerate_bycond - baseline;

%% get tuning curves and OSI and DSI
for lc = 1:length(lightconds)
    tuning_curve(lc,:) = spikerate_bycond(oriinds,lc,runconds==0);
    tuning_curve_norm(lc,:) = spikerate_bycond_norm(oriinds,lc,runconds==0);          % baseline subtracted
    tuning_curveSE(lc,:) = spikerateSE_bycond(oriinds,lc,runconds==0);
    if sum(oriconds<=360) > 1
        [OSI(lc),OSI_CV(lc),DSI(lc),DSI_CV(lc)] = calcOSIDSI(tuning_curve(lc,:),oris');       % using only STATIONARY trials (NOT baseline-subtracted, but check whether this is correct!!!!)
    end
end

%% get PSTHs
% useful trial indicators:
vis_trials = find(trial_type(:,1)==1);
blank_trials = find(trial_type(:,1)==0);
binsize = .025;         % 25 ms
edges = [0:binsize:total_time_ms/1000];      % in sec
[~,psthV] = make_psth_v2(binsize,edges,ismember(1:num_trials,vis_trials),spike_raster,all_light);
[~,psthB] = make_psth_v2(binsize,edges,ismember(1:num_trials,blank_trials),spike_raster,all_light);

%% get waveform information
disp('Getting waveform information...')
exp_path = cd;          % assuming CD was set in analysis_master
if ~isempty(rez)
    if ~exist('rez.mat','file')     % if not in current dir (for when concatenating files)
        exp_path = fileparts(exp_path);
    end
    [waveforms_microV,max_ch,shank] = readWaveformsFromRez_K2(unit,exp_path,rez);      % default to read waveforms from Rez
elseif exist(sprintf('/%s/Cluster_%s_waveforms.mat',exp_path,num2str(unit)),'file')
    load(sprintf('/%s/Cluster_%s_waveforms.mat',exp_path,num2str(unit)));
else
    % still need to make compatible with old phy-sorted data
end

% determine layer
if exist('layers.mat','file')           % in current path
    load('layers.mat');
else
    define_layers(25,32,exp_path,1);        %** currently hard-coded for 32ch NN probes - need to change!!
    load('layers.mat');
end
if exist('layers_shank','var') && ~isempty(shank)
    layers = layers(layers_shank==shank);
end
layer = layers(max_ch);           %TEMP +1


%% save results
disp('Saving....')
% firing rate stuff
FRs.psthVisual = psthV;
FRs.psthBlank = psthB;
FRs.prestim = [spikerate_prestim; spikerateSE_prestim];
FRs.visual.ev = [spikerate_visual; spikerateSE_visual];
FRs.visual.evstart = [spikerate_visual_early; spikerateSE_visual_early];
FRs.visual.evlightonset = [spikerate_visual_lighton; spikerateSE_visual_lighton];
FRs.visual.evlate = [spikerate_visual_half; spikerateSE_visual_half];
FRs.blank.ev = [spikerate_blank; spikerateSE_blank];
FRs.blank.evstart = [spikerate_blank_early; spikerateSE_blank_early];
FRs.blank.evlightonset = [spikerate_blank_lighton; spikerateSE_blank_lighton];
FRs.blank.evlate = [spikerate_blank_half; spikerateSE_blank_half];
FRs.onset = [spikerate_onset; spikerateSE_onset];
FRs.baseline = baseline;
FRs.baselinedef = normalizing;

% tuning stuff
tuning.curve = tuning_curve;
tuning.normcurve = tuning_curve_norm;
if sum(oriconds<=360) > 1
    tuning.OSI = OSI;
    tuning.OSI_CV = OSI_CV;
    tuning.DSI = DSI;
    tuning.DSI_CV = DSI_CV;
end
if find(strcmp(IVs,'s_freq'))
    tuning.SF = [SF_FR; SF_FR_SE];
end
if find(strcmp(IVs,'t_period'))
    tuning.TF = [TF_FR; TF_FR_SE];
end

% general unit info
% unitinfo.name = sprintf('%s_%d',exp_name,unit);
unitinfo.name = unit;
unitinfo.layer = layer;
unitinfo.numspikes = length(unit_times);
unitinfo.rast = spike_raster;

% waveforms
waveforms.microV = waveforms_microV;
waveforms.shank = shank;
waveforms.max_ch = max_ch;
    

%% plot stuff!
if makeplots
    disp('Plotting...')
    
    % for keeping track of which plots to make
    nplots = sum([plot_rast,plot_psthV,plot_psthB,plot_rasterZ,plot_psthPref,plot_run,plot_FRs,plot_ori,plot_SF,plot_TF,plot_waveforms]);
    if nplots <= 8
        yy = 2;
        xx = ceil(nplots/yy);
    else
        yy = 3;
        xx = ceil(nplots/yy);
    end
    count_plots = 1;
    
    % make appropriate figure legends
    if strcmp(exp_type,'trains')
        for i = 2:length(lightconds)        % for graphing purposes
            legend_labels{i-1} = sprintf('%dHz',lightconds(i));
        end
    elseif strcmp(exp_type,'intensities')
        legend_labels = {'Low light','Medium light','High light'};
        full_legend_labels = {'Light OFF','Low light','Medium light','High light'};
    else
        legend_labels = {'Light ON'};
    end
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2;0 .8 1; 0 0 1]; % for graphing purposes (first is black, last is green)
    
    % make 2-by-4 figure array containing raster plot, PSTHs for visual and
    % blank trials, zoom-in of PSTH, running v. stationary barplot, prestim
    % v evoked v blanks barplot, orientation tuning curve, and waveform
    fig_title = ['Cluster_' num2str(unit)];
    clust_fig = figure('name', fig_title); 
    [type,idx] = sort(all_light);    % sort trials by light conditions
    
    % Raster plot
    if plot_rast
        subplot(yy,xx,count_plots)
        make_raster_plot_v2(spike_raster,params.prestim,total_time_ms/1000,all_light,params.av_light_start,params.pulse_dur)
        set(gca,'FontSize',18);
        title('Raster plot')
        count_plots = count_plots+1;
    end
    
    % PSTH - visual trials
    if plot_psthV
        subplot(yy,xx,count_plots)
        binsize = .025;         % 25 ms
    %     binsize = .05;
        make_psth_plot_v2(psthV,binsize,params.prestim,params.stimtime,total_time_ms/1000,trial_type(:,lightvar),params.av_light_start,params.light_dur);
        title('PSTH - visual trials','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end
    
    % PSTH - blank trials
    if plot_psthB
        subplot(yy,xx,count_plots)
        make_psth_plot_v2(psthB,binsize,params.prestim,params.stimtime,total_time_ms/1000,trial_type(:,lightvar),params.av_light_start,params.light_dur);
        title('PSTH - blank trials','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end
    
    % raster - zoom
    if plot_rasterZ
        subplot(yy,xx,count_plots)
        make_raster_plot_v2(spike_raster(:,round(min(av_light_start))-100:round(min(av_light_start))+299),-(((round(min(av_light_start))-100))/1000-params.prestim),length(round(min(av_light_start))-100:round(min(av_light_start))+299)/1000,all_light,.1*ones(length(av_light_start)),params.pulse_dur)
        title('Raster plot (zoom)','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end

    % PSTH - preferred orientation trials
    if plot_psthPref
        subplot(yy,xx,count_plots)
        [~,pref_ori] = max(abs(tuning_curve_norm(1,:)));    % 'preferred orientation' = maximum change from baseline (i.e., not just ori with largest FR)
        [~,psthP] = make_psth_v2(binsize,edges,ismember(1:num_trials,find(trial_type(:,orivar)==oriconds(pref_ori))),spike_raster,all_light);
        make_psth_plot_v2(psthP,binsize,params.prestim,params.stimtime,total_time_ms/1000,trial_type(:,lightvar),params.av_light_start,params.light_dur); % make PSTH from trials with preferred orientation only
        title('PSTH - Preferred Orientation','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end
    
    % PSTH - preferred SF and TF trials
    if plot_psthPrefSFTF
        subplot(yy,xx,count_plots)
        [~,pref_SF] = max(abs(SF_FR(:,1) - repmat(baseline,size(SF_FR,1),1))); 
        [~,pref_TF] = max(abs(TF_FR(:,1) - repmat(baseline,size(SF_FR,1),1))); 
        [~,psthPsftf] = make_psth_v2(binsize,edges,ismember(1:num_trials,find((trial_type(:,SFvar)==SFs(pref_SF))&(trial_type(:,TFvar)==TFs(pref_TF)))),spike_raster,all_light);
        make_psth_plot_v2(psthPsftf,binsize,params.prestim,params.stimtime,total_time_ms/1000,trial_type(:,lightvar),params.av_light_start,params.light_dur); % make PSTH from trials with preferred orientation only
        title('PSTH - Preferred SF&TF','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end
    
    % PSTH - standard SF and TF trials
    if plot_psthstdSFTF
        subplot(yy,xx,count_plots)
        [~,psthstd] = make_psth_v2(binsize,edges,ismember(1:num_trials,find((trial_type(:,SFvar)==.04)&(trial_type(:,TFvar)==30))),spike_raster,all_light);
        make_psth_plot_v2(psthstd,binsize,params.prestim,params.stimtime,total_time_ms/1000,trial_type(:,lightvar),params.av_light_start,params.light_dur); % make PSTH from trials with preferred orientation only
        title('PSTH - Standard SF&TF','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end
    
    % running v stationary firing rates
    if plot_run
        subplot(yy,xx,count_plots)
        for lc = 1:length(lightconds)
            [spikerate_run(lc) spikerateSE_run(lc)] = calc_firing_rates(spikes_ev(:,1),find(trial_type(:,runvar)==1&trial_type(:,lightvar)==lightconds(lc)),ev_lighttime);     % TEMP - 2nd column for iC++ exps
            [spikerate_stat(lc) spikerateSE_stat(lc)] = calc_firing_rates(spikes_ev(:,1),find(trial_type(:,runvar)==0&trial_type(:,lightvar)==lightconds(lc)),ev_lighttime);
        end
        runbar = bargraph([spikerate_run; spikerate_stat],...
            [spikerateSE_run; spikerateSE_stat]);
        set(get(gca,'YLabel'),'String','Mean FR (spikes/s)','FontSize',18)
        set(gca,'XTicklabel',sprintf('Running\n Stationary'))
        for i = 1:length(lightconds)
            set(runbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
        end
        title('Firing rate - running vs. stationary','FontSize',24)
        legend off
        set(gca,'FontSize',18);
        xax = get(gca,'xtick');
        xlim([xax(1)-.5 xax(end)+.5])
        count_plots = count_plots+1;
    end
    
    % prestim vs evoked vs blank FRs
    if plot_FRs
        subplot(yy,xx,count_plots)
        FR = [spikerate_prestim; spikerate_onset; spikerate_visual];
        SE = [spikerateSE_prestim; spikerateSE_onset; spikerateSE_visual];
%         if sum(FR(end,:)~=FR(end-1,:))==0       % cheat if using same evoked window for diff starttimes
%             FR = FR(1:end-(length(lightconds)-2),:);
%             SE = SE(1:end-(length(lightconds)-2),:);
%         end
        if length(unique(trial_type(:,1))) > 1          % if blank trials
    %         for lc = 1:length(lightconds)
    %             if lightconds(lc)
    %                 [spikerate_blank(lc) spikerateSE_blank(lc)] = calc_firing_rates(spikes_ev(:,lc-1),intersect(blank_trials,light_trials{lc-1}),lighttime);    % TEMP for iC++ exp - column indicator
    %             else
    %                 [spikerate_blank(lc) spikerateSE_blank(lc)] = calc_firing_rates(spikes_ev(:,1),intersect(blank_trials,nolight_trials),lighttime);
    %             end
    %         end
%             FR = [FR(1:end-1,:); spikerate_visual_lighton; spikerate_visual_early; spikerate_visual_half; spikerate_blank(1,:)];
%             SE = [SE(1:end-1,:); spikerateSE_visual_lighton; spikerateSE_visual_early; spikerateSE_visual_half; spikerateSE_blank(1,:)];
%             xcondslabel = sprintf('Prestim\n Onset\n LightonS\n Early\n Late\n Blanks');
            FR = [FR(1:end-1,:); spikerate_visual_lighton; spikerate_visual; spikerate_blank(1,:)];
            SE = [SE(1:end-1,:); spikerateSE_visual_lighton; spikerateSE_visual; spikerateSE_blank(1,:)];
            xcondslabel = sprintf('Prestim\n Onset\n LightonS\n Evoked\n Blanks');
        else
            xcondslabel = sprintf('Prestim\n Onset\n Evoked');
        end
        evokedbar = bargraph(FR,SE);
        set(get(gca,'YLabel'),'String','Mean FR (spikes/sec)','FontSize',24)
        set(gca,'XTicklabel',xcondslabel,'FontSize',18)
        for i = 1:length(lightconds)
            set(evokedbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
        end
        title('Firing rate','FontSize',24)
        legend off
         xax = get(gca,'xtick');
        xlim([xax(1)-.5 xax(end)+.5])
        count_plots = count_plots+1;
    end
    
    % Tuning curves (stationary trials, baseline subtracted)
    if plot_ori
        subplot(yy,xx,count_plots)
        vals = [];
        maxoris = [];
        sigconds = [];
        for lc = 1:length(lightconds)
            shadedErrorBar(oris,tuning_curve(lc,:),spikerateSE_bycond(oriinds,lc,runconds==0),{'Color',color_mat(lc,:),'linewidth',2},1);  % plot ori tuning in NO RUN trials (NOT baseline-subtracted)
            hold on
%             plot(oris,repmat(spikerate_bycond(oriinds==0,lc,runconds==0),length(oris),1),'Color',color_mat(lc,:),'linewidth',1,'linestyle','--');  % plot blank trial FRs for comparison
%             if tuned_sig(lc) < .05
%                 [val,maxori] = max(tuning_curve_norm(lc,:));
%                 vals = [vals val];
%                 maxoris = [maxoris maxori];
%                 sigconds = [sigconds lc];
%             end
        end
        yax = get(gca,'YLim');
        for ii = 1:length(vals)
            plot(oris(maxoris(ii)), yax(2)-.5*ii, '*','Color',color_mat(sigconds(ii),:))
        end
        ylabel('Firing rate (Hz)','FontSize',24)
        xlabel('Orientation (degrees)','FontSize',24)
        xlim([0 max(oris)])
        title('Orientation tuning during light period','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end

    % SF
    if plot_SF
        subplot(yy,xx,count_plots)
        
         for lc = 1:length(lightconds)
            shadedErrorBar(SFs,SF_FR(:,lc)'-repmat(baseline,1,size(SF_FR,1)),SF_FR_SE(:,lc)',{'Color',color_mat(lc,:),'linewidth',2},1);  % plot ori tuning in NO RUN trials (NOT baseline-subtracted)
            hold on
         end
         ylabel('Firing rate (Hz)','FontSize',24)
        xlabel('SF (cpd)','FontSize',24)
        
%         evokedbar = bargraph(SF_FR,SF_FR_SE);
%         set(get(gca,'YLabel'),'String','Mean FR (spikes/sec)','Fontsize',10)
%         set(gca,'XTicklabel',num2str(SFs),'Fontsize',10)
%         for i = 1:length(lightconds)
%                 set(evokedbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
%         end
        title('SF tuning','FontSize',24)
        legend off
        count_plots = count_plots+1;
    end
    
        
    % TF
    if plot_TF
        subplot(yy,xx,count_plots)
        
        for lc = 1:length(lightconds)
            shadedErrorBar(60./TFs,TF_FR(:,lc)'-repmat(baseline,1,size(TF_FR,1)),TF_FR_SE(:,lc)',{'Color',color_mat(lc,:),'linewidth',2},1);  % plot ori tuning in NO RUN trials (NOT baseline-subtracted)
            hold on
         end
         ylabel('Firing rate (Hz)','FontSize',24)
        xlabel('TF (Hz)','FontSize',24)
        
%         evokedbar = bargraph(TF_FR,TF_FR_SE);
%         set(get(gca,'YLabel'),'String','Mean FR (spikes/sec)','Fontsize',10)
%         set(gca,'XTicklabel',num2str(TFs),'Fontsize',10)
%         for i = 1:length(lightconds)
%                 set(evokedbar(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
%         end
        title('TF tuning','FontSize',24)
        legend off
        count_plots = count_plots+1;
    end
    
    % Plot waveforms
    if plot_waveforms
        subplot(yy,xx,count_plots)
        t = linspace(0,(size(waveforms_microV,1)-1)/20,size(waveforms_microV,1));    % convert to ms
        plot(t,waveforms_microV,'LineWidth',2);
        xlim([0 max(t)])
        hold on 
        title(sprintf('Average waveform of spikes (Max ch: %d)',max_ch),'FontSize',24)
        ylabel('Amplitude (uV)','FontSize',24)
        xlabel('Time (ms)','FontSize',24)
        set(gca,'FontSize',18);
        count_plots = count_plots+1;
    end

     %% save figs
    if ~exist(fig_dir)
        mkdir(fig_dir);
    end
    xSize = 24; ySize = 11;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[0 0 xSize*50 ySize*50])
    if layer == 5
        layer_name = 'L5A';
    elseif layer == 5.5
        layer_name = 'L5B';
    elseif layer == 2.5
        layer_name = 'L23';
    elseif strfind(exp_name,'LP')
        layer_name = num2str(max_ch);      % in absence of cortical layers (i.e. in LP), just use channel # (numbered from top to bottom0
    else
        layer_name = sprintf('L%s',num2str(layer));
    end
    if ~isempty(shank)
        layer_name = strcat(sprintf('shank%d_',shank),layer_name);
    end
    save_clust_name= sprintf('%s\\%s_%s_Cluster%d',fig_dir,exp_name,layer_name,unit);
    print(clust_fig,'-dpng',strcat(save_clust_name,'.png'))
%     print2eps(save_clust_name,clust_fig)
    close all
end

end

function rast = make_raster(unit_times,field_trials,Fs,total_time)

% unit_times = spike times for current cluster (in 1000Hz sampling rate),
    % starting from 1
% field_trials - trial num x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% Fs = sampling rate (in hz)
% total_time = total time (in sec)

num_trials = size(field_trials,1);
dt = 1/Fs;
time = 0:dt:total_time-dt;
rast = zeros(num_trials,length(time));
for t = 1:num_trials
    first_spike = find(unit_times >= field_trials(t,1), 1, 'first');
    last_spike = find(unit_times < (field_trials(t,1)+total_time*1000), 1, 'last');
    rast(t,unit_times(first_spike:last_spike) - field_trials(t,1)+1) = 1;  % spike times (in secs) from trial start
end

end


function [spikerate, spikerate_SE] = calc_firing_rates(spikes,...
    which_trials,period)

% spikes = column vector of number of spikes by trial for unit of interest
% which_trials = vector of which trials you want to include
% period = length of period (in seconds) that number of spikes was counted 
    % from (e.g. evoked: 1.8s)

spikerate = mean(spikes(which_trials,:)./period);
spikerate_SE = std(spikes(which_trials,:)./period)./sqrt(length(which_trials));

end

