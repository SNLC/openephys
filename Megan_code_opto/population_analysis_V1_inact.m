function population_analysis_V1_inact(pop,mani)
% new version created 10/21/20 to better accomodate inactivation experiments with variable different opto conditions experiments (MAK)
% checked throuroughly 1/9/21
% inputs:
% pop = population (e.g., 'L6')
% mani = manipulation (e.g., 'ChR2')

% If downloading data from Mendeley, start from line 286 ('Identify
% experiment parameters for multi-experiment analyses') (but will need to
% set your own main_dir and fig_dir

vissigdef = 'blanks';  % defined coparison for determining visual responsivity: 'prepost' for comparing pre- to post-stimulus in preferred visual trials, or 'blanks' to compare preferred visual to blank trials

%% set up directories
main_dir = 'H:\LPproject\LPresults';
if strcmpi(pop,'L6') && strcmpi(mani,'halo_laser')
    exp_paths = {'H:\LPproject\MH19\V1\2018-12-04_19-59-16_V1_Halo',...
        'H:\LPproject\MH20\V1\2018-12-04_12-26-57_V1_Halo',...
        'H:\LPproject\MH25\V1\2019-04-10_13-20-15_V1_Halo',...
        'H:\LPproject\MH23\V1\2019-04-11_12-33-16_V1_Halo_4.5mW',...
        'H:\LPproject\MH24\V1\2019-04-12_11-57-09_V1_Halo_4.5mw',...
        'H:\LPproject\MH22\V1\2019-04-26_12-33-34_V1_Halo_4.5mW'};
     shanks = {0,0,[0,1],[0,1],[0,1],[0,1]};
elseif strcmpi(pop,'L6') && strcmpi(mani,'halo_led')
        exp_paths = {'Z:\L6\MH27\V1\2019-08-28_13-09-09_V1_Halo',...
        'Z:\L6\MH28\V1\2019-08-27_17-38-29_V1_halo',...
        'Z:\L6\MH31\V1\2020-01-02_16-27-43_Halo',...
        'Z:\L6\MH32\V1\2020-01-03_17-07-12_Halo'};
    shanks = {0,0,0,0};
elseif strcmpi(pop,'L6') && strcmpi(mani,'halo_all')
    exp_paths = {'H:\LPproject\MH19\V1\2018-12-04_19-59-16_V1_Halo',...
        'H:\LPproject\MH25\V1\2019-04-10_13-20-15_V1_Halo',...
        'H:\LPproject\MH23\V1\2019-04-11_13-00-08_V1_Halo_5.5mW',...
        'H:\LPproject\MH24\V1\2019-04-12_11-27-33_V1_Halo',...
        'Z:\L6\MH27\V1\2019-08-28_13-09-09_V1_Halo',...
        'Z:\L6\MH28\V1\2019-08-27_17-38-29_V1_halo',...
        'Z:\L6\MH30\V1\2020-01-03_12-56-56_Halo',...
        'Z:\L6\MH31\V1\2020-01-02_16-27-43_Halo',...
        'Z:\L6\MH32\V1\2020-01-03_17-07-12_Halo',...
        'Z:\L6\MH33\V1\2020-04-09_18-50-28_halo',...
        'Z:\L6\MH34\V1\2020-04-10_12-35-28_halo',...
        'Z:\L6\MH35\V1\2020-04-10_16-53-34_halo',...
        'Z:\L5&L6\MD9\V1pen1\2020-11-25_16-16-51_driftgrat_2LEDs',...
        'Z:\L5&L6\MD9\V1pen2\2020-11-26_11-20-32_driftgrat_2LEDs'};
    shanks = {0,[0,1],[0,1],[0,1],0,0,0,0,0,0,0,0,0,0};
    lightconds = {2,2,2,2,2,2,1,1,1,1,1,1,1,1};
    lightcond = 1;
elseif strcmpi(pop,'L6') && strcmpi(mani,'halo_wRF')
    exp_paths = {'Z:\L6\MH27\V1\2019-08-28_13-09-09_V1_Halo',...
        'Z:\L6\MH28\V1\2019-08-27_17-38-29_V1_halo',...
        'Z:\L6\MH30\V1\2020-01-03_12-56-56_Halo',...
        'Z:\L6\MH31\V1\2020-01-02_16-27-43_Halo',...
        'Z:\L6\MH33\V1\2020-04-09_18-50-28_halo',...
        'Z:\L6\MH34\V1\2020-04-10_12-35-28_halo',...
        'Z:\L6\MH35\V1\2020-04-10_16-53-34_halo',...
        'Z:\L5&L6\MD9\V1pen1\2020-11-25_16-16-51_driftgrat_2LEDs',...
        'Z:\L5&L6\MD9\V1pen2\2020-11-26_11-20-32_driftgrat_2LEDs'};
    shanks = {0,0,0,0,0,0,0,0,0};
    lightconds = {2,2,1,1,1,1,1,1,1};
    lightcond = 1;
elseif strcmpi(pop,'L6') && strcmpi(mani,'dta')
    exp_paths = {'Z:\dTa\MDTA2\V1\2020-03-12_12-32-48_driftgrat_DTA',...
        'Z:\dTa\MDTA3\V1\2020-03-13_17-32-38_driftgrat_DTA'};
    shanks={0,0};
    lightconds = {[],[]};
    lightcond = 1;
elseif strcmpi(pop,'L6') && strcmpi(mani,'ChR2')
    exp_paths = {'H:\LPproject\M38\2019-05-24_12-46-43_V1_diffintensities',...
        'J:\LPproject\M39\2019-05-23_12-05-08_V1_diffintensities',...
        'J:\LPproject\M40\2019-05-25_14-09-42_V1_diffintensities',...
        'J:\LPproject\M32\V1\2019-02-14_16-40-51_V1_diffintensities',...
        'J:\LPproject\M36\V1\2019-02-21_15-36-40_V1_diffintensities'};
    shanks = {[0,1],[0,1],[0,1],[0],[0]};
elseif strcmpi(pop,'L6') && strcmpi(mani,'ChR2_trains')
    exp_paths = {'H:\LPproject\M38\2019-05-24_13-52-50_V1_trains',...
        'J:\LPproject\M39\2019-05-23_13-12-13_V1_trains',...
        'J:\LPproject\M40\2019-05-25_15-07-54_V1_trains'};
    shanks = {[0,1],[0,1],[0,1]};
elseif strcmpi(pop,'L5') && strcmpi(mani,'ChR2_dilute')
    exp_paths = {'J:\LPproject\NPRN4\2019-04-03_14-21-04_V1_diffintensities'};
    shanks = {0};
elseif strcmpi(pop,'L5') && strcmpi(mani,'ChR2')
    exp_paths = {'J:\LPproject\NPRN5\2019-04-03_16-57-37_V1_diffintensities'};
    shanks = {0};
elseif strcmpi(pop,'L5') && strcmpi(mani,'halo')
    exp_paths = {'H:\LPproject\NPRH5\V1\2019-08-16_11-41-27_V1_Halo',...
        'J:\LPproject\NPRH18\V1\2019-11-12_14-36-49_V1_halo',...
        'J:\LPproject\NPRH19\V1\2019-11-13_13-07-41_V1_halo',...
        'J:\LPproject\NPRH21\V1\2019-12-13_17-37-53_Halo',...
        'J:\LPproject\NPRH23\V1\2019-12-19_12-09-54_Halo',...
        'J:\LPproject\NPRH24\V1\2019-12-18_12-16-55_Halo',...
        'J:\LPproject\NPRH27\V1\2020-01-28_13-01-16_Halo_REAL',...
        'J:\LPproject\NPRH37\V1\2020-06-19_17-08-03_Halo'};
    shanks = {0,0,0,0,0,0,0,0};
    lightconds = {2,1,1,1,1,1,1,1};
    lightcond = 1;
elseif strcmpi(pop,'L5') && strcmpi(mani,'GtACR')
    exp_paths = {'Z:\L5\stGtACR\NPRSC5\V1\2020-07-22_13-49-52_driftgrat_2LEDs',...
        'Z:\L5\stGtACR\NPRSC4\V1\2020-07-21_15-45-08_driftgrat_2LEDs',...
        'Z:\L5\stGtACR\NPRSC6\V1\2020-09-02_14-37-11_driftgrat_2LEDs',...
        'Z:\L5\stGtACR\NPRSC8\V1\2020-08-31_16-14-35_driftgrat_2LEDs'};
    shanks = {0,0,0,0};
    lightconds = {[1],[1],[2],[2]}; % in case of multiple light conds, which do you want to analyze(1 means first lightcond -excluding no light cond)
    lightcond = 1;
elseif strcmpi(pop,'L5') && strcmpi(mani,'gtacr_2LEDs')
    exp_paths = {'Z:\L5&L6\DFlp8\V1\2020-11-06_14-11-36_driftgrat_2LEDs'};
    shanks = {[0,1]};
    lightconds = {[2 1 3]}; 
    lightcond = 3;
elseif strcmpi(pop,'L6') && strcmpi(mani,'halo_2LEDs')
    exp_paths = {'Z:\L5&L6\MD9\V1pen1\2020-11-25_16-16-51_driftgrat_2LEDs',...
        'Z:\L5&L6\MD9\V1pen2\2020-11-26_11-20-32_driftgrat_2LEDs',...
        'Z:\L5&L6\MDCTRL1\V1\2020-12-19_18-32-08_halo_2LEDs'};
    shanks = {0,0,[0 1]};
    lightconds = {[2 1 3],[2 1 3],[2 1 3]}; 
    lightcond = 3;
elseif strcmpi(pop,'L5ET') 
    exp_paths = {'Z:\L5&L6\MD5\V1-pen1\2020-11-07_15-25-10_driftgrat_2LEDS_REALREAL',...
        'Z:\L5&L6\MD5\V1-pen2\2020-11-07_19-45-43_driftgrat_2LEDs',...
        'Z:\L5&L6\MD6\V1\2020-11-23_13-45-33_driftgrat_2LEDs',...
        'Z:\L5&L6\DFlp8\V1\2020-11-06_14-11-36_driftgrat_2LEDs',...
        'Z:\L5\fDIO-Halo\DFlp7\V1\2020-10-02_16-56-05_driftgrat_2LEDs',...
        'Z:\L5&L6\DM2\V1\2020-10-12_17-27-26_driftgrat_2LEDs',...
        'Z:\L5&L6\DM3\V1-pen1\2020-10-14_15-16-53_driftgrat_2LEDs',...
        'Z:\L5&L6\DM3\V1-pen2\2020-10-14_20-12-02_driftgrat_2LEDs'};
    shanks = {[0 1],[0 1],[0],[0 1],[0],[0 1], [0 1], [0 1]};
    lightconds = {2,2,2,2,2,1,1,1};
    lightcond = 1;
elseif strcmpi(pop,'L5') && strcmpi(mani,'dta')
    exp_paths = {'Z:\dTa\DDTA3\V1\2020-01-31_15-27-25_driftgrat',...
        'Z:\dTa\DDTA5\V1\2020-03-09_15-24-23_driftgrat_DTA',...
        'Z:\dTa\DDTA6\V1-pen1\2020-08-27_12-22-45_driftgrat',...
        'Z:\dTa\DDTA6\V1-pen2\2020-08-27_15-04-45_driftgrat'};
    shanks={0,0,[0;1],[0:1]};
    lightconds = {[],[],[],[]};
    lightcond = 1;
elseif strcmpi(pop,'L6') && strcmpi(mani,'GtACR')
    exp_paths = {'Z:\L6\stGtACR\MG2\V1\2020-08-07_13-59-20_driftgrat_2lightpwrs',...
        'Z:\L6\stGtACR\MG6\V1\2020-12-20_19-42-41_haloREAL'};
    shanks = {0,0};
    lightconds = {2,1};
    lightcond = 1;
elseif strcmpi(pop,'L5&6') && strcmpi(mani,'2LEDs')  % these were the L5-halo L6-gtacr experiments
    exp_paths = {'Z:\L5&L6\DM2\V1\2020-10-12_17-27-26_driftgrat_2LEDs',...
        'Z:\L5&L6\DM3\V1-pen1\2020-10-14_15-16-53_driftgrat_2LEDs',...
        'Z:\L5&L6\DM3\V1-pen2\2020-10-14_20-12-02_driftgrat_2LEDs'};
    shanks = {[0 1], [0 1], [0 1]};
    lightconds = {[1:3], [1:3], [1:3]}; 
    lightcond = 3;
elseif strcmpi(pop,'L5&6') && contains(mani,'2LEDs_flip')
    exp_paths = {'Z:\L5&L6\MD5\V1-pen1\2020-11-07_15-25-10_driftgrat_2LEDS_REALREAL',...
        'Z:\L5&L6\MD5\V1-pen2\2020-11-07_19-45-43_driftgrat_2LEDs',...
        'Z:\L5&L6\MD6\V1\2020-11-23_13-45-33_driftgrat_2LEDs',...
        'Z:\L5&L6\MD7\V1\2020-11-24_14-24-32_driftgrat_2LEDs',...
        'Z:\L5&L6\MD10\V1\2020-12-21_15-56-41_driftgrat_2LEDsREAL'};
    shanks = {[0 1], [0 1],[0 1], [0 1], [0 1]};
    lightconds = {[2 1 3], [2 1 3],[2 1 3], [2 1 3], [2 1 3]}; % in order: L5, L6, L5&6
    lightcond = 3;
elseif strcmpi(pop,'V1') && strcmpi(mani,'PVChR2')
    exp_paths = {'Z:\V1Inactivation\VH3\V1\2018-05-17_13-45-22_V1_PVChR2_take2'};
    shanks = {[0]};
    lightconds = {[2]};
    lightcond = 1;
elseif strcmpi(pop,'V1') && strcmpi(mani,'DlxChR2')
    exp_paths = {'Z:\V1Inactivation\VSC3\V1\2020-06-04_12-22-21_HaloChR2',...
        'Z:\V1Inactivation\VSC4\V1\2020-06-04_20-16-50_HaloChR2'};
    shanks = {0,0};
    lightconds = {1,1};
    lightcond = 1;
elseif strcmpi(pop,'V1') && strcmpi(mani,'AAV')
    exp_paths = {'Z:\V1Inactivation\VH3\V1\2018-05-17_13-45-22_V1_PVChR2_take2',...
        'Z:\V1Inactivation\VSC3\V1\2020-06-04_12-22-21_HaloChR2',...
        'Z:\V1Inactivation\VSC4\V1\2020-06-04_20-16-50_HaloChR2'};
    shanks = {0,0,0};
    lightconds = {2,1,1};
    lightcond = 1;
elseif strcmpi(pop,'V1') && strcmpi(mani,'Ai32')
    exp_paths = {'Z:\V1Inactivation\VSC5\V1\2020-07-28_14-55-27_driftgrat_2LEDs',...
        'Z:\V1Inactivation\VSC6\V1\2020-07-29_18-25-56_driftgrat_2LEDs',...
        'Z:\V1Inactivation\VES3\V1\2020-05-07_15-47-43_ChR2_2LEDs'};
    shanks = {[0],0,0};
    lightconds = {[1],[1],[1]};
    lightcond = 1;
end


if ~exist(main_dir,'dir')
    mkdir(main_dir)
end
cd(main_dir)

fig_dir =  sprintf('%s\\%s\\%s_%s_figs',main_dir,'InactFigures_V1',pop,mani);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

%% load experiment data
params = [];
unitinfo = [];
FRs = [];
tuning = [];
waveforms = [];
refr_idx = [];
cR = [];
uQ = [];
refV = [];
exp_num = [];

for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
    % get experiment name
    out = regexp(exp_path,'\\','split');
    an_name{i} = out{end-1};       % animal name
    if contains(an_name{i},'LP','ignorecase',1) || contains(an_name{i},'V1','ignorecase',1) || contains(an_name{i},'TRN','ignorecase',1)
        an_name{i} = strcat(out{end-2},'_',an_name{i});
    elseif contains(an_name{i},'day','ignorecase',1) 
        an_name{i} = out{end-2};
    end
    if ~isempty(strfind(out{end},'shank'))      % if exp_path is a shank subdirectory w/in experimental folder
        inds = regexp(out{end-1},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
        exp_name = strcat(out{end-1}(1:inds(1)-1),sprintf('_%s',out{end}));     % experiment name includes shank
    else
        inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
        exp_name = out{end}(1:inds(1)-1);
        if sum(isstrprop(exp_name,'digit'))>4       % this means it's an OpenEphys file (because it names files starting with the date)
            inds = regexp(out{end},'_\D','start');
            if strfind(out{end}(inds(1):end),an_name{i})
                exp_name = out{end}(inds(1)+1:end);
            else
                exp_name = strcat(an_name{i},out{end}(inds(1):end));
            end
        end
    end
    exp_dir = strcat(main_dir,'\',exp_name);
    if ~exist(exp_dir,'dir') && exist(sprintf('%s\\%s\\%s',main_dir,an_name{i},exp_name))       % need better solution for this
        exp_dir = sprintf('%s\\%s\\%s',main_dir,an_name{i},exp_name);
    end
    fprintf(sprintf('Loading experiment %s\n',exp_name))
    
    % load results
    cd(exp_dir)
    s = dir; 
    for ii=1:length(s)
        if strfind(s(ii).name,'_results') 
            results_file = s(ii).name;
        elseif strfind(s(ii).name,'cluster')
            cluster_file = s(ii).name;
        end
    end
    if ~exist('cluster_file','var')      % need better solution for this
        cluster_file = sprintf('%s\\%s\\%s_cluster_quality.mat',main_dir,an_name{i},an_name{i});
    end
    exp = importdata(results_file,'-mat');     % load results mat
    clust = importdata(cluster_file,'-mat');
    clear cluster_file results_file

    exp_num = [exp_num i*ones(1,length(exp.unitinfo))];
    params = [params exp.params];
    unitinfo = [unitinfo exp.unitinfo];
    FRs = [FRs exp.FRs];
    if isfield(exp.tuning,'SF')
        exp.tuning = rmfield(exp.tuning,'SF');      % remove SF and TF fields, if present (for now...just so I can concatenate results from prior experiments)
        exp.tuning = rmfield(exp.tuning,'TF');
    end
    tuning = [tuning exp.tuning];
    waveforms = [waveforms exp.waveforms];
    uQ = [uQ clust.uQ'];
    cR = [cR clust.cR'];
    refV = [refV clust.refV];
    
    load(fullfile(exp_path,'layers.mat'))
    CSDlayers{i} = layers;       % V1 specific - get layer info from CSD
    if exist('layers_shank','var')
        CSDshank{i} = layers_shank;
    else
        if length(layers)>64
            error(strcat('Need to rerun CSD for exp ',an_name))
        else
            CSDshank{i} = zeros(length(layers),1);
        end
    end
        
end

%% Identify experiment parameters for multi-experiment analyses

% first, set up analysis window (esp. important for halo experiments where light
% starts before vis stim)
for n = 1:length(params)
    if ~isempty(lightconds{1})      % opto experiments
        if round(max(params(n).av_light_start(lightconds{n})),2) < params(n).prestim % if LED started BEFORE vis stim
            startT(n) = 1500*params(n).prestim+1;        % manually set start time to prestim amount of time after vis stim start (e.g.startT = 751, 500ms after vis stim start)
        else
            startT(n) = round(max(params(n).av_light_start(lightconds{n}))*1000);
        end
        endT(n) = round(min(params(n).av_light_start(lightconds{n})+params(n).light_dur(lightconds{n}+1))*1000);
    else            % non-opto experiments
        startT(n) = 1500*params(n).prestim+1;        % manually set start time to prestim amount of time after vis stim start (e.g.startT = 1001, 500ms after vis stim start)
        endT(n) = 4500*params(n).prestim+1;  % manually set end time to 2250
    end
end
window = [max(startT) min(endT)]; % use common time window across experiments (in ms)
ev_lighttime = (diff(window)+1)/1000; % in sec

% identify different trial types
for n = 1:length(params)
    if strcmp(params(n).IVs,'s_freq')
        vis_trials{n} = find((params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1)&(params(n).trial_type(:,strcmpi(params(n).IVs,'s_freq'))==.04)&(params(n).trial_type(:,strcmpi(params(n).IVs,'t_period'))==30));  % in case of multiple SFs and TFs, only compare regular 2Hz and .04cpd trials (for now)
    else
        vis_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1);
    end
    blank_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==0);
    nolight_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==0);
    conds{n} = unique(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit')));  % the actual values assigned to light conditions 
    for lc=1:length(lightconds{n})  
        light_trials{n}{lc} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1));
        vislight_trials{n}{lc} = intersect(vis_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1)));
        blanklight_trials{n}{lc} = intersect(blank_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1)));
    end
    run_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==1);
    stat_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==0);
    fprintf('%s: %d total trials, %6.2f%% stationary \n',params(n).exp_name,size(params(n).trial_type,1),(length(stat_trials{n})/size(params(n).trial_type,1))*100)
end

%% identify which units should be included for analyses, pt. 1 (based on shank location and cluster quality)
incl_units = zeros(1,length(unitinfo)); 
if isfield(waveforms,'shank')
    shk = [waveforms.shank];
    if isempty(shk); shk = zeros(1,length(unitinfo)); end    % if single shank probe, shank field may be empty
else
    shk = zeros(1,length(unitinfo));
end

for n = 1:length(params)
    which_units = find(exp_num==n);
    incl_units(which_units(ismember(shk(exp_num==n),shanks{n}))) = 1;  % only include units from designated shank(s)
end
clean_units = find((uQ>16)&(refV<.5)&(incl_units));    % only include units with <.5% refractory period violations, unit quality > 16, and those from designated shanks 
incl_units = find(incl_units);      % necessary for later..use incl_units instead of clean_units to define LP boundaries

% also exclude low-FR cells (<.25Hz in blank and visual trials)
FRb = [FRs(incl_units).baseline];  %  baseline (from blank trials w/ no running, light or vis stim)
FRv = arrayfun(@(x) x.visual.ev(1,1),FRs(incl_units)); % visually-evoked FR from visual, stationary, no light trials
quiet_cells =max([FRv; FRb])<.25;
incl_units(quiet_cells) = [];
FRb(quiet_cells) = [];
FRv(quiet_cells) = [];

%% Calculate visual significance prior to getting distances from first
% and last ch and deciding which units are "clean". Use first and last
% visually significant channels to determine borders of LP

onsetT=.2;  

for i = 1:length(incl_units)      % for each unit
    nn = incl_units(i);
    tuning_curve{i} = tuning(nn).curve(:,:);
    oris = unique(params(exp_num(nn)).trial_type(:,strcmp(params(exp_num(nn)).IVs,'ori')));
        oris(oris>=999) = [];
    % wilcoxen rank-sum test to test for significant and visual- and light-modulation
    % for visual modulation, find preferred direction trials - only use THESE
    % trials to test for significant visual modulation (in case of extremely
    % tuned cells)
    [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-FRb(i)));   % direction w/ biggest change from baseline (defined from no-light condition)
    prefori_trials{i} = find(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori'))==oris(prefdir_deg(i)));
    if contains(mani,'dta','IgnoreCase',1) || strcmpi(vissigdef,'prepost') % if dta exp without blank trials, compare pre- and post-FRs (OR if i specifically indicate I wanna look at visual responsivity this way)
        vis_sig(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),window(1):window(1)+1000*params(exp_num(nn)).prestim-1),2),... % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
            sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1:1000*params(exp_num(nn)).prestim),2));
        vis_sig_ons(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2),...  % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
            sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1:onsetT*1000),2));
    else % if comparing visual to blank trials, use UNPAIRED parametric text (ie wilcoxen ranksum)
        vis_sig(i) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),... % significance of visual response (evoked periods of 1000ms in blank vs. visual trials) % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2));
        vis_sig_ons(i) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2),...  % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials) % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2));   % changed from params(exp_num(nn)).onset to .2s 6/8/20 because 100ms may be too fast for vis-evoked response in LP
    end
    
    for lc = 1:length(lightconds{exp_num(nn)})  % changed 1/7/21 to specifically use preferred ori trials and stationary trials for assessing sig light effect on visual responses
        light_sig(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),...       % NEW 6/25/20 - separately assessing light significant for visual vs. blank trials  % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2));% significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
        light_sig_bl(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),...       % NEW 6/25/20 - separately assessing light significant for visual vs. blank trials  % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2));% significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
        if params(exp_num(nn)).av_light_start(lc) < params(exp_num(nn)).prestim    % if light started prestim, you can use visual+blank trials to test significance of light-onset response
            light_sig_ons(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(nolight_trials{exp_num(nn)},stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2),...
            sum(unitinfo(nn).rast(intersect(light_trials{exp_num(nn)}{lc},stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2)); % significance of light-evoked response at light onset (200ms after light onset in nolight vs light trials  - changed from params(exp_num(nn)).onset to .2s 12/28/20 % Stationary trials ONLY (new addition 1/7/21)
        else % otherwise, just use blank trials
            light_sig_ons(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2),...
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2)); % significance of light-evoked response at light onset (200ms after light onset in nolight vs light trials  - changed from params(exp_num(nn)).onset to .2s 12/28/20 % Stationary trials ONLY (new addition 1/7/21)
        end
    end
%     if lightcond > 1   % otherwise, do kruskalwallis (because >2 groups of trial types)
%         rast{1} = sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
%         tri_group{1} = zeros(length(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)})),1)';
%         rast_bl{1} = sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
%         tri_group_bl{1} = zeros(length(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)})),1)';
%         for lc = 1:length(lightconds{exp_num(nn)})
%             rast{lc+1} = sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
%             tri_group{lc+1} = lc*ones(length(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)})),1)';
%             rast_bl{lc+1} = sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
%             tri_group_bl{lc+1} = lc*ones(length(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)})),1)';
%         end
%         [light_sig(i,lightcond),~,unit_stats(i)] = kruskalwallis([rast{:}],[tri_group{:}],'off');
%         mc(:,:,i) = multcompare(unit_stats(i),'ctype','dunn-sidak','display','off');
%         [light_sig_bl(i,lightcond),~,unit_stats_bl(i)] = kruskalwallis([rast_bl{:}],[tri_group_bl{:}],'off');
%         mc_bl(:,:,i) = multcompare(unit_stats(i),'ctype','dunn-sidak','display','off');
% %         light_sig_ons(i,1:length(lightconds{exp_num(nn)})) = nan(1,length(lightconds{exp_num(nn)})); % JUST A PLACEHOLDER  - not actually using for these experiments
%     end
%     
   % next, check tuning significance and get tuning curves
     if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
        % evaluate significance or orientation tuning
        for lc = 1:length(lightconds{exp_num(nn)}) + 1  % FIXED 1/5/21 - was excluding no-light trials!
            ori_trials = cell(1, length(oris)/2);
            for o = 1:length(oris)/2
                if lc == 1
                    ori_trials{o} = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==conds{exp_num(nn)}(1)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in no-light trials
                else
                    ori_trials{o} = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==conds{exp_num(nn)}(lightconds{exp_num(nn)}(lc-1)+1)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
                end
            end
            n_oritrials = cellfun(@(x) length(x),ori_trials);
            for o = 1:length(oris)/2
                if length(unique(n_oritrials)) > 1    % in cases of unequal # of trials per ori
                    ori_trials{o}(randperm(n_oritrials(o),n_oritrials(o)-min(n_oritrials))) = []; % if unequal # of trials per ori, randomly choose trial(s) to exclude
                end
                tuning_trials{exp_num(nn)}(:,o,lc) = sum(unitinfo(nn).rast(ori_trials{o},window(1):window(2)),2); % numbers of spikes across trials of given light condition for each orientation (by column)
            end
            tuned_sig(i,lc) = T2Hot1(tuning_trials{exp_num(nn)}(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
        end
    end
end
%% Get depth information (**this section is V1-specific: get layer and also look at depth relative to L4-L5 border)

% get layer info
layer = [unitinfo(incl_units).layer];

% now, use visual significance to approximate dorsal and ventral boundaries of LP
distfromlastch = nan(1,length(incl_units));
distfromfirstch = distfromlastch;
distfromgran = distfromlastch;

% get probe info 
p = eval('probemap_64D_bottom_func');   % assumes all V1 recordings were with either 64D or 128AN probes
Zchan = flipud(sort(p.z(p.shaft==1))); % from top to bottom (bottom=0)

count=1;
for n = 1:length(params)        % for each exp
    for sh = 1:length(shanks{n})    % for each shank in exp
        shk_units{count} = find((exp_num(incl_units)==n) & (shk(incl_units)==shanks{n}(sh)));
        allchs = sort([waveforms(incl_units(shk_units{count})).max_ch]);
        firstch(count) = min(allchs);
        granch(count) = find(CSDlayers{n}(CSDshank{n}==shanks{n}(sh))<=4 & CSDlayers{n}(CSDshank{n}==shanks{n}(sh))~=0,1,'last');
        % previously, using layer assignment outputs from unit analysis
        % (and thus can be inaccurate in experiments with few units
        % recorded)
%         if sum(layer(shk_units{count})<=4) % layer 4 and/or L2/3 detected
%             granch(count) = allchs(find(layer(shk_units{count}(inds))<=4,1,'last')); % gets the channel number (1:64)
%         else
%             granch(count) = allchs(find(sort(layer(shk_units{count}))>=5,1,'first'))-1;
%         end
        lastch(count) = max(allchs);
        distfromlastch(shk_units{count}) = Zchan(lastch(count))-Zchan([waveforms(incl_units(shk_units{count})).max_ch]);   % negative values are actually good (above last ch)
        distfromfirstch(shk_units{count}) = Zchan(firstch(count))-Zchan([waveforms(incl_units(shk_units{count})).max_ch]); % positive values are below first ch
        distfromgran(shk_units{count}) = Zchan([waveforms(incl_units(shk_units{count})).max_ch])-Zchan(granch(count));
        count = count+1;
    end
end

% adjust for new clean_units
clean_units = intersect(clean_units,incl_units);
distfromlastch = distfromlastch(ismember(incl_units,clean_units));          % only includes GOOD units
distfromfirstch = distfromfirstch(ismember(incl_units,clean_units));          % only includes GOOD units
distfromgran = distfromgran(ismember(incl_units,clean_units));
tuning_curve=tuning_curve(ismember(incl_units,clean_units));
prefdir_deg = prefdir_deg(ismember(incl_units,clean_units));
prefori_trials = prefori_trials(ismember(incl_units,clean_units));
vis_sig = vis_sig(ismember(incl_units,clean_units));
vis_sig_ons = vis_sig_ons(ismember(incl_units,clean_units));
if exist('light_sig','var')  % if opto
    light_sig = light_sig(ismember(incl_units,clean_units),:);
    light_sig_bl = light_sig_bl(ismember(incl_units,clean_units),:);
    light_sig_ons = light_sig_ons(ismember(incl_units,clean_units),:);
end
tuned_sig = tuned_sig(ismember(incl_units,clean_units),:);
FRb = FRb(ismember(incl_units,clean_units));
exp_num = exp_num(clean_units);
layer = layer(ismember(incl_units,clean_units));

%% test for light modulation

% first, need to load important FR, tuning and psth data
num_lcs = length(lightconds{1});  % currently, must choose to look at same number of lightconds from every experiment
FRev = nan(length(clean_units),num_lcs+1); 
FRearly = FRev;
FRlate = FRev;
FRonset = FRev;
FRbl = FRev;
FRvison = FRev;
OSI_CV = FRev;
OSI = FRev;
DSI_CV = FRev;
DSI = FRev;
psthV = nan(num_lcs+1,size(FRs(1).psthVisual,2),length(clean_units)); % #conds x #timebins x #units 
psthBl = psthV;
renormcurve = nan(num_lcs+1,size(tuning(1).curve,2),length(clean_units));
for n = 1:length(params)
    FRev(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.ev(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';   %vis-evoked
    FRearly(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evstart(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';   % early evoked period
    FRlate(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evlate(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % late evoked period
    FRonset(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evlightonset(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % vis-evoked light onset
    FRbl(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.blank.ev(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % blank, evoked period
    FRvison(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.onset(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % vis stim onset
    % also get orientation values for later use
    OSI_CV(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.OSI_CV(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    OSI(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.OSI(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    DSI_CV(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.DSI_CV(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    DSI(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.DSI(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    renormcurve(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.curve([1 lightconds{n}+1],:), tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(tuning(1).normcurve,2),sum(exp_num==n));
    % and get PSTHs for later use
    psthV(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.psthVisual([1 lightconds{n}+1],:), FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(FRs(1).psthVisual,2),sum(exp_num==n));
    psthBl(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.psthBlank([1 lightconds{n}+1],:), FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(FRs(1).psthBlank,2),sum(exp_num==n));
end

% recalculate tuning using normalized tuning curves (additively or
% multiplicatively scaled to account for baseline FR change due to light)
renormcurve2 = renormcurve;     % multiplicatively rather than additively scaling
for i=2:num_lcs+1 % for each lightcond
    renormcurve(i,:,:) = squeeze(renormcurve(i,:,:)) - repmat(diff(FRev(:,[1 i]),[],2)',size(renormcurve,2),1); % using FRev instead of FRbl since halo is more effective in blank than visual trials?
    scalef = FRev(:,1)./FRev(:,i);
    renormcurve2(i,:,:) = squeeze(renormcurve2(i,:,:)).* repmat(scalef',size(renormcurve2,2),1);
end
for n=1:length(clean_units)
    for lc = 1:num_lcs+1
        [newOSI(n,lc),newOSI_CV(n,lc),newDSI(n,lc),newDSI_CV(n,lc)] = calcOSIDSI(squeeze(renormcurve(lc,:,n)),oris');       % baseline-subtracted (baseline from each light condition separately)
%         [newOSI2(n,lc),newOSI_CV2(n,lc),newDSI2(n,lc),newDSI_CV2(n,lc)] = calcOSIDSI(squeeze(renormcurve2(lc,:,n)),oris');       % normalized to pref deg (i.e. multiplicative scaling) -  **same as orig OSI vals - OSI calcs not affected by scaling!
    end
end

% calculate light modulation indexes
for ii = 2:size(FRev,2)
    lightmod(:,ii-1) = (diff(FRev(:,[1 ii]),[],2))./sum(FRev(:,[1 ii]),2);
    lightmod_early(:,ii-1) = (diff(FRearly(:,[1 ii]),[],2))./sum(FRearly(:,[1 ii]),2);
    lightmod_late(:,ii-1) = (diff(FRlate(:,[1 ii]),[],2))./sum(FRlate(:,[1 ii]),2);
    lightmod_onset(:,ii-1) = (diff(FRonset(:,[1 ii]),[],2))./sum(FRonset(:,[1 ii]),2);
    lightmod_bl(:,ii-1) = (diff(FRbl(:,[1 ii]),[],2))./sum(FRbl(:,[1 ii]),2);
end

%% get preferred stimulus FR for each condition (but preferred stimulus defined in no light condition)
FRpref = nan(size(FRev));
for i = 1:length(clean_units)
    for ii = 1:num_lcs+1
        if ii == 1
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(i)}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime; % only stationary trials (as for other FRs)
        else
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(i)}{ii-1}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;  % only stationary trials (as for other FRs)
        end
    end
end

% and visual modulation (changed 1/28/21) - calculated from preferred vis
% trials (but onset is still from ALL vis trials)
vismod = (FRpref(:,1) - FRb')./(FRpref(:,1) + FRb');        % using baseline (from blank trials w/ no running, light or vis stim)
vismod_on = (FRvison(:,1)-FRb')./(FRvison(:,1)+FRb');

%% get running modulation index (just for V1)
FRev_run = nan(size(FRev));
FRev_stat = FRev_run;
FRbl_run = FRev_run;
FRbl_stat = FRev_run;
for i=1:length(clean_units)
    if length(run_trials{exp_num(i)})/size(unitinfo(clean_units(i)).rast,1) >=.05 % if at least 5% of trials were running trials
        FRev_run(i,1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(vis_trials{exp_num(i)},nolight_trials{exp_num(i)}),run_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
        FRev_stat(i,1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(vis_trials{exp_num(i)},nolight_trials{exp_num(i)}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
        FRbl_run(i,1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},nolight_trials{exp_num(i)}),run_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
        FRbl_stat(i,1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},nolight_trials{exp_num(i)}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
        for ii = 1:num_lcs
            FRev_run(i,ii+1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(vis_trials{exp_num(i)},light_trials{exp_num(i)}{ii}),run_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
            FRev_stat(i,ii+1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(vis_trials{exp_num(i)},light_trials{exp_num(i)}{ii}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
            FRbl_run(i,ii+1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},light_trials{exp_num(i)}{ii}),run_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
            FRbl_stat(i,ii+1) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},light_trials{exp_num(i)}{ii}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;
        end
    else
        FRev_run(i,:) = nan(1,num_lcs+1);
        FRev_stat(i,:) = nan(1,num_lcs+1);
        FRbl_run(i,:) = nan(1,num_lcs+1);
        FRbl_stat(i,:) = nan(1,num_lcs+1);
    end
end
        
% calculate "running modulation index" (FRrun-FRstat)/(FRrun+FRstat)
runmod = reshape((FRev_run(:)-FRev_stat(:))./(FRev_run(:)+FRev_stat(:)),length(clean_units),num_lcs+1); % checked - this reshaping avoids another for loop!
runmod_bl = reshape((FRbl_run(:)-FRbl_stat(:))./(FRbl_run(:)+FRbl_stat(:)),length(clean_units),num_lcs+1);


%% get waveform props
for i = 1:length(clean_units)
    [t2p_t(i),t2p_r(i),fwhm(i)] = get_waveform_props(waveforms(clean_units(i)).microV,params(exp_num(i)).amp_sr);
end

reg_cells = find(t2p_t>=.475);
FS_cells = find(t2p_t<.475);

% plot trough-to-peak times
figure; 
hold on;
plot(t2p_t,t2p_r,'o')
xlabel('trough-to-peak time (ms)','Fontsize',18)
ylabel('trough-to-peak-ratio','Fontsize',18)
yax=get(gca,'ylim');
line([.475 .475],yax,'linestyle','--','color','k')
print(gcf, '-dpng','FSvsRS')
print(gcf, '-painters','-depsc','FSvsRS')

%% Figure time! (**this section is V1-specific)

cd(fig_dir)
% get different cell types
visual_cells = find((vis_sig < .025)|(vis_sig_ons < .025));
nonvisual_cells = find(~ismember(1:size(vis_sig,2),visual_cells));
tuned_cells = intersect(find(tuned_sig(:,1) < .05),visual_cells);
if contains(pop,'L6','ignorecase',1)        % L6 == green
    mani_layer = 6;
    IG_layer = 5;       % other infragranular layer
elseif contains(pop,'L5&6','ignorecase',1)
    mani_layer = [5 6];
    IG_layer = [5 6];
elseif contains(pop,'L5','ignorecase',1) 
    mani_layer = 5;
    IG_layer = 6;       % other infragranular layer
elseif contains(pop,'V1','ignorecase',1)
    mani_layer = unique(floor(layer));
else
    mani_layer = 5; % TEMP - for now, focus on L5 for V1 inactivation experiments
end
RSmani_cells = reg_cells(ismember(floor(layer(reg_cells)'),mani_layer));
if exist('light_sig','var')
    [~, pthresh, ~, ~] = fdr_bh(light_sig(:,lightcond),.1); % 10% FDR
    pthresh = max(pthresh,.01); % don't allow pthreshBl to be <.01
    fprintf('new corrected pvalue for 10percent FDR - blank lightsig: %6.4f \n', pthresh)
    supp_cells = reg_cells(find(ismember(floor(layer(reg_cells)'),mani_layer) & lightmod(reg_cells,lightcond)<-.33 & light_sig(reg_cells,lightcond)<=pthresh)); 
    other_cells = RSmani_cells(~ismember(RSmani_cells,supp_cells));
else
    supp_cells = [];
    other_cells = [];
end

%% set up coloring
% if contains(mani,'halo','ignorecase',1)
%     color_mat = [0 0 0; .9 0 .3; 0.6350, 0.0780, 0.1840]; % for graphing purposes (first is black, last is green)
%     nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];      % make red for halo
if contains(pop,'L6','ignorecase',1)        % L6 == green
    color_mat = [.166 .674 .188; .166 .674 .188;.166 .674 .188];
    nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];
elseif contains(pop,'L5&6','ignorecase',1)
    color_mat = [0 .2 .9; .166 .674 .188; .083 .567 .600]; % black, blue, green, teal
elseif contains(pop,'L5','ignorecase',1)    % L5 == blue
    color_mat = [0 .2 .9; 0 .2 .9; 0 .2 .9; 0 .2 .9];
    nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];
elseif contains(pop,'V1','ignorecase',1)
    color_mat = [.083 .567 .6];     % teal
else
    color_mat = [0 .8 1; 0 0 1; 0 0.5 .4]; % % for lighter-shade dots
    nonvis_color_mat = [.5 .5 .5; .75 .8 1; .7 .8 .7];  % for lighter-shade dots
end

%% barplots - by layer

types = zeros(1,size(FRev,1));
if exist('light_sig','var')
types(supp_cells) = 1;
    
    boxplot_by_layer(lightmod(reg_cells,:)',layer(reg_cells),types(reg_cells),'Light modulation index','Lightmodulation (RS evoked - all conditions)',color_mat)
    boxplot_by_layer(lightmod_bl(reg_cells,:)',layer(reg_cells),types(reg_cells),'Light modulation index','Lightmodulation_blanks (RS evoked - all conditions)',color_mat)
    boxplot_by_layer(lightmod(FS_cells,:)',layer(FS_cells),types(FS_cells),'Light modulation index','Lightmodulation (FS evoked - all conditions)',color_mat)
    boxplot_by_layer(lightmod_bl(FS_cells,:)',layer(FS_cells),types(FS_cells),'Light modulation index','Lightmodulation_blanks (FS evoked - all conditions)',color_mat)
end

boxplot_by_layer(runmod(reg_cells,1)',layer(reg_cells),types(reg_cells),'Runmod','Runmod (RS evoked - all conditions)',color_mat)
boxplot_by_layer(abs(runmod(reg_cells,1))',layer(reg_cells),types(reg_cells),'abs(Runmod)','Runmod(abs) (RS evoked - all conditions)',color_mat)
boxplot_by_layer(vismod(reg_cells,1)',layer(reg_cells),types(reg_cells),'Vismod','Vismod (RS evoked - all conditions)',color_mat)
boxplot_by_layer(abs(vismod(reg_cells,1))',layer(reg_cells),types(reg_cells),'abs(Vismod)','Vismod(abs) (RS evoked - all conditions)',color_mat)
boxplot_by_layer(OSI_CV(reg_cells,1)',layer(reg_cells),types(reg_cells),'OSI-CV','OSI-CV (RS evoked - all conditions)',color_mat)
boxplot_by_layer(DSI_CV(reg_cells,1)',layer(reg_cells),types(reg_cells),'DSI-CV','DSI-CV (RS evoked - all conditions)',color_mat)
boxplot_by_layer(FRb(reg_cells),layer(reg_cells),types(reg_cells),'Spont. FR','Spont FR (RS evoked - all conditions)',color_mat)

boxplot_by_layer(runmod(FS_cells,1)',layer(FS_cells),types(FS_cells),'Runmod','Runmod (FS evoked - all conditions)',color_mat)
boxplot_by_layer(abs(runmod(FS_cells,1))',layer(FS_cells),types(FS_cells),'abs(Runmod)','Runmod(abs) (FS evoked - all conditions)',color_mat)
boxplot_by_layer(vismod(FS_cells,1)',layer(FS_cells),types(FS_cells),'Vismod','Vismod (FS evoked - all conditions)',color_mat)
boxplot_by_layer(abs(vismod(FS_cells,1))',layer(FS_cells),types(FS_cells),'abs(Vismod)','Vismod(abs) (FS evoked - all conditions)',color_mat)
boxplot_by_layer(OSI_CV(FS_cells,1)',layer(FS_cells),types(FS_cells),'OSI-CV','OSI-CV (FS evoked - all conditions)',color_mat)
boxplot_by_layer(DSI_CV(FS_cells,1)',layer(FS_cells),types(FS_cells),'DSI-CV','DSI-CV (FS evoked - all conditions)',color_mat)
boxplot_by_layer(FRb(FS_cells),layer(FS_cells),types(FS_cells),'Spont. FR','Spont FR (FS evoked - all conditions)',color_mat)

%% plot distance from bottom of LP by lightmod
if exist('light_sig','var')
    figure;
    subplot(111)
    plot(lightmod(:,lightcond),distfromlastch','.','MarkerSize',24,'color',color_mat(lightcond,:))
    % hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
    h = get(gca,'ytick');
    set(gca,'yticklabel',h*25);
    xlim([-1 1])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--')
    % legend('low','high')
    xlabel('LMI_v_i_s','Fontsize',16)
    ylabel(strcat('Distance from bottom of LP (µm)'),'Fontsize',16)
    print(gcf, '-dpng','lightmodbydepth_bottom')
    print(gcf,'-painters','-depsc','lightmodbydepth_bottom')

    figure;
    subplot(111)
    plot(lightmod(:,lightcond),abs(distfromfirstch)','.','MarkerSize',24,'color',color_mat(lightcond,:))
    % hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
    h = get(gca,'ytick');
    % set(gca,'yticklabel',h*25);
    view(0,270)
    xlim([-1 1])
    ylim([0 max(abs(distfromfirstch))])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    % legend('low','high')
    xlabel('LMI_v_i_s','Fontsize',24)
    ylabel(strcat('Depth (µm)'),'Fontsize',24)
    set(gca,'fontsize',18,'linewidth',2)
    print(gcf, '-dpng','lightmodbydepth_top')
    print(gcf,'-painters','-depsc','lightmodbydepth_top')

    figure;
    subplot(111)
    plot(lightmod(:,lightcond),distfromgran','.','MarkerSize',24,'color',color_mat(lightcond,:))
    h = get(gca,'ytick');
    % set(gca,'yticklabel',h*25);
    xlim([-1 1])
    % ylim([min(distfromgran),max(distfromgran)])
    ylim([-700 525])        % currently manually set so that they're identical for L5 and L6
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    % legend('low','high')
    xlabel('LMI_v_i_s','Fontsize',24)
    ylabel(strcat('Depth from L4-5 border (µm)'),'Fontsize',18)
    set(gca,'fontsize',18,'linewidth',2)
    print(gcf, '-dpng','lightmodbydepth_gran')
    print(gcf,'-painters','-depsc','lightmodbydepth_gran')
end

if contains(mani,'2LEDs','ignorecase',1)
    LED_titles = {'L5&L6CTs inactivated','L5CTs inactivated','L6CTs inactivated'};
    plotord = [3 1 2];
    figure; 
    for ii = 1:size(lightmod,2)
        i=plotord(ii);
        subplot(1,size(lightmod,2),ii)
        scatter(lightmod(:,i),distfromgran',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0],'markerfacealpha',1)
        hold on;
%         scatter(lightmod(intersect(supp_cells,find(floor(layer)==5)),i),distfromgran(intersect(supp_cells,find(floor(layer)==5)))',75,'filled','markerfacecolor',color_mat(1,:),'markerfacealpha',1)  % color by layer assignment
%         scatter(lightmod(intersect(supp_cells,find(floor(layer)==6)),i),distfromgran(intersect(supp_cells,find(floor(layer)==6)))',75,'filled','markerfacecolor',color_mat(2,:),'markerfacealpha',1)
        scatter(lightmod(supp_cells,i),distfromgran(supp_cells)',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',1)  % color by lightcond
        h = get(gca,'ytick');
        % set(gca,'yticklabel',h);
        xlim([-1 1])
        ylim([min(distfromgran),max(distfromgran)])
        yax = get(gca,'YLim');
        line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
        % legend('low','high')
        xlabel('LMI_v_i_s','Fontsize',24)
        ylabel(strcat('Depth from L4-5 border (µm)'),'Fontsize',18)
        title(LED_titles{ii})
        set(gca,'fontsize',18,'linewidth',2)
    end
    set(gcf,'Paperposition',[0 0 20 6])
    print(gcf, '-dpng','lightmodbydepth_fromgran_bytype')
    print(gcf,'-painters','-depsc','lightmodbydepth_fromgran_bytype')
end

%% scatter plots

for i = 1:num_lcs
    % FR light vs no light - visual vs nonvisual units, visual trials
    plot_scatter(FRev(:,[1 i+1]), types, {[.85 .85 .85],color_mat(i,:)}, [1 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis',i), {'Other','Inact cells'}, 2)     

    % FR light vs no light - visual vs nonvisual units, blank trials
    plot_scatter(FRbl(:,[1 i+1]), types, {[.85 .85 .85],color_mat(i,:)}, [1 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis_bl',i), {'Other','Inact cells'}, 2)    

    % FR light vs no light - visual vs nonvisual units, preferred trials
    plot_scatter(FRpref(:,[1 i+1]), types,  {[.85 .85 .85],color_mat(i,:)}, [1 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis_pref',i), {'Other','Inact cells'}, 2)    

    % change in evoked FR light vs no light - visual vs nonvisual units
    plot_scatter(FRev(:,[1 i+1])-FRbl(:,[1 i+1]), types,  {[.85 .85 .85],color_mat(i,:)}, [1 1], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON)', sprintf('%s_%d','FR_vischange',i), {'Other','Inact cells'}, 2)     

    % compare visual vs. blank lightmod 
    plot_scatter([lightmod(supp_cells,i) lightmod_bl(supp_cells,i)], ones(1,length(supp_cells)), {color_mat(i,:)}, 1, 'LMI_v_i_s', 'LMI_s_p_o_n_t', sprintf('%s_%d','Lightmod_visVSbl',i), {'Inact cells'}, 1)    

%     % compare effects on other cells within and outside of inactivated
%     % layer
    if strcmpi(pop,'V1')
        plot_scatter(FRev(:,[1,i+1]),t2p_t<.475,{[.5 .5 .5],color_mat(i,:)},[1 1],'FRev (spks/s - light OFF)','FRev (spks/s - light ON)',sprintf('%s_%d','FR_FScells',i), {'RS','FS'}, 2)
    else
        plot_scatter(FRev(floor(layer)==mani_layer&t2p_t>=.475,[1 i+1]),types(floor(layer)==mani_layer&t2p_t>=.475), {[.5 .5 .5],color_mat(i,:)}, [1 1], 'Vis-evoked FR (Spks/s-light OFF)', 'Vis-evoked FR (Spks/s-light ON)', sprintf('%s_%d','FR_supplayer',i), {'Other','Inact cells'}, 2)
%         plot_scatter(FRev(floor(layer)~=mani_layer&t2p_t>=.475,[1 i+1]),types(floor(layer)~=mani_layer&t2p_t>=.475), {[.5 .5 .5],color_mat(i,:)}, [1 1], 'Vis-evoked FR (Spks/s-light OFF)', 'Vis-evoked FR (Spks/s-light ON)', sprintf('%s_%d','FR_otherlayers',i), {}, 1)
        plot_scatter(FRev(floor(layer)~=mani_layer&t2p_t>=.475,[1 i+1]),floor(layer(floor(layer)~=mani_layer&t2p_t>=.475))~=IG_layer, {[.75 .75 .75],[.25 .25 .25]}, [1 1], 'Vis-evoked FR (Spks/s-light OFF)', 'Vis-evoked FR (Spks/s-light ON)', sprintf('%s_%d','FR_otherlayers',i), {strcat('L',num2str(IG_layer)),'Other layers'}, 2)
        plot_scatter(FRev(FS_cells,[1 i+1]),floor(layer(FS_cells))==mani_layer, {[.5 .5 .5],color_mat(i,:)}, [1 1], 'Vis-evoked FR (Spks/s-light OFF)', 'Vis-evoked FR (Spks/s-light ON)', sprintf('%s_%d','FR_FScells',i), {'Other-layer FS','Within-layer FS'}, 2)
    end
    
    % change in tuning (for visual cells only)
    plot_scatter(OSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1], 'OSI(CV) (light OFF)', 'OSI(CV) (light ON)', sprintf('%s_%d','OSI(CV)_bylightmod',i), {'Other','Inact cells'}, 2)    
    plot_scatter(OSI(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1],'OSI (light OFF)', 'OSI (light ON)', sprintf('%s_%d','OSI_bylightmod',i), {'Other','Inact cells'}, 2)     
    plot_scatter(DSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1],'DSI(CV) (light OFF)', 'DSI(CV) (light ON)', sprintf('%s_%d','DSI(CV)_bylightmod',i), {'Other','Inact cells'}, 2)     
    plot_scatter(DSI(visual_cells,[1 i+1]), types(visual_cells),{[.85 .85 .85],color_mat(i,:)}, [1 1],'DSI (light OFF)', 'DSI (light ON)', sprintf('%s_%d','DSI_bylightmod',i), {'Other','Light-suppressed'}, 2)     
     plot_scatter(newOSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1], 'OSI(CV) (light OFF)', 'OSI(CV) (light ON)', sprintf('%s_%d','OSI(CV)_bylightmod_baselinenorm',i), {'Other','Inact cells'}, 2)     % first lightcond pwr
    plot_scatter(newOSI(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1],'OSI (light OFF)', 'OSI (light ON)', sprintf('%s_%d','OSI_bylightmod_baselinenorm',i), {'Other','Light-activated','Inact cells'}, 2)     % first lightcond pwr
    plot_scatter(newDSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1],'DSI(CV) (light OFF)', 'DSI(CV) (light ON)', sprintf('%s_%d','DSI(CV)_bylightmod_baselinenorm',i), {'Other','Inact cells'}, 2)     % first lightcond pwr
    plot_scatter(newDSI(visual_cells,[1 i+1]), types(visual_cells), {[.85 .85 .85],color_mat(i,:)}, [1 1] ,'DSI (light OFF)', 'DSI (light ON)', sprintf('%s_%d','DSI_bylightmod_baselinenorm',i), {'Other','Inact cells'}, 2)     % first lightcond pwr

    close all
end

%%  make population-level average PSTHs (normalized to prestim baseline FR)
% thus, baseline = 1 
binsize = .025;  %% currently hardcoded!!! 25ms bins

% determine baseline FRs (for visual and blank trials separately)
mean_bs = repmat(mean(psthV(:,1:200/(binsize*1000),:),2),1,size(psthV,2),1);    % define baseline from first 200ms
mean_bs_bl = repmat(mean(psthBl(:,1:200/(binsize*1000),:),2),1,size(psthBl,2),1);
mean_bs(mean_bs==0) = nan;     % because otherwise could get infinity when normalizing by baseline
mean_bs_bl(mean_bs_bl==0) = nan;

norm_psth = psthV./mean_bs; % normalizes to prestim baseline, per condition (so that prestim=1)
norm_psth_bl = psthBl./mean_bs_bl;

norm_mean = nanmean(norm_psth,3); % average across units
norm_se = nanstd(norm_psth,0,3)./sqrt(size(norm_psth,3));
norm_mean_bl = nanmean(norm_psth_bl,3);
norm_se_bl = nanstd(norm_psth_bl,0,3)./sqrt(size(norm_psth_bl,3));
% and separately for supp and enh cells
norm_supp_mean = nanmean(norm_psth(:,:,supp_cells),3); 
norm_supp_se = nanstd(norm_psth(:,:,supp_cells),0,3)./sqrt(length(supp_cells));
norm_enh_mean = nanmean(norm_psth(:,:,other_cells),3);  % **this is only difference from LP script- "enh-cells" are actually non-supp cells
norm_enh_se = nanstd(norm_psth(:,:,other_cells),0,3)./sqrt(length(other_cells));

% set up parameters for plotting
for n = 1:length(params)
    tri_st(n) = params(n).prestim;
    tri_visdur(n) = params(n).stimtime;
    if exist('light_sig','var')
        start_times(n,:) = round(params(n).av_light_start(lightconds{n}),2)-params(n).prestim;  % find LED start times in diff exps, and round to nearest hundreth
        stim_durs(n,:) = round(params(n).light_dur(lightconds{n}+1),1);
    else
        start_times(n,:) = nan;
    end
end

if length(unique(tri_st))>1 || length(unique(tri_visdur))>1
    error('Disparate trial times across experiments')
else
    tri_st = tri_st(1);
    tri_visdur = tri_visdur(1);
end

if contains(mani,'halo','ignorecase',1)
    patch_color = [.9 .1 .1];  % red/orange patch for halo
else
    patch_color = [0 .1 1];  % default patch color is blue
end

% VISUAL trials
pop_fig= figure;
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean(1,:), norm_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean(i+1,:), norm_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_norm';
print(pop_fig,'-dpng',save_fig_name)
print(pop_fig,'-painters','-depsc',save_fig_name)

% BLANK trials
pop_fig2= figure;
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean_bl(1,:), norm_se_bl(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean_bl(i+1,:), norm_se_bl(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_blank_norm';
print(pop_fig2,'-dpng',save_fig_name)
print(pop_fig2,'-painters','-depsc',save_fig_name)

% VISUAL trials - by type
pop_fig3= figure;
subplot(121)
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_supp_mean(1,:), norm_supp_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_supp_mean(i+1,:), norm_supp_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
title(sprintf('Suppressed RS L%s units (n=%d)',num2str(mani_layer),length(supp_cells)),'fontsize',14);
set(gca,'fontsize',16,'linewidth',2);

subplot(122)
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_enh_mean(1,:), norm_enh_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_enh_mean(i+1,:), norm_enh_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
title(sprintf('Other RS L%s units (n=%d)',num2str(mani_layer),length(other_cells)),'fontsize',14);
set(gca,'fontsize',16,'linewidth',2);

set(gcf, 'Position', [100, 100, 1000, 420])
save_fig_name = 'PopulationPSTH_bytype_norm';
print(pop_fig3,'-dpng',save_fig_name)
print(pop_fig3,'-painters','-depsc',save_fig_name)

%% pie graph of number of suppressed (putative opto-tagged) vs. other RSunits
if exist('light_sig','var')
    figure;
    pie_data = [length(supp_cells) length(other_cells)];
    piegraph = pie(pie_data);
    piegraph_labels = {'Inact','Other'};
    piegraph_labels = piegraph_labels(pie_data>0);
    colormap([color_mat(lightcond,:); color_mat(lightcond,:); 1 1 1])
    % set(piegraph([1:2:end]),'edgecolor',color_mat(1,:))
    set(piegraph([2:2:end]),'fontsize',16)
    for i=2:2:length(piegraph)
        set(piegraph(i),'string',sprintf('%s (%s)',piegraph_labels{i/2},piegraph(i).String))
    end
    title(sprintf('%d RS L%s units',length(RSmani_cells),num2str(mani_layer)),'fontsize',16)
    print(gcf,'-dpng','Piegraph_lighteffects')
    print(gcf,'-painters','-depsc','Piegraph_lighteffects')
end

%% save stuff

num_units = length(clean_units);
fileID = fopen('results.txt','w');
fprintf(fileID,'Number of visually responsive cells: %d of %d\r\n',length(visual_cells),num_units);
fprintf(fileID,'Percent visually responsive: %.2f\r\n',100*length(visual_cells)/num_units);
fprintf(fileID,'Number of tuned cells: %d of %d\r\n',length(tuned_cells),num_units);
fprintf(fileID,'Percent significantly tuned: %.2f\r\n', 100*length(tuned_cells)/num_units);
fprintf(fileID,'Number of regular-spiking cells: %d of %d\r\n',length(reg_cells),num_units);
fprintf(fileID,'Percent regular-spiking: %.2f\r\n', 100*length(reg_cells)/num_units);
fprintf(fileID,'Number of fast-spiking cells: %d of %d\r\n',length(FS_cells),num_units);
fprintf(fileID,'Percent fast-spiking: %.2f\r\n', 100*length(FS_cells)/num_units);
fclose(fileID);

% save some more stuff for L5 vs L6 experiment comparisons
good_units = [unitinfo(clean_units).name];
exps = {params.exp_name};
if exist('light_sig','var')
    save('results.mat','lightmod','lightmod_bl','OSI_CV','DSI_CV','runmod','vismod','FRev','FRbl','supp_cells','other_cells','good_units','exps','exp_num','layer')
else
    save('results.mat','OSI_CV','DSI_CV','runmod','vismod','FRev','FRb','t2p_t','good_units','exps','exp_num','layer')
end

end


function plot_scatter(data, color_var, colors, alpha, xlab, ylab, title, leg, lobf)
% if lobf = 1, make one lobf regardless of variables; if lobf = 2, make
% separate lobfs for each variable
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
    scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75, 'filled','markerfaceColor', colors{i},'markerfacealpha',alpha(i));
end
set(gca,'Fontsize',18,'linewidth',2)
xlabel(xlab,'Fontsize',18)
ylabel(ylab,'Fontsize',18)
xax = get(gca,'XLim');
yax = get(gca,'YLim');

if prod(xax)<0 && abs(prod(xax))<1 && prod(yax)<0 && abs(prod(yax))<1
    xlim([-1 1])
    ylim([-1 1])
    xax = get(gca,'XLim');
    yax = get(gca,'YLim');
else
    xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
    ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
end
x = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)),100);
y=x;
plot(x,y,'k--','color','k');
line([min(xax(1),yax(1)) max(xax(2),yax(2))],[0 0],'color','k')
line([0 0], [min(xax(1),yax(1)) max(xax(2),yax(2))],'color','k')
if lobf == 1
    fittedX = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)), 200);
    coeffs = polyfitZero(data(:,1), data(:,2), 1);
    fittedY = polyval([0 coeffs], fittedX);
    plot(fittedX,fittedY,'color', colors{1},'linewidth',2)
elseif lobf == 2
    fittedX = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)), 200);
    for i = 1:length(vars)
        if sum(color_var==vars(i)) > 1    % comment out if you don't want to separately calculate lobfs
            coeffs(i,:) = polyfitZero(data(color_var==vars(i),1), data(color_var==vars(i),2), 1);     
            fittedY(i,:) = polyval([0 coeffs(i,:)], fittedX);
            plot(fittedX,fittedY(i,:),'color', colors{i},'linewidth',2)
%             scatter(fittedX,fittedY(i,:),20,'o','filled','markerfacecolor',colors{i},'markerfacealpha',alpha(i))
        end
    end

else    % plot median
    for i = 1:length(vars)
        plot(median(data(color_var==vars(i),1)),median(data(color_var==vars(i),2)), 'marker','+', 'Color', colors{i},'markersize',28,'linewidth',4);
    end
end
if ~isempty(leg)
    l=legend(leg,'location','best');
    set(l,'fontsize',12)
end
print(f, '-dpng',title)
print(f,'-painters','-depsc',title)
end
