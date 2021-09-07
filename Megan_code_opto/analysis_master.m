function analysis_master(exp_path,exp_system,exp_type,plots)
% Master code for loading raw data into matlab, extracting experimental
% conditions and trial types, calling 'unit_analysis.m' and plotting
% figures for each unit (optional), and saving all units into a single data
% structure. Compatible with either Intan or OpenEphys data.
% MAK 8/2/17 (modified from intan_analysis_master.m)

% Inputs:
% exp_path = path to experiment folder
    % (e.g.'J:\LPproject\D4\2017-07-07_14-43-31_LP_diffintensities')
% exp_system = 'Intan' or 'OpenEphys' (case insensitive)
% exp_type = 'ramp', 'trains', 'intensities' or 'size' 
% plots = 1 if you want figures, otherwise 0

% Output: *_results.mat file saved in exp_path


%% get experiment name and load raw data
% get name
out = regexp(exp_path,'\\','split');
an_name = out{end-1};       % animal name
if contains(an_name,'LP','ignorecase',1) || contains(an_name,'V1','ignorecase',1) || contains(an_name,'TRN','ignorecase',1) || contains(an_name,'SC','ignorecase',1)
    an_name = strcat(out{end-2},'_',an_name);
elseif contains(an_name,'day','ignorecase',1) 
    an_name = out{end-2};
end
if ~isempty(strfind(out{end},'shank'))      % if exp_path is a shank subdirectory w/in experimental folder
    inds = regexp(out{end-1},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
    exp_name = strcat(out{end-1}(1:inds(1)-1),sprintf('_%s',out{end}));     % experiment name includes shank
else
    inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
    exp_name = out{end}(1:inds(1)-1);
    if sum(isstrprop(exp_name,'digit'))>4       % this means it's an OpenEphys file (because it names files starting with the date)
        inds = regexp(out{end},'_\D','start');
        if strfind(out{end}(inds(1):end),an_name)
            exp_name = out{end}(inds(1)+1:end);
        else
            exp_name = strcat(an_name,out{end}(inds(1):end));
        end
    end
end

% load data
disp(strcat('Loading raw data from experiment ',exp_name))
cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      
else
    if strcmpi(exp_system,'intan')
        intan2matlab(exp_path);      
    elseif strcmpi(exp_system,'openephys')
        openEphys2matlab(exp_path);
    else
        error('Unrecognized acquisition system')
    end
    load(sprintf('%s/data.mat',exp_path))      % data from intan2matlab.m
end

% data.mat includes the following variables:
% trials - num_trials x 2 matrix of trial start and end times (in seconds)
% field_trials - trial num x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% time_index - vector of time stamps of each sample (in seconds), with 1000
    % samples per second - starts from 0!
% amp_sr - original amplifier sampling rate
% epoc - vector of 1s and 0s indicating whether each sample was part of a 
    % trial (1) or not (0). Still in original sampling rate (amp_sr)
% mmvt -  vector of 1s and 0s indicating whether mouse was running (1) or  
    % not (0) at time of sample. Still in original sampling rate (amp_sr)
% encd - vector of analog output from movement encoder. Still in original 
    % sampling rate (amp_sr)
% light - is 0 if not an optogenetics experiment, but otherwise, vectors of
    % whether LED was on (1) or off (0)
% re - vector for the timestamps (rising phase of photodiode signal) in 
    % sec. you need to remove the first 2 and the last 2 timestamps. 
    % Timestamps signify the end of a four second black to gray period.
% photo - vector of analog output from photodiode. In original sampling rate

% get number of channels
if strcmpi(exp_system, 'intan')
    [amplifier_channels,~,~,...
    ~,~,~,...
    ~] = read_Intan_RHD2000_file;
    chansInDat = size(amplifier_channels,2);
    if strfind(exp_path,'shank')
        chansInDat = chansInDat/4;      % cheat for LP experiments in which I separately analyzed 4 shanks of 128DN probe
    end
elseif strcmpi(exp_system,'openephys')
    filenames = dir;
    filenames = {filenames.name};
    chansInDat = length(cell2mat(strfind(filenames,'CH')));       % number of files with 'CH' in file name (i.e. continuous channel files)
end
    
% set up where to store results
% if contains(exp_name,'M') || contains(exp_name,'D') || contains(exp_name,'VH')
    full_dir = 'H:\LPproject\LPresults';        % where you want to save data
% elseif contains(exp_name,'TH')
%     full_dir = 'H:\Tlx3project\Haloresults';
% else
%     full_dir = 'H:\Tlx3project\ChR2results';
% end

main_dir = strcat(full_dir, '\',an_name);
exp_dir = strcat(main_dir,'\',exp_name);
if ~exist(exp_dir,'dir')
    mkdir(exp_dir);     % will save everything in experiment subfolder
end
all_fig_dir = sprintf('%s\\Figures',exp_dir);
if plots
    if exist(all_fig_dir,'dir')
        delete(sprintf('%s/*',all_fig_dir));        % delete and remake all figures
    end
end

%% load clustering info
disp('Loading clustering data...')
% if clustering done with kilosort
if exist('spike_times.npy','file')       % one of the outputs of kilosort
    % load kilosort data
    if exist('rez.mat','file')
        load('rez.mat');          % load rez.mat output from kilo (for waveform extraction)
        % get cluster group
        if exist('cluster_group.tsv','file')
            fid=fopen('cluster_group.tsv');     % new phy2 output
        else
            fid = fopen('cluster_groups.csv');
        end
        cat_data = 0;
    elseif exist('..\rez.mat','file')       % this means files were concatenated
        cat_data = 1;
        load('..\rez.mat');
        % get cluster group
        if exist('..\cluster_group.tsv','file')
            fid=fopen('..\cluster_group.tsv');     % new phy2 output
        else
            fid = fopen('..\cluster_groups.csv');
        end
    else
        warning('missing rez.mat\n')
    end
        
    spike_times = readNPY('spike_times.npy');           % here, spike times from ALL clusters
    clusters = readNPY('spike_clusters.npy');
%     if ~exist(fullfile(exp_path,'dropped_spikes.mat'),'file')
%         [spike_times,clusters] = clean_light_artifacts(field_trials,LED,spike_times,clusters,amp_sr);      % NEW MAK addition 2/26/17 - get rid of light artifacts mistaken for spikes - modified and confirmed works 8/9/17
%     else
%         load(fullfile(exp_path,'dropped_spikes.mat'));
%         spike_times(dropped_spikes) = [];
%         clusters(dropped_spikes) = [];
%     end
    
    % read column headers
    C_text = textscan(fid, '%s', 1, 'delimiter', ',');
    % read group info
    grp_dat = textscan(fid, '%f %s', 'delimiter', ',');
    fclose(fid);
    units = grp_dat{1};             % cluster numbers
    cluster_group = grp_dat{2};     % 'good', 'noise', etc.
    good_units = units(strcmp(cluster_group,'good'));
else        % load phy data
    kwik_file = sprintf('%s/amplifier.kwik',exp_path);   
    info = hdf5info(kwik_file);
    spike_times = hdf5read(kwik_file, '/channel_groups/0/spikes/time_samples');
    clusters = hdf5read(kwik_file, '/channel_groups/0/spikes/clusters/main');   % currently only from 1st shank!
    rez = [];
end

if exist('LED','var')       % if optogenetics experiment
%     led = light(izx);
    disp('Loading light parameters...')
    if ~exist('light_params.mat','file')
        [params.all_light,params.pulse_dur,params.lighttime,params.av_light_start] = get_lightstim_v3(exp_path,exp_type);   % changed from v2 9/2/20 MAK
        params.light_dur = params.pulse_dur;
    else
        load('light_params.mat');
        params.all_light = all_light;
        params.light_dur = pulse_dur;
        params.pulse_dur = params.light_dur;
        params.lighttime = floor(min(pulse_dur(pulse_dur>0))*1000)/1000;
        if strcmp(exp_type,'trains')
            params.lighttime = round(max(params.light_dur(params.light_dur>0)));
            params.light_dur(params.light_dur>0) = params.lighttime;
            params.pulse_dur(params.pulse_dur>0) = min(params.pulse_dur(params.pulse_dur>0));
        end
        params.av_light_start = av_light_start;
        clear all_light pulse_dur av_light_start
    end
end

%% get the details of the experiment (trial types, prestim time, etc.)
disp('Loading trial variables and conditions...')
[params.prestim,params.poststim,params.stimtime,params.trial_type,params.IVs] = get_exp_params(exp_path,exp_type);
params.onset = .1;         % hardcoded (in sec)
params.exp_name = exp_name;
params.exp_type = exp_type;
params.amp_sr = amp_sr;
params.nchs = chansInDat;

%% run unit analysis for each unit in this experiment
% if exist(sprintf('%s/good_units.txt',exp_path),'file')
%     good_units = load(sprintf('%s/good_units.txt',exp_path));   % if you predetermined which clusters to look at
% end

num_units = length(good_units);
for i = 1:num_units
    fprintf(sprintf('Performing analysis on %s cluster %s \n',exp_name,num2str(good_units(i))))
    if exist('LED','var') 
        [unitinfo(i),FRs(i),tuning(i),waveforms(i)] = unit_analysis_opto(good_units(i),field_trials,spike_times(clusters==good_units(i)),params,rez,'blanks',plots,all_fig_dir);
    else
        [unitinfo(i),FRs(i),tuning(i),waveforms(i)] = unit_analysis_noopto(good_units(i),field_trials,spike_times(clusters==good_units(i)),params,rez,'prestim',plots,all_fig_dir);
    end
end

%% evaluate cluster quality
if ~exist(sprintf('%s\\%s_cluster_quality.mat',main_dir,an_name),'file') && ~exist(sprintf('%s\\%s_cluster_quality.mat',exp_dir,exp_name),'file')
    disp('Running cluster quality analysis')
%     [refr_idx, SNR] = eval_clusters(good_units,spike_times,clusters,chansInDat,amp_sr,exp_system,exp_path);
%     [refr_idx, ~] = eval_clusters(good_units,spike_times,clusters,chansInDat,amp_sr,exp_system,exp_path);
    if exist('rez.mat','file')
        [cids, uQ, cR, isiV, refV] = sqKilosort.computeAllMeasures(exp_path);
    elseif exist('..\rez.mat','file')
        [cids, uQ, cR, isiV, refV] = sqKilosort.computeAllMeasures(fileparts(exp_path));
    end
    if cat_data
%         save(sprintf('%s\\%s_cluster_quality.mat',main_dir,an_name),'SNR','refr_idx','uQ','cR','isiV')
        save(sprintf('%s\\%s_cluster_quality.mat',main_dir,an_name),'cids','isiV','refV','uQ','cR')
    else
%         save(sprintf('%s\\%s_cluster_quality.mat',exp_dir,exp_name),'SNR','refr_idx','uQ','cR','isiV')
        save(sprintf('%s\\%s_cluster_quality.mat',exp_dir,exp_name),'cids','isiV','refV','uQ','cR')
    end
end

%% save full experimental results
disp(strcat('Saving results.mat for experiment ',exp_name))

%  save to experiment folder within full_dir
% try
chk = whos('unitinfo');
if chk.bytes>2e9    % if greater than 2GB
    disp('unitinfo too big to save with default version 6 - using ver. 7.3 (will take longer)')
    save(sprintf('%s\\%s_results.mat',exp_dir,exp_name),'params','exp_type','unitinfo','FRs','tuning','waveforms','-nocompression','-v7.3')
else
    save(sprintf('%s\\%s_results.mat',exp_dir,exp_name),'params','exp_type','unitinfo','FRs','tuning','waveforms')
end
    % catch
%     rast = arrayfun(@(x) x.rast ,unitinfo,'uniformoutput',0);   % if unitinfo variable is too big to save
%     unitinfo = rmfield(unitinfo,'rast');
%     save(sprintf('%s\\%s_results.mat',exp_dir,exp_name),'params','exp_type','unitinfo','FRs','tuning','waveforms','rast','-v7.3')
% end

end