function RFanalysis_master(exp_path)
% inputs:
% exp_path - path to the sparsenoise experiment, as a string (e.g.,
    % 'Y:\L5\stGtACR\NPRSC5\LPl\2020-07-22_14-45-12_sparsenoise_2LEDs')
% maybe try to add in later...
% good_units(OPTIONAL) - specify which units (by their cluster name) you
    % want to calc RFs for
    
% for openEphys data
binsize = 10;      % in ms

%% get experiment name and load raw data
% get name
out = regexp(exp_path,'\\','split');
an_name = out{end-1};       % animal name
if contains(an_name,'LP') || contains(an_name,'V1') || contains(an_name,'SC')
    an_name = strcat(out{end-2},'_',an_name);
elseif contains(an_name,'day','ignorecase',1) 
    an_name = out{end-2};
end
inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
exp_name = out{end}(1:inds(1)-1);
if sum(isstrprop(exp_name,'digit'))>4       % this means it's an OpenEphys file (because it names files starting with the date)
    inds = regexp(out{end},'_\D','start');
    exp_name = strcat(out{end-1},out{end}(inds(1):end));  
    if strfind(out{end}(inds(1):end),an_name)
        exp_name = out{end}(inds(1)+1:end);
    else
        exp_name = strcat(an_name,out{end}(inds(1):end));
    end
    exp_name = exp_name(1:strfind(exp_name,'sparsenoise')-2); % remove "sparsenoise" from end of exp_name
    
end

disp(strcat('Loading raw data from experiment ',exp_name))
cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      
else
    openEphys2matlab_sparsenoise(exp_path);
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
% stim_times - 
% start_times - 

% set up where to store results
full_dir = 'H:\LPproject\LPresults';        % where you want to save data
main_dir = strcat(full_dir, '\',an_name);
all_fig_dir = sprintf('%s\\RFfigures',main_dir);
if ~exist(all_fig_dir)
    mkdir(all_fig_dir)
end
% if plots
%     if exist(all_fig_dir,'dir')
%         delete(sprintf('%s/*',all_fig_dir));        % delete and remake all figures
%     end
% end

%% load clustering info (if units to analyze weren't specified in good_units entry)
% if nargin(good_units) < 2  % if units to analyze weren't pre-specified 
    disp('Loading clustering data...')
    % load kilosort data
    spike_times = readNPY('spike_times.npy');           % here, spike times from ALL clusters
    clusters = readNPY('spike_clusters.npy');
    if exist('..\cluster_group.tsv','file')
        fid=fopen('..\cluster_group.tsv');     % new phy2 output
    else
        fid = fopen('..\cluster_groups.csv');
    end
    % if exist('cluster_groups.csv','file')
    %     fid = fopen('cluster_groups.csv');
    % elseif exist('..\cluster_groups.csv','file')
    %     % get cluster group
    %     fid = fopen('..\cluster_groups.csv');
    % end

    % read column headers
    C_text = textscan(fid, '%s', 1, 'delimiter', ',');
    % read group info
    grp_dat = textscan(fid, '%f %s', 'delimiter', ',');
    fclose(fid);
    units = grp_dat{1};             % cluster numbers
    cluster_group = grp_dat{2};     % 'good', 'noise', etc.
    good_units = units(strcmp(cluster_group,'good'));
% end

%% get stimulus info
stim_file = dir(fullfile(exp_path,'_*_*_*.mat'));       % output format for stimulus info .mat file
fprintf(sprintf('Loading stimulus file %s\n',stim_file.name'))
stimdat = load(fullfile(exp_path,stim_file.name),'-mat');
analyze_file = dir(fullfile(exp_path,'*.analyzer'));
fprintf(sprintf('Loading analyzer file %s\n',analyze_file.name'))
load(fullfile(exp_path,analyze_file.name),'-mat');  % under variable name 'Analyzer'
nX = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'Nx'), Analyzer.P.param)}{3};
nY = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'Ny'), Analyzer.P.param)}{3};
nBW = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'bw_bit'), Analyzer.P.param)}{3};
% num_reps = length(fieldnames(stimdat))-1;       % assumes one of the fields is for monitor refresh rate ('frate')
num_reps = sum(cellfun(@(x) contains(x,'randlog'),fieldnames(stimdat)));
num_stim = size(stimdat.randlog_T1.seqs.xseq,1);
num_sqs = size(stimdat.randlog_T1.seqs.xseq,2); % when multiple squares are presented at once

conds = cellfun(@(x) x.val{strcmp(Analyzer.loops.conds{1}.symbol,'light_bit')}, Analyzer.loops.conds);
diff_conds = unique(conds);
trialcond = cellfun(@(x) x.repeats{1}.trialno,Analyzer.loops.conds);
lightcond = zeros(1,num_reps);
for i = 1:length(diff_conds)
    lightcond(ismember(trialcond,find(conds==diff_conds(i)))) = diff_conds(i);
end
if length(diff_conds)>1
    lightexp = 1;
else
    lightexp = 0;
end

% downsample LED to 1000Hz
div             = amp_sr/1000;
zx              = 1:div:length(photo);
izx             = floor(zx);
if lightexp
    light = LED(:,izx);
end
% pho = photo(izx);
   
T = 1000/(60/Analyzer.P.param{14}{3});    % length of each stim presentation (in ms). 14th entry is h_per
light_T = Analyzer.LP.pulse_dur(1);       % length of light pulse duration (in ms)
if ~exist(fullfile(exp_path,'stim_trials.mat'),'file') % NEW 1/16/21 - save stim_trials because it takes so long to load...
    % make matrix of "trials" for each sparse noise stimulus presentation
    stim_trials = zeros(num_reps*num_stim,5);       % 5 columns for x pos, y pos, b/w, and lightcond, and trialnum
    for n = 1:num_reps
        xx = eval(sprintf('stimdat.randlog_T%d.seqs.xseq',n));
        yy = eval(sprintf('stimdat.randlog_T%d.seqs.yseq',n));
        bw = eval(sprintf('stimdat.randlog_T%d.seqs.bwseq',n));
        bw(bw==1) = -1; % 1 means BLACK
        bw(bw==2) = 1;  % 2 means WHITE
        for nn = 1:num_stim
            stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),1:3) = [xx(nn,:)' yy(nn,:)' bw(nn)*ones(num_sqs,1)]; % WILL NEED TO CHANGE WHEN LOGFILE IS FIXED TO SAVE BW FOR EXTRA SQUARES
            stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),5) = (n-1)*num_stim+nn;
            %         stim_trials((n-1)*num_stim+nn,1:3) = [xx(nn) yy(nn) bw(nn)]; 
    %         thresh = unique([round(max(light,[],2),2)*min(T,light_T)*.95; sum(round(max(light,[],2),2)*min(T,light_T)*.95)]); % this and next two lines are new 9/10/20 to account for multiple LEDs and LED conditions
    %         if lightexp && sum(sum(light(:,stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T),2)) > min(thresh)
    %             stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),4) = find(thresh < sum(sum(light(:,stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T),2)),1,'last');  
            if lightexp && sum(light(1,stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T)) > round(max(light(1,:))-std(light(1,:)),2)*min(T,light_T)     
                if size(light,1)>1 && lightexp && sum(light(2,stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T)) > round(max(light(2,:))-std(light(2,:)),2)*min(T,light_T)
                    stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),4) = 3;
                else
                    stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),4) = 1;
                end
            elseif lightexp && size(light,1)>1 && sum(light(2,stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T)) > round(max(light(2,:))-std(light(2,:)),2)*min(T,light_T)  
                %         if find(light(stim_times((n-1)*num_stim+1):stim_times((n)*num_stim))>4) 	% TEMP
    %             stim_trials((n-su1)*num_stim+nn,4) = 1;
                stim_trials((n-1)*(num_stim*num_sqs)+((nn-1)*num_sqs+1:nn*num_sqs),4) = 2;
            end
        end
    end
    save(fullfile(exp_path,'stim_trials.mat'),'stim_trials')
else
    load(fullfile(exp_path,'stim_trials.mat'))
end

% % much faster way of deriving stim_trials(:,4) (but potentially less
% % accurate if light control wasn't 100% precise...)
% lightcond = find(strcmp(Analyzer.loops.conds{1}.symbol,'light_bit'));
% conds = cellfun(@(x) x.val{lightcond}, Analyzer.loops.conds);
% diff_conds = unique(conds);
% trial_type = zeros(1,size(field_trials,1));
% for c = 1:length(conds)
%     trials{:,c} = cellfun(@(x) x.trialno, Analyzer.loops.conds{c}.repeats);     % store trials from different conditions in separate columns
%     trial_type(trials{:,c}) = c;        
% end
% tcond = conds(trial_type);
% stimL = (60/Analyzer.P.param{14}{3})*light_T/1000;
% new_stim_trials = nan(num_stim,1);
% for t=1:length(tcond)
%     if tcond(t) == 1
%         new_stim_trials(:,t) = repmat([repmat(1,stimL,1); repmat(2,stimL,1)],num_stim/stimL/2,1);
%     elseif tcond(t) == 2
%         new_stim_trials(:,t) = repmat([repmat(2,stimL,1); repmat(1,stimL,1)],num_stim/stimL/2,1);
%     elseif tcond(t) == 3
%         new_stim_trials(:,t) = repmat([repmat(3,stimL,1); repmat(0,stimL,1)],num_stim/stimL/2,1);
%     else
%         new_stim_trials(:,t) = repmat([repmat(0,stimL,1); repmat(3,stimL,1)],num_stim/stimL/2,1);
%     end
% end
% test = new_stim_trials(:);


lconds = unique(stim_trials(:,4));
bwconds = unique(stim_trials(:,3));
if lightexp && sum(stim_trials(:,4)>=1) ~= sum(lconds>0)*size(stim_trials,1)/length(lconds)
    warning('Unequal numbers of light and no-light trials detected')
end

% make "key" matrix for all possible sparse noise stimuli
count = 1;
stim_key = nan(nBW*nY*nX,4);
for lc = 1:length(lconds)
    for b = 1:nBW
        for y = 1:nY
            for x = 1:nX
                stim_key(count,:) = [x y bwconds(b) lconds(lc)];
                count = count+1;
            end
        end
    end
end

%% find running trials
if exist('encdA','var')
    move_trials = Intan_digital_movement([stim_times stim_times+250],encdA,encdB,0)';
elseif exist('mouse','var')    % using optical mouse for old recordings - might want to change to encoder?
    move_trials = zeros(1,size(field_trials,1));
    moveVec = move_trials;
    for i = 1:size(field_trials,1)
        moveVec(i) = sum(mouse(field_trials(i,1)*20:field_trials(i,2)*20));         % *20 because mouse is still in original amp_sr
    end
    move_trials(find(moveVec)) = 1;
else
    move_trials = zeros(1,size(field_trials,1));
end

%% run RF analysis for each unit in this experiment

num_units = length(good_units);
% get units' shanks
full_exp_path = fileparts(exp_path);
load(fullfile(full_exp_path,'rez.mat'));          % load rez.mat output from kilo (for waveform extraction)

for i = 1:num_units
    fprintf(sprintf('Performing RF analysis on %s cluster %s \n',exp_name,num2str(good_units(i))))
    unit_times = spike_times(clusters==good_units(i));
    unit_times_ds = floor(unit_times./(amp_sr/1000));   % change sample #s to account for downsampling
    unit_times_ds = unit_times_ds + 1; % has to be +1 because spike_times starts at 0, but the min possible field_trials value could be 1
    [~,max_ch(i),shank(i)] = readWaveformsFromRez_K2(good_units(i),full_exp_path,rez);      % default to read waveforms from Rez
    [psth_norm{i}, rfMap{i}, stats(i)] = sparsenoiseRF(unit_times_ds,stim_times,binsize,T,stim_trials,move_trials,stim_key);
    % save figures
    save_name = sprintf('%s\\%s_%s_%d_Cluster%d',all_fig_dir,exp_name,strcat('shank',num2str(shank(i))),max_ch(i),good_units(i));
    set(gcf, 'PaperPosition', [0 0 24 12]);
    print(gcf,'-dpng',save_name)
    close all
end

%% save results
save(fullfile(main_dir,'RFresults'),'good_units','shank','max_ch','psth_norm','rfMap','stats')

end
