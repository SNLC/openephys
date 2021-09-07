function experiment_analysis_dTA(exp_path,probe,area)


params = [];
% exp_type = [];
unitinfo = [];
FRs = [];
tuning = [];
waveforms = [];
refr_idx = [];
% SNR = [];
% isiV = [];
cR = [];
uQ = [];
refV = [];
exp_num = [];

full_dir = 'H:\LPproject\LPresults';

% get experiment name
out = regexp(exp_path,'\\','split');
an_name = out{end-1};       % animal name
if contains(an_name,'LP','ignorecase',1) || contains(an_name,'V1','ignorecase',1) || contains(an_name,'TRN','ignorecase',1)
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

main_dir = strcat(full_dir, '\',an_name);
exp_dir = strcat(main_dir,'\',exp_name);
if ~exist(exp_dir,'dir') && exist(sprintf('%s\\%s\\%s',main_dir,an_name,exp_name))       % need better solution for this
    exp_dir = sprintf('%s\\%s\\%s',main_dir,an_name,exp_name);
end
fprintf(sprintf('Processing experiment %s\n',exp_name))
% if ~exist(exp_dir,'dir')
%     analysis_master(exp_path,'OpenEphys','step',0)
% end


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
    cluster_file = sprintf('%s\\%s\\%s_cluster_quality.mat',full_dir,an_name,an_name);
end
exp = importdata(results_file,'-mat');     % load results mat
clust = importdata(cluster_file,'-mat');
clear cluster_file results_file

if isfield(exp.waveforms,'shank')
    shk = [exp.waveforms.shank];
else
    shk = ones(1,length(exp.FRs));
end
shanks = unique(shk);
for sh = 1:length(shanks)    
    if exist(sprintf('%s/shank%d_LPchannels.txt',exp_path,shanks(sh)),'file')
        channels{sh} = load(sprintf('%s/shank%d_LPchannels.txt',exp_path,shanks(sh)));   % if you predetermined which clusters to look at
    elseif exist(sprintf('%s/good_channels.txt',exp_path),'file')
        channels = load(sprintf('%s/good_channels.txt',exp_path));   % if you predetermined which clusters to look at
    end
end

params = [params exp.params];
%     exp_type = [exp_type exp.exp_type];
unitinfo = [unitinfo exp.unitinfo];
FRs = [FRs exp.FRs];
if isfield(exp.tuning,'SF')
    exp.tuning = rmfield(exp.tuning,'SF');      % for M12, remove SF and TF fields (for now...just so I can concatenate results from prior experiments)
    exp.tuning = rmfield(exp.tuning,'TF');
end
tuning = [tuning exp.tuning];
waveforms = [waveforms exp.waveforms];
%     refr_idx = [refr_idx clust.refr_idx];
%     SNR = [SNR clust.SNR];
%     isiV = [isiV clust.isiV];
uQ = [uQ clust.uQ'];
cR = [cR clust.cR'];
refV = [refV clust.refV];


%% calculate relevant significance values
% using kruskal-wallis (non-parametric, indep samples) for testing vis and
% light modulation, using hotellings for significant orientation tuning.
% Are these the right tests to use??
if strcmp(params.IVs,'s_freq')
    vis_trials = find((params.trial_type(:,strcmpi(params.IVs,'visual'))==1)&(params.trial_type(:,strcmpi(params.IVs,'s_freq'))==.04)&(params.trial_type(:,strcmpi(params.IVs,'t_period'))==30));  % in case of multiple SFs and TFs, only compare regular 2Hz and .04cpd trials (for now)
else
    vis_trials = find(params.trial_type(:,strcmpi(params.IVs,'visual'))==1);
end
blank_trials = find(params.trial_type(:,strcmpi(params.IVs,'visual'))==0);
nolight_trials = find(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==0);
lightconds = unique(params.trial_type(:,strcmpi(params.IVs,'light_bit')));
% if length(lightconds) > 3
%     lightconds= lightconds([1 2 end],:);    % if experiment had low, medium and high intensity light conditions, drop the medium condition (b/c M12 only has low and high)
% end
for lc=2:length(lightconds)  % assumes first condition is no-light condition
    light_trials{lc-1} = find(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==lightconds(lc));
    vislight_trials{lc-1} = intersect(vis_trials,find(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==lightconds(lc)));
    blanklight_trials{lc-1} = intersect(blank_trials,find(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==lightconds(lc)));
end
run_trials= find(params.trial_type(:,strcmpi(params.IVs,'running'))==1);
stat_trials = find(params.trial_type(:,strcmpi(params.IVs,'running'))==0);

% all_units = 1:length(unitinfo);
% units = all_units;
% distfromlastch = nan(1,length(units));
% distfromfirstch = distfromlastch;
% for n = 1:length(params)
%     if iscell(channels{n})
%         shk = [waveforms(exp_num==n).shank];
%         shks = shanks;
%         exp_trials = find(exp_num==n);
%         for ss = 1:length(shks)
%             distfromlastch(exp_trials(shk==shks(ss))) = max(channels{n}{ss})-[waveforms(exp_trials(shk==shks(ss))).max_ch];  
%             distfromfirstch(exp_trials(shk==shks(ss))) = min(channels{n}{ss})-[waveforms(exp_trials(shk==shks(ss))).max_ch];
%         end
%     else
%         distfromlastch(exp_num==n) = max(channels{n})-[waveforms(exp_num==n).max_ch]; 
%         distfromfirstch(exp_num==n) = min(channels{n})-[waveforms(exp_num==n).max_ch]; 
%     end
% end
% units(distfromlastch<0|distfromfirstch>0|isnan(distfromlastch)) = [];
% clean_units = units(SNR(units)>=1.5&refr_idx(units)<1);      % only include units that pass SNR and refractory period thresholds
% distfromlastch = distfromlastch(clean_units);          % only includes GOOD units
% distfromfirstch = distfromfirstch(clean_units);          % only includes GOOD units
% 
% FRb = [FRs(clean_units).baseline];  %  baseline (from blank trials w/ no running, light or vis stim)
% 
% for i = 1:length(clean_units)      % for each unit
%     nn = clean_units(i);
%     tuning_curve{i} = tuning(nn).curve(:,:);
%     oris = unique(params.trial_type(:,strcmp(params.IVs,'ori')));
%         oris(oris>=999) = [];
%     % kruskal-wallis test to test for significant and visual- and light-modulation
%     % for visual modulation, find preferred direction trials - only use THESE
%     % trials to test for significant visual modulation (in case of extremely
%     % tuned cells)
%     [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-repmat(FRb(i),1,size(tuning_curve,3))));   % direction w/ biggest change from baseline
%     prefori_trials(i,:) = find(params.trial_type(:,strcmpi(params.IVs,'ori'))==oris(prefdir_deg(i)));
%     vis_sig(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials(i,:),nolight_trials),round(1000*(params.av_light_start(1)))+1:round(1000*(params.av_light_start(1)))+params.lighttime*1000),2)'...
%         sum(unitinfo(nn).rast(intersect(blank_trials,nolight_trials),round(1000*(params.av_light_start(1)))+1:round(1000*(params.av_light_start(1)))+params.lighttime*1000),2)'],...
%         [ones(1,length(intersect(prefori_trials(i,:),nolight_trials))) 2*ones(1,length(intersect(blank_trials,nolight_trials)))],'off');    % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
%     vis_sig_ons(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials(i,:),nolight_trials),1000*params.prestim+1:1000*(params.prestim+params.onset)),2)'... 
%         sum(unitinfo(nn).rast(intersect(blank_trials,nolight_trials),1000*params.prestim+1:1000*(params.prestim+params.onset)),2)'],...
%         [ones(1,length(intersect(prefori_trials(i,:),nolight_trials))) 2*ones(1,length(intersect(blank_trials,nolight_trials)))],'off');         % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
%     for lc = 1:length(lightconds)-1
%         light_sig(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials,round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)))+params.lighttime*1000),2)'...       % currently using ALL light trials to evaluate light significance (visual+blank)
%             sum(unitinfo(nn).rast(light_trials{lc},round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)))+params.lighttime*1000),2)'],...
%             [ones(1,length(nolight_trials)) 2*ones(1,length(light_trials{lc}))],'off');    % significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
%         light_sig_ons(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials,round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)+params.onset))),2)'...
%             sum(unitinfo(nn).rast(light_trials{lc},round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)+params.onset))),2)'],...
%             [ones(1,length(nolight_trials)) 2*ones(1,length(light_trials{lc}))],'off');         % significance of light-evoked response at light onset (100ms after light onset in nolight vs light trials (each condition separately))
%     end
%     
%    % next, check tuning significance and get tuning curves
%      if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
%         % evaluate significance or orientation tuning
%         for o = 1:length(oris)/2
%             for lc = 1:length(lightconds)
%                 ori_trials = find((ismember(params.trial_type(:,strcmpi(params.IVs,'ori')),oris([o o+length(oris)/2])))&(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==lightconds(lc)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
%                 tuning_trials(:,o,lc) = sum(unitinfo(nn).rast(ori_trials,round(1000*(params.av_light_start(1)))+1:round(1000*(params.av_light_start(1)))+params.lighttime*1000),2); % numbers of spikes across trials of given light condition for each orientation (by column)
% 
%             end
%             tuning_curve_collapse(i,:,o)   = mean([tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)],2);         % average evoked FR of orientations collapsed across directions
% %             tuning_curve(i,:,[o o+length(oris)/2]) = [tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)];
%         end
%         for lc = 1:length(lightconds)      % currently, NOT separating running and stationary trials
%             tuned_sig(i,lc) = T2Hot1(tuning_trials(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
%         end
%     end
% end

% NEW (6/15/18) - calculate visual significance prior to getting distances from first
% and last ch and deciding which units are "clean". Use first and last
% visually significant channels to determine borders of LP


% good_SNR = find(SNR>=1.5&refr_idx<.1); % only include units that pass SNR and refractory period thresholds
 good_SNR = find(refV<.5);
% good_isi = find(isiV<.1);       % only include units with <10% ISI violations
% clean_units = intersect(good_SNR,good_isi(cR<.3 | isnan(cR)));   % and exclude units with 30% or more contamination rate of other good_isi units (include NANs because doesn't necessarily mean they're bad)
% clean_units = intersect(good_SNR,good_isi);    
good_uQ = find(uQ>16);       % only include units with <25% contamination "false positive" rate
clean_units = intersect(good_SNR,good_uQ); 

FRb = [FRs(clean_units).baseline];  %  baseline (from blank trials w/ no running, light or vis stim)
for i = 1:length(clean_units)      % for each unit
    nn = clean_units(i);
    tuning_curve{i} = tuning(nn).curve(:,:);
    oris = unique(params.trial_type(:,strcmp(params.IVs,'ori')));
        oris(oris>=999) = [];
    % kruskal-wallis test to test for significant and visual- and light-modulation
    % for visual modulation, find preferred direction trials - only use THESE
    % trials to test for significant visual modulation (in case of extremely
    % tuned cells)
    [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-repmat(FRb(i),1,size(tuning_curve,3))));   % direction w/ biggest change from baseline
    prefori_trials{i} = find(params.trial_type(:,strcmpi(params.IVs,'ori'))==oris(prefdir_deg(i)));
    vis_sig(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials),1000*params.prestim*2+1:3*1000*params.prestim),2),... % significance of visual response (compared to prestim - MAY NEED TO CHANGE)
        sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials),1:1000*params.prestim),2));
    vis_sig_ons(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials),1000*params.prestim+1:1000*2*params.prestim),2),...  % significance of visual response (compared to prestim - MAY NEED TO CHANGE)
        sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials),1:1000*params.prestim),2));
%     for lc = 1:length(lightconds)-1
%         light_sig(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials,round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)))+1000),2)'...       % currently using ALL light trials to evaluate light significance (visual+blank)
%             sum(unitinfo(nn).rast(light_trials{lc},round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)))+1000),2)'],...
%             [ones(1,length(nolight_trials)) 2*ones(1,length(light_trials{lc}))],'off');    % significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
%         light_sig_ons(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials,round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)+params.onset))),2)'...
%             sum(unitinfo(nn).rast(light_trials{lc},round(1000*(params.av_light_start(lc)))+1:round(1000*(params.av_light_start(lc)+params.onset))),2)'],...
%             [ones(1,length(nolight_trials)) 2*ones(1,length(light_trials{lc}))],'off');         % significance of light-evoked response at light onset (100ms after light onset in nolight vs light trials (each condition separately))
%     end
    
   % next, check tuning significance and get tuning curves
     if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
        % evaluate significance or orientation tuning
        for o = 1:length(oris)/2
            for lc = 1:length(lightconds)
                ori_trials = find((ismember(params.trial_type(:,strcmpi(params.IVs,'ori')),oris([o o+length(oris)/2])))&(params.trial_type(:,strcmpi(params.IVs,'light_bit'))==lightconds(lc)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
                if length(ori_trials)>40
                    ori_trials(randi(length(ori_trials),1)) = [];
                end
                tuning_trials(:,o,lc) = sum(unitinfo(nn).rast(ori_trials,1000+1:1000+1000),2); % numbers of spikes across trials of given light condition for each orientation (by column)
                
            end
            tuning_curve_collapse(i,:,o)   = mean([tuning(nn).curve(:,o) tuning(nn).curve(:,o+length(oris)/2)],2);         % average evoked FR of orientations collapsed across directions
%             tuning_curve(i,:,[o o+length(oris)/2]) = [tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)];
        end
        for lc = 1:length(lightconds)      % currently, NOT separating running and stationary trials
            tuned_sig(i,lc) = T2Hot1(tuning_trials(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
        end
    end
end

distfromlastch = nan(1,length(clean_units));
distfromfirstch = distfromlastch;
vis_units = find(vis_sig<.025|vis_sig_ons<.025);      % CHANGED 7/21/19 (.025 b/c bonferroni correction for 2 tests)

if isfield(waveforms,'shank')
    shk = [waveforms.shank];
    if isempty(shk); shk = zeros(1,length(unitinfo)); end    % if single shank probe, shank field may be empty
else
    shk = zeros(1,length(unitinfo));
end

% get probe info (NEW 1/22/19)
p = eval(sprintf('probemap_%s_func',probe));
Zchan = flipud(sort(p.z(p.shaft==1))); % from top to bottom (bottom=0)
shanks = unique(p.shaft)-1; % -1 to set first shank to 0
for sh = 1:length(shanks)    % for each shank in exp
    shk_units{sh} = find(shk(clean_units)==shanks(sh));
%     firstch(count) = min([waveforms(clean_units(intersect(shk_units{count},vis_units))).max_ch]);
    vischs = sort([waveforms(clean_units(intersect(shk_units{sh},vis_units))).max_ch]);    
    firstch(sh) = min(vischs); % assumes first channel with good unit must be in thalamus
    if firstch > 1
        if Zchan(firstch(sh))==Zchan(firstch(sh)-1) % for probes in hexagonal orientation, might leave out channel that is actually same height as "firstch"
            firstch(sh) = firstch(sh)-1;
        end
    end
    lastch(sh) = max(vischs);    % need better solution for when shank passes through not very visual area
    if strcmpi(area,'trn') && (sum(abs(diff(Zchan(vischs))))) > 100      % if more than 50um separates consecutively located visually-responsive units in TRN experiment
        lastinTRN = find(abs(diff(Zchan(vischs)))>100,1,'first');
        lastch(sh) = vischs(lastinTRN);
    end
    if lastch < length(Zchan)
        if Zchan(lastch(sh))==Zchan(lastch(sh)+1)
            lastch(sh) = lastch(sh)+1;
        end
    end
%     distfromlastch(shk_units{sh}) = lastch(sh)-[waveforms(clean_units(shk_units{sh})).max_ch];
%     distfromfirstch(shk_units{sh}) = firstch(sh)-[waveforms(clean_units(shk_units{sh})).max_ch];
    distfromlastch(shk_units{sh}) = Zchan([waveforms(clean_units(shk_units{sh})).max_ch]) - Zchan(lastch(sh)); % get actual distance, in um
    distfromfirstch(shk_units{sh}) = Zchan(firstch(sh)) - Zchan([waveforms(clean_units(shk_units{sh})).max_ch]);
end
unit_chk = 1:length(clean_units);
% unit_chk(distfromlastch<0|distfromfirstch>0|isnan(distfromlastch)) = [];
unit_chk(distfromlastch<0|distfromfirstch<0|isnan(distfromlastch)) = [];
clean_units = clean_units(unit_chk);
distfromlastch = distfromlastch(unit_chk);          % only includes GOOD units
distfromfirstch = distfromfirstch(unit_chk);          % only includes GOOD units
tuning_curve=tuning_curve(unit_chk);
prefdir_deg = prefdir_deg(unit_chk);
prefori_trials = prefori_trials(unit_chk);
vis_sig = vis_sig(unit_chk);
vis_sig_ons = vis_sig_ons(unit_chk);
% light_sig = light_sig(unit_chk,:);
% light_sig_ons = light_sig_ons(unit_chk,:);
tuning_curve_collapse = tuning_curve_collapse(unit_chk,:,:);
tuned_sig = tuned_sig(unit_chk,:);
FRb = FRb(unit_chk);

%% test for light modulation
% verified this combo of reshape, cell2mat and arrayfun yields accurate
% results! MAK 1/31/18
% numconds = arrayfun(@(x) length(unique(x.all_light)),params,'uniformoutput',1);     % number of light conds in each experiment
numconds = 1;
conds = 1:min(numconds)-1;
FRev = reshape(cell2mat(arrayfun(@(x) x.visual.ev(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';    %vis-evoked
FRearly = reshape(cell2mat(arrayfun(@(x) x.visual.evstart(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';   % early evoked period
FRlate = reshape(cell2mat(arrayfun(@(x) x.visual.evlate(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';  % late evoked period
FRonset = reshape(cell2mat(arrayfun(@(x) x.visual.evlightonset(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % vis-evoked light onset
FRbl = reshape(cell2mat(arrayfun(@(x) x.blank.ev(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % blank, evoked period
FRvison = reshape(cell2mat(arrayfun(@(x) x.onset(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % vis stim onset
% for ii = 2:size(FRev,2)
%     lightmod(:,ii-1) = (diff(FRev(:,[1 ii]),[],2))./sum(FRev(:,[1 ii]),2);
%     lightmod_early(:,ii-1) = (diff(FRearly(:,[1 ii]),[],2))./sum(FRearly(:,[1 ii]),2);
%     lightmod_late(:,ii-1) = (diff(FRlate(:,[1 ii]),[],2))./sum(FRlate(:,[1 ii]),2);
%     lightmod_onset(:,ii-1) = (diff(FRonset(:,[1 ii]),[],2))./sum(FRonset(:,[1 ii]),2);
% end

% and visual modulation
vismod = (FRev(:,1) - FRb')./(FRev(:,1) + FRb');        % using baseline (from blank trials w/ no running, light or vis stim)
vismod_on = (FRvison(:,1)-FRb')./(FRvison(:,1)+FRb');

% get orientation values for later use
OSI_CV = reshape(cell2mat(arrayfun(@(x) x.OSI_CV(1,[conds end]), tuning(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; 
OSI = reshape(cell2mat(arrayfun(@(x) x.OSI(1,[conds end]), tuning(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; 
DSI_CV = reshape(cell2mat(arrayfun(@(x) x.DSI_CV(1,[conds end]), tuning(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; 
DSI = reshape(cell2mat(arrayfun(@(x) x.DSI(1,[conds end]), tuning(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; 


%% get waveform props
for i = 1:length(clean_units)
    [t2p_t(i),t2p_r(i),fwhm(i)] = get_waveform_props(waveforms(clean_units(i)).microV,params.amp_sr);
end

% %% NEW test: lightmod for trains experiments
if strfind(exp_path,'trains')
    for ii=1:2      % for two light conditions >1Hz
        resp_t_v{ii} = nan(1,length(clean_units));
        ppr_v{ii} = nan(length(clean_units),lightconds(ii+2)-1);    
        resp_dur_v{ii} = nan(length(clean_units),lightconds(ii+2));
        resp_t_bl{ii} = resp_t_v{ii};
        ppr_bl{ii} = ppr_v{ii};
        resp_dur_bl{ii} = resp_dur_v{ii};
        for i = 1:length(clean_units)
            vis_inds = (params.trial_type(:,1)==1);
            [resp_t_v{ii}(i), resp_dur_v{ii}(i,:), ppr_v{ii}(i,:)]  = ppanalysis(params.prestim, params.stimtime, round(mean(params.av_light_start))-params.prestim, max(params.light_dur), lightconds(ii+2), unitinfo(clean_units(i)).rast(vis_inds,:), params.trial_type(vis_inds,strcmpi(params.IVs,'light_bit'))) ;
           [resp_t_bl{ii}(i), resp_dur_bl{ii}(i,:), ppr_bl{ii}(i,:)]  = ppanalysis(params.prestim, params.stimtime, round(mean(params.av_light_start))-params.prestim, max(params.light_dur), lightconds(ii+2), unitinfo(clean_units(i)).rast(~vis_inds,:), params.trial_type(~vis_inds,strcmpi(params.IVs,'light_bit'))) ;
        end
        nonans_bl = ~isnan(ppr_bl{ii}(:,1));  
        nonans_v = ~isnan(ppr_v{ii}(:,1));

        figure;
        h=histogram(ppr_bl{ii}(:,1),'binwidth',.25,'facecolor','b');
        ylim([0 max(histcounts(ppr_bl{ii}(:,1),'binwidth',.25))+1])
        xlim([1-ceil(max(ppr_bl{ii}(:,1)))+1 ceil(max(ppr_bl{ii}(:,1)))])
        yax = get(gca,'ylim');
        line([1 1],yax,'linestyle','--','color','k','linewidth',2)
        ylabel('Number of units','fontsize',24)
        xlabel('Spike count ratio: 2nd/1st pulse','fontsize',24)
        set(gca,'fontsize',18,'linewidth',2)
        print(gcf, '-dpng',sprintf('Paired pulse histogram %dHz (blanks)',lightconds(ii+2)))

        figure;
        h=histogram(ppr_v{ii}(:,1),'binwidth',.25,'facecolor','b');
        ylim([0 max(histcounts(ppr_v{ii}(:,1),'binwidth',.25))+1])
        xlim([1-ceil(max(ppr_v{ii}(:,1)))+1 ceil(max(ppr_v{ii}(:,1)))])
        yax = get(gca,'ylim');
        line([1 1],yax,'linestyle','--','color','k','linewidth',2)
        ylabel('Number of units','fontsize',24)
        xlabel('Spike count ratio: 2nd/1st pulse','fontsize',24)
        set(gca,'fontsize',18,'linewidth',2)
        print(gcf, '-dpng',sprintf('Paired pulse histogram %dHz (visual)',lightconds(ii+2)))

        figure;
        b=bar([1:lightconds(ii+2)-1],nanmedian(ppr_bl{ii}));
        ppiqr_bl = iqr(ppr_bl{ii}(nonans_bl,:))/2;
        set(b,'XData',[2:lightconds(ii+2)])
        xax = get(gca,'xlim');
        line(xax,[1 1],'linestyle','--','color','k','linewidth',2)
        hold on;
        for x=1:lightconds(ii+2)-1
    %         plot((x+1)*ones(1,sum(~isnan(ppr_bl(:,x)))),ppr_bl(~isnan(ppr_bl(:,x)),x),'k.')
            line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+ppiqr_bl(x)],'color','k')
        end
        ylabel('Spike count ratio (pulse x/pulse 1)','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        print(gcf, '-dpng',sprintf('Spike count ratio %dHz (blanks)',lightconds(ii+2)))

        figure;
        b=bar([1:lightconds(ii+2)-1],nanmedian(ppr_v{ii}));
        ppiqr_v = iqr(ppr_v{ii}(nonans_v,:))/2;
        set(b,'XData',[2:lightconds(ii+2)])
        xax = get(gca,'xlim');
        line(xax,[1 1],'linestyle','--','color','k','linewidth',2)
        hold on;
        for x=1:lightconds(ii+2)-1
    %         plot((x+1)*ones(1,sum(~isnan(ppr_v(:,x)))),ppr_v(~isnan(ppr_v(:,x)),x),'k.')
            line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+ppiqr_v(x)],'color','k')
        end
        ylabel('Spike count ratio (pulse x/pulse 1)','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        print(gcf, '-dpng',sprintf('Spike count ratio %dHz (visual)',lightconds(ii+2)))

        figure;
        b=bar([1:lightconds(ii+2)],nanmedian(resp_dur_bl{ii}));
        iqr_bl= iqr(resp_dur_bl{ii}(nonans_bl,:))/2;
        set(b,'XData',[1:lightconds(ii+2)])
        hold on;
        for x=1:lightconds(ii+2)
            line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+iqr_bl(x)],'color','k')
        end
        ylabel('Median # of significant bins','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        print(gcf, '-dpng',sprintf('Pulse response duration %dHz (blanks)',lightconds(ii+2)))

        figure;
        b=bar([1:lightconds(ii+2)],nanmedian(resp_dur_v{ii}));
        iqr_v = iqr(resp_dur_v{ii}(nonans_v,:))/2;
        set(b,'XData',[1:lightconds(ii+2)])
        hold on;
        for x=1:lightconds(ii+2)
            line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+iqr_v(x)],'color','k')
        end
        ylabel('Median # of significant bins','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        print(gcf, '-dpng',sprintf('Pulse response duration %dHz (visual)',lightconds(ii+2)))
    end
end
%% F1/F0 response analysis
binsize = .025;
pref_psth = nan(length(0:binsize:params.stimtime-binsize),length(conds)+1,length(clean_units));     % n timebins x num conds x  num units
Fratio = nan(length(clean_units),length(conds)+1);
for i = 1:length(clean_units)
%     all_light = params.all_light;
    vis_start = params.prestim*1000;       % in ms
    vis_end = params.poststim*1000;      % in ms
    which_trials = ismember(1:size(unitinfo(clean_units(i)).rast,1),prefori_trials{i});
    [~,tmp_psth] = make_psth_v2(binsize,0:binsize:(size(unitinfo(clean_units(i)).rast,2)-vis_start)/1000,which_trials,unitinfo(clean_units(i)).rast(:,vis_start+1:end-vis_end),zeros(size(unitinfo(clean_units(i)).rast,1)));
    pref_psth(:,:,i) = tmp_psth'-FRs(clean_units(i)).psthBlank(:,21:end)';       % subtract baseline!! (baseline from each light cond in order to look at whether light impacts F1/F0 independent of any gain change. does this make sense??)
    % but how should I handle psths with negative values?? (i.e. units
    % suppressed by vis stim)
    Fratio(i,:) = calc_F1F0(pref_psth(:,:,i),binsize,2);        % **currently hardcoded for 2Hz tfreq - will need to change!!
end

%% get preferred stimulus FR for each condition (but preferred stimulus defined in no light condition)
% durs = [params(:).light_dur];
% if contains(exp_path,'trains','ignorecase',1)
%     dur = max(durs(durs>0));
% else
%     dur = min(durs(durs>0));
% end
% light_times = round([max([params(:).av_light_start]) max([params(:).av_light_start])+dur].*1000);   % start with latest light start time, end after minimum light duration that isn't 0
light_times = [1000 2000];
FRpref = zeros(length(clean_units),length(lightconds));
for i = 1:length(clean_units)
    for ii = 1:length(lightconds)
        if ii == 1
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(prefori_trials{i},nolight_trials),light_times(1)+1:light_times(2)),2))/(diff(light_times)/1000);
        else
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(prefori_trials{i},light_trials{ii-1}),light_times(1)+1:light_times(2)),2))/(diff(light_times)/1000);
        end
    end
end

%%
cd(exp_dir)
% get different cell types
visual_cells = find((vis_sig < .025)|(vis_sig_ons < .025));
nonvisual_cells = find(~ismember(1:length(vis_sig),visual_cells));
% light_cells= find(min(light_sig,[],2)<.01);   % find cells with significant effect in any light condition
tuned_cells = find(tuned_sig(:,1) < .05);

onset_cells = find(vis_sig_ons<.05 & vis_sig>=.05);
sust_cells = find(vis_sig_ons<.05&vismod_on'<0&vis_sig<.05&vismod'<0 | vis_sig_ons<.05&vismod_on'>0&vis_sig<.05&vismod'>0);
sust_act_cells =  find(vis_sig_ons<.05&vismod_on'>0&vis_sig<.05&vismod'>0);
sust_sup_cells = find(vis_sig_ons<.05&vismod_on'<0&vis_sig<.05&vismod'<0);
rev_cells = find(vis_sig_ons<.05&vismod_on'<0&vis_sig<.05&vismod'>0 | vis_sig_ons<.05&vismod_on'>0&vis_sig<.05&vismod'<0);
delay_act_cells = find(vis_sig_ons>=.05 & vis_sig < .05 & vismod'>0);
delay_sup_cells = find(vis_sig_ons>=.05 & vis_sig < .05 & vismod'<0);
trans_cells = find(vis_sig_ons<.05 & vismod_on'>0 & vis_sig>=.05 | vis_sig_ons<.05 & vismod_on'>0 & vis_sig < .05 & vismod'<0);

% orthFR_delta(isnan(orthFR_delta)) = 0;
% prefFR_delta(isnan(prefFR_delta)) = 0;
% 
reg_cells = find(t2p_t>=.4);
FS_cells = find(t2p_t<.4);

% % find diff cell types
% supp_cells = intersect(find((light_sig(1,:)<.05)&(light_sig(2,:)<.05)&(light_sig(3,:)<.05)),find(sum(sign(lightmod),2)==-3));
% % significantly light suppressed in all conditions
% enh_cells = intersect(find((light_sig(1,:)<.05)&(light_sig(2,:)<.05)&(light_sig(3,:)<.05)),find(sum(sign(lightmod),2)==3));
% % significantly light activated in all conditions
% complex_cells = intersect(find(min(light_sig)<.05),find(abs(sum(sign(lightmod),2))< 3));
% % significantly affected by light in at least one condition, but direction
% % of light modulation depends on light intensity
% all_supp = intersect(find((light_sig(1,:)<.05)&(light_sig(2,:)<.05)&(light_sig(3,:)<.05)),intersect(find(sum(sign(lightmod_early),2)==-3),find(sum(sign(lightmod_late),2)==-3)));
% % significantly light modulated in all light conditions, and late and early
% % periods are suppressed in all conditions
% all_enh = intersect(find((light_sig(1,:)<.05)&(light_sig(2,:)<.05)&(light_sig(3,:)<.05)),intersect(find(sum(sign(lightmod_early),2)==3),find(sum(sign(lightmod_late),2)==3)));
% % significantly light modulated in all light conditions, and late and early
% % periods are both enhanced in all conditions
% quick_enh = intersect(find(light_sig(3,:)<.05),find((lightmod_onset(:,3)>0)));
% % significantly light modulated in high light condition; onset is enhanced
% onthensupp = intersect(quick_enh,find((lightmod_late(:,3)<0)));
% % significantly light modulated in high light condition; onset is enhanced
% % and late period is suppressed (in high light condition)
% onthenenh = intersect(quick_enh,find((lightmod_late(:,3)>0)&(lightmod_early(:,3)<0)));
% % significantly light modulated in high light condition; onset is enhanced,
% % early period is suppressed but late period is enhanced (in high light condition)
% lowenh_highsupp = intersect(find((light_sig(1,:)<.05)&(light_sig(3,:)<.05)),find((sign(lightmod_late(:,3))==-1)&(sign(lightmod_late(:,1))==1)));
% % significantly light modulated in low and high light conditions, but
% % enhanced in low light while suppressed in high light (considering late
% % period)
% lowsupp_highenh = intersect(find((light_sig(1,:)<.05)&(light_sig(3,:)<.05)),find((sign(lightmod_late(:,3))+sign(lightmod_early(:,3))==2)&(sign(lightmod_late(:,1))+sign(lightmod_early(:,1))==-2)));
% % significantly light modulated in low and high light conditions, but
% % suppressed in low light while enhanced in high light (considering early and late
% % periods)
% delayact = intersect(find(min(light_sig)<.05),find((lightmod_early(:,3)>0)&(sign(lightmod_onset(:,2))+sign(lightmod_onset(:,3))<0)));
% % significantly light modulated in any light cond; suppressed at onset (med & high conds) and then
% % activated during early period

%% test significance of overall light modulation and OSI/DSI change using Wilcoxen signed-rank
% lightsig_all_low = signrank(FRev(:,1),FRev(:,2));
% lightsig_all_low_dir = sign(nanmedian(FRev(:,2))-nanmedian(FRev(:,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_all_high = signrank(FRev(:,1),FRev(:,end));
% lightsig_all_high_dir = sign(nanmedian(FRev(:,end))-nanmedian(FRev(:,1)));
% lightsig_bl_low = signrank(FRbl(:,1),FRbl(:,2));
% lightsig_bl_low_dir = sign(nanmedian(FRbl(:,2))-nanmedian(FRbl(:,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_bl_high = signrank(FRbl(:,1),FRbl(:,end));
% lightsig_bl_high_dir = sign(nanmedian(FRbl(:,end))-nanmedian(FRbl(:,1)));
% lightsig_vis_low = signrank(FRev(visual_cells,1),FRev(visual_cells,2));
% lightsig_vis_low_dir = sign(nanmedian(FRev(visual_cells,2))-nanmedian(FRev(visual_cells,1)));
% lightsig_vis_high = signrank(FRev(visual_cells,1),FRev(visual_cells,end));
% lightsig_vis_high_dir = sign(nanmedian(FRev(visual_cells,end))-nanmedian(FRev(visual_cells,1)));
% lightsig_blvis_low = signrank(FRbl(visual_cells,1),FRbl(visual_cells,2));
% lightsig_blvis_low_dir = sign(nanmedian(FRbl(visual_cells,2))-nanmedian(FRbl(visual_cells,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_blvis_high = signrank(FRbl(visual_cells,1),FRbl(visual_cells,end));
% lightsig_blvis_high_dir = sign(nanmedian(FRbl(visual_cells,end))-nanmedian(FRbl(visual_cells,1)));
% lightsig_nonvis_low = signrank(FRev(nonvisual_cells,1),FRev(nonvisual_cells,2));
% lightsig_nonvis_low_dir = sign(nanmedian(FRev(nonvisual_cells,2))-nanmedian(FRev(nonvisual_cells,1)));
% lightsig_nonvis_high = signrank(FRev(nonvisual_cells,1),FRev(nonvisual_cells,end));
% lightsig_nonvis_high_dir = sign(nanmedian(FRev(nonvisual_cells,end))-nanmedian(FRev(nonvisual_cells,1)));
% lightsig_blnonvis_low = signrank(FRbl(nonvisual_cells,1),FRbl(nonvisual_cells,2));
% lightsig_blnonvis_low_dir = sign(nanmedian(FRbl(nonvisual_cells,2))-nanmedian(FRbl(nonvisual_cells,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_blnonvis_high = signrank(FRbl(nonvisual_cells,1),FRbl(nonvisual_cells,end));
% lightsig_blnonvis_high_dir = sign(nanmedian(FRbl(nonvisual_cells,end))-nanmedian(FRbl(nonvisual_cells,1)));
% lightsig_vispref_low = signrank(FRpref(:,1),FRpref(:,2));
% lightsig_vispref_low_dir = sign(nanmedian(FRpref(:,2))-nanmedian(FRpref(:,1)));
% lightsig_vispref_high = signrank(FRpref(:,1),FRpref(:,end));
% lightsig_vispref_high_dir = sign(nanmedian(FRpref(:,end))-nanmedian(FRpref(:,1)));
% lightsig_visFRdelt_low = signrank(FRev(:,1)-FRbl(:,1),FRev(:,2)-FRbl(:,2));
% lightsig_visFRdelt_low_dir = sign(nanmedian(FRev(:,2)-FRbl(:,2))-nanmedian(FRev(:,1)-FRbl(:,1)));
% lightsig_visFRdelt_high = signrank(FRev(:,1)-FRbl(:,1),FRev(:,end)-FRbl(:,end));
% lightsig_visFRdelt_high_dir = sign(nanmedian(FRev(:,end)-FRbl(:,end))-nanmedian(FRev(:,1)-FRbl(:,1)));
% lightsig_prefFRdelt_low_pref = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,2)-FRbl(:,2));
% lightsig_prefFRdelt_low_pref_dir = sign(nanmedian(FRpref(:,2)-FRbl(:,2))-nanmedian(FRpref(:,1)-FRbl(:,1)));
% lightsig_prefFRdelt_high = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,end)-FRbl(:,end));
% lightsig_prefFRdelt_high_dir = sign(nanmedian(FRpref(:,end)-FRbl(:,end))-nanmedian(FRpref(:,1)-FRbl(:,1)));
% 
% lightsig_vals = [lightsig_all_low lightsig_all_low_dir; lightsig_all_high lightsig_all_high_dir; lightsig_bl_low lightsig_bl_low_dir; lightsig_bl_high lightsig_bl_high_dir;...
%     lightsig_vis_low lightsig_vis_low_dir; lightsig_vis_high lightsig_vis_high_dir; lightsig_blvis_low lightsig_blvis_low_dir; lightsig_blvis_high lightsig_blvis_high_dir;...
%     lightsig_nonvis_low lightsig_nonvis_low_dir; lightsig_nonvis_high lightsig_nonvis_high_dir; lightsig_blnonvis_low lightsig_blnonvis_low_dir; lightsig_blnonvis_high lightsig_blnonvis_high_dir;...
%     lightsig_vispref_low lightsig_vispref_low_dir; lightsig_vispref_high lightsig_vispref_high_dir; lightsig_visFRdelt_low lightsig_visFRdelt_low_dir; lightsig_visFRdelt_high lightsig_visFRdelt_high_dir;...
%     lightsig_prefFRdelt_low_pref lightsig_prefFRdelt_low_pref_dir; lightsig_prefFRdelt_high lightsig_prefFRdelt_high_dir];
% 
% osiCVsig_tuned_low = signrank(OSI_CV(tuned_cells,1),OSI_CV(tuned_cells,2));
% osiCVsig_tuned_low_dir = sign(nanmedian(OSI_CV(tuned_cells,2))-nanmedian(OSI_CV(tuned_cells,1)));
% osiCVsig_tuned_high = signrank(OSI_CV(tuned_cells,1),OSI_CV(tuned_cells,end));
% osiCVsig_tuned_high_dir = sign(nanmedian(OSI_CV(tuned_cells,end))-nanmedian(OSI_CV(tuned_cells,1)));
% osisig_tuned_low = signrank(OSI(tuned_cells,1),OSI(tuned_cells,2));
% osisig_tuned_low_dir = sign(nanmedian(OSI(tuned_cells,2))-nanmedian(OSI(tuned_cells,1)));
% osisig_tuned_high = signrank(OSI(tuned_cells,1),OSI(tuned_cells,end));
% osisig_tuned_high_dir = sign(nanmedian(OSI(tuned_cells,end))-nanmedian(OSI(tuned_cells,1)));
% dsiCVsig_tuned_low = signrank(DSI_CV(tuned_cells,1),DSI_CV(tuned_cells,2));
% dsiCVsig_tuned_low_dir = sign(nanmedian(DSI_CV(tuned_cells,2))-nanmedian(DSI_CV(tuned_cells,1)));
% dsiCVsig_tuned_high = signrank(DSI_CV(tuned_cells,1),DSI_CV(tuned_cells,end));
% dsiCVsig_tuned_high_dir = sign(nanmedian(DSI_CV(tuned_cells,end))-nanmedian(DSI_CV(tuned_cells,1)));
% dsisig_tuned_low = signrank(DSI(tuned_cells,1),DSI(tuned_cells,2));
% dsisig_tuned_low_dir = sign(nanmedian(DSI(tuned_cells,2))-nanmedian(DSI(tuned_cells,1)));
% dsisig_tuned_high = signrank(DSI(tuned_cells,1),DSI(tuned_cells,end));
% dsisig_tuned_high_dir = sign(nanmedian(DSI(tuned_cells,end))-nanmedian(DSI(tuned_cells,1)));
% osiCVsig_vis_low = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,2));
% osiCVsig_vis_low_dir = sign(nanmedian(OSI_CV(visual_cells,2))-nanmedian(OSI_CV(visual_cells,1)));
% osiCVsig_vis_high = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,end));
% osiCVsig_vis_high_dir = sign(nanmedian(OSI_CV(visual_cells,end))-nanmedian(OSI_CV(visual_cells,1)));
% osisig_vis_low = signrank(OSI(visual_cells,1),OSI(visual_cells,2));
% osisig_vis_low_dir = sign(nanmedian(OSI(visual_cells,2))-nanmedian(OSI(visual_cells,1)));
% osisig_vis_high = signrank(OSI(visual_cells,1),OSI(visual_cells,end));
% osisig_vis_high_dir = sign(nanmedian(OSI(visual_cells,end))-nanmedian(OSI(visual_cells,1)));
% dsiCVsig_vis_low = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,2));
% dsiCVsig_vis_low_dir = sign(nanmedian(DSI_CV(visual_cells,2))-nanmedian(DSI_CV(visual_cells,1)));
% dsiCVsig_vis_high = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,end));
% dsiCVsig_vis_high_dir = sign(nanmedian(DSI_CV(visual_cells,end))-nanmedian(DSI_CV(visual_cells,1)));
% dsisig_vis_low = signrank(DSI(visual_cells,1),DSI(visual_cells,2));
% dsisig_vis_low_dir = sign(nanmedian(DSI(visual_cells,2))-nanmedian(DSI(visual_cells,1)));
% dsisig_vis_high = signrank(DSI(visual_cells,1),DSI(visual_cells,end));
% dsisig_vis_high_dir = sign(nanmedian(DSI(visual_cells,end))-nanmedian(DSI(visual_cells,1)));
% tuningsig_vals = [osiCVsig_tuned_low osiCVsig_tuned_low_dir; osiCVsig_tuned_high osiCVsig_tuned_high_dir; osisig_tuned_low osisig_tuned_low_dir;...
%     osisig_tuned_high osisig_tuned_high_dir; dsiCVsig_tuned_low dsiCVsig_tuned_low_dir; dsiCVsig_tuned_high dsiCVsig_tuned_high_dir;...
%     dsisig_tuned_low dsisig_tuned_low_dir; dsisig_tuned_high dsisig_tuned_high_dir; osiCVsig_vis_low osiCVsig_vis_low_dir;...
%     osiCVsig_vis_high osiCVsig_vis_high_dir; osisig_vis_low osisig_vis_low_dir; osisig_vis_high osisig_vis_high_dir;...
%     dsiCVsig_vis_low dsiCVsig_vis_low_dir; dsiCVsig_vis_high dsiCVsig_vis_high_dir; dsisig_vis_low dsisig_vis_low_dir; dsisig_vis_high dsisig_vis_high_dir];

%%
% vis_type = nan(1,length(clean_units));
% vis_type(trans_cells) = 1;
% vis_type(sust_act_cells) = 2;
% vis_type(delay_act_cells) = 3;
% vis_type(delay_sup_cells) = 4;
% % barplot_by_layer(lightmod(visual_cells,end)',vis_type(visual_cells),ones(1,length(visual_cells)),'Light modulation index','Lightmodulation (evoked)')
% vis_types = unique(vis_type(~isnan(vis_type)));
% for n = 1:length(vis_types)
%     mean_data(n,:) = nanmean(lightmod(vis_type==vis_types(n),:),1);
%     se_data(n,:) = nanstd(lightmod(vis_type==vis_types(n),:),[],1)./sqrt(size(lightmod(vis_type==vis_types(n),:),1));
% end
% fig = figure;
% h = bar(mean_data);
% hold on;
% 
% numbars = size(mean_data, 1);
% numconds = size(mean_data, 2);
% groupwidth = min(0.8, numconds/(numconds+1.5));
% for i = 1:numconds
% % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% x(i,:) = (1:numbars) - groupwidth/2 + (2*i-1) * groupwidth / (2*numconds); % Aligning error bar with individual bar
% errorbar(x(i,:), mean_data(:,i),se_data(:,i), se_data(:,i), 'k', 'linestyle', 'none');
% end
% set(get(gca,'YLabel'),'String','Light modulation index','Fontsize',24)
% set(gca,'XTicklabel','Transient| Sustained| Delay-act| Delay-supp','Fontsize',18)
% color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2];% for graphing purposes (first is black, last is green)
% for i = 1:length(h)
%     set(h(i),'FaceColor',color_mat(i+1,:),'EdgeColor',color_mat(i+1,:));
% end
% type_ind = zeros(numconds,sum(~isnan(vis_type)));
% vis_type_clean = vis_type(~isnan(vis_type));
% for i = 1:numconds
%     type_ind(i,vis_type_clean==1) = x(i,1);
%     type_ind(i,vis_type_clean==2) = x(i,2);
%     type_ind(i,vis_type_clean==3) = x(i,3);
%     type_ind(i,vis_type_clean==4) = x(i,4);
% end
% plot(type_ind,lightmod(~isnan(vis_type),:)','k.','MarkerSize',18)
% print(gcf, '-dpng','lightmodbyvisresp')

%% plot distance from bottom of LP by lightmod
figure;
subplot(111)
plot(vismod(:,end),distfromlastch','.','MarkerSize',24)
% hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
% h = get(gca,'ytick');
% set(gca,'yticklabel',h*25);
xlim([-1 1])
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--')
% legend('low','high')
xlabel('Visual modulation index ','Fontsize',16)
ylabel(strcat('Distance from bottom of LP (in um)'),'Fontsize',16)
print(gcf, '-dpng','lightmodbydepth_bottom')

figure;
subplot(111)
plot(vismod(:,end),abs(distfromfirstch)','.','MarkerSize',24)
% hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
% h = get(gca,'ytick');
% set(gca,'yticklabel',h*25);
view(0,270)
xlim([-1 1])
ylim([0 max(abs(distfromfirstch))])
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
% legend('low','high')
xlabel('Visual modulation index ','Fontsize',24)
ylabel(strcat('Distance from top of LP (by ch)'),'Fontsize',24)
set(gca,'fontsize',18,'linewidth',2)
print(gcf, '-dpng','lightmodbydepth_top')

% trying by shank...
shank = unique(shk(clean_units));    % this is INCORRECT - shk is only from last experiment, not clean units
dist = unique(distfromfirstch);
for i = 1:length(dist)
    for sh = 1:length(shank)
        mean_lm(i,sh) = nanmean(vismod((distfromfirstch==dist(i)&shk(clean_units)==shank(sh)),end));
    end
end
mean_lm(isnan(mean_lm)) = 0;
figure;
for i = 1:size(mean_lm,2)
    subplot(1,size(mean_lm,2),i)
    bar(unique(distfromfirstch),mean_lm(:,i))
    hold on
    plot(distfromfirstch(shk(clean_units)==shank(i))',vismod(shk(clean_units)==shank(i),end),'.','color',[0 .8 .7])
    view(90,90)
    ylim([-1 1])
% %     xlim([0 max(abs(distfromfirstch))])
%     h = get(gca,'xtick');
%     set(gca,'xticklabel',h*25);
    title(sprintf('shank%d',i))
    ylabel('Visual modulation index ','Fontsize',12)
    xlabel(strcat('Depth in LP (in um)'),'Fontsize',12)
end
    print(gcf, '-dpng','lightmodbyshank')


% %%
if contains(exp_path,'halo','ignorecase',1)
    color_mat = [0 0 0; .9 0 .3; 0.6350, 0.0780, 0.1840]; % for graphing purposes (first is black, last is green)
    nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];      % make red for halo
else
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4]; % % for lighter-shade dots
    nonvis_color_mat = [.5 .5 .5; .75 .8 1; .7 .8 .7];  % for lighter-shade dots
end
% % % FR light vs no light, by shank
% % plot_scatter(FRev(:,[1 2]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byshk', {'Most medial','','','Most lateral'}, 1)     % first lightcond pwr
% % plot_scatter(FRev(:,[1 end]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byshk', {'Most medial','','','Most lateral'}, 1)   % last lightcond pwr
% % % FR light vs no light - HIGH pwr, light onset, by shank
% % plot_scatter(FRonset(:,[1 end]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FRonset_high_byshk', {'Most medial','','','Most lateral'}, 1)
% 
% % FR light vs no light - visual vs nonvisual units
% plot_scatter(FRev(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRev(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byvis', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% % if strcmpi(exp_type,'trains')
% %     plot_scatter(FRev(:,[1 end-1]), (vis_sig < .05)|(vis_sig_ons < .05), {[.7 .8 .7],[.7 0 1]}, 'Spks/s (light OFF)', 'Spks/s (light ON - med)', 'FR_med_byvis', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% % end
% 
% % FR light vs no light - visual vs nonvisual units, blank trials
% plot_scatter(FRbl(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_bl', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRbl(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byvis_bl', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% 
% plot_scatter(FRvison(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_onset', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRvison(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byvis_onset', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% 
% % FR light vs no light - visual vs nonvisual units, preferred trials
% plot_scatter(FRpref(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRpref(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byvis_pref', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% 
% % change in preferred FR light vs no light - visual vs nonvisual units
% plot_scatter(FRpref(:,[1 2])-FRbl(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-low)', 'FR_low_vischange_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRpref(:,[1 end])-FRbl(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked change in FR (Spks/s-light ON-high)', 'FR_high_vischange_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% 
% % change in evoked FR light vs no light - visual vs nonvisual units
% plot_scatter(FRev(:,[1 2])-FRbl(:,[1 2]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(2,:),color_mat(2,:)}, 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-low)', 'FR_low_vischange', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRev(:,[1 end])-FRbl(:,[1 end]), (vis_sig < .05)|(vis_sig_ons < .05), {nonvis_color_mat(3,:),color_mat(3,:)}, 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-high)', 'FR_high_vischange', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% 

%% average PSTHs
% color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4]; % for graphing purposes (first is black, last is green)
psthV(:,:,:) = reshape(cell2mat(arrayfun(@(x) x.psthVisual([conds end],:), FRs(clean_units),'uniformoutput',0)),length(conds)+1,size(FRs(1).psthVisual,2),length(clean_units));
visUp = clean_units(vismod(visual_cells)>0);
visDown = clean_units(vismod(visual_cells)<0);

% test_psth = psthV(:,:,visUp)./repmat(repmat(max(max(psthV(:,:,visUp))),size(psthV,1),1),1,size(psthV,2));
% test_mean = mean(test_psth,3);
% test_se = std(test_psth,0,3)./sqrt(size(test_psth,3));
% pop_fig = figure;
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([0 1])
% shadedErrorBar([-.475:.025:2],test_mean(1,:), test_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
% hold on;
% shadedErrorBar([-.475:.025:2],test_mean(4,:), test_se(4,:), {'Color', [0 .8 1],'linewidth',2},1);
% line([0 0],[0 1],'Color','r','LineStyle','--')
% 
% test_psth = psthV(:,:,sust_sup_cells)./repmat(repmat(max(max(psthV(:,:,sust_sup_cells))),size(psthV,1),1),1,size(psthV,2));
% test_mean = mean(test_psth,3);
% test_se = std(test_psth,0,3)./sqrt(size(test_psth,3));
% pop_fig2= figure;
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([0 1])
% shadedErrorBar([-.475:.025:2],test_mean(1,:), test_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
% hold on;
% shadedErrorBar([-.475:.025:2],test_mean(4,:), test_se(4,:), {'Color', [0 .8 1],'linewidth',2},1);
% line([0 0],[0 1],'Color','r','LineStyle','--')

rep_baseline = reshape(repmat(repmat(FRb,size(psthV,2),1),size(psthV,1),1),size(psthV,1),size(psthV,2),size(psthV,3));
bs_psth = psthV-rep_baseline;       % baseline-subtracted psth
% test_psth = bs_psth./repmat(repmat(max(max(bs_psth)),size(bs_psth,1),1),1,size(bs_psth,2));     % normalize
psthZ = zscore(bs_psth,0,2);        % zscore across timepoints
mean_bs = repmat(mean(bs_psth(:,1:20,:),2),1,100,1);    % prestim currently hardcoded! 20 time bins x 25ms each = 500ms
std_bs = repmat(std(bs_psth(:,1:20,:),[],2),1,100,1);
psthZ = (bs_psth-mean_bs)./std_bs;      % setting 0 to the mean during prestim period only

% mean_visUp = mean(psthZ(:,:,clean_units(visUp)),3);
test_mean = mean(psthZ,3);
test_se = std(psthZ,0,3)./sqrt(size(psthZ,3));
pop_fig3= figure;
xlim([-.475 2])   % ticks mark the END of 25ms bins
hold on;
for i = 1:size(psthZ,1)
    shadedErrorBar([-.5:.025:1.975],test_mean(i,:), test_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
% ylim([-.5 .5])
yax = get(gca,'YLim');
% yax = [-min(abs(yax)) min(abs(yax))];
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
% start_times = round([params(:).av_light_start],2)-unique([params(:).prestim]);  % round to nearest hundredth
% stim_durs = round([params(:).light_dur],2);
% stim_durs = stim_durs(stim_durs>0);
% if length(unique(start_times))>1
%     for ii = 1:length(unique(start_times))
%         line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
%         line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
%     end
% else
%     patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% end
ylim(yax)
xlabel('Time from visual stimulus onset (sec)','fontsize',24)
ylabel('Z-scored firing rate (Spks/s)','fontsize',24)
set(gca,'fontsize',18,'linewidth',2);
save_fig_name = 'PopulationPSTH_zscore';
print(pop_fig3,'-dpng',save_fig_name)
print2eps(save_fig_name,pop_fig3)

% % population psths: onset, sustained, rev, delay activated, delay suppressed
% pop_psth = figure;
% subplot(141)
% ons_mean = mean(test_psth(:,:,trans_cells),3);
% ons_se = std(test_psth(:,:,trans_cells),0,3)./sqrt(length(trans_cells));
% hold on;
% for i = 1:size(test_psth,1)
%     shadedErrorBar([-.475:.025:2],ons_mean(i,:), ons_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
% end
% title('Transient cells','fontsize',14);
% xlabel('Time from visual stim (sec)','fontsize',14)
% ylabel('Normalized firing rate (Spks/s)','fontsize',14)
% set(gca,'fontsize',14);
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([-1 1])
% yax = get(gca,'YLim');
% patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% line([0 0],yax,'Color','r','LineStyle','--')
% 
% subplot(142)
% sust_mean = mean(test_psth(:,:,sust_cells),3);
% sust_se = std(test_psth(:,:,sust_cells),0,3)./sqrt(length(sust_cells));
% hold on;
% for i = 1:size(test_psth,1)
%     shadedErrorBar([-.475:.025:2],sust_mean(i,:), sust_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
% end
% title('Sustained response cells','fontsize',14);
% xlabel('Time from visual stim (sec)','fontsize',14)
% ylabel('Normalized firing rate (Spks/s)','fontsize',14)
% set(gca,'fontsize',14);
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([-1 1])
% yax = get(gca,'YLim');
% patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% line([0 0],yax,'Color','r','LineStyle','--')
% 
% % subplot(153)
% % rev_mean = mean(test_psth(:,:,rev_cells),3);
% % rev_se = std(test_psth(:,:,rev_cells),0,3)./sqrt(length(rev_cells));
% % shadedErrorBar([-.475:.025:2],rev_mean(1,:), rev_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
% % hold on;
% % shadedErrorBar([-.475:.025:2],rev_mean(4,:), rev_se(4,:), {'Color', [0 .8 1],'linewidth',2},1);
% % title('Reverse response cells','fontsize',14);
% % xlabel('Time from visual stim (sec)','fontsize',14)
% % ylabel('Normalized firing rate (Spks/s)','fontsize',14)
% % set(gca,'fontsize',14);
% % xlim([-.475 2])   % ticks mark the END of 25ms bins
% % ylim([0 1])
% % yax = get(gca,'YLim');
% % patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% % line([0 0],[0 1],'Color','r','LineStyle','--')
% 
% subplot(143)
% delay_act_mean = mean(test_psth(:,:,delay_act_cells),3);
% delay_act_se = std(test_psth(:,:,delay_act_cells),0,3)./sqrt(length(delay_act_cells));
% hold on;
% for i = 1:size(test_psth,1)
%     shadedErrorBar([-.475:.025:2],delay_act_mean(i,:), delay_act_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
% end
% title('Delay activated cells','fontsize',14);
% xlabel('Time from visual stim (sec)','fontsize',14)
% ylabel('Normalized firing rate (Spks/s)','fontsize',14)
% set(gca,'fontsize',14);
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([-1 1])
% yax = get(gca,'YLim');
% patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% line([0 0],yax,'Color','r','LineStyle','--')
% 
% subplot(144)
% delay_sup_mean = mean(test_psth(:,:,delay_sup_cells),3);
% delay_sup_se = std(test_psth(:,:,delay_sup_cells),0,3)./sqrt(length(delay_sup_cells));
% hold on;
% for i = 1:size(test_psth,1)
%     shadedErrorBar([-.475:.025:2],delay_sup_mean(i,:), delay_sup_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
% end
% title('Delay suppressed cells','fontsize',14);
% xlabel('Time from visual stim (sec)','fontsize',14)
% ylabel('Normalized firing rate (Spks/s)','fontsize',14)
% set(gca,'fontsize',14);
% xlim([-.475 2])   % ticks mark the END of 25ms bins
% ylim([-1 1])
% yax = get(gca,'YLim');
% patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
% xSize =30; ySize = 11;
%     xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%     set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
%     set(gcf,'Position',[0 0 xSize*50 ySize*50])
% line([0 0],yax,'Color','r','LineStyle','--')
% print(pop_psth, '-dpng','PopulationPSTH_subtypes')
% print2eps('PopulationPSTH_subtypes',pop_psth)
% 
% test_v = psthV(:,21:40,:)./repmat(repmat(max(max(psthV(:,21:40,:))),size(psthV,1),1),1,length(21:40));
% test_l = psthV(3:4,40:59,:)./repmat(repmat(max(max(psthV(3:4,40:59,:))),2,1),1,length(40:59));

% %% boxplot
% figure;
% subplot(141)
% boxplot(lightmod_onset,'Labels',{'low','med','high'})
% hold on
% for i = 1:size(lightmod_onset,2)
%     plot(i*ones(1,length(reg_cells)),lightmod_onset(reg_cells,i),'.')
%     plot(i*ones(1,length(FS_cells)),lightmod_onset(FS_cells,i),'r.')
% end
% ylim([-1 1])
% set(gca,'fontsize',16)
% set(findobj(gca,'Type','text'),'FontSize',14)
% xlabel('Light power','Fontsize',18)
% ylabel('Light modulation index','Fontsize',18)
% title('0-100ms post-light onset','fontsize',18)
% subplot(142)
% boxplot(lightmod_early,'Labels',{'low','med','high'})
% hold on
% for i = 1:size(lightmod_early,2)
%     plot(i*ones(1,length(reg_cells)),lightmod_onset(reg_cells,i),'.')
%     plot(i*ones(1,length(FS_cells)),lightmod_onset(FS_cells,i),'r.')
% end
% ylim([-1 1])
% set(gca,'fontsize',16)
% set(findobj(gca,'Type','text'),'FontSize',14)
% xlabel('Light power','Fontsize',18)
% % ylabel('Light modulation index','Fontsize',18)
% title('100-500ms post-light onset','fontsize',18)
% subplot(143)
% boxplot(lightmod_late,'Labels',{'low','med','high'})
% hold on
% for i = 1:size(lightmod_late,2)
%     plot(i*ones(1,length(reg_cells)),lightmod_onset(reg_cells,i),'.')
%     plot(i*ones(1,length(FS_cells)),lightmod_onset(FS_cells,i),'r.')
% end
% ylim([-1 1])
% set(gca,'fontsize',16)
% set(findobj(gca,'Type','text'),'FontSize',14)
% xlabel('Light power','Fontsize',18)
% % ylabel('Light modulation index','Fontsize',18)
% title('500-1000ms post-light onset','fontsize',18)
% subplot(144)
% boxplot(lightmod,'Labels',{'low','med','high'})
% hold on
% for i = 1:size(lightmod,2)
%     plot(i*ones(1,length(reg_cells)),lightmod_onset(reg_cells,i),'.')
%     plot(i*ones(1,length(FS_cells)),lightmod_onset(FS_cells,i),'r.')
% end
% ylim([-1 1])
% set(gca,'fontsize',16)
% set(findobj(gca,'Type','text'),'FontSize',14)
% xlabel('Light power','Fontsize',18)
% % ylabel('Light modulation index','Fontsize',18)
% title('Overall','fontsize',18)
% 
% set(gcf,'Position', [100, 100, 1500, 500]);
% drawnow
% export_fig ('lightmod_boxplots', '-png','-r600','-zbuffer');
% 
% %% OSI 
% plot_scatter(OSI_CV(tuned_cells,[1 end]), ones(1,length(tuned_cells)), {'k'}, 'OSI(CV) (light OFF)', 'OSI(CV) (light ON - high)', 'OSI(CV)', {'Tuned cells'}, 1)     % first lightcond pwr
% plot_scatter(OSI(tuned_cells,[1 end]), ones(1,length(tuned_cells)), {'k'}, 'OSI (light OFF)', 'OSI (light ON - high)', 'OSI', {'Tuned cells'}, 1)     % first lightcond pwr
% plot_scatter(DSI_CV(tuned_cells,[1 end]), ones(1,length(tuned_cells)), {'k'}, 'DSI(CV) (light OFF)', 'DSI(CV) (light ON - high)', 'DSI(CV)', {'Tuned cells'}, 1)     % first lightcond pwr
% plot_scatter(DSI(tuned_cells,[1 end]), ones(1,length(tuned_cells)), {'k'}, 'DSI (light OFF)', 'DSI (light ON - high)', 'DSI', {'Tuned cells'}, 1)     % first lightcond pwr
% 
% plot_scatter(OSI_CV(tuned_cells,[1 end]), lightmod(tuned_cells,end)>0, {'k','b'}, 'OSI(CV) (light OFF)', 'OSI(CV) (light ON - high)', 'OSI(CV)_bylightmod', {'Light-suppressed','Light-enhanced'}, 1)     % first lightcond pwr
% plot_scatter(OSI(tuned_cells,[1 end]), lightmod(tuned_cells,end)>0, {'k','b'}, 'OSI (light OFF)', 'OSI (light ON - high)', 'OSI_bylightmod', {'Light-suppressed','Light-enhanced'}, 1)     % first lightcond pwr
% plot_scatter(DSI_CV(tuned_cells,[1 end]), lightmod(tuned_cells,end)>0, {'k','b'}, 'DSI(CV) (light OFF)', 'DSI(CV) (light ON - high)', 'DSI(CV)_bylightmod', {'Light-suppressed','Light-enhanced'}, 1)     % first lightcond pwr
% plot_scatter(DSI(tuned_cells,[1 end]), lightmod(tuned_cells,end)>0 , {'k','b'}, 'DSI (light OFF)', 'DSI (light ON - high)', 'DSI_bylightmod', {'Light-suppressed','Light-enhanced'}, 1)     % first lightcond pwr

% %%
% osi_fig = figure('name','OSI Change (light-nolight) vs. Light modulation');
% plot(diff(OSI(:,[1 4]),[],2),lightmod(:,3),'r.','MarkerSize',24)
% hold on
% xlabel('change in OSI CV (light - no light)','Fontsize',24)
% ylabel('Light modulation index','Fontsize',24)
% line([0 0],[-1 1],'Color','k')
% line([-1 1],[0 0],'Color','k')
% set(gca,'fontsize',18)
% set(l,'fontsize',18)
% 
% %% additive vs multiplicative changes
% 
% ptrg
% figure;
% plot(orthFR_delta(tuned_cells),prefFR_delta(tuned_cells),'.','MarkerSize',24)
% hold on
% xmax = max(max(abs(orthFR_delta(tuned_cells))),5);
% ymax = max(max(abs(prefFR_delta(tuned_cells))),5);
% xlim([-xmax xmax])
% ylim([-ymax ymax])
% xax = get(gca,'XLim');
% yax = get(gca,'YLim');
% x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
% y=x;
% plot(x,y,'k--');
% line(xax,[0 0],'Color','k')
% line([0,0],yax,'Color','k')
% xlabel('% change in orthogonal FR','fontsize',24)
% ylabel('% change in preferred FR','fontsize',24)
% h = get(gca,'ytick');
% f = get(gca,'xtick');
% set(gca,'yticklabel',h*100);
% set(gca,'xticklabel',f*100);
% plot(nanmedian(orthFR_delta(tuned_cells)),nanmedian(prefFR_delta(tuned_cells)),'k+','MarkerSize',24)
% print(gcf, '-dpng','FRchange_orthbypref_tuned')
% 
% figure;
% plot(orthFR_delta,prefFR_delta,'.','MarkerSize',24)
% hold on
% xmax = max(max(abs(orthFR_delta)),5);
% ymax = max(max(abs(prefFR_delta)),5);
% xlim([-xmax xmax])
% ylim([-ymax ymax])
% xax = get(gca,'XLim');
% yax = get(gca,'YLim');
% x = [min(xax(1),yax(1)):1:max(xax(2),yax(2))];
% y=x;
% plot(x,y,'k--');
% line(xax,[0 0],'Color','k')
% line([0,0],yax,'Color','k')
% xlabel('% change in orthogonal FR','fontsize',24)
% ylabel('% change in preferred FR','fontsize',24)
% h = get(gca,'ytick');
% f = get(gca,'xtick');
% set(gca,'yticklabel',h*100);
% set(gca,'xticklabel',f*100);
% plot(nanmedian(orthFR_delta),nanmedian(prefFR_delta),'k+','MarkerSize',24)
% print(gcf, '-dpng','FRchange_orthbypref_all')

%%
num_units = length(clean_units);
fileID = fopen('experiment_summary.txt','w');
fprintf(fileID,'Number of visually responsive cells: %d of %d\r\n',length(visual_cells),num_units);
fprintf(fileID,'Percent visually responsive: %.2f\r\n',100*length(visual_cells)/num_units);
% fprintf(fileID,'Number of light-modulated cells: %d of %d\r\n',length(light_cells),num_units);
% fprintf(fileID,'Percent light-modulated: %.2f\r\n', 100*length(light_cells)/num_units);
fprintf(fileID,'Number of tuned cells: %d of %d\r\n',length(tuned_cells),num_units);
fprintf(fileID,'Percent significantly tuned: %.2f\r\n', 100*length(tuned_cells)/num_units);
fprintf(fileID,'Number of regular-spiking cells: %d of %d\r\n',length(reg_cells),num_units);
fprintf(fileID,'Percent regular-spiking: %.2f\r\n', 100*length(reg_cells)/num_units);
fprintf(fileID,'Number of fast-spiking cells: %d of %d\r\n',length(FS_cells),num_units);
fprintf(fileID,'Percent fast-spiking: %.2f\r\n', 100*length(FS_cells)/num_units);
fprintf(fileID,'Number of linear cells: %.2f\r\n', sum(Fratio(:,1)>1));
fprintf(fileID,'Percent linear cells: %.2f\r\n', 100*sum(Fratio(:,1)>1)/size(Fratio,1));
% fprintf(fileID,'Number of units suppressed in all conditions: %d of %d\r\n',length(supp_cells),num_units);
% fprintf(fileID,'Percent suppressed: %.2f\r\n', 100*length(supp_cells)/num_units);
% fprintf(fileID,'Number of units activated in all conditions: %d of %d\r\n',length(enh_cells),num_units);
% fprintf(fileID,'Percent activated: %.2f\r\n', 100*length(enh_cells)/num_units);
% fprintf(fileID,'Number of units visually responsive at onset only: %d of %d\r\n',length(onset_cells),num_units);
% fprintf(fileID,'Percent visual onset-responsive: %.2f\r\n', 100*length(onset_cells)/num_units);
% fprintf(fileID,'Number of transiently-activated units: %d of %d\r\n',length(trans_cells),num_units);
% fprintf(fileID,'Percent transiently-activated: %.2f\r\n', 100*length(trans_cells)/num_units);
% fprintf(fileID,'Number of units visually responsive throughout full stim duration: %d of %d\r\n',length(sust_cells),num_units);
% fprintf(fileID,'Percent with sustained visual response: %.2f\r\n', 100*length(sust_cells)/num_units);
% fprintf(fileID,'Number of units whose onset response is in opposite direction of sustained evoked response: %d of %d\r\n',length(rev_cells),num_units);
% fprintf(fileID,'Percent of units with reversing visual response: %.2f\r\n', 100*length(rev_cells)/num_units);
% fprintf(fileID,'Number of units visually activated after delay: %d of %d\r\n',length(delay_act_cells),num_units);
% fprintf(fileID,'Percent delay-activated cells: %.2f\r\n', 100*length(delay_act_cells)/num_units);
% fprintf(fileID,'Number of units visually suppressed after delay: %d of %d\r\n',length(delay_sup_cells),num_units);
% fprintf(fileID,'Percent delay-suppressed cells: %.2f\r\n', 100*length(delay_sup_cells)/num_units);
% fprintf(fileID,'Number of units whos direction of light modulation depends on power: %d of %d\r\n',length(complex_cells),num_units);
% fprintf(fileID,'Percent with power-dependent light modulation: %.2f\r\n', 100*length(complex_cells)/num_units);
% fprintf(fileID,'Number of units suppressed in all conditions: %d of %d\r\n',length(all_supp),num_units);
% fprintf(fileID,'Percent suppressed: %.2f\r\n', 100*length(all_supp)/num_units);
% fprintf(fileID,'Number of units enhanced in all conditions: %d of %d\r\n',length(all_enh),num_units);
% fprintf(fileID,'Percent enhanced: %.2f\r\n', 100*length(all_enh)/num_units);
% fprintf(fileID,'Number of units with quick onset response: %d of %d\r\n',length(quick_enh),num_units);
% fprintf(fileID,'Percent with quick onset response: %.2f\r\n', 100*length(quick_enh)/num_units);
% fprintf(fileID,'Number of units with delayed activation: %d of %d\r\n',length(delayact),num_units);
% fprintf(fileID,'Percent with delayed activation response: %.2f\r\n', 100*length(delayact)/num_units);
% fprintf(fileID,'Number of units with quick onset but then suppressed: %d of %d\r\n',length(onthensupp),num_units);
% fprintf(fileID,'Percent with onset then suppression: %.2f\r\n', 100*length(onthensupp)/num_units);
% fprintf(fileID,'Number of units with quick onset, then suppressed then enhanced: %d of %d\r\n',length(onthenenh),num_units);
% fprintf(fileID,'Percent on-off-on: %.2f\r\n', 100*length(onthenenh)/num_units);
% fprintf(fileID,'Number of units suppressed with high power but enhanced with low power: %d of %d\r\n',length(lowenh_highsupp),num_units);
% fprintf(fileID,'Percent low-enhanced and high-suppressed: %.2f\r\n', 100*length(lowenh_highsupp)/num_units);
% fprintf(fileID,'Number of units enhanced with high power but suppressed with low power: %d of %d\r\n',length(lowsupp_highenh),num_units);
% fprintf(fileID,'Percent low-suppressed and high-enhanced: %.2f\r\n', 100*length(lowsupp_highenh)/num_units);
fclose(fileID);

% fileID2 = fopen('summary_stats.txt','w');
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, VISUAL trials, ALL cells: %6.3f %d\r\n',lightsig_vals(1,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, VISUAL trials, ALL cells: %6.3f %d\r\n',lightsig_vals(2,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, BLANK trials, ALL cells: %6.3f %d\r\n',lightsig_vals(3,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, BLANK trials, ALL cells: %6.3f %d\r\n',lightsig_vals(4,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, VISUAL trials, VISUAL cells: %6.3f %d\r\n',lightsig_vals(5,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, VISUAL trials, VISUAL cells: %6.3f %d\r\n',lightsig_vals(6,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, BLANK trials, VISUAL cells: %6.3f %d\r\n',lightsig_vals(7,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, BLANK trials, VISUAL cells: %6.3f %d\r\n',lightsig_vals(8,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, VISUAL trials, NONVISUAL cells: %6.3f %d\r\n',lightsig_vals(9,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, VISUAL trials, NONVISUAL cells: %6.3f %d\r\n',lightsig_vals(10,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in LOW condition, BLANK trials, NONVISUAL cells: %6.3f %d\r\n',lightsig_vals(11,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on FR in HIGH condition, BLANK trials, NONVISUAL cells: %6.3f %d\r\n',lightsig_vals(12,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on pref FR in LOW condition: %6.3f %d\r\n',lightsig_vals(13,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on pref FR in HIGH condition: %6.3f %d\r\n',lightsig_vals(14,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on visual FR change in LOW condition: %6.3f %d\r\n',lightsig_vals(15,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on visual FR change in HIGH condition: %6.3f %d\r\n',lightsig_vals(16,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on pref visual FR change in LOW condition: %6.3f %d\r\n',lightsig_vals(17,:));
% fprintf(fileID2,'Signed-rank test of sig light effect on pref visual FR change in HIGH condition: %6.3f %d\r\n',lightsig_vals(18,:));
% 
% fprintf(fileID2,'\r\n');
% fprintf(fileID2,'Signed-rank test of sig OSI_CV change in LOW condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(1,:));
% fprintf(fileID2,'Signed-rank test of sig OSI_CV change in HIGH condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(2,:));
% fprintf(fileID2,'Signed-rank test of sig OSI change in LOW condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(3,:));
% fprintf(fileID2,'Signed-rank test of sig OSI change in HIGH condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(4,:));
% fprintf(fileID2,'Signed-rank test of sig DSI_CV change in LOW condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(5,:));
% fprintf(fileID2,'Signed-rank test of sig DSI_CV change in HIGH condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(6,:));
% fprintf(fileID2,'Signed-rank test of sig DSI change in LOW condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(7,:));
% fprintf(fileID2,'Signed-rank test of sig DSI change in HIGH condition, TUNED cells: %6.3f %d\r\n',tuningsig_vals(8,:));
% fprintf(fileID2,'Signed-rank test of sig OSI_CV change in LOW condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(9,:));
% fprintf(fileID2,'Signed-rank test of sig OSI_CV change in HIGH condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(10,:));
% fprintf(fileID2,'Signed-rank test of sig OSI change in LOW condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(11,:));
% fprintf(fileID2,'Signed-rank test of sig OSI change in HIGH condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(12,:));
% fprintf(fileID2,'Signed-rank test of sig DSI_CV change in LOW condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(13,:));
% fprintf(fileID2,'Signed-rank test of sig DSI_CV change in HIGH condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(14,:));
% fprintf(fileID2,'Signed-rank test of sig DSI change in LOW condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(15,:));
% fprintf(fileID2,'Signed-rank test of sig DSI change in HIGH condition, ALL cells: %6.3f %d\r\n',tuningsig_vals(16,:));
% 
% fclose(fileID2);

fileID3 = fopen('cleanunits.txt','w');
fprintf(fileID3,'%d\r\n',[unitinfo(clean_units).name]);
fclose(fileID3);
end

function plot_scatter(data, color_var, colors, xlab, ylab, title, leg, lobf)
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
    plot(data(color_var==vars(i),1), data(color_var==vars(i),2),'.', 'Color', colors{i},'MarkerSize',28);
end
set(gca,'Fontsize',18,'linewidth',2)
xlabel(xlab,'Fontsize',18)
ylabel(ylab,'Fontsize',18)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)),100);
y=x;
plot(x,y,'k--');
line([min(xax(1),yax(1)) max(xax(2),yax(2))],[0 0],'color','k')
line([0 0], [min(xax(1),yax(1)) max(xax(2),yax(2))],'color','k')
xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
if lobf
    fittedX = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)), 200);
    for i = 1:length(vars)
        if sum(color_var==vars(i)) > 1
            coeffs(i,:) = polyfitZero(data(color_var==vars(i),1), data(color_var==vars(i),2), 1);
            fittedY(i,:) = polyval([0 coeffs(i,:)], fittedX);
            plot(fittedX,fittedY(i,:),'color', colors{i},'linewidth',2)
        end
    end
end
l=legend(leg,'location','best');
set(l,'fontsize',18)
print(f, '-dpng',title)
print2eps(title,f)        % doesn't seem to work with new matlab...
end

