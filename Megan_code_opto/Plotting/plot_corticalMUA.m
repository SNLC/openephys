function plot_corticalMUA(exp_path,probe,lighttrial) 
% probe: e.g., '64D_botton'
% lighttrial: 1 (light trial) or 0 (no-light trial)

cd(exp_path)

% get animal and exp_name so that you can load experiment results 
out = regexp(exp_path,'\\','split');
an_name = out{end-1};       % animal name
if contains(an_name,'LP','ignorecase',1) || contains(an_name,'V1','ignorecase',1) || contains(an_name,'TRN','ignorecase',1)
    an_name = strcat(out{end-2},'_',an_name);
elseif contains(an_name,'day','ignorecase',1) 
    an_name = out{end-2};
end
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
results_folder = fullfile('H:\LPproject\LPresults',an_name,exp_name);

% determine which channels to load
filenames = dir;
filenames = {filenames.name};
% nchs = length(cell2mat(strfind(filenames,'CH')));       % number of files with 'CH' in file name (i.e. continuous channel files)
first_half = exist(fullfile(exp_path,'100_CH1.continuous'),'file');
second_hs = exist(fullfile(exp_path,'100_CH193.continuous'),'file');    % assumes using chs 64-127 on second headstage
if first_half
    firstch = 1;
elseif second_hs 
    if exist(fullfile(exp_path,'100_CH129.continuous'),'file')
        firstch = 129;
    else
        firstch = 193;
    end
else
    firstch = 65;
end

% identify example trial
load('data.mat');
load('layers.mat');
[params.prestim,params.poststim,params.stimtime,params.trial_type,params.IVs] = get_exp_params(exp_path,'step');
if lighttrial
    extrials = find(params.trial_type(:,1) & params.trial_type(:,3)>0 & params.trial_type(:,4)==0,5,'first');   % first visual trial
else
    extrials = find(params.trial_type(:,1) & params.trial_type(:,3)==0 & params.trial_type(:,4)==0,5,'first');   % first visual trial
end

if lighttrial
    load('light_params.mat');
    nolightsamps = av_light_start*amp_sr+200:(av_light_start+light_dur)*amp_sr-200;
end 
   
% get probemap
p = eval(sprintf('probemap_%s_func',probe));
shks = unique(p.shaft);
for s=1:length(shks)
    chs = find(p.shaft==shks(s));
    [~,inds] = sort(p.z(chs));  % sorts from bottom to top
    chan_order = p.channels(chs(inds));
    chan_order = flipud(chan_order);    % change from top to bottom
    chan_order = chan_order+1; % from 1-64
    nchs = length(inds);

    run = 1;
    while run
        extrial = extrials(run);
        stimtime = field_trials(extrial,1).*(amp_sr/1000); % back into original sampling rate
        exsamps = [stimtime : stimtime + params.stimtime*amp_sr + 2*params.prestim*amp_sr-1];
        timevec = [-params.prestim:1/amp_sr:params.stimtime+params.prestim];
        timevec = timevec(1:end-1);

        chMat = nan(nchs,length(exsamps));
        [bf,af] = butter(2,150/(20000/2),'high');    % build low-pass filter (<150Hz)
        for ch=1:nchs
            [datin, ~, ~] =  load_open_ephys_data_faster(sprintf('%s\\100_CH%d.continuous',exp_path,chan_order(ch)+firstch-1));
            chMat(ch,:) = filtfilt(bf,af,datin(exsamps));
        end

        % scaleF = max(max(abs(chMat),[],2))/10;
        scaleF = 30;        % so that all graphs I make have the same scale
        % determine which channels to plot
        if isfile(sprintf('%s\\%s_results.mat',results_folder,exp_name)) % if already run through analysis_master, use max_ch values to determine which chs to include
            load(sprintf('%s\\%s_results.mat',results_folder,exp_name),'waveforms');
            max_ch = [waveforms.max_ch];
            startch = min(max_ch);
            lastch = max(max_ch);
        else        % otherwise, just find channels with voltage deflections large enough to feasibly be spikes
            if lighttrial   % need to exclude light artifacts
                nolight_chMat = chMat(:,nolightsamps);
                startch = find(max(abs(nolight_chMat),[],2)>mean(abs(nolight_chMat(:)))+6*std(abs(nolight_chMat(:))),1,'first');
                lastch = find(max(abs(nolight_chMat),[],2)>mean(abs(nolight_chMat(:)))+6*std(abs(nolight_chMat(:))),1,'last');
            else
                startch = find(max(abs(chMat),[],2)>mean(abs(chMat(:)))+6*std(abs(chMat(:))),1,'first');
                lastch = find(max(abs(chMat),[],2)>mean(abs(chMat(:)))+6*std(abs(chMat(:))),1,'last');
            end
        end
        colors = {[.75 .75 .75],[.85 .325 .098],[.635 .078 .184],[0 .447 .741],[.166 .674 .188]}; % <L2/3 = grey, L2/3=orange, L4=red, L5=blue, L6=green 
        layers_graph = floor(layers); % keep layers 5A&B together
        layers_graph(1:startch-1)=0;
        difflayers=unique(layers_graph);
        figure; hold on;
        for ch=startch:lastch
            plot(timevec,chMat(ch,:)./scaleF - (ch-startch)*12,'color',colors{layers_graph(ch)==difflayers})
        end
        ylim([(lastch-startch+1)*-12 12])
        yax = get(gca,'ylim');
        line([0 0],yax,'Color','k','LineStyle','-','linewidth',.5)
        line([params.stimtime params.stimtime],yax,'Color','k','LineStyle','-','linewidth',.5)
        yticklabels('')
        yticks([])
        if lighttrial
            patch([-params.prestim+av_light_start -params.prestim+av_light_start -params.prestim+av_light_start+light_dur -params.prestim+av_light_start+light_dur -params.prestim+av_light_start],[yax(1) yax(2) yax(2) yax(1) yax(1)], [.8 0 0], 'LineStyle', 'none', 'FaceAlpha',.07 );
        end
        xlabel('Time from vis stim onset (s)')
        ylabel('Channels (inf. <  > sup.)')
        set(gca,'linewidth',2,'fontsize',16)
        box(gca,'on')
        title(strcat(exp_name, ' Shank',num2str(s),' Trial#',num2str(extrial)),'fontsize',16,'Interpreter', 'none')
        rerun = input('Do you want to try the next trial? y/n ','s');
        if strcmp(rerun,'y')
            run=run+1;
        else 
            run=0;
        end
    end
    print(gcf,'-dpng',fullfile(results_folder,'Figures','LaminarMUA'))
    print(gcf,'-painters','-depsc',fullfile(results_folder,'Figures','LaminarMUA'))
end


return
