function make_artifacttemplate
% took 14min to make templates
exp_paths = {'F:\LPproject-L5\NPRH24\LP\2019-12-18_12-16-55_Halo',...
    'F:\LPproject-L5\NPRH24\LP\2019-12-18_12-53-20_sparsenoise'};
adcLight = 7;

tic
for n=1:length(exp_paths)   % make separate artifact templates for each experiment (in case light conds were different)
    exp_path = exp_paths{n};
    cd(exp_path)
    filenames = dir;
    filenames = {filenames.name};
    where_adcs = cellfun(@(x) ~isempty(x), strfind(filenames,'ADC'));
    adcfiles = filenames(where_adcs);       % in case you aren't just using ADC inputs 1 and 2
    nADC = length(adcfiles); 
    nchs = length(cell2mat(strfind(filenames,'CH')));       % number of files with 'CH' in file name (i.e. continuous channel files)
    first_half = exist(fullfile(exp_path,'100_CH1.continuous'),'file');
    second_hs = exist(fullfile(exp_path,'100_CH193.continuous'),'file');    % assumes using chs 64-127 on second headstage
    if first_half
        firstch = 1;
    elseif second_hs 
        firstch = 193;
    else
        firstch = 65;
    end

    [ADCin,~,~] = load_open_ephys_data_faster(sprintf('%s\\100_ADC%d.continuous',exp_path,adcLight));
    thresh=std(diff(ADCin))*20;    % alt: median(ADCin)*2;  
    ampthresh = median(ADCin)*2;
    lightON{n} = find(diff(ADCin)>thresh&abs(ADCin(2:end))>ampthresh);  % find time points at which light turned ON (passes both slope and amplitude thresholds)
    lightON{n}(diff(lightON{n})==1) = []; %use second time points (i.e., NOT +1) because those will be the acual peaks
    lightOFF{n} = find(diff(ADCin)<-thresh&abs(ADCin(1:end-1))>ampthresh);
    lightOFF{n}(diff(lightOFF{n})==1) = [];
    datlength(n) = length(ADCin);
    if length(lightOFF) ~= length(lightON)
        error('Error: %d onsets but %d light offsets detected',length(lightON{n}),length(lightOFF{n}))
    end
    
    artsON{n} = zeros(nchs,length(lightON{n}),4001);   
    artsOFF{n} = zeros(nchs,length(lightOFF{n}),4001);   
    for ch=1:nchs       % do separately for each channel
        [datin, ~, ~] =  load_open_ephys_data_faster(sprintf('%s\\100_CH%d.continuous',exp_path,firstch+ch-1));
        if length(datin) ~= datlength(n)
            error('ADC and DAT files are unequal lengths')
        end
%         artsON = zeros(length(lightON{n}),4001);     % first, collect each artifact before creating an average template
%         artsOFF = zeros(length(lightOFF{n}),4001);
        for i=1:length(lightON{n})     % for each light onset
            [~,maxtON{n}(ch,i)] = min(datin(lightON{n}(i):lightON{n}(i)+19));  % look for peak time within first 1ms (20samps)
            [~,maxtOFF{n}(ch,i)] = min(datin(lightOFF{n}(i):lightOFF{n}(i)+19));
            blON = mean(datin(lightON{n}(i)+maxtON{n}(ch,i)-50:lightON{n}(i)+maxtON{n}(ch,i)-11));
            blOFF = mean(datin(lightOFF{n}(i)+maxtOFF{n}(ch,i)-50:lightOFF{n}(i)+maxtOFF{n}(ch,i)-11));
            artsON{n}(ch,i,:) = datin(lightON{n}(i)+maxtON{n}(ch,i)-50:lightON{n}(i)+maxtON{n}(ch,i)+3950)-blON;
            artsOFF{n}(ch,i,:) = datin(lightOFF{n}(i)+maxtOFF{n}(ch,i)-50:lightOFF{n}(i)+maxtOFF{n}(ch,i)+3950)-blOFF;
        end
        ontimes{ch} = unique(maxtON{n}(ch,:));
        offtimes{ch} = unique(maxtOFF{n}(ch,:));    
    end
    
    ontimes_all = unique(mode(maxtON{n}));  % most common timepoints across channels
    offtimes_all = unique(mode(maxtOFF{n}));
    maxtON{n}(maxtON{n}<min(ontimes_all)) = min(ontimes_all);
    maxtON{n}(maxtON{n}>max(ontimes_all)) = max(ontimes_all);
    maxtOFF{n}(maxtOFF{n}<min(offtimes_all)) = min(offtimes_all);
    maxtOFF{n}(maxtOFF{n}>max(offtimes_all)) = max(offtimes_all);
    for ii=1:length(ontimes_all)  %max(cellfun(@(x) length(x),ontimes,'uniformoutput',1))
        artON_all{n}{ii} = zeros(nchs,4001);    % collect 2.5ms before (50samps) and 197.5ms after (3950samps) light onset
    end
    for nn=1:length(offtimes_all)  %max(cellfun(@(x) length(x),offtimes,'uniformoutput',1))
        artOFF_all{n}{nn} = zeros(nchs,4001);       % same for offset
    end
    for ch=1:nchs
        for ii = 1:length(ontimes_all)
            artON_all{n}{ii}(ch,:) = mean(artsON{n}(ch,maxtON{n}(ch,:)==ontimes_all(ii),:))./.195; % scale to match binary file
        end
        for nn = 1:length(offtimes_all)
            artOFF_all{n}{nn}(ch,:) = mean(artsOFF{n}(ch,maxtOFF{n}(ch,:)==offtimes_all(nn),:))./.195;
        end
    %     artON_all{n} = artON_all{n}./.195;      % scale to match binary file
    %     artOFF_all{n} = artOFF_all{n}./.195;
    end
end
toc

cd ..       % cd into main folder (where exp_paths are subfolders)
save('artifact_templates.mat','artON_all', 'artOFF_all');
lightON_all = lightON{1}'+mode(maxtON{1}(:));
lightOFF_all = lightOFF{1}'+mode(maxtOFF{1}(:));
for n=2:length(exp_paths)   % combine onset and offset times for concatenated bin files
    lightON_all = [lightON_all lightON{n}'+mode(maxtON{n})+datlength(n-1)];
    lightOFF_all = [lightOFF_all lightOFF{n}'+mode(maxtOFF{n})+datlength(n-1)];
end
maxtON = cell2mat(maxtON);
maxtOFF = cell2mat(maxtOFF);
% ontimes = cell2mat(ontimes);
% offtimes = cell2mat(offtimes);
    
%%make new binary file
fid=fopen('F:\LPproject-L5\NPRH24\LP\NPRH24_LP_all_artremoval.bin','r+');
% fid=fopen('F:\LPproject-L5\NPRH24\NPRH24_LP_all.bin','r');
for i=1:length(lightON_all)
    fseek(fid,2*(nchs*(lightON_all(i)-1)-50*nchs),-1);
    raw = fread(fid,4001*nchs,'int16');
    rawmat = reshape(raw,nchs,4001);
    exp = find(i<=cumsum(cellfun(@(x) length(x),lightON,'uniformoutput',1)),1,'first');  % get which exp num you're on
    baseline = mean(rawmat(:,1:40),2);  % get baseline (before light onset )
    artsub = raw-(artON_all{exp}{mode(maxtON(:,i))==ontimes_all}(:)+repmat(baseline,4001,1))+repmat(mean(rawmat(:,1:40),2),4001,1);      % assumes same ontimes for all channels??
    fseek(fid,2*(nchs*(lightON_all(i)-1)-50*nchs),-1);
    fwrite(fid,artsub,'int16');
end

for i=1:length(lightOFF_all)
    fseek(fid,2*(nchs*(lightOFF_all(i)-1)-50*nchs),-1);
    raw = fread(fid,4001*nchs,'int16');
    rawmat = reshape(raw,nchs,4001);
    baseline = mean(rawmat(:,1:40),2);  % get baseline (before light onset)
    exp = find(i<=cumsum(cellfun(@(x) length(x),lightOFF,'uniformoutput',1)),1,'first');  % get which exp num you're on
    artsub = raw-(artOFF_all{exp}{maxtOFF(i)==offtimes}(:)+repmat(baseline,4001,1))-repmat(baseline,4001,1);
    fseek(fid,2*(nchs*(lightOFF_all(i)-1)-50*nchs),-1);
    fwrite(fid,artsub,'int16');
end
fclose(fid);