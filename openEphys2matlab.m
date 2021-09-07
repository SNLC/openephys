function openEphys2matlab(exp_path)
%extracts info from open Ephys files
%to do: write code to detect all
%continuous files in exp_path, separate aux and electrodes and load all
%channels

% files = dir(exp_path);
% nfiles = length(files)-2;
% eventfile = {};
% allFiles = cell(1,nfiles);
% for i = 1:nfiles
%     allFiles{i} = files(i+2).name;
%     if strcomp(allFiles{i},'all_channels.events')
%         eventfile = fullfile(exp_path,'all_channels.events');
%     elseif i==nfiles & isempty(eventfile)
%         error('all_channels.events is missing')
%     end
% end
% 
% 
% nADC = 0;

    

%todo - initialize channels w/ size info
cd(exp_path)
filenames = dir;
filenames = {filenames.name};
nchs = length(cell2mat(strfind(filenames,'CH')));       % number of files with 'CH' in file name (i.e. continuous channel files)
where_adcs = cellfun(@(x) ~isempty(x), strfind(filenames,'ADC'));
adcfiles = filenames(where_adcs);       % in case you aren't just using ADC inputs 1 and 2
nADC = length(adcfiles);      % number of files with 'ADC' in file name (i.e. analog input files)

%load electrode channels
first_half = exist(fullfile(exp_path,'100_CH1.continuous'),'file');
second_hs = exist(fullfile(exp_path,'100_CH193.continuous'),'file');    % assumes using chs 64-127 on second headstage
if first_half
    contfile = fullfile(exp_path,'100_CH1.continuous');
elseif second_hs
    contfile = fullfile(exp_path,'100_CH193.continuous');
else
    contfile = fullfile(exp_path,'100_CH65.continuous');
end
[~, dataTime, dataInfo] = load_open_ephys_data(contfile);       % changed from using 'faster' because I think times are less accurate?? - MAK 4/13/18
% dataTime = dataTime./dataInfo(1).header.sampleRate;       % uncomment if using load_open_ephys_data_faster
nsamples = length(dataTime);
% dataTime is in SECONDS

%load ADC inputs
for i = 1:nADC
    [ADCin(i,:),ADCTime(i,:),ADCinfo(i,:)] = load_open_ephys_data_faster(fullfile(exp_path,adcfiles{i}));
    [ADCin_2(i,:),ADCTime_2(i,:),ADCinfo_2(i,:)] = load_open_ephys_data_faster(fullfile(exp_path,adcfiles{i}));
end
ADCTime = ADCTime./ADCinfo(1).header.sampleRate;

%load all_channels.events
eventfile = fullfile(exp_path,'all_channels.events');
[events,eventTime,info] = load_open_ephys_data(eventfile);
 amp_sr = info.header.sampleRate;
% for i = 1:length(eventTime)                           % only use when timing between events and continuous files are off
%     [~,eventIdx(i)] = min(abs(dataTime-eventTime(i))); 
% end

%  % check:
% [dataTime(eventIdx(1:20)) eventTime(1:20)]
% dataTime(eventIdx(nsamples-10:nsamples))
% eventTime(nsamples-10:nsamples)

% check for gaps in dataTime and adjust dataTime and eventTimes accordingly
if sum(diff(dataTime).*amp_sr >= 2)   % if theres a gap of 2 or more samples between data points
    warning(sprintf('Warning: there are %d gaps in the data files',sum(diff(dataTime).*amp_sr >= 2)))
    gap_ts = find(diff(dataTime).*amp_sr>=2);
    for tt = 1:length(gap_ts)
        gap(tt) = diff([dataTime(gap_ts(tt)),dataTime(gap_ts(tt)+1)]);
        dataTime(gap_ts(tt)+1:end) = dataTime(gap_ts(tt)+1:end) - gap(tt) + 1/amp_sr;
        eventTime(eventTime*amp_sr>=gap_ts(tt)) = eventTime(eventTime*amp_sr>=gap_ts(tt)) - gap(tt) + 1/amp_sr;
    end
end
eventIdx = floor((eventTime-dataTime(1))*amp_sr+1);

% in case of weird event where analog and data file lengths don't match
if size(ADCin,2)~=nsamples
    warning('Amplifier and analog data file lengths dont match.')
    shorter = min(nsamples,size(ADCin,2));
    err = (dataTime(shorter)-ADCTime(1,shorter)')*1000;
    disp(strcat('Difference (in ms) from end of data and analog time vectors:',num2str(err)))
    aa = input('Force them to be same length? Y/N: ','s');
    if strcmp(aa,'Y')
        ADCin = ADCin(:,1:shorter);
        nsamples = shorter;
        dataTime = dataTime(1:shorter);
        disp('Forcing analog and data vectors to be same length')
    else
        error('Amplifier and analog data file lengths dont match.')
    end
end
clear ADCTime

%define digital channels (subtract 1 from label on i/o boards)
epocCH = 0;
encdACH = 1;
encdBCH = 2;

%define analog channels
photo = ADCin(1,:);
if size(ADCin,1)>1  % if opto
    LED = ADCin(2:end,:);   % changed 9/2/20 - before was just ADCin(2,:) (MAK)
end
clear ADCin

%digital events
epoc = zeros(size(dataTime));
encdA = zeros(size(dataTime));
encdB = zeros(size(dataTime));

%find epoc on and off times and fill out binary vector same size as data
epocOn = eventIdx(events==epocCH&info.eventId&info.eventType==3);
epocOff = eventIdx(events==epocCH&~info.eventId&info.eventType==3);

%in case of corrupted files:
epocOn = epocOn(ismember(epocOn,1:nsamples));
epocOff = epocOff(ismember(epocOff,1:nsamples));

% downsample to 1000Hz to save memory
LN              = nsamples;
div             = amp_sr/1000;
zx              = 1:div:LN;
izx             = floor(zx);
time_index      = dataTime(izx)-dataTime(1);    % downsample from 20000 hz to 1000 hz
clear izx zx
% epocOn = epocOn(izx);
% epocOff = epocOff(izx);

% if epocOn(1)<epocOff(1) 
    if size(epocOn,1 )> size(epocOff,1)   % if experiment got cut off in middle of trial...
        epocOn = epocOn(1:size(epocOff,1));         % ...drop the last trial
    elseif length(epocOn)==(length(epocOff)-1)
        if epocOn(1)<epocOff(2)
            epocOff(1) = [];
        end
    end
    if size(epocOn) == size(epocOff)
        field_trials = [epocOn epocOff];
        trials = [dataTime(epocOn) dataTime(epocOff)];    % find and ismember in case of corrupted files
    else
        error('Number of epoc start and ends time dont match.')
    end
    for i = 1:length(epocOn)
        if i<=length(epocOff)
            epoc(epocOn(i):epocOff(i)) = ones(size(epocOn(i):epocOff(i)));
        else
            epoc(epocOn(i):end) = ones(1,nsamples-epocOn(i));
        end
    end
% end

% load analyzer - check for postdelay
s = dir; 
for i=1:length(s)
    if strfind(s(i).name,'.analyzer') 
        analyze_file = s(i).name;
        load(sprintf('%s/%s',exp_path,analyze_file),'-mat') ;    % load analyzer file with stimulus info
        if size(field_trials,1) > Analyzer.M.NumTrials          % TEMP - when files were concatenated, cut off extra trials from end
            field_trials = field_trials(1:Analyzer.M.NumTrials,:);
        end
        postdelay = Analyzer.P.param{2}{3};
        postdelay_samps = postdelay*amp_sr;
        field_trials(:,2) = field_trials(:,2)-postdelay_samps;
    end
end

%downsample field_trials for Megan's code:
field_trials =  floor(field_trials./div)+1;

%find encdA on and off times and fill out binary vector same size as data
encdAOn = eventIdx(events==encdACH&info.eventId&info.eventType==3);
encdAOff = eventIdx(events==encdACH&~info.eventId&info.eventType==3);

if encdAOn(1)>encdAOff(1)
    encdA(1:encdAOff(1)) = ones(encdAOff(1),1);
    for i = 1:length(encdAOn)
        if i+1<=length(encdAOff)
            encdA(encdAOn(i):encdAOff(i+1)) = ones(size(encdAOn(i):encdAOff(i+1)));
        else
            encdA(encdAOn(i):end) = ones(nsamples-encdAOn(i)+1,1);
        end
    end
else
    for i = 1:length(encdAOn)
        if i<=length(encdAOff)
            encdA(encdAOn(i):encdAOff(i)) = ones(size(encdAOn(i):encdAOff(i)));
        else
            encdA(encdAOn(i):end) = ones(nsamples-encdAOn(i)+1,1);
        end
    end
end


%find encdB on and off times and fill out binary vector same size as data
encdBOn = eventIdx(events==encdBCH&info.eventId&info.eventType==3);
encdBOff = eventIdx(events==encdBCH&~info.eventId&info.eventType==3);

if encdBOn(1)>encdBOff(1)
    encdB(1:encdBOff(1)) = ones(encdBOff(1),1);
    for i = 1:length(encdBOn)
        if i+1<=length(encdBOff)
            encdB(encdBOn(i):encdBOff(i+1)) = ones(size(encdBOn(i):encdBOff(i+1)));
        else
            encdB(encdBOn(i):end) = ones(nsamples-encdBOn(i)+1,1);
        end
    end
else
    for i = 1:length(encdBOn)
        if i<=length(encdBOff)
            encdB(encdBOn(i):encdBOff(i)) = ones(size(encdBOn(i):encdBOff(i)));
        else
            encdB(encdBOn(i):end) = ones(nsamples-encdBOn(i)+1,1);
        end
    end
end

clear dataTime dataInfo ADCinfo adcfiles encdAOff encdAOn encdBOff encdBOn epocOff epocOn events eventTime LN s
%getSyncTimes for "re"
[re]= getSyncTimesRevCorr_AC(photo',(1/amp_sr));

%todo - write code for when recording starts during a trial ->
%epocOff(1)<epocOn(1)

%save variables in data.mat
cd(exp_path)
if exist('LED','var')
    save('data.mat', 'trials','field_trials','amp_sr','photo','LED','epoc','encdA','encdB','re','time_index')
else
    save('data.mat', 'trials','field_trials','amp_sr','photo','epoc','encdA','encdB','re','time_index')
end

end

function [re]= getSyncTimesRevCorr_AC(x,sp)   %%x is Phot, sp is sampling period (1/24414.1)

sigma = 15/4881.82/sp;          %%I am iffy on why the sigma is this   NEED TO CHECK ON THIS (MAK)
size = length(x);
z = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-z .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
gaussFilterNorm = gaussFilter';
%sig = 15/4882.81/sp;

%h = fspecial('gaussian', [length(x) 1], sig);

x = ifft(fft(x).*abs(fft(gaussFilterNorm)));

thresh =(max(x)+min(x))*.75;        % changed from .5 2/26/19 MAK 

x = (sign(x-thresh) + 1)/2;

x = diff(x);
id = find(x<0); %Get rid of falling edges
x(id) = 0;

re_i = find(x);  %Index values of rising edges
re = re_i*sp;           %%this will be the vector for the timestamps (rising phase of photodiode signal). you need to remove the first 2 and the last 2 timestamps. Timestamps signify the end of a four second balack to gray period.
%re_r = round(re);
end