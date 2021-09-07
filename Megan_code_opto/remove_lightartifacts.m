function remove_lightartifacts(datfile,adcLight,exp_paths)
% updated 4/22/20 to improve the artifact template that is removed (MAK)
% took 14min to make templates
% exp_paths = {'F:\LPproject-L5\NPRH24\LP\2019-12-18_12-16-55_Halo',...
%     'F:\LPproject-L5\NPRH24\LP\2019-12-18_12-53-20_sparsenoise'};
% adcLight = 7;

tic
for n=1:length(exp_paths)   % make separate artifact templates for each experiment (in case light conds were different)
    exp_path = exp_paths{n};
%     cd(exp_path)
    filenames = dir(exp_path);
    filenames = {filenames.name};
    nchs = length(cell2mat(strfind(filenames,'CH')));       % number of files with 'CH' in file name (i.e. continuous channel files)
    for i = 1:length(adcLight)
        [ADCin,~,~] = load_open_ephys_data_faster(sprintf('%s\\100_ADC%d.continuous',exp_path,adcLight(i)));
        thresh=std(diff(ADCin))*20;    % alt: median(ADCin)*2;  
        ampthresh = median(ADCin)*2;
        lightON{n}{i} = find(diff(ADCin)>thresh&abs(ADCin(2:end))>ampthresh);  % find time points at which light turned ON (passes both slope and amplitude thresholds)
        lightON{n}{i}(diff(lightON{n}{i})==1) = []; %use second time points (i.e., NOT +1) because those will be the acual peaks
        lightOFF{n}{i} = find(diff(ADCin)<-thresh&abs(ADCin(1:end-1))>ampthresh);
        lightOFF{n}{i}(diff(lightOFF{n}{i})==1) = [];
    end
    if length(adcLight)>1
        newlightON{n} = [union(lightON{n}{:})];
        newlightON{n}(diff(newlightON{n})==1) = []; % again, go with later timepoints
        newlightOFF{n} = [union(lightOFF{n}{:})];
        newlightOFF{n}(diff(newlightOFF{n})==1) = []; % again, go with later timepoints
%         chkON = length(newlightON{n})==numel([lightON{n}{:}])-.5*length(lightON{n}{1});     % check that theres the right number of trials in exps with 2 LEDs and 4 lightconds (1ON, 2ON, 1&2ON, OFF)
%         chkOFF = length(newlightOFF{n})==numel([lightOFF{n}{:}])-.5*length(lightOFF{n}{1});     
%         if ~chkON || ~chkOFF
%             error('Error: incorrect numbers of light trials in dual-LED exp were detected')
%         end
    else
        newlightON{n} = lightON{n}{1};
        newlightOFF{n} = lightOFF{n}{1};
    end
    datlength(n) = length(ADCin);
    if length(newlightOFF) ~= length(newlightON)
        error('Error: %d onsets but %d light offsets detected',length(newlightON{n}),length(newlightOFF{n}))
    end 
    fprintf('Exp %d: %d onsets and %d light offsets detected\n',n,length(newlightON{n}),length(newlightOFF{n}))
    chk = input('Are there any "trials" that need to be removed? y/n ', 's');
    if strcmp(chk,'y')
        badtris = input('Which trials do you need to remove? ');
        for c=1:length(badtris)
            newlightON{n}(badtris(c)) = [];
            newlightOFF{n}(badtris(c)) = [];
        end
    fprintf('Exp %d: %d onsets and %d light offsets detected\n',n,length(newlightON{n}),length(newlightOFF{n}))
    end
end
   
% cd ..       % cd into main folder (where exp_paths are subfolders)
lightON_all = newlightON{1}';
lightOFF_all = newlightOFF{1}';
for n=2:length(exp_paths)   % combine onset and offset times for concatenated bin files
    lightON_all = [lightON_all newlightON{n}'+datlength(n-1)];
    lightOFF_all = [lightOFF_all newlightOFF{n}'+datlength(n-1)];
end
    
%%revise binary file
fid=fopen(datfile,'r+');
for i=1:length(lightON_all)
    fseek(fid,2*(nchs*(lightON_all(i)-1)-50*nchs),-1);
    raw = fread(fid,4001*nchs,'int16');
    rawmat = reshape(raw,nchs,4001);
    [bf,af] = butter(2,200/(20000/2),'low');    % build low-pass filter (<200Hz)
    [bf2,af2] = butter(2,[20/(20000/2) 200/(20000/2)]);    % build bandpass filter (20-200Hz)
    sig_artremov = zeros(size(rawmat));
    for ch = 1:nchs
        baseline = mean(rawmat(ch,1:50));   
        tStart = find(abs(rawmat(ch,:)-baseline)>=5*std(rawmat(ch,1:50)),1,'first'); % first timepoint >5stds beyond baseline
        [~,tPeak] = min(rawmat(ch,tStart:tStart+20));   % minimum within 20ms of tStart (changed from abs(max) to min to be exclusive to light onsets) 
        tPeak = tPeak + tStart-1;
        if isempty(tPeak)
            error('no peak timepoint detected (ch#%d, lightonset#%d)',ch,i)
        end
%         tArtEnd = find(abs(chFilt(tPeak:end))>=abs(rawmat(ch,tPeak:end)),1,'first') + tPeak-1;
        tArtEnd = tPeak+20; % take out first ms from light onset  (**assumes 20khz sampling rate)
        chFilt = filtfilt(bf, af, rawmat(ch,tArtEnd:end)); % low-pass filtered signal (for constructing artifact template)
        chFilt2 = filtfilt(bf2, af2, rawmat(ch,tArtEnd+100:end)); % band-pass filtered signal (for replacing low-freq oscillations after artifact removal)
%         art = [rawmat(ch,tStart:tArtEnd-1) chFilt(tArtEnd:end)];    % artifact: peak voltage deflection from raw trace + low-pass-filtered raw signal
        art = [rawmat(ch,tStart:tArtEnd-1) chFilt];    % artifact: peak voltage deflection from raw trace + low-pass-filtered raw signal
        artremov= [rawmat(ch,1:tStart-1) rawmat(ch,tStart:end)-art+baseline];
        tTrans = tArtEnd + 100+find((artremov(tArtEnd+100:end-1)>chFilt2(1:end-1)+rawmat(ch,end))&(artremov(tArtEnd+101:end)<chFilt2(2:end)+rawmat(ch,end)),1,'first'); % find timepoint at which to start adding chFilt2
        sig_artremov(ch,:) = artremov;
        if isempty(tTrans)  % if chFilt2+rawmat(ch,end) doesn't overlap with artremov
%             [~,tTrans_alt] = min(abs(chFilt2(101:end)+rawmat(ch,end)-artremov(tArtEnd+100:end))); % point at which chFilt2+rawmat(ch,end) is closest to artremov (arbitrarily start looking from 5ms after tArtEnd)
%             tTrans_alt = tTrans_alt+tArtEnd+100-1;
%             ramp_st = find(artremov(tArtEnd:end)-baseline<2*std(artremov(tArtEnd:end)),1,'first')+tArtEnd-1; % point at which to start the ramp (first point <2stds away from baseline)
%             ramp_st = find(artremov(tArtEnd:end)-baseline<2*std(artremov(tArtEnd:end)),1,'first')+tArtEnd-1; % point at which to start the ramp (first point <2stds away from baseline)
                tTrans_alt = tArtEnd+100+find((artremov(tArtEnd+100:end-1)>chFilt2(1:end-1)+baseline)&(artremov(tArtEnd+101:end)<chFilt2(2:end)+baseline),1,'first'); % find timepoint at which to start adding chFilt2
                sig_artremov(ch,tArtEnd:end) = artremov(tArtEnd:end)+linspace(0,rawmat(ch,end)-artremov(end),length(tArtEnd:length(artremov))); % % make linear rampup to rawmat(ch,end) level
                sig_artremov(ch,tTrans_alt:end) = sig_artremov(ch,tTrans_alt:end)+chFilt2(tTrans_alt-(tArtEnd+100)+1:end); % and then from there, add chFilt2 
%             sig_artremov(ch,ramp_st:tTrans_alt-1) = artremov(ramp_st:tTrans_alt-1)+linspace(0,chFilt2(tTrans_alt)+rawmat(ch,end)-baseline,length(ramp_st:tTrans_alt-1)); % make linear rampup to chFilt2 level
%             sig_artremov(ch,tTrans_alt:end) = artremov(tTrans_alt:end)+rawmat(ch,end)-artremov(end)+chFilt2(tTrans_alt:end); % and then from there, add chFilt2 (rawmat(ch,end)-artremov(end) so that rawmat(ch,end)==sig_artremov(ch,end))
%             sig_artremov(ch,tArtEnd:end) = artremov(tArtEnd:end)+linspace(0,rawmat(ch,end)-artremov(tArtEnd-1),length(chFilt(tArtEnd:end)));
        else
            sig_artremov(ch,tTrans:end) = artremov(tTrans:end)+rawmat(ch,end)-artremov(end)+chFilt2(tTrans-(tArtEnd+100)+1:end); % rawmat(ch,end)-artremov(end) so that rawmat(ch,end)==sig_artremov(ch,end)
        end
    end
    sig_artremov = sig_artremov(:);
    fseek(fid,2*(nchs*(lightON_all(i)-1)-50*nchs),-1);
    fwrite(fid,sig_artremov,'int16');
end

for i=1:length(lightOFF_all)
    fseek(fid,2*(nchs*(lightOFF_all(i)-1)-50*nchs),-1);
    raw = fread(fid,4001*nchs,'int16');
    rawmat = reshape(raw,nchs,4001);
    [bf,af] = butter(2,200/(20000/2),'low');    % build low-pass filter (<200Hz)
    [bf2,af2] = butter(2,[20/(20000/2) 200/(20000/2)]);    % build bandpass filter (20-200Hz)
    sig_artremov = zeros(size(rawmat));
    for ch = 1:nchs
          baseline = mean(rawmat(ch,1:50));
        tStart = find(abs(rawmat(ch,:)-baseline)>=5*std(rawmat(ch,1:50)),1,'first'); % first timepoint >5stds beyond baseline
        [~,tPeak] = max(rawmat(ch,tStart:tStart+40));  % maximum within 2ms of tStart (changed from abs(max) to max to be exclusive to light offsets ) ** 50 instead of 20 because offset takes longer to peak and is more variable than onset for some reason
        tPeak = tPeak + tStart-1;
        if isempty(tPeak)
            error('no peak timepoint detected (ch#%d, lightoffset#%d)',ch,i)
        end
%         tArtEnd = find(abs(chFilt(tPeak:end))>=abs(rawmat(ch,tPeak:end)),1,'first') + tPeak-1;
        tArtEnd = tPeak+20; % take out first ms from light onset (**assumes 20khz sampling rate)
        chFilt = filtfilt(bf, af, rawmat(ch,tArtEnd:end)); % low-pass filtered signal (for constructing artifact template)
        chFilt2 = filtfilt(bf2, af2, rawmat(ch,tArtEnd+100:end)); % band-pass filtered signal (for replacing low-freq oscillations after artifact removal)
%         art = [rawmat(ch,tStart:tArtEnd-1) chFilt(tArtEnd:end)];    % artifact: peak voltage deflection from raw trace + low-pass-filtered raw signal
        art = [rawmat(ch,tStart:tArtEnd-1) chFilt];    % artifact: peak voltage deflection from raw trace + low-pass-filtered raw signal
        artremov= [rawmat(ch,1:tStart-1) rawmat(ch,tStart:end)-art+baseline];
        tTrans = tArtEnd + 100+find((artremov(tArtEnd+100:end-1)>chFilt2(1:end-1)+rawmat(ch,end))&(artremov(tArtEnd+101:end)<chFilt2(2:end)+rawmat(ch,end)),1,'first'); % find timepoint at which to start adding chFilt2
        sig_artremov(ch,:) = artremov;
        if isempty(tTrans)  % if chFilt2+rawmat(ch,end) doesn't overlap with artremov
%             [~,tTrans_alt] = min(abs(chFilt2(101:end)+rawmat(ch,end)-artremov(tArtEnd+100:end))); % point at which chFilt2+rawmat(ch,end) is closest to artremov (arbitrarily start looking from 5ms after tArtEnd)
%             tTrans_alt = tTrans_alt+tArtEnd+100-1;
%             ramp_st = find(artremov(tArtEnd:end)-baseline<2*std(artremov(tArtEnd:end)),1,'first')+tArtEnd-1; % point at which to start the ramp (first point <2stds away from baseline)
%             ramp_st = find(artremov(tArtEnd:end)-baseline<2*std(artremov(tArtEnd:end)),1,'first')+tArtEnd-1; % point at which to start the ramp (first point <2stds away from baseline)
                tTrans_alt = tArtEnd+100+find((artremov(tArtEnd+100:end-1)>chFilt2(1:end-1)+baseline)&(artremov(tArtEnd+101:end)<chFilt2(2:end)+baseline),1,'first'); % find timepoint at which to start adding chFilt2
                sig_artremov(ch,tArtEnd:end) = artremov(tArtEnd:end)+linspace(0,rawmat(ch,end)-artremov(end),length(tArtEnd:length(artremov))); % make linear rampup to rawmat(ch,end) level
                sig_artremov(ch,tTrans_alt:end) = sig_artremov(ch,tTrans_alt:end)+chFilt2(tTrans_alt-(tArtEnd+100)+1:end); % and then from there, add chFilt2 
%             sig_artremov(ch,ramp_st:tTrans_alt-1) = artremov(ramp_st:tTrans_alt-1)+linspace(0,chFilt2(tTrans_alt)+rawmat(ch,end)-baseline,length(ramp_st:tTrans_alt-1)); % make linear rampup to chFilt2 level
%             sig_artremov(ch,tTrans_alt:end) = artremov(tTrans_alt:end)+rawmat(ch,end)-artremov(end)+chFilt2(tTrans_alt:end); % and then from there, add chFilt2 (rawmat(ch,end)-artremov(end) so that rawmat(ch,end)==sig_artremov(ch,end))
%             sig_artremov(ch,tArtEnd:end) = artremov(tArtEnd:end)+linspace(0,rawmat(ch,end)-artremov(tArtEnd-1),length(chFilt(tArtEnd:end)));
        else
            sig_artremov(ch,tTrans:end) = artremov(tTrans:end)+rawmat(ch,end)-artremov(end)+chFilt2(tTrans-(tArtEnd+100)+1:end); % rawmat(ch,end)-artremov(end) so that rawmat(ch,end)==sig_artremov(ch,end)
        end
    end
    sig_artremov = sig_artremov(:);
    fseek(fid,2*(nchs*(lightOFF_all(i)-1)-50*nchs),-1);
    fwrite(fid,sig_artremov,'int16');
end
fclose(fid);