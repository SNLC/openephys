function split_exp_sts(system,path,nCh)

% if kilosort was run on concatenated data, create separate spike_times and
% spike_clusters files for each separate experiment
%
% path - main direction where experiment subdirectories are stored
% EDITED 1/15/21 - add 'nCh' as input (# of channels) to check for same #
% of samples among diff experiments

cd(path)
dirs = regexp(genpath(path),['[^;]*'],'match');
exp_dirs = dirs(cellfun(@(x) ~strcmp(x,path),dirs,'uniformoutput',1));
exp_dirs = exp_dirs(~contains(exp_dirs,'csd','ignorecase',1));  % exclude csd files
exp_dirs = exp_dirs(~contains(exp_dirs,'phy','ignorecase',1));  % exclude .phy files
exp_dirs = exp_dirs(~contains(exp_dirs,'archive','ignorecase',1));  % exclude archive files
exp_dirs = exp_dirs(~contains(exp_dirs,'sync','ignorecase',1));  % exclude sync files
exp_dirs = exp_dirs(~contains(exp_dirs,'fields','ignorecase',1));  % exclude fields files
exp_dirs = exp_dirs(~contains(exp_dirs,'shank','ignorecase',1));  % exclude shank folders
exp_dirs = exp_dirs(~contains(exp_dirs,'sort','ignorecase',1));  % exclude folders with other sorting version

% load kilosort files for concatenated data (they will automatically be in
% chronological order)
all_sts = readNPY(fullfile(path,'spike_times.npy'));           % here, spike times from ALL clusters
all_clusts = readNPY(fullfile(path,'spike_clusters.npy'));


for n = 1:length(exp_dirs)
    if strcmpi(system,'openephys')
        first_half = exist(fullfile(exp_dirs{n},'100_CH1.continuous'),'file');
        second_hs = exist(fullfile(exp_dirs{n},'100_CH193.continuous'),'file');    % assumes using chs 64-127 on second headstage
        if first_half
            contfile = fullfile(exp_dirs{n},'100_CH1.continuous');
            contfile2 = sprintf('%s\\100_CH%d.continuous',exp_dirs{n},nCh);
        elseif second_hs
            contfile = fullfile(exp_dirs{n},'100_CH193.continuous');
            contfile2 = sprintf('%s\\100_CH%d.continuous',exp_dirs{n},256);
        else
            contfile = fullfile(exp_dirs{n},'100_CH65.continuous');
            contfile2 = sprintf('%s\\100_CH%d.continuous',exp_dirs{n},64+nCh);
        end
        [~, dataTime, dataInfo] = load_open_ephys_data(contfile);
        [~, dataTime2, dataInfo2] = load_open_ephys_data(contfile2);
        if length(dataTime) ~= length(dataTime2)
            if ~rem(abs(length(dataTime2)-length(dataTime)),1024)    
                dataTime = dataTime2;   % use last one (because that's how the openephystorawbinary script handles it)
                dataInfo = dataInfo2;
                fprintf('Continuous files have unequal number of batches. Using last continous file. \n')
            else
                error('Unequal number of samples among continuous files')
            end
        end
    %     starttime(n) = dataTime(1)*dataInfo.header.sampleRate-dataInfo.header.blockLength;   % changed 6/4/18
        if dataTime(1) ~= 0
            startblock = dataInfo.nsamples(1);
    %         endblock = 0;
        else
            startblock = 0;     % not sure why sometimes starttime is 0 vs 1024
    %         endblock = dataInfo.nsamples(1);
        end
        starttime(n) = dataTime(1)*dataInfo.header.sampleRate-startblock;       % actually, shouldn't this just be 0? (mak 4/9/19)
        endtime(n) = starttime(n) + length(dataTime) - 1;
    %     endtime(n) = floor(dataTime(end)*dataInfo.header.sampleRate)+endblock;     
    %     endtime(n) = dataInfo.header.blockLength * length(dataInfo.ts)+ starttime(n) - 1;
        if n>1
            %spike_times files per experiment starts with at the time of the experiment (reset for each exp for time 0), but may go longer than experiment (if there were subsequent experiments)
    %         writeNPY(all_sts(all_sts>=starttime(n))-starttime(n),fullfile(exp_dirs{order==n},'spike_times.npy'));       
    %         writeNPY(all_clusts(all_sts>=starttime(n)),fullfile(exp_dirs{order==n},'spike_clusters.npy'));   
    %         if starttime(n) <= starttime(n-1)            % starts over if press play button to stop between recordings
    %             nsamps = endtime(n)-starttime(n);
                nsamps = length(dataTime);
                starttime(n) = endtime(n-1)+1;      
                endtime(n) = starttime(n)+nsamps-1;
    %         else                % need to check this!!!
    %             err=starttime(n)-startblock -endtime(n-1);
    %             starttime(n) = starttime(n) - err + 1;
    %             endtime(n) = endtime(n)-err +1;
    %         end
            writeNPY(all_sts(all_sts>=starttime(n))-starttime(n),fullfile(exp_dirs{n},'spike_times.npy'));     % changed to -starttime MAK 6/12  
            writeNPY(all_clusts(all_sts>=starttime(n)),fullfile(exp_dirs{n},'spike_clusters.npy'));
            % need to add saving cluster_groups.csv and rez and spike_templates.npy and channel_postions.npy in each experiment
            % subfolder
        end
    elseif strcmpi(system,'intan')
        fileinfo = dir(fullfile(exp_dirs{n},'time.dat'));
        num_samples = fileinfo.bytes/4; % int32 = 4 bytes
        fid = fopen(fullfile(exp_dirs{n},'time.dat'), 'r');
        time = fread(fid, num_samples, 'int32'); 
        
        starttime(n) = time(1);
        endtime(n) = time(end);
        if n>1
            starttime(n) = endtime(n-1)+1;
            endtime(n) = starttime(n)+num_samples-1;
            writeNPY(all_sts(all_sts>=endtime(n-1))-starttime(n),fullfile(exp_dirs{n},'spike_times.npy'));     % changed to -starttime MAK 6/12  
            writeNPY(all_clusts(all_sts>=endtime(n-1)),fullfile(exp_dirs{n},'spike_clusters.npy'));
        end
%         [amplifier_channels,aux_input_channels,board_adc_channels,...
%         board_dig_in_channels,spike_triggers,supply_voltage_channels,...
%         frequency_parameters] = read_Intan_RHD2000_file;    % script from Intan

    end
end

writeNPY(all_sts,fullfile(exp_dirs{1},'spike_times.npy'));       %also write full spike times and cluster files into first experiment directory
writeNPY(all_clusts,fullfile(exp_dirs{1},'spike_clusters.npy')); 

end