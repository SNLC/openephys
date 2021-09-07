function [waveforms,max_ch,shank] = readWaveformsFromRez(unit,exp_path,rez)
% unit = cluster identifier (e.g. 47)
% exp_path = experiment directory
% rez = rez.mat output from kilosort
    % load(sprintf('%s\\rez.mat',exp_path))

clusters = readNPY(sprintf('%s\\spike_clusters.npy',exp_path));
templates = readNPY(sprintf('%s\\spike_templates.npy',exp_path));
chans = readNPY(sprintf('%s\\channel_positions.npy',exp_path));

n = hist(templates(clusters==unit),0:double(max(templates)));
[~,m] = max(n);

waveforms(:,:) = rez.Wraw(:,:,m)';
waveforms = waveforms(22:62,:).*.195;     % not at all sure about .195...
nchs = size(chans,1);
if range(chans(:,1))<50
    shankchs = nchs;
else
    shankchs = nchs/ceil(range(chans(:,1))/max(diff(unique(chans(:,1)))));     % chs per shank 
end

% figure out max_ch            
% index = rez.ypos(m);      % CHANGED 2/7/20 - found small % of units for
% whom this index doesnt match max(waveforms), in which case the latter is
% correct
[~,index] =  max(range(waveforms));
shank = ceil(index/shankchs)-1;       
[~,pos] = sort(chans(shank*shankchs+1:(shank+1)*shankchs,2));
max_ch = find(flipud(pos)==(index-(shank*shankchs)));

if shankchs == nchs
    shank = [];
end

end