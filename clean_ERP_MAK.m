function [Final_ERP probemap_chrem] = clean_ERP_MAK (ERP,timing,probemap_chrem,shank)
% modified by MAK 6/6/17 - type in ch to remove rather than selecting

% EGM edited 5/18/2016 to suit 64D probes & track which channels were
% deleted & do normalization

norm_ERP = ERP;
for i=1:length(ERP(:,1))
    if ~isnan(ERP(i,1))
        norm_ERP(i,:) = ERP(i,:)-nanmean(ERP,1);
    end
end

figure;
plot(timing,ERP','DisplayName',sprintf('ERP_shk%d',shank),'YDataSource',sprintf('ERP_shk%d',shank));figure(gcf)
xlabel({'Time from stimulus onset (s)'});
ylabel({'Voltage (uV)'})
xlim([-0.05 0.4])
fig_title=sprintf('%s','Raw ERP ');
title(strrep(fig_title,'_','\_'));

for i=1:size(ERP,1)
    legend_cell{1,i} = num2str(i);
end
legend(legend_cell);

figure;
plot(timing,norm_ERP','DisplayName',sprintf('ERP_shk%d',shank),'YDataSource',sprintf('ERP_shk%d',shank));figure(gcf)
xlabel({'Time from stimulus onset (s)'});
ylabel({'Voltage (uV)'})
xlim([-0.05 0.4])
fig_title=sprintf('%s','Raw ERP - Normalized to mean');
title(strrep(fig_title,'_','\_'));

for i=1:size(ERP,1)
    legend_cell{1,i} = num2str(i);
end
legend(legend_cell);

reply = input('Do you want to remove any electodes? Y/N: ','s');
if reply == 'Y' || reply == 'y';
    electrode_position = input('Which electrode would you like to remove? ', 's');
    i = str2num(electrode_position);
    ERP(i,:)=nan;
    probemap_chrem(i,shank) = nan;
end

while reply ~= 'N' && reply ~= 'n'
    close all
    %display corrected
    norm_ERP = ERP;
    for i=1:length(ERP(:,1))
        if ~isnan(ERP(i,1))
            norm_ERP(i,:) = ERP(i,:)-nanmean(ERP,1);
        end
    end
    figure;plot (timing,ERP', 'DisplayName', sprintf('ERP_shk%d',shank), 'YDataSource', sprintf('ERP_shk%d',shank)); figure(gcf);
    xlabel({'Time from stimulus onset (s)'});
    ylabel({'Amplitude(uV)'});
    xlim([-0.05 0.4])
    fig_title=sprintf('%s','Electrode removed ERP - RAW');
    title(strrep(fig_title,'_','\_'));
    for i=1:size(ERP,1)
        legend_cell{1,i} = num2str(i);
    end
    legend(legend_cell);
    figure;plot (timing,norm_ERP', 'DisplayName', sprintf('ERP_shk%d',shank), 'YDataSource', sprintf('ERP_shk%d',shank)); figure(gcf);
    xlabel({'Time from stimulus onset (s)'});
    ylabel({'Amplitude(uV)'});
    xlim([-0.05 0.4])
    fig_title=sprintf('%s','Electrode removed ERP - Normalized to mean');
    title(strrep(fig_title,'_','\_'));
    for i=1:size(ERP,1)
        legend_cell{1,i} = num2str(i);
    end
    legend(legend_cell);
    % remove another electrode
    reply = input('Do you want to remove another electode? Y/N: ','s');
    if reply == 'Y' || reply == 'y'
        electrode_position = input('Which electrode would you like to remove? ', 's');
        i = str2num(electrode_position);
        ERP(i,:)=nan;
        probemap_chrem(i,shank) = nan;
    end
end

%display figure
close all
disp('Plot new and improved ERP and save mat file for CSD analysis');
figure;plot (timing,ERP', 'DisplayName', sprintf('ERP_shk%d',shank), 'YDataSource', sprintf('ERP_shk%d',shank)); figure(gcf);
xlabel({'Time from stimulus onset (s)'});
ylabel({'Amplitude(mV)'});
xlim([-0.05 0.4])
fig_title=sprintf('%s','Final ERP - RAW ');
title(strrep(fig_title,'_','\_'));
for i=1:size(ERP,1)
    legend_cell{1,i} = num2str(i);
end
legend(legend_cell);
figure;plot (timing,norm_ERP', 'DisplayName', sprintf('ERP_shk%d',shank), 'YDataSource', sprintf('ERP_shk%d',shank)); figure(gcf);
xlabel({'Time from stimulus onset (s)'});
ylabel({'Amplitude(mV)'});
xlim([-0.05 0.4])
fig_title=sprintf('%s','Final ERP - NORMALIZED ');
title(strrep(fig_title,'_','\_'));
for i=1:size(ERP,1)
    legend_cell{1,i} = num2str(i);
end
legend(legend_cell);
% end
Final_ERP=norm_ERP;
end