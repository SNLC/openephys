function population_analysis(proj, area, pop, exp_type)

% proj = string indicating which project (e.g., 'LP', or 'Tlx')
% area = e.g., 'LPlateral', 'LPmedial', 'LGN'
% pop = string indicating which population of interest (e.g., 'driver',
% 'modulator', 'ChR2', 'Halo', etc.)
% exp_type = 'ramp', 'trains', 'step', or 'size' 

%% set up directories, identify experiments to analyze
if strcmpi(proj,'LP')
    main_dir = 'H:\LPproject\LPresults';
    if strcmpi(pop,'driver')
        if strcmpi(exp_type,'step')&& strcmpi(area, 'LPlateral')
            exp_paths = {'H:\LPproject\D3\2017-05-25_17-18-54_LP_diffintensities',...
                'J:\LPproject\D22\2018-02-05_10-23-13_LP_diffintensities',...
                'J:\LPproject\D27\2018-04-27_16-00-16_LP_diffintensities'};
%                 'J:\LPproject\D28\2018-04-26_13-18-40_LP_diffintensities'};
%             shanks= {[0:3],[2 3],[0:3],[0]};        % hardcode here which shanks to look at for which experiments
                shanks= {[0:3],[2 3],[0:3]};        % hardcode here which shanks to look at for which experiments
                probe = '128D_bottom';   % doesn't matter if it was DN or D
                lightcond = 2;
        elseif strcmpi(exp_type,'step')&& strcmpi(area, 'LPmedial')
            exp_paths = {'J:\LPproject\D22\2018-02-05_10-23-13_LP_diffintensities'};
            shanks= {[0 1]};        % hardcode here which shanks to look at for which experiments
        elseif strcmpi(exp_type,'trains')&& strcmpi(area, 'LPlateral')
            % olde experiments
            exp_paths = {'H:\LPproject\D3\2017-05-25_18-23-19_LP_trains',...
                'J:\LPproject\D27\2018-04-27_17-04-08_LP_trains'};
            shanks={[0:3],[0:3]};
            probe = '128D_bottom';   % doesn't matter if it was DN or D
                lightcond = 2;
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area, 'LPlateral')
%             exp_paths = {'J:\LPproject\D26\2018-04-03_15-58-56_LP_diffintensities_inLP'};
            exp_paths = {'H:\LPproject\RBP10\2019-01-29_17-57-51_LP_diffintensitiesFULL_inLP',...
                'H:\LPproject\RBP11\2019-01-29_13-14-24_LP_diffintensities_inLP',...
                'H:\LPproject\RBP12\2019-01-28_11-56-30_LP_diffintensities_inLP',...
                'H:\LPproject\RBP4\2018-12-05_12-27-05_LP_diffintensities'};
            shanks = {[1],[1],[0,1],[0]};
            probe = '64G_bottom';  % doesn't matter if it was DN or D
            lightcond = 3; 
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area, 'LPrm')
%             exp_paths = {'J:\LPproject\D26\2018-04-03_15-58-56_LP_diffintensities_inLP'};
            exp_paths = {'H:\LPproject\RBP10\2019-01-29_17-57-51_LP_diffintensitiesFULL_inLP',...
                'H:\LPproject\RBP11\2019-01-29_13-14-24_LP_diffintensities_inLP',...
                'H:\LPproject\RBP12\2019-01-28_11-56-30_LP_diffintensities_inLP',...
                'H:\LPproject\RBP4\2018-12-05_12-27-05_LP_diffintensities'};
            shanks = {[0],[0],[0],[0]};
            probe = '64G_bottom';  % doesn't matter if it was DN or D
            lightcond = 3;
        elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
            exp_paths = {'J:\LPproject\NPRH18\LP\2019-11-12_18-06-34_LP_halo',...
                'J:\LPproject\NPRH19\LP\2019-11-13_17-32-45_LP_halo'};
            shanks = {[2:3],[2]};
            probe = '128D_bottom';
            lightcond = 1;
        end
    elseif strcmpi(pop,'modulator')
        if strcmpi(exp_type,'step') && strcmpi(area, 'LPlateral')
               %old experiments:
%             exp_paths = {'H:\LPproject\M7\M7_LP_3diffintensities_170203_155737\shank4',...
%     'H:\LPproject\M7\M7_LP_3diffintensities_170203_155737\shank3',...
%     'H:\LPproject\M10\M10_LP_3diffintensities_170211_144026\shank3',...
%         'H:\LPproject\M10\M10_LP_3diffintensities_170211_144026\shank2'};
% %         'J:\LPproject\M12\2017-06-09_16-57-03_LP_SFTF'};    % excluding for now because too much is different  
%         shanks = {[0], [0], [0],[0]};
                % new experiments:
                exp_paths = {'H:\LPproject\M17\2018-07-23_17-20-47_LP_diffintensities',...
                    'H:\LPproject\M18\2018-07-23_13-12-11_LP_diffintensities',...
                    'H:\LPproject\M20\2018-09-25_12-44-06_LP_diffintensities_real',...
                    'H:\LPproject\M21\2018-09-25_17-53-08_LP_diffintensities',...
                    'H:\LPproject\M7\M7_LP_diffintensities_170203_155737',...
                    'H:\LPproject\M10\M10_LP_3diffintensities_170211_144026'};
                shanks = {[3],[2,3],[1],[0,1],[2:3],[1:2]};
                probe = '128D_bottom';  % doesn't matter if it was DN or D
                lightcond = 2;  % which lightcond to analyze (excluding no light condition; eg., of 1hz, 10hz and 20hz conditions lightcond=2 will analyze 10hz condition
        elseif strcmpi(exp_type,'trains')  && strcmpi(area,'LPlateral')
            exp_paths = {'H:\LPproject\M17\2018-07-23_18-25-12_LP_trains',...
                'H:\LPproject\M18\2018-07-23_14-17-35_LP_trains',...
                'H:\LPproject\M20\2018-09-25_13-48-42_LP_trains',...
                'H:\LPproject\M21\2018-09-25_18-56-22_LP_trains',...
                'H:\LPproject\M7\M7_LP_trains_170203_172900',...
                'H:\LPproject\M10\M10_LP_trains_170211_154907'};
            shanks = {[3],[2,3],[1],[0,1],[2:3],[1:2]};
                probe = '128D_bottom';   % doesn't matter if it was DN or D
                lightcond = 2;
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area, 'LPlateral')
            % old experiments
%             exp_paths = {'J:\LPproject\M16\2018-02-09_19-07-52_LP_diffintensities_inLP',...
%                 'J:\LPproject\M14\2018-02-10_14-15-29_LP_diffintensities_inLP'};
%             shanks = {[0],[0]};
            exp_paths = {'H:\LPproject\M22\2018-10-09_18-06-10_LP_diffintensities',...
                'H:\LPproject\M23\2018-10-09_13-09-21_LP_diffintensities',...
                'H:\LPproject\M25\2018-12-07_12-28-16_LP_diffintensities_inLP'};
%             'H:\LPproject\M24\2018-10-10_12-24-59_LP_diffintensities',...      % LP looks damaged in this experiment
            shanks = {[0],[0],[0]};
            probe = '64G_bottom';  
            lightcond = 3;
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area, 'LPlateral_only')
            % older experiments when I was FULLY in LP and using higher
            % light powers
            exp_paths = {'H:\LPproject\M14\2018-02-10_14-15-29_LP_diffintensities_inLP',...
                'H:\LPproject\M15\2018-02-10_17-31-47_LP_diffintensities_inLP',...
                'H:\LPproject\M16\2018-02-09_19-07-52_LP_diffintensities_inLP'};
            shanks = {[0:1],[0:1],[0:1]};
            probe = '64G_bottom';  
            lightcond = 2;
        elseif strcmpi(exp_type,'trains_inLP') && strcmpi(area,'LPlateral')
            exp_paths = {'H:\LPproject\M22\2018-10-09_19-10-35_LP_trains',...
                'H:\LPproject\M23\2018-10-09_14-13-12_LP_trains',...
                'H:\LPproject\M25\2018-12-07_13-37-55_LP_trains_inLP'}; % **so far not including M24
            shanks = {[0],[0],[0]};
            probe = '64G_bottom';  
            lightcond = 2;
        elseif strcmpi(exp_type,'trains_inLP') && strcmpi(area,'LPlateral_only')
            exp_paths = {'H:\LPproject\M14\2018-02-10_15-21-20_LP_trains_inLP',...
                'H:\LPproject\M15\2018-02-10_18-34-57_LP_trains_inLP',...
                'H:\LPproject\M16\2018-02-09_20-10-37_LP_trains_inLP'}; % 
            shanks = {[0:1],[0:1],[0:1]};
            probe = '64G_bottom';  
            lightcond = 2;
        elseif strcmpi(exp_type,'step_halo') && strcmpi(area, 'LPlateral')
            % these are the ones with too high light power...
%             exp_paths = {'J:\LPproject\MH15\2018-05-23_16-36-18_LP_Halo',...      
%                 'J:\LPproject\MH12\2018-05-09_18-07-28_Halo_LP',...
%                 'J:\LPproject\MH11\2018-05-09_12-08-51_Halo_LP',...
%                 'H:\LPproject\MH9\2018-03-27_19-39-26_Halo_LP',...
%                 'H:\LPproject\MH7\LP\2018-03-28_12-36-06_Halo_LP'};
%                 shanks = {[1 2],[2 3],[1],[1],[3]};

        exp_paths = {...
        'H:\LPproject\MH18\LP\2018-08-24_14-49-10_LP_Halo',...   
            'J:\LPproject\MH19\LP\2018-12-04_21-43-54_LP_Halo',...
                'J:\LPproject\MH20\LP\2018-12-04_14-05-51_LP_Halo',...
                'J:\LPproject\MH25\LP\2019-04-10_17-05-08_LP_Halo',...
                'J:\LPproject\MH24\LP\2019-04-12_13-33-18_LP_halo_real'};
            shanks = {[0],[0],[0],[1],[3]};
            probe = '128D_bottom';  % doesn't matter if it was DN or D
            lightcond = 2;
        elseif strcmpi(exp_type,'step') && strcmpi(area,'LPmedial')
           exp_paths =  {'H:\LPproject\M10\M10_LP_3diffintensities_170211_144026',...
               'H:\LPproject\M17\2018-07-23_17-20-47_LP_diffintensities',...
                    'H:\LPproject\M18\2018-07-23_13-12-11_LP_diffintensities',...
                    'H:\LPproject\M20\2018-09-25_12-44-06_LP_diffintensities_real'};
        shanks = {3, 0:2, 0:1,0};
        probe = '128D_bottom';
        lightcond = 2;
        elseif strcmpi(exp_type,'trains') && strcmpi(area,'LPmedial')
           exp_paths =  {'H:\LPproject\M10\M10_LP_trains_170211_154907',...
               'H:\LPproject\M17\2018-07-23_18-25-12_LP_trains',...
                'H:\LPproject\M18\2018-07-23_14-17-35_LP_trains',...
                'H:\LPproject\M20\2018-09-25_13-48-42_LP_trains'};
        shanks = {3, 0:2, 0:1,0};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightcond = 2;
        elseif strcmpi(exp_type,'step') && strcmpi(area, 'LGN')
                % OLD experiments:
%             exp_paths = {'H:\LPproject\M7\M7_LP_3diffintensities_170203_155737\shank1',...
%     'H:\LPproject\M10\M10_LP_3diffintensities_170211_144026\shank1'};
%             shanks = {[0], [0]};
                % NEW experiments:
             exp_paths = {'H:\LPproject\M20\2018-09-25_12-44-06_LP_diffintensities_real',...
                'H:\LPproject\M21\2018-09-25_17-53-08_LP_diffintensities',...
                'H:\LPproject\M36\LP\2019-02-21_12-39-57_LP_diffintensities',...
                'H:\LPproject\M7\M7_LP_diffintensities_170203_155737'};
            shanks = {[3],[2,3],[0:2],[0]};     % 9/13/19 - changed to EXCLUDE last M36 shank which was outside terminal field
            probe = '128D_bottom';  % doesn't matter if it was DN or D
            lightcond = 2;
        elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LGN')
                % these are the ones with too high light power...
%             exp_paths = {'J:\LPproject\MH15\2018-05-23_16-36-18_LP_Halo',...
%                 'J:\LPproject\MH11\2018-05-09_12-08-51_Halo_LP',...
%                 'H:\LPproject\MH9\2018-03-27_19-39-26_Halo_LP'};
%                 shanks = {[3],[2 3],[2 3]};
            exp_paths = {...
                'H:\LPproject\MH18\LP\2018-08-24_14-49-10_LP_Halo',... 
                'J:\LPproject\MH19\LP\2018-12-04_21-43-54_LP_Halo',...
                'J:\LPproject\MH20\LP\2018-12-04_14-05-51_LP_Halo',...
                'J:\LPproject\MH21\LP\2018-12-03_17-56-28_LP_Halo',...
                'J:\LPproject\MH25\LP\2019-04-10_17-05-08_LP_Halo',...
                'J:\LPproject\MH22\LP\2019-04-26_14-46-21_LP_Halo',...
                'H:\LPproject\MH23\LP\2019-04-11_14-40-06_LP_Halo'};
            shanks = {[1],[1,2],[1,2],[2],[2],[0:3],[0:1] };
            probe = '128D_bottom';  % doesn't matter if it was DN or D
            lightcond = 2;
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area,'LGN')
            exp_paths = {'H:\LPproject\M22\2018-10-09_18-06-10_LP_diffintensities',...
                'H:\LPproject\M23\2018-10-09_13-09-21_LP_diffintensities',...
                'H:\LPproject\M25\2018-12-07_12-28-16_LP_diffintensities_inLP',...
                'H:\LPproject\M26\2018-12-06_16-06-11_LP_diffintensities_real_inLP'};
            shanks = {[1],[1],[1],[0:1]};
            probe = '64G_bottom';  
            lightcond = 3;
        elseif strcmpi(exp_type,'trains_inLP') && strcmpi(area,'LGN')
            exp_paths = {'H:\LPproject\M22\2018-10-09_19-10-35_LP_trains',...
                'H:\LPproject\M23\2018-10-09_14-13-12_LP_trains',...
                'H:\LPproject\M25\2018-12-07_13-37-55_LP_trains_inLP',...
                'H:\LPproject\M26\2018-12-06_17-15-45_LP_trains_inLP'};
            shanks = {[1],[1],[1],[0:1]};
            probe = '64G_bottom';  
            lightcond = 2;
        elseif strcmpi(exp_type,'trains')  && strcmpi(area,'LGN')
            % old...
%             exp_paths = {'H:\LPproject\M7\M7_LP_trains_170203_172900'};         % add more?!
%             shanks = {[2 3]};
            exp_paths = {'H:\LPproject\M20\2018-09-25_13-48-42_LP_trains',...
                'H:\LPproject\M21\2018-09-25_18-56-22_LP_trains',...
                'H:\LPproject\M36\LP\2019-02-21_13-42-50_LP_trains',...
                'H:\LPproject\M7\M7_LP_trains_170203_172900'};
            shanks = {[3],[2,3],[0:2],[0]}; % 9/13/19 - changed to EXCLUDE last M36 shank which was outside terminal field
            probe = '128D_bottom';  % doesn't matter if it was DN or D
            lightcond = 2;
        
        elseif strcmpi(exp_type,'step') && strcmpi(area,'TRN')
            exp_paths = {'J:\LPproject\M27\2018-12-17_12-32-20_TRN_diffintensities',...
                'J:\LPproject\M32\TRN_pen2\2019-02-13_19-35-45_TRN_diffintensities_pen2',...
                'J:\LPproject\M33\2019-02-13_12-58-05_TRN_diffintensities',...
                'J:\LPproject\M32\TRN_pen1\2019-02-13_16-54-13_TRN_diffintensitiesFULL'};
%                 'J:\LPproject\M29\2018-12-18_12-21-53_TRN_diffintensities'};

            shanks = {[0],[0],[0],[0]};
            probe = '64D_bottom';
            lightcond = 2;
        elseif strcmpi(exp_type,'trains') && strcmpi(area,'TRN')
            exp_paths = {'H:\LPproject\M27\2018-12-17_13-35-15_TRN_trains',...
                'J:\LPproject\M32\TRN_pen2\2019-02-13_20-32-58_TRN_trains_pen2',...
                'J:\LPproject\M33\2019-02-13_14-01-38_TRN_trains',...
                'H:\LPproject\M32\TRN_pen1\2019-02-13_17-50-52_TRN_trains'};    % trials cut off early, but still looks good (but requires commenting out hotellings test for tuning b/c uneven number of orientation trials)
            % 'J:\LPproject\M29\2018-12-18_13-24-55_TRN_trains',...
            shanks = {[0],[0],[0],[0]};
            probe = '64D_bottom';
            lightcond = 2;
        end
    elseif strcmpi(pop,'CTRL')
        if strcmpi(exp_type,'step_inLP') && strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\MCTRL1\2018-10-08_14-31-08_LP_diffintensities'};  % i think I was actually in LD? 
            shanks = {[0]};
        elseif strcmpi(exp_type,'step_inLP') && strcmpi(area,'LGN')
%             exp_paths =
%             {'J:\LPproject\MCTRL1\2018-10-08_14-31-08_LP_diffintensities'};
%             % this experiment was crappy
            exp_paths = {'H:\LPproject\MCTRL4\2019-02-22_18-20-07_LP_diffintensities_inLP'};
            shanks = {[0,1]};
            probe = '64G_bottom';
            lightcond = 3;
        elseif strcmpi(exp_type,'trains_inLP') && strcmpi(area,'LGN')
            exp_paths = {'H:\LPproject\MCTRL4\2019-02-22_19-23-28_LP_trains_inLP'};
            shanks = {[0,1]};
            probe = '64G_bottom';
            lightcond = 3;
        elseif strcmpi(exp_type,'step') && strcmpi(area,'LPlateral')
            exp_paths = {'H:\LPproject\MCTRL2\2019-01-21_14-58-23_LP_diffintensities',...
                'H:\LPproject\MCTRL3\LP\2019-02-22_14-37-55_LP_diffintensities'};
            shanks = {[0],[2:3]};
            probe = '128DN_bottom';
            lightcond = 2;
        elseif strcmpi(exp_type,'trains') && strcmpi(area,'LPlateral')
            exp_paths = {'H:\LPproject\MCTRL2\2019-01-21_16-02-55_LP_trains',...
                'H:\LPproject\MCTRL3\LP\2019-02-22_15-40-44_LP_trains'};
            shanks = {[0],[2:3]};
            probe = '128DN_bottom';
            lightcond = 2;
        elseif strcmpi(exp_type,'step') && strcmpi(area,'LGN')
            exp_paths = {'H:\LPproject\MCTRL2\2019-01-21_14-58-23_LP_diffintensities'};
            shanks = {[1:3]};
            probe = '128DN_bottom';
            lightcond = 2;
        elseif strcmpi(exp_type,'trains') && strcmpi(area,'LGN')
            exp_paths = {'H:\LPproject\MCTRL2\2019-01-21_16-02-55_LP_trains'};
            shanks = {[1:3]};
            probe = '128DN_bottom';
            lightcond = 2;
        end
    elseif strcmpi(pop,'ChR2')
    elseif strcmpi(pop,'Halo')
    elseif strcmpi(pop,'PVChR2')
%         if strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\VH2\2018-05-16_19-05-31_LP_PVChR2',...
%                 'J:\LPproject\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
%             shanks = {[1],[2 3]};
        if strcmpi(area,'LPlateral')
            exp_paths = {'J:\LPproject\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
            shanks = {[2]};
        elseif strcmpi(area,'LPlateral_caudal')
            exp_paths = {'J:\LPproject\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
            shanks = {[3]};
        elseif strcmpi(area,'LPrm')
            exp_paths = {'J:\LPproject\VH2\2018-05-16_19-05-31_LP_PVChR2'};
            shanks = {[0:1]};
        elseif strcmpi(area,'LGN')
            exp_paths = {'J:\LPproject\VH2\2018-05-16_19-05-31_LP_PVChR2'};
            shanks = {[3]};
        elseif strcmpi(area,'LD')
            exp_paths = {'J:\LPproject\VH2\2018-05-16_19-05-31_LP_PVChR2',...
                'J:\LPproject\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
            shanks = {[2],[0 1]};
        elseif strcmpi(exp_type,'step') && strcmpi(area,'LPmedial')
            exp_paths = {'J:\LPproject\VH2\2018-05-16_19-05-31_LP_PVChR2'};
            shanks = {[0]};
        end
    end
elseif strcmpi(proj, 'Tlx')
    main_dir = sprintf('H:\\Tlx3project\\%sresults',pop);
    if strcmpi(pop,'halo')
            exp_paths = {'J:\TH8\Day2\2018-05-08_11-55-38_V1_Halo_pen2', 'H:\Tlx3project\TH10\2018-05-08_15-25-21_V1_Halo'};
            shanks = {[0],[0]};
    end
end

% set up colors
if contains(area,'lgn','ignorecase',1)
    area_color = [.85 .325 .098];   % LGN=orange
elseif contains(area,'LP','ignorecase',1)
    area_color = [.494 .184 .556];   % LP = purple
elseif contains(area,'TRN','ignorecase',1)
    area_color = [.5 .5 .5];
end

if ~exist(main_dir,'dir')
    mkdir(main_dir)
end
cd(main_dir)

fig_dir =  sprintf('%s\\%s_%s_%s_figs',main_dir,area,pop,exp_type);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

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

for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
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
    exp_dir = strcat(main_dir,'\',exp_name);
    if ~exist(exp_dir,'dir') && exist(sprintf('%s\\%s\\%s',main_dir,an_name,exp_name))       % need better solution for this
        exp_dir = sprintf('%s\\%s\\%s',main_dir,an_name,exp_name);
    end
    fprintf(sprintf('Processing experiment %s\n',exp_name))
%     if ~exist(exp_dir,'dir')
%     if i > 4
% %         get_lightstim_v2(exp_path,exp_type)
%         analysis_master(exp_path,'OpenEphys','step',0)
% % if i > 10
% %         get_lightstim_v2(exp_path,'step')
%         intan_analysis_master(exp_path,'step',1);
% end
% 
%     end
    
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
        cluster_file = sprintf('%s\\%s\\%s_cluster_quality.mat',main_dir,an_name,an_name);
    end
    exp = importdata(results_file,'-mat');     % load results mat
    clust = importdata(cluster_file,'-mat');
    clear cluster_file results_file
    
%     for sh = 1:length(shanks{i})    
%         if exist(sprintf('%s/shank%d_LPchannels.txt',exp_path,shanks{i}(sh)),'file')
%             channels{i}{sh} = load(sprintf('%s/shank%d_LPchannels.txt',exp_path,shanks{i}(sh)));   % if you predetermined which clusters to look at
%         elseif exist(sprintf('%s/good_channels.txt',exp_path),'file')
%             channels{i} = load(sprintf('%s/good_channels.txt',exp_path));   % if you predetermined which clusters to look at
%         end
%     end

    exp_num = [exp_num i*ones(1,length(exp.unitinfo))];
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
end


%% calculate relevant significance values
% using kruskal-wallis (non-parametric, indep samples) for testing vis and
% light modulation, using hotellings for significant orientation tuning.
% Are these the right tests to use??
for n = 1:length(params)
    if strcmp(params(n).IVs,'s_freq')
        vis_trials{n} = find((params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1)&(params(n).trial_type(:,strcmpi(params(n).IVs,'s_freq'))==.04)&(params(n).trial_type(:,strcmpi(params(n).IVs,'t_period'))==30));  % in case of multiple SFs and TFs, only compare regular 2Hz and .04cpd trials (for now)
    else
        vis_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1);
    end
    blank_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==0);
    nolight_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==0);
    lightconds{n} = unique(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit')));
%     if length(lightconds{n}) > 3
%         lightconds{n}= lightconds{n}([1 2 end],:);    % if experiment had low, medium and high intensity light conditions, drop the medium condition (b/c M12 only has low and high)
%     end
    for lc=2:length(lightconds{n})  % assumes first condition is no-light condition
        light_trials{n}{lc-1} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==lightconds{n}(lc));
        vislight_trials{n}{lc-1} = intersect(vis_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==lightconds{n}(lc)));
        blanklight_trials{n}{lc-1} = intersect(blank_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==lightconds{n}(lc)));
    end
    run_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==1);
    stat_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==0);
end

% all_units = 1:length(unitinfo);
% units = all_units;
% distfromlastch = nan(1,length(units));
% distfromfirstch = distfromlastch;
% for n = 1:length(params)
%     if iscell(channels{n})
%         shk = [waveforms(exp_num==n).shank];
%         shks = shanks{n};
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
%     oris = unique(params(exp_num(nn)).trial_type(:,strcmp(params(exp_num(nn)).IVs,'ori')));
%         oris(oris>=999) = [];
%     % kruskal-wallis test to test for significant and visual- and light-modulation
%     % for visual modulation, find preferred direction trials - only use THESE
%     % trials to test for significant visual modulation (in case of extremely
%     % tuned cells)
%     [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-repmat(FRb(i),1,size(tuning_curve,3))));   % direction w/ biggest change from baseline
%     prefori_trials{exp_num(nn)}(i,:) = find(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori'))==oris(prefdir_deg(i)));
%     vis_sig(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials{exp_num(nn)}(i,:),nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2)'...
%         sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2)'],...
%         [ones(1,length(intersect(prefori_trials{exp_num(nn)}(i,:),nolight_trials{exp_num(nn)}))) 2*ones(1,length(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)})))],'off');    % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
%     vis_sig_ons(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials{exp_num(nn)}(i,:),nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2)'... 
%         sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2)'],...
%         [ones(1,length(intersect(prefori_trials{exp_num(nn)}(i,:),nolight_trials{exp_num(nn)}))) 2*ones(1,length(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)})))],'off');         % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
%     for lc = 1:length(lightconds{exp_num(nn)})-1
%         light_sig(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2)'...       % currently using ALL light trials to evaluate light significance (visual+blank)
%             sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2)'],...
%             [ones(1,length(nolight_trials{exp_num(nn)})) 2*ones(1,length(light_trials{exp_num(nn)}{lc}))],'off');    % significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
%         light_sig_ons(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2)'...
%             sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2)'],...
%             [ones(1,length(nolight_trials{exp_num(nn)})) 2*ones(1,length(light_trials{exp_num(nn)}{lc}))],'off');         % significance of light-evoked response at light onset (100ms after light onset in nolight vs light trials (each condition separately))
%     end
%     
%    % next, check tuning significance and get tuning curves
%      if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
%         % evaluate significance or orientation tuning
%         for o = 1:length(oris)/2
%             for lc = 1:length(lightconds{exp_num(nn)})
%                 ori_trials = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==lightconds{exp_num(nn)}(lc)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
%                 tuning_trials{exp_num(nn)}(:,o,lc) = sum(unitinfo(nn).rast(ori_trials,round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2); % numbers of spikes across trials of given light condition for each orientation (by column)
% 
%             end
%             tuning_curve_collapse(i,:,o)   = mean([tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)],2);         % average evoked FR of orientations collapsed across directions
% %             tuning_curve(i,:,[o o+length(oris)/2]) = [tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)];
%         end
%         for lc = 1:length(lightconds{exp_num(nn)})      % currently, NOT separating running and stationary trials
%             tuned_sig(i,lc) = T2Hot1(tuning_trials{exp_num(nn)}(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
%         end
%     end
% end

% NEW (6/15/18) - calculate visual significance prior to getting distances from first
% and last ch and deciding which units are "clean". Use first and last
% visually significant channels to determine borders of LP

incl_units = zeros(1,length(unitinfo));
if isfield(waveforms,'shank')
    shk = [waveforms.shank];
    if isempty(shk); shk = zeros(1,length(unitinfo)); end    % if single shank probe, shank field may be empty
else
    shk = zeros(1,length(unitinfo));
end

for n = 1:length(params)
    which_units = find(exp_num==n);
    incl_units(which_units(ismember(shk(exp_num==n),shanks{n}))) = 1;
end
% good_SNR = find((incl_units)&(SNR>=1.5&refr_idx<.1)); % only include units that pass SNR and refractory period thresholds
 good_SNR = find((incl_units)&(refV<.5));
% for n=1:length(params)
% %     isiV{n} = sqKilosort.isiViolations(fileparts(exp_paths{n}));
%     [~, uQ{n}, cR{n}, isiV{n}] = sqKilosort.computeAllMeasures(fileparts(exp_paths{n}));
%     uQ{n} = uQ{n}';
%     cR{n} = cR{n}';
% end
% isiV = [isiV{:}];
% uQ = [uQ{:}];
% cR = [cR{:}];
% good_isi = find(isiV<.1);       % only include units with <10% ISI violations
good_uQ = find(uQ>16);       % only include units with <25% contamination "false positive" rate
% clean_units = intersect(good_SNR,good_isi(cR<.3 | isnan(cR)));   % and exclude units with 30% or more contamination rate of other good_isi units (include NANs because doesn't necessarily mean they're bad)
clean_units = intersect(good_SNR,good_uQ);    
incl_units = find(incl_units);      % NEW MAK 5/6/19

FRb = [FRs(incl_units).baseline];  %  baseline (from blank trials w/ no running, light or vis stim)
for i = 1:length(incl_units)      % for each unit
    nn = incl_units(i);
    tuning_curve{i} = tuning(nn).curve(:,:);
    oris = unique(params(exp_num(nn)).trial_type(:,strcmp(params(exp_num(nn)).IVs,'ori')));
        oris(oris>=999) = [];
    % wilcoxen rank-sum test to test for significant and visual- and light-modulation
    % for visual modulation, find preferred direction trials - only use THESE
    % trials to test for significant visual modulation (in case of extremely
    % tuned cells)
    [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-repmat(FRb(i),1,size(tuning_curve,3))));   % direction w/ biggest change from baseline
    prefori_trials{i} = find(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori'))==oris(prefdir_deg(i)));
    vis_sig(i) = ranksum(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2),... % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
        sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2));
    vis_sig_ons(i) = ranksum(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2),...  % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
        sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2));
    for lc = 1:length(lightconds{exp_num(nn)})-1
        light_sig(i,lc) = ranksum(sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2),...       % currently using ALL light trials to evaluate light significance (visual+blank)
            sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2));% significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
        light_sig_ons(i,lc) = ranksum(sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2),...
            sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2)); % significance of light-evoked response at light onset (100ms after light onset in nolight vs light trials (each condition separately))
    end
%   % previously using kruskalwallis - changed 7/21/19 MAK
%     vis_sig(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2)'...
%         sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2)'],...
%         [ones(1,length(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}))) 2*ones(1,length(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)})))],'off');    % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
%     vis_sig_ons(i) = kruskalwallis([sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2)'... 
%         sum(unitinfo(nn).rast(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+params(exp_num(nn)).onset)),2)'],...
%         [ones(1,length(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}))) 2*ones(1,length(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)})))],'off');         % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
%     for lc = 1:length(lightconds{exp_num(nn)})-1
%         light_sig(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2)'...       % currently using ALL light trials to evaluate light significance (visual+blank)
%             sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)))+params(exp_num(nn)).lighttime*1000),2)'],...
%             [ones(1,length(nolight_trials{exp_num(nn)})) 2*ones(1,length(light_trials{exp_num(nn)}{lc}))],'off');    % significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
%         light_sig_ons(i,lc) = kruskalwallis([sum(unitinfo(nn).rast(nolight_trials{exp_num(nn)},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2)'...
%             sum(unitinfo(nn).rast(light_trials{exp_num(nn)}{lc},round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+params(exp_num(nn)).onset))),2)'],...
%             [ones(1,length(nolight_trials{exp_num(nn)})) 2*ones(1,length(light_trials{exp_num(nn)}{lc}))],'off');         % significance of light-evoked response at light onset (100ms after light onset in nolight vs light trials (each condition separately))
%     end
    
   % next, check tuning significance and get tuning curves
     if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
        % evaluate significance or orientation tuning
        for o = 1:length(oris)/2
            for lc = 1:length(lightconds{exp_num(nn)})
                ori_trials = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==lightconds{exp_num(nn)}(lc)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
%                 if length(ori_trials)>40      % not sure what's the purpose of this?? (MAK 5/8/19)
%                     ori_trials(randi(length(ori_trials),1)) = [];
%                 end
                tuning_trials{exp_num(nn)}(:,o,lc) = sum(unitinfo(nn).rast(ori_trials,round(1000*(params(exp_num(nn)).av_light_start(1)))+1:round(1000*(params(exp_num(nn)).av_light_start(1)))+params(exp_num(nn)).lighttime*1000),2); % numbers of spikes across trials of given light condition for each orientation (by column)
                
            end
            tuning_curve_collapse(i,:,o)   = mean([tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)],2);         % average evoked FR of orientations collapsed across directions
%             tuning_curve(i,:,[o o+length(oris)/2]) = [tuning(nn).curve([1 2 end],o) tuning(nn).curve([1 2 end],o+length(oris)/2)];
        end
        for lc = 1:length(lightconds{exp_num(nn)})      % currently, NOT separating running and stationary trials
            tuned_sig(i,lc) = T2Hot1(tuning_trials{exp_num(nn)}(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
        end
%         tuned_sig = zeros(length(incl_units),4);    % hack for when ori_trials are uneven
    end
end

distfromlastch = nan(1,length(incl_units));
distfromfirstch = distfromlastch;
count=1;
vis_units = find(vis_sig<.025|vis_sig_ons<.025);      % CHANGED 7/21/19 (.025 b/c bonferroni correction for 2 tests)
% vis_units_strict = find(vis_sig<.01|vis_sig_ons<.01);

% get probe info (NEW 3/26/19)
p = eval(sprintf('probemap_%s_func',probe));
Zchan = flipud(sort(p.z(p.shaft==1))); % from top to bottom (bottom=0)
for n = 1:length(params)        % for each exp
    for sh = 1:length(shanks{n})    % for each shank in exp
        shk_units{count} = find((exp_num(incl_units)==n) & (shk(incl_units)==shanks{n}(sh)));
        vischs = sort([waveforms(incl_units(intersect(shk_units{count},vis_units))).max_ch]);
        firstch(count) = min(vischs);
        if firstch > 1
            if Zchan(firstch(count))==Zchan(firstch(count)-1) % for probes in hexagonal orientation, might leave out channel that is actually same height as "firstch"
                firstch(count) = firstch(count)-1;
            end
        end
        lastch(count) = max(vischs);
        if strcmpi(area,'trn') && sum(abs(diff(Zchan(vischs))) > 100)      % if more than 100um separates consecutively located visually-responsive units in TRN experiment 
            lastinTRN = find(abs(diff(Zchan(vischs)))> 100,1,'first');
            lastch(count) = vischs(lastinTRN);
        elseif strcmpi(area,'trn') && sum(Zchan(vischs)-Zchan(vischs(1))<=-400)     % don't include more than 400um of units for TRN
            lastinTRN = find(Zchan(vischs)-Zchan(vischs(1)) <= -400,1,'first');
            lastch(count) = vischs(lastinTRN-1);
        end
        if lastch < length(Zchan)
            if Zchan(lastch(count))==Zchan(lastch(count)+1)
                lastch(count) = lastch(count)+1;
            end
        end
        distfromlastch(shk_units{count}) = Zchan(lastch(count))-Zchan([waveforms(incl_units(shk_units{count})).max_ch]);   % negative values are actually good (above last ch)
        distfromfirstch(shk_units{count}) = Zchan(firstch(count))-Zchan([waveforms(incl_units(shk_units{count})).max_ch]); % positive values are below first ch
        count = count+1;
    end
end
unit_chk = 1:length(incl_units);
unit_chk(distfromlastch>0|distfromfirstch<0|isnan(distfromlastch)) = [];
clean_units = intersect(clean_units,incl_units(unit_chk));

%% NEW 7/21/19 - check "clean_units" across experiments
files = dir(main_dir);
files = {files.name};
if contains(exp_type,'inLP')
    othermat = cellfun(@(x) contains(x,strcat(area,'_',pop))&&contains(x,'inLP')&&contains(x,'.mat'),files);
else
    othermat = cellfun(@(x) contains(x,strcat(area,'_',pop))&&~contains(x,'inLP')&&~contains(x,'halo')&&contains(x,'.mat'),files);
end
if contains(exp_type,'halo')    % temp hack - currently don't need to corroborate across experiments for any halo experiments?
    othermat = [];
end
clean_units_expnum = exp_num(clean_units);
save(sprintf('%s/%s_%s_%s_cleanunits.mat',main_dir,area,pop,exp_type),'clean_units','clean_units_expnum')    % save after locating other mats and before changing clean_units
others = find(othermat);
for f=1:sum(othermat)
    other_clean{f} = load(fullfile(main_dir,files{others(f)})); % doesn't matter if current .mat is included since it's intersecting
    if length(unique(clean_units_expnum))~= length(unique(other_clean{f}.clean_units_expnum))   % TEMP - currently assumes that if there are inequal experiment numbers, the last experiment is the one that doesn't have both experimental runs
        incl_expnums = min(length(unique(clean_units_expnum)),length(unique(other_clean{f}.clean_units_expnum)));
        clean_units = [intersect(other_clean{f}.clean_units(other_clean{f}.clean_units_expnum<=incl_expnums),clean_units(clean_units_expnum<=incl_expnums)) clean_units(clean_units_expnum>incl_expnums)];
    else
        clean_units = intersect(other_clean{f}.clean_units,clean_units);
    end
end
% adjust for new clean_units
distfromlastch = distfromlastch(ismember(incl_units,clean_units));          % only includes GOOD units
distfromfirstch = distfromfirstch(ismember(incl_units,clean_units));          % only includes GOOD units
tuning_curve=tuning_curve(ismember(incl_units,clean_units));
prefdir_deg = prefdir_deg(ismember(incl_units,clean_units));
prefori_trials = prefori_trials(ismember(incl_units,clean_units));
vis_sig = vis_sig(ismember(incl_units,clean_units));
vis_sig_ons = vis_sig_ons(ismember(incl_units,clean_units));
light_sig = light_sig(ismember(incl_units,clean_units),:);
light_sig_ons = light_sig_ons(ismember(incl_units,clean_units),:);
tuning_curve_collapse = tuning_curve_collapse(ismember(incl_units,clean_units),:,:);
tuned_sig = tuned_sig(ismember(incl_units,clean_units),:);
FRb = FRb(ismember(incl_units,clean_units));

%% test for light modulation
% verified this combo of reshape, cell2mat and arrayfun yields accurate
% results! MAK 1/31/18
numconds = arrayfun(@(x) length(unique(x.all_light)),params,'uniformoutput',1);     % number of light conds in each experiment
conds = 1:min(numconds)-1;
FRev = reshape(cell2mat(arrayfun(@(x) x.visual.ev(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';    %vis-evoked
FRearly = reshape(cell2mat(arrayfun(@(x) x.visual.evstart(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';   % early evoked period
FRlate = reshape(cell2mat(arrayfun(@(x) x.visual.evlate(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))';  % late evoked period
FRonset = reshape(cell2mat(arrayfun(@(x) x.visual.evlightonset(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % vis-evoked light onset
FRbl = reshape(cell2mat(arrayfun(@(x) x.blank.ev(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % blank, evoked period
FRvison = reshape(cell2mat(arrayfun(@(x) x.onset(1,[conds end]), FRs(clean_units),'uniformoutput',0)),length(conds)+1,length(clean_units))'; % vis stim onset
for ii = 2:size(FRev,2)
    lightmod(:,ii-1) = (diff(FRev(:,[1 ii]),[],2))./sum(FRev(:,[1 ii]),2);
    lightmod_early(:,ii-1) = (diff(FRearly(:,[1 ii]),[],2))./sum(FRearly(:,[1 ii]),2);
    lightmod_late(:,ii-1) = (diff(FRlate(:,[1 ii]),[],2))./sum(FRlate(:,[1 ii]),2);
    lightmod_onset(:,ii-1) = (diff(FRonset(:,[1 ii]),[],2))./sum(FRonset(:,[1 ii]),2);
    lightmod_bl(:,ii-1) = (diff(FRbl(:,[1 ii]),[],2))./sum(FRbl(:,[1 ii]),2);
end

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
    [t2p_t(i),t2p_r(i),fwhm(i)] = get_waveform_props(waveforms(clean_units(i)).microV,params(exp_num(clean_units(i))).amp_sr);
end

% %% NEW test: lightmod for trains experiments
cd(fig_dir)
all_lightconds = [lightconds{:}];
if contains(exp_type,'trains')    
    % check if different experiments used same light trains conditions
%     if length(unique(max(all_lightconds)))==1    % if max lightcond was uniform across experiments (this is kinda a hack for M7 experiment...)
%        highf_cond = 2;
%     else
%        highf_cond = 1;
%     end
    for ii=1:sum(mean(all_lightconds,2)>1)      % for  light conditions >1Hz
        cons_exps = find(all_lightconds(ii+2,:)==mode(all_lightconds(ii+2,:)));  % which experiments used consistent trains frequencies per condition
        clean_units_trains = clean_units(ismember(exp_num(clean_units),cons_exps));
        resp_t_v{ii} = nan(length(clean_units_trains),lightconds{1}(ii+2));
        ppr_v{ii} = nan(length(clean_units_trains),lightconds{1}(ii+2)-1);    
        resp_dur_v{ii} = nan(length(clean_units_trains),lightconds{1}(ii+2));
        resp_t_bl{ii} = resp_t_v{ii};
        ppr_bl{ii} = ppr_v{ii};
        ppr2_v{ii} = ppr_v{ii};
        ppr2_bl{ii} = ppr_v{ii};
        resp_dur_bl{ii} = resp_dur_v{ii};
        for i = 1:length(clean_units_trains)
             vis_inds = (params(exp_num(clean_units_trains(i))).trial_type(:,1)==1);
            % fourier analysis
            binsize = .005;     % 5ms
            light_bins = (round(mean(params(exp_num(clean_units_trains(i))).av_light_start)))/binsize+1:(round(mean(params(exp_num(clean_units_trains(i))).av_light_start))+max(params(exp_num(clean_units_trains(i))).light_dur))/binsize;   % bins from which to extract threshold (light period)
            [~,psth] = make_psth_v2(binsize,0:binsize:params(exp_num(clean_units_trains(i))).stimtime+params(exp_num(clean_units_trains(i))).prestim,1:sum(vis_inds),unitinfo(clean_units_trains(i)).rast(vis_inds,:),params(exp_num(clean_units_trains(i))).all_light(vis_inds));  
            [~,psth_bl] = make_psth_v2(binsize,0:binsize:params(exp_num(clean_units_trains(i))).stimtime+params(exp_num(clean_units_trains(i))).prestim,1:sum(~vis_inds),unitinfo(clean_units_trains(i)).rast(~vis_inds,:),params(exp_num(clean_units_trains(i))).all_light(~vis_inds));  
            F_Z(ii,i) = calc_Ftrains(psth(ii+2,light_bins)',.005,lightconds{1}(ii+2));
            F_Zbl(ii,i) = calc_Ftrains(psth_bl(ii+2,light_bins)',.005,lightconds{1}(ii+2));
            
%             if F_Z(ii,i) > 3      % tried using this to determine
%             "hz-activated" units, but it captured a lot of units without
%             clear, quantifiable spike outputs
%                 [resp_t_v{ii}(i,:), resp_dur_v{ii}(i,:), ppr_v{ii}(i,:)]  =  ppanalysis_v2(psth(ii+2,light_bins),binsize, max(params(exp_num(clean_units_trains(i))).light_dur), lightconds{1}(ii+2), params(exp_num(clean_units_trains(i))).all_light(vis_inds)) ;
%                [resp_t_bl{ii}(i,:), resp_dur_bl{ii}(i,:), ppr_bl{ii}(i,:)]  = ppanalysis_v2(psth_bl(ii+2,light_bins),binsize, max(params(exp_num(clean_units_trains(i))).light_dur), lightconds{1}(ii+2), params(exp_num(clean_units_trains(i))).all_light(~vis_inds)) ;
%           
        [resp_t_v{ii}(i,:), resp_dur_v{ii}(i,:), ppr_v{ii}(i,:), ppr2_v{ii}(i,:)]  = ppanalysis(params(exp_num(clean_units_trains(i))).prestim, params(exp_num(clean_units_trains(i))).stimtime, round(mean(params(exp_num(clean_units_trains(i))).av_light_start))-params(exp_num(clean_units_trains(i))).prestim, max(params(exp_num(clean_units_trains(i))).light_dur), lightconds{1}(ii+2), unitinfo(clean_units_trains(i)).rast(vis_inds,:), params(exp_num(clean_units_trains(i))).trial_type(vis_inds,strcmpi(params(exp_num(clean_units_trains(i))).IVs,'light_bit'))) ;
           [resp_t_bl{ii}(i,:), resp_dur_bl{ii}(i,:), ppr_bl{ii}(i,:), ppr2_bl{ii}(i,:)]  = ppanalysis(params(exp_num(clean_units_trains(i))).prestim, params(exp_num(clean_units_trains(i))).stimtime, round(mean(params(exp_num(clean_units_trains(i))).av_light_start))-params(exp_num(clean_units_trains(i))).prestim, max(params(exp_num(clean_units_trains(i))).light_dur), lightconds{1}(ii+2), unitinfo(clean_units_trains(i)).rast(~vis_inds,:), params(exp_num(clean_units_trains(i))).trial_type(~vis_inds,strcmpi(params(exp_num(clean_units_trains(i))).IVs,'light_bit'))) ;
           
        end

       nonans_bl{ii} = ~isnan(ppr_bl{ii}(:,1));  
        nonans_v{ii} = ~isnan(ppr_v{ii}(:,1));
%         hz_actV{ii} = F_Z>3;        % include those greater than 3stds 
%         hz_actBl{ii} = F_Zbl>3;
        
        % check if paired pulse ratios are significantly different from 1
        for x=1:size(ppr_v{ii},2)
            [p_v{ii}(x)] = signtest(ppr_v{ii}(:,x),1);      % check if paired pulse ratios are significantly different from 1
            [p_bl{ii}(x)] = signtest(ppr_bl{ii}(:,x),1);    % this and above changed 7/21/19 from signrank to signtest (makes less assumptions about data distribution - MAK)
        end
        
        figure;
        h=histogram(ppr_bl{ii}(:,1),'binwidth',.25,'facecolor','b');
        ylim([0 max(histcounts(ppr_bl{ii}(:,1),'binwidth',.25))+1])
        xlim([max(-5,-1*ceil(max(ppr_bl{ii}(:,1)))+2) ceil(max(ppr_bl{ii}(:,1)))])
        yax = get(gca,'ylim');
        line([1 1],yax,'linestyle','--','color','k','linewidth',2)
        ylabel('Number of units','fontsize',24)
        xlabel('Spike count ratio: 2nd/1st pulse','fontsize',24)
        set(gca,'fontsize',18,'linewidth',2)
        set(h,'facecolor',area_color)
        print(gcf, '-dpng',sprintf('Paired pulse histogram %dHz (blanks)',lightconds{1}(ii+2)))
        print(gcf, '-painters','-depsc',sprintf('Paired pulse histogram %dHz (blanks)',lightconds{1}(ii+2)))

        figure;
        h=histogram(ppr_v{ii}(:,1),'binwidth',.25,'facecolor','b');
        ylim([0 max(histcounts(ppr_v{ii}(:,1),'binwidth',.25))+1])
        xlim([max(-5,-1*ceil(max(ppr_v{ii}(:,1)))+2) ceil(max(ppr_v{ii}(:,1)))])
        yax = get(gca,'ylim');
        line([1 1],yax,'linestyle','--','color','k','linewidth',2)
        ylabel('Number of units','fontsize',24)
        xlabel('Spike count ratio: 2nd/1st pulse','fontsize',24)
        set(gca,'fontsize',18,'linewidth',2)
        set(h,'facecolor',area_color)
        print(gcf, '-dpng',sprintf('Paired pulse histogram %dHz (visual)',lightconds{1}(ii+2)))
        print(gcf, '-painters','-depsc',sprintf('Paired pulse histogram %dHz (visual)',lightconds{1}(ii+2)))
       

        figure;
        b_bl=bar([1:lightconds{1}(ii+2)-1],nanmedian(ppr_bl{ii}));
        ppiqr_bl = iqr(ppr_bl{ii}(nonans_bl{ii},:))/2;
        set(b_bl,'XData',[2:lightconds{1}(ii+2)],'facecolor',area_color)
        xax = get(gca,'xlim');
        line(xax,[1 1],'linestyle','--','color','k','linewidth',2)
        hold on;
        for x=1:lightconds{1}(ii+2)-1
    %         plot((x+1)*ones(1,sum(~isnan(ppr_bl(:,x)))),ppr_bl(~isnan(ppr_bl(:,x)),x),'k.')
            line([b_bl.XData(x) b_bl.XData(x)],[b_bl.YData(x) b_bl.YData(x)+ppiqr_bl(x)],'color','k')
        end
        set(gca,'fontsize',18,'linewidth',2)
        ylabel('Spike count ratio (pulse x/pulse 1)','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        yax_bl = get(gca,'ylim');
        plot(b_bl.XData(p_bl{ii}<.05),(yax_bl(2)-.1*yax_bl(2))*ones(1,sum(p_bl{ii}<.05)),'k*')   % black asterisk above all sig vals in ppr figs
        print(gcf, '-dpng',sprintf('Spike count ratio %dHz (blanks)',lightconds{1}(ii+2)))
        print(gcf, '-painters','-depsc',sprintf('Spike count ratio %dHz (blanks)',lightconds{1}(ii+2)))

        figure;
        b_v=bar([1:lightconds{1}(ii+2)-1],nanmedian(ppr_v{ii}));
        ppiqr_v = iqr(ppr_v{ii}(nonans_v{ii},:))/2;
        set(b_v,'XData',[2:lightconds{1}(ii+2)],'facecolor',area_color)
        xax = get(gca,'xlim');
        line(xax,[1 1],'linestyle','--','color','k','linewidth',2)
        hold on;
        for x=1:lightconds{1}(ii+2)-1
    %         plot((x+1)*ones(1,sum(~isnan(ppr_v(:,x)))),ppr_v(~isnan(ppr_v(:,x)),x),'k.')
            line([b_v.XData(x) b_v.XData(x)],[b_v.YData(x) b_v.YData(x)+ppiqr_v(x)],'color','k')
        end
        set(gca,'fontsize',18,'linewidth',2)
        ylabel('Spike count ratio (pulse x/pulse 1)','fontsize',16)
        xlabel('Pulse number','fontsize',16)
        yax_v = get(gca,'ylim');
        plot(b_v.XData(p_v{ii}<.05),(yax_v(2)-.1*yax_v(2))*ones(1,sum(p_v{ii}<.05)),'k*')   % black asterisk above all sig vals in ppr figs
        print(gcf, '-dpng',sprintf('Spike count ratio %dHz (visual)',lightconds{1}(ii+2)))
        print(gcf, '-painters','-depsc',sprintf('Spike count ratio %dHz (visual)',lightconds{1}(ii+2)))

%         figure;
%         b=bar([1:lightconds{1}(ii+2)],nanmedian(resp_dur_bl{ii}));
%         iqr_bl= iqr(resp_dur_bl{ii}(nonans_bl{ii},:))/2;
%         set(b,'XData',[1:lightconds{1}(ii+2)],'facecolor',area_color)
%         hold on;
%         for x=1:lightconds{1}(ii+2)
%             line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+iqr_bl(x)],'color','k')
%         end
%         set(gca,'fontsize',18,'linewidth',2)
%         ylabel('Median # of significant bins','fontsize',16)
%         xlabel('Pulse number','fontsize',16)
%         print(gcf, '-dpng',sprintf('Pulse response duration %dHz (blanks)',lightconds{1}(ii+2)))
%         print(gcf, '-painters','-depsc',sprintf('Pulse response duration %dHz (blanks)',lightconds{1}(ii+2)))
% 
%         figure;
%         b=bar([1:lightconds{1}(ii+2)],nanmedian(resp_dur_v{ii}));
%         iqr_v = iqr(resp_dur_v{ii}(nonans_v{ii},:))/2;
%         set(b,'XData',[1:lightconds{1}(ii+2)],'facecolor',area_color)
%         hold on;
%         for x=1:lightconds{1}(ii+2)
%             line([b.XData(x) b.XData(x)],[b.YData(x) b.YData(x)+iqr_v(x)],'color','k')
%         end
%         set(gca,'fontsize',18,'linewidth',2)
%         ylabel('Median # of significant bins','fontsize',16)
%         xlabel('Pulse number','fontsize',16)
%         print(gcf, '-dpng',sprintf('Pulse response duration %dHz (visual)',lightconds{1}(ii+2)))
%         print(gcf,'-painters','-depsc',sprintf('Pulse response duration %dHz (visual)',lightconds{1}(ii+2)))

    % save stats
        if ii==1
            fileID = fopen('trains_PPRstats.txt','w');
        else
            fileID = fopen('trains_PPRstats.txt','a');
        end
        fprintf(fileID,'Signed-rank test of paired-pulse ratios, %d Hz, visual trials: %s \r\n',lightconds{1}(ii+2),num2str(round(p_v{ii},4)));
        fprintf(fileID,'Signed-rank test of paired-pulse ratios, %d Hz, blank trials: %s \r\n',lightconds{1}(ii+2),num2str(round(p_bl{ii},4)));
        fclose(fileID);
    end
    
end

%% F1/F0 response analysis

binsize = .025;
% pref_psth = nan(length(0:binsize:params(exp_num(clean_units(i))).stimtime-binsize),length(conds)+1,length(clean_units)); 
psthV(:,:,:) = reshape(cell2mat(arrayfun(@(x) x.psthVisual([conds end],:), FRs(clean_units),'uniformoutput',0)),length(conds)+1,size(FRs(1).psthVisual,2),length(clean_units));
% n timebins x num conds x  num units
Fratio = nan(length(clean_units),length(conds)+1);
for i = 1:length(clean_units)
    all_light = params(exp_num(clean_units(i))).all_light;
    vis_start = params(exp_num(clean_units(i))).prestim*1000;       % in ms
    vis_end = params(exp_num(clean_units(i))).poststim*1000;      % in ms
%     % only for preferred visual stim trials
%     which_trials = ismember(1:size(unitinfo(clean_units(i)).rast,1),prefori_trials{i});
%     [~,tmp_psth] = make_psth_v2(binsize,0:binsize:(size(unitinfo(clean_units(i)).rast,2)-vis_start)/1000,which_trials,unitinfo(clean_units(i)).rast(:,vis_start+1:end-vis_end),all_light);
%     pref_psth(:,:,i) = tmp_psth'-FRs(clean_units(i)).psthBlank(:,21:end)';       % subtract baseline!! (baseline from each light cond in order to look at whether light impacts F1/F0 independent of any gain change. does this make sense??)
%     Fratio(i,:) = calc_F1F0(pref_psth(:,:,i),binsize,2);        % **currently hardcoded for 2Hz tfreq - will need to change!!    
%
    % for all visual stim trials:
    % but how should I handle psths with negative values?? (i.e. units
    % suppressed by vis stim)
    Fratio(i,:) = calc_F1F0(psthV(:,(vis_start/1000)/binsize+1:end,i)',binsize,2);        % **currently hardcoded for 2Hz tfreq - will need to change!!
end

% for trains experiments, also find power at 10hz frequency in no light and
% 10hz conds
Fratio_10hz = nan(length(clean_units),2);
if contains(exp_type,'trains')
    for i = 1:length(clean_units)
        for ii=[1,3]    % manually set for first (no light) and third (10Hz) light conds
            light_start = round(params(exp_num(clean_units(i))).av_light_start(lightcond));
            light_dur = params(exp_num(clean_units(i))).light_dur(lightcond+1);
            Fratio_10hz(i,round(ii/2)) = calc_F1F0(psthV(ii,light_start/binsize+1:(light_start+light_dur)/binsize,i)',binsize,10);
        end
    end
end

%% NEW (10/28/19) - burst vs tonic firing stuff
burstrate = zeros(length(clean_units),1);
burstrate_light = zeros(length(clean_units),size(all_lightconds,1)-1);
burstnum = burstrate;
burstnumLight = burstrate_light;
burstresprate = burstrate;
burstresprate_light = burstrate_light;
for i = 1:length(clean_units)
%     timeinds = round(1000*(params(exp_num(clean_units(i))).av_light_start(1)))+1-100:round(1000*(params(exp_num(clean_units(i))).av_light_start(1)))+params(exp_num(clean_units(i))).lighttime*1000+4;      % CHANGED 3/2/20 - -100 AND +4 b/c otherwise first and last burst spikes in trials would be excluded
    timeinds = round(1000*(params(exp_num(clean_units(i))).av_light_start(1)))+1:round(1000*(params(exp_num(clean_units(i))).av_light_start(1)))+params(exp_num(clean_units(i))).lighttime*1000;      
    visnolight_trials = intersect(stat_trials{exp_num(clean_units(i))},intersect(nolight_trials{exp_num(clean_units(i))},vis_trials{exp_num(clean_units(i))}));  % visual and STATIONARY trials
%     visRast = unitinfo(clean_units(i)).rast(visnolight_trials,timeinds);
    [burstrate(i),burstnum(i),burstresprate(i),avtonicrate(i),avburstrate(i)] = calc_bursting(unitinfo(clean_units(i)).rast(visnolight_trials,:),timeinds,0);
%     timeinds_lighttime = [101:length(timeinds)-4];  % bounds of timeinds indicating when light was ACTUALLY on
%     [rows, cols] = find(visRast);
%     [trialn,ord] = sort(rows);   %sort by trials
%     isis = [cols(ord(1)); diff(cols(ord))]; % interspike intervals by trial (use ms index of first spike for it's isi, since we don't actually know when prior spike occurred
%     % an isi for every spike - e.g., isis(5) indicates how many ms PRECEDED the 5th spike
% %     trialn = trialn(isis>=0);
% %     isis = isis(isis>=0);   % isis>0 are transitions between trials
%     isis(isis<=0) = cols(ord(isis<=0)); % again, use ms index of first spike in a trial for it's isi
%     burst_mid = find(isis(1:end-1)<=4); % spikes in middle and end of bursts (preceeded by <=4ms ISIs)
%     burst_st = find(isis(1:end-1)>=100 & isis(2:end)<=4);   % spikes at start of bursts (followed by <=4ms ISI but preceded by >=100ms)
%     bursts = union(burst_mid,burst_st); % all burst indices
%     while sum(~ismember(burst_mid-1,bursts))    % fixed 3/2/20 - previously only did this once, which still left in very small number incorrect burst spikes
%         premid = burst_mid(~ismember(burst_mid-1,bursts));  % burst_mid spikes whose preceding spikes were NOT burst spikes
%         bursts(ismember(bursts,premid)) = [];
%         burst_mid(ismember(burst_mid,premid)) = [];
%     end
%     tonics = trialn(~ismember(1:length(isis),bursts));
%     for n=1:length(visnolight_trials)
%         tonic_tri(n) = sum(tonics==n);
%         burst_tri(n) = sum(trialn(bursts)==n);
%     end
%     tonic_rate(i) = mean(tonic_tri);
%     burst_rate(i) = mean(burst_tri);
    
%     % calculate number of spikes in bursts
%     burst_st(burst_st>length(isis)-4) = [];  % if start of burst is too close to end, won't be able to count # spikes in that burst
%     burstmat = cumsum([burst_st ones(length(burst_st),25)],2);
%     burstmat(burstmat>length(isis)) = length(isis);
%     burstlength = nan(1,length(burst_st));
%     for b=1:size(burstmat,1)    % for each start of a burst
%         if ~isempty(find(isis(burstmat(b,2:end))>4,1,'first'))
%             burstlength(b) = find(isis(burstmat(b,2:end))>4,1,'first');   % how many spikes before >4ms ISI
%         end
%     end
%     burstnum(i) = nanmean(burstlength);    % save unit's mean burstlength
%     
% %     burstrate(i) = bursts/sum(visRast(:)); % how many of all spikes in lightstim period were part of bursts
% %         burstrate(i) = length(bursts)/(length(isis)-1); % I am not sure why this was -1 before?? (changed 3/2/20 MAK)
%     burstrate(i) = length(bursts)/(sum(sum(visRast(:,timeinds_lighttime)))); % how many of all spikes in lightstim period were part of bursts
% %     burstresprate(i) = length(burst_st)/(length(isis)-1-length(burst_mid));   % how many of all RESPONSES (as opposed to individual spikes) were the starts of bursts
%     burstresprate(i) = length(burst_st)/(sum(sum(visRast(:,timeinds_lighttime)))-length(burst_mid));    

    for ii = 1:size(all_lightconds,1)-1
        vislight_trials = intersect(stat_trials{exp_num(clean_units(i))},intersect(light_trials{exp_num(clean_units(i))}{ii},vis_trials{exp_num(clean_units(i))})); % visual and STATIONARY trials
        vislightRast = unitinfo(clean_units(i)).rast(vislight_trials,timeinds);
        [burstrate_light(i,ii),burstnumLight(i,ii),burstresprate_light(i,ii),avtonicrate_light(i,ii),avburstrate_light(i,ii)] = calc_bursting(unitinfo(clean_units(i)).rast(vislight_trials,:),timeinds,0);
        
%         [rowsLight, colsLight] = find(vislightRast(:,timeinds_lighttime));
%         [trialnLight,ordLight] = sort(rowsLight);
%         isisLight = diff(colsLight(ordLight)); 
%         trialnLight = trialnLight(isisLight>0);
%         isisLight = isisLight(isisLight>0);
%         burst_mid_light = find(isisLight(1:end-1)<=4);
%         burst_st_light = find(isisLight(1:end-1)>=100 & isisLight(2:end)<=4);
%         bursts_light = union(burst_mid_light,burst_st_light); % all burst indices
%         while sum(~ismember(burst_mid_light-1,bursts_light))
%             premid_light = burst_mid_light(~ismember(burst_mid_light-1,bursts_light));  % burst_mid spikes whose preceding spikes were NOT burst spikes
%             bursts_light(ismember(bursts_light,premid_light)) = [];
%             burst_mid_light(ismember(burst_mid_light,premid_light)) = [];
%         end
%         tonics_light = trialnLight(~ismember(1:length(isisLight),bursts_light));
%         tonic_tri_light = zeros(1,length(vislight_trials));
%         burst_tri_light = tonic_tri_light;
%         for n=1:length(vislight_trials)
%             tonic_tri_light(n) = sum(tonics_light==n);
%             burst_tri_light(n) = sum(trialnLight(bursts_light)==n);
%         end
%         tonic_rate_light(i,ii) = mean(tonic_tri_light);
%         burst_rate_light(i,ii) = mean(burst_tri_light); % mean across trials
%     
%         % calculate number of spikes in bursts
%         burst_st_light(burst_st_light>length(isisLight)-4) = [];  % if start of burst is too close to end, won't be able to count # spikes in that burst
%         burstmatLight = cumsum([burst_st_light ones(length(burst_st_light),25)],2);
%         burstmatLight(burstmatLight>length(isisLight)) = length(isisLight);
%         burstlengthLight = nan(1,length(burst_st_light));
%         for b=1:size(burstmatLight,1)    % for each start of a burst
%             if ~isempty(find(isisLight(burstmatLight(b,2:end))>4,1,'first'))
%                 burstlengthLight(b) = find(isisLight(burstmatLight(b,2:end))>4,1,'first');   % how many spikes before >4ms ISI
%             end
%         end
%         burstnumLight(i,ii) = nanmean(burstlengthLight);    % save unit's mean burstlength
%         
% %         burstrate_light(i,ii) = length(bursts_light)/(length(isisLight)-1);% I am not sure why this was -1 before?? (changed 3/2/20 MAK)
%         burstrate_light(i,ii) = length(bursts_light)/(sum(sum(vislightRast(:,timeinds_lighttime)))); % how many of all spikes in lightstim period were part of bursts (collapsed across all trials!)
% %         burstresprate_light(i,ii) = length(burst_st_light)/(length(isisLight)-1-length(burst_mid_light));   % how many of all RESPONSES (as opposed to individual spikes) were the starts of bursts
%         burstresprate_light(i,ii) = length(burst_st_light)/(sum(sum(vislightRast(:,timeinds_lighttime)))-length(burst_mid_light));   % how many of all RESPONSES (as opposed to individual spikes) were the starts of bursts
    end
end

%% get preferred stimulus FR for each condition (but preferred stimulus defined in no light condition)
durs = [params(:).light_dur];
if contains(exp_type,'trains','ignorecase',1)
    dur = max(durs(durs>0));
elseif contains(exp_type,'step','ignorecase',1)
    dur = min(durs(durs>0));
end
light_times = round([max([params(:).av_light_start]) max([params(:).av_light_start])+dur].*1000);   % start with latest light start time, end after minimum light duration that isn't 0
FRpref = zeros(length(clean_units),length(lightconds{1}));
for i = 1:length(clean_units)
    for ii = 1:length(lightconds{exp_num(clean_units(i))})
        if ii == 1
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(clean_units(i))}),stat_trials{exp_num(clean_units(i))}),light_times(1)+1:light_times(2)),2))/(diff(light_times)/1000); % edited to include only stationary trials (MAK - 5/8/19)
        else
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(clean_units(i))}{ii-1}),stat_trials{exp_num(clean_units(i))}),light_times(1)+1:light_times(2)),2))/(diff(light_times)/1000);  % edited to include only stationary trials (MAK - 5/8/19)
        end
    end
end

%%
cd(fig_dir)
% get different cell types
visual_cells = find((vis_sig < .025)|(vis_sig_ons < .025));  % CHANGED 7/21/19 (.025 b/c bonferroni correction for 2 tests)
nonvisual_cells = find(~ismember(1:length(vis_sig),visual_cells));
% light_cells= find(min(light_sig,[],2)<.05);   % find cells with significant effect in any light condition
if contains(exp_type,'trains','ignorecase',1) || contains(exp_type,'inLP','ignorecase',1)
    light_cells = find(light_sig(:,lightcond)<.05);
    supp_cells = find((light_sig(:,lightcond)<.05) & sign(lightmod(:,lightcond))==-1);
    enh_cells = find((light_sig(:,lightcond)<.05) & sign(lightmod(:,lightcond))==1);
else
    light_cells = find(sum(light_sig<.05,2)>1);     % significant in 2 or more light conditions
    % supp_cells = find((light_sig(:,1)<.05)&(light_sig(:,2)<.05)&(light_sig(:,3)<.05)&(sum(sign(lightmod),2)<0));
    supp_cells = find((sum(light_sig<.05,2)>1)&(sum(sign(lightmod),2)<0)); % significant in 2 or more light conditions, and direction of lightmod was - in 2+ conditions
    % significantly light suppressed in 2/3 conditions
    % enh_cells = find((light_sig(:,1)<.05)&(light_sig(:,2)<.05)&(light_sig(:,3)<.05)&(sum(sign(lightmod),2)>0));
    enh_cells = find((sum(light_sig<.05,2)>1)&(sum(sign(lightmod),2)>0));% significant in 2 or more light conditions, and direction of lightmod was + in 2+ conditions
end
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


% significantly light activated in 2/3 conditions
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
% lightsig_all_high = signrank(FRev(:,1),FRev(:,lightcond+1));
% lightsig_all_high_dir = sign(nanmedian(FRev(:,lightcond+1))-nanmedian(FRev(:,1)));
% lightsig_bl_low = signrank(FRbl(:,1),FRbl(:,2));
% lightsig_bl_low_dir = sign(nanmedian(FRbl(:,2))-nanmedian(FRbl(:,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_bl_high = signrank(FRbl(:,1),FRbl(:,lightcond+1));
% lightsig_bl_high_dir = sign(nanmedian(FRbl(:,lightcond+1))-nanmedian(FRbl(:,1)));
% lightsig_vis_low = signrank(FRev(visual_cells,1),FRev(visual_cells,2));
% lightsig_vis_low_dir = sign(nanmedian(FRev(visual_cells,2))-nanmedian(FRev(visual_cells,1)));
% lightsig_vis_high = signrank(FRev(visual_cells,1),FRev(visual_cells,lightcond+1));
% lightsig_vis_high_dir = sign(nanmedian(FRev(visual_cells,lightcond+1))-nanmedian(FRev(visual_cells,1)));
% lightsig_blvis_low = signrank(FRbl(visual_cells,1),FRbl(visual_cells,2));
% lightsig_blvis_low_dir = sign(nanmedian(FRbl(visual_cells,2))-nanmedian(FRbl(visual_cells,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_blvis_high = signrank(FRbl(visual_cells,1),FRbl(visual_cells,lightcond+1));
% lightsig_blvis_high_dir = sign(nanmedian(FRbl(visual_cells,lightcond+1))-nanmedian(FRbl(visual_cells,1)));
% lightsig_nonvis_low = signrank(FRev(nonvisual_cells,1),FRev(nonvisual_cells,2));
% lightsig_nonvis_low_dir = sign(nanmedian(FRev(nonvisual_cells,2))-nanmedian(FRev(nonvisual_cells,1)));
% lightsig_nonvis_high = signrank(FRev(nonvisual_cells,1),FRev(nonvisual_cells,lightcond+1));
% lightsig_nonvis_high_dir = sign(nanmedian(FRev(nonvisual_cells,lightcond+1))-nanmedian(FRev(nonvisual_cells,1)));
% lightsig_blnonvis_low = signrank(FRbl(nonvisual_cells,1),FRbl(nonvisual_cells,2));
% lightsig_blnonvis_low_dir = sign(nanmedian(FRbl(nonvisual_cells,2))-nanmedian(FRbl(nonvisual_cells,1))); % 1 if light increased FR; -1 if it decreased FR
% lightsig_blnonvis_high = signrank(FRbl(nonvisual_cells,1),FRbl(nonvisual_cells,lightcond+1));
% lightsig_blnonvis_high_dir = sign(nanmedian(FRbl(nonvisual_cells,lightcond+1))-nanmedian(FRbl(nonvisual_cells,1)));
% lightsig_vispref_low = signrank(FRpref(:,1),FRpref(:,2));
% lightsig_vispref_low_dir = sign(nanmedian(FRpref(:,2))-nanmedian(FRpref(:,1)));
% lightsig_vispref_high = signrank(FRpref(:,1),FRpref(:,lightcond+1));
% lightsig_vispref_high_dir = sign(nanmedian(FRpref(:,lightcond+1))-nanmedian(FRpref(:,1)));
% lightsig_visFRdelt_low = signrank(FRev(:,1)-FRbl(:,1),FRev(:,2)-FRbl(:,2));
% lightsig_visFRdelt_low_dir = sign(nanmedian(FRev(:,2)-FRbl(:,2))-nanmedian(FRev(:,1)-FRbl(:,1)));
% lightsig_visFRdelt_high = signrank(FRev(:,1)-FRbl(:,1),FRev(:,lightcond+1)-FRbl(:,lightcond+1));
% lightsig_visFRdelt_high_dir = sign(nanmedian(FRev(:,lightcond+1)-FRbl(:,lightcond+1))-nanmedian(FRev(:,1)-FRbl(:,1)));
% lightsig_prefFRdelt_low_pref = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,2)-FRbl(:,2));
% lightsig_prefFRdelt_low_pref_dir = sign(nanmedian(FRpref(:,2)-FRbl(:,2))-nanmedian(FRpref(:,1)-FRbl(:,1)));
% lightsig_prefFRdelt_high = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,lightcond+1)-FRbl(:,lightcond+1));
% lightsig_prefFRdelt_high_dir = sign(nanmedian(FRpref(:,lightcond+1)-FRbl(:,lightcond+1))-nanmedian(FRpref(:,1)-FRbl(:,1)));
% 
% lightsig_vals = [lightsig_all_low lightsig_all_low_dir; lightsig_all_high lightsig_all_high_dir; lightsig_bl_low lightsig_bl_low_dir; lightsig_bl_high lightsig_bl_high_dir;...
%     lightsig_vis_low lightsig_vis_low_dir; lightsig_vis_high lightsig_vis_high_dir; lightsig_blvis_low lightsig_blvis_low_dir; lightsig_blvis_high lightsig_blvis_high_dir;...
%     lightsig_nonvis_low lightsig_nonvis_low_dir; lightsig_nonvis_high lightsig_nonvis_high_dir; lightsig_blnonvis_low lightsig_blnonvis_low_dir; lightsig_blnonvis_high lightsig_blnonvis_high_dir;...
%     lightsig_vispref_low lightsig_vispref_low_dir; lightsig_vispref_high lightsig_vispref_high_dir; lightsig_visFRdelt_low lightsig_visFRdelt_low_dir; lightsig_visFRdelt_high lightsig_visFRdelt_high_dir;...
%     lightsig_prefFRdelt_low_pref lightsig_prefFRdelt_low_pref_dir; lightsig_prefFRdelt_high lightsig_prefFRdelt_high_dir];
% 
% osiCVsig_tuned_low = signrank(OSI_CV(tuned_cells,1),OSI_CV(tuned_cells,2));
% osiCVsig_tuned_low_dir = sign(nanmedian(OSI_CV(tuned_cells,2))-nanmedian(OSI_CV(tuned_cells,1)));
% osiCVsig_tuned_high = signrank(OSI_CV(tuned_cells,1),OSI_CV(tuned_cells,lightcond+1));
% osiCVsig_tuned_high_dir = sign(nanmedian(OSI_CV(tuned_cells,lightcond+1))-nanmedian(OSI_CV(tuned_cells,1)));
% osisig_tuned_low = signrank(OSI(tuned_cells,1),OSI(tuned_cells,2));
% osisig_tuned_low_dir = sign(nanmedian(OSI(tuned_cells,2))-nanmedian(OSI(tuned_cells,1)));
% osisig_tuned_high = signrank(OSI(tuned_cells,1),OSI(tuned_cells,lightcond+1));
% osisig_tuned_high_dir = sign(nanmedian(OSI(tuned_cells,lightcond+1))-nanmedian(OSI(tuned_cells,1)));
% dsiCVsig_tuned_low = signrank(DSI_CV(tuned_cells,1),DSI_CV(tuned_cells,2));
% dsiCVsig_tuned_low_dir = sign(nanmedian(DSI_CV(tuned_cells,2))-nanmedian(DSI_CV(tuned_cells,1)));
% dsiCVsig_tuned_high = signrank(DSI_CV(tuned_cells,1),DSI_CV(tuned_cells,lightcond+1));
% dsiCVsig_tuned_high_dir = sign(nanmedian(DSI_CV(tuned_cells,lightcond+1))-nanmedian(DSI_CV(tuned_cells,1)));
% dsisig_tuned_low = signrank(DSI(tuned_cells,1),DSI(tuned_cells,2));
% dsisig_tuned_low_dir = sign(nanmedian(DSI(tuned_cells,2))-nanmedian(DSI(tuned_cells,1)));
% dsisig_tuned_high = signrank(DSI(tuned_cells,1),DSI(tuned_cells,lightcond+1));
% dsisig_tuned_high_dir = sign(nanmedian(DSI(tuned_cells,lightcond+1))-nanmedian(DSI(tuned_cells,1)));
% osiCVsig_vis_low = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,2));
% osiCVsig_vis_low_dir = sign(nanmedian(OSI_CV(visual_cells,2))-nanmedian(OSI_CV(visual_cells,1)));
% osiCVsig_vis_high = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,lightcond+1));
% osiCVsig_vis_high_dir = sign(nanmedian(OSI_CV(visual_cells,lightcond+1))-nanmedian(OSI_CV(visual_cells,1)));
% osisig_vis_low = signrank(OSI(visual_cells,1),OSI(visual_cells,2));
% osisig_vis_low_dir = sign(nanmedian(OSI(visual_cells,2))-nanmedian(OSI(visual_cells,1)));
% osisig_vis_high = signrank(OSI(visual_cells,1),OSI(visual_cells,lightcond+1));
% osisig_vis_high_dir = sign(nanmedian(OSI(visual_cells,lightcond+1))-nanmedian(OSI(visual_cells,1)));
% dsiCVsig_vis_low = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,2));
% dsiCVsig_vis_low_dir = sign(nanmedian(DSI_CV(visual_cells,2))-nanmedian(DSI_CV(visual_cells,1)));
% dsiCVsig_vis_high = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,lightcond+1));
% dsiCVsig_vis_high_dir = sign(nanmedian(DSI_CV(visual_cells,lightcond+1))-nanmedian(DSI_CV(visual_cells,1)));
% dsisig_vis_low = signrank(DSI(visual_cells,1),DSI(visual_cells,2));
% dsisig_vis_low_dir = sign(nanmedian(DSI(visual_cells,2))-nanmedian(DSI(visual_cells,1)));
% dsisig_vis_high = signrank(DSI(visual_cells,1),DSI(visual_cells,lightcond+1));
% dsisig_vis_high_dir = sign(nanmedian(DSI(visual_cells,lightcond+1))-nanmedian(DSI(visual_cells,1)));
% tuningsig_vals = [osiCVsig_tuned_low osiCVsig_tuned_low_dir; osiCVsig_tuned_high osiCVsig_tuned_high_dir; osisig_tuned_low osisig_tuned_low_dir;...
%     osisig_tuned_high osisig_tuned_high_dir; dsiCVsig_tuned_low dsiCVsig_tuned_low_dir; dsiCVsig_tuned_high dsiCVsig_tuned_high_dir;...
%     dsisig_tuned_low dsisig_tuned_low_dir; dsisig_tuned_high dsisig_tuned_high_dir; osiCVsig_vis_low osiCVsig_vis_low_dir;...
%     osiCVsig_vis_high osiCVsig_vis_high_dir; osisig_vis_low osisig_vis_low_dir; osisig_vis_high osisig_vis_high_dir;...
%     dsiCVsig_vis_low dsiCVsig_vis_low_dir; dsiCVsig_vis_high dsiCVsig_vis_high_dir; dsisig_vis_low dsisig_vis_low_dir; dsisig_vis_high dsisig_vis_high_dir];


%% %% test significance of overall light modulation and OSI/DSI change using Wilcoxen signed-rank
% 7/27/19 - changed "_dir"s to reflect sign of median of differences,
% rather than sign of difference of medians (MAK)
for i=1:length(conds)
    lightsig_all(i) = signrank(FRev(:,1),FRev(:,i+1));
    lightsig_all_dir(i) = sign(nanmedian(FRev(:,i+1)-(FRev(:,1)))); % 1 if light increased FR; -1 if it decreased FR
    lightsig_bl(i) = signrank(FRbl(:,1),FRbl(:,i+1));
    lightsig_bl_dir(i) = sign(nanmedian(FRbl(:,i+1)-(FRbl(:,1)))); % 1 if light increased FR; -1 if it decreased FR
    lightsig_vis(i) = signrank(FRev(visual_cells,1),FRev(visual_cells,i+1));
    lightsig_vis_dir(i) = sign(nanmedian(FRev(visual_cells,i+1)-(FRev(visual_cells,1))));
    lightsig_blvis(i) = signrank(FRbl(visual_cells,1),FRbl(visual_cells,i+1));
    lightsig_blvis_dir(i) = sign(nanmedian(FRbl(visual_cells,i+1)-(FRbl(visual_cells,1)))); % 1 if light increased FR; -1 if it decreased FR
    lightsig_nonvis(i) = signrank(FRev(nonvisual_cells,1),FRev(nonvisual_cells,i+1));
    lightsig_nonvis_dir(i) = sign(nanmedian(FRev(nonvisual_cells,i+1)-(FRev(nonvisual_cells,1))));
    lightsig_blnonvis(i) = signrank(FRbl(nonvisual_cells,1),FRbl(nonvisual_cells,i+1));
    lightsig_blnonvis_dir(i) = sign(nanmedian(FRbl(nonvisual_cells,i+1)-(FRbl(nonvisual_cells,1)))); % 1 if light increased FR; -1 if it decreased FR
    lightsig_vispref(i) = signrank(FRpref(:,1),FRpref(:,i+1));
    lightsig_vispref_dir(i) = sign(nanmedian(FRpref(:,i+1)-(FRpref(:,1))));
    lightsig_visFRdelt(i) = signrank(FRev(:,1)-FRbl(:,1),FRev(:,i+1)-FRbl(:,i+1));
    lightsig_visFRdelt_dir(i) = sign(nanmedian(FRev(:,i+1)-FRbl(:,i+1)-(FRev(:,1)-FRbl(:,1))));
    lightsig_prefFRdelt(i) = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,i+1)-FRbl(:,i+1));
    lightsig_prefFRdelt_dir(i) = sign(nanmedian(FRpref(:,i+1)-FRbl(:,i+1)-(FRpref(:,1)-FRbl(:,1))));

    osiCVsig_tuned(i) = signrank(OSI_CV(tuned_cells,1),OSI_CV(tuned_cells,i+1));
    osiCVsig_tuned_dir(i) = sign(nanmedian(OSI_CV(tuned_cells,i+1)-(OSI_CV(tuned_cells,1))));
    osisig_tuned(i) = signrank(OSI(tuned_cells,1),OSI(tuned_cells,i+1));
    osisig_tuned_dir(i) = sign(nanmedian(OSI(tuned_cells,i+1)-(OSI(tuned_cells,1))));
    dsiCVsig_tuned(i) = signrank(DSI_CV(tuned_cells,1),DSI_CV(tuned_cells,i+1));
    dsiCVsig_tuned_dir(i) = sign(nanmedian(DSI_CV(tuned_cells,i+1)-(DSI_CV(tuned_cells,1))));
    dsisig_tuned(i) = signrank(DSI(tuned_cells,1),DSI(tuned_cells,i+1));
    dsisig_tuned_dir(i) = sign(nanmedian(DSI(tuned_cells,i+1)-(DSI(tuned_cells,1))));
    osiCVsig_vis(i) = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,i+1));
    osiCVsig_vis_dir(i) = sign(nanmedian(OSI_CV(visual_cells,i+1)-(OSI_CV(visual_cells,1))));
    osisig_vis(i) = signrank(OSI(visual_cells,1),OSI(visual_cells,i+1));
    osisig_vis_dir(i) = sign(nanmedian(OSI(visual_cells,i+1)-(OSI(visual_cells,1))));
    dsiCVsig_vis(i) = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,i+1));
    dsiCVsig_vis_dir(i) = sign(nanmedian(DSI_CV(visual_cells,i+1)-(DSI_CV(visual_cells,1))));
    dsisig_vis(i) = signrank(DSI(visual_cells,1),DSI(visual_cells,i+1));
    dsisig_vis_dir(i) = sign(nanmedian(DSI(visual_cells,i+1)-(DSI(visual_cells,1))));
    
    burstrate_sig(i) = signrank(burstrate,burstrate_light(:,i));
    burstrate_dir(i) = sign(nanmedian(burstrate_light(:,i)-burstrate));
    burstresprate_sig(i) = signrank(burstresprate,burstresprate_light(:,i));
    burstresprate_dir(i) = sign(nanmedian(burstresprate_light(:,i)-burstresprate));
end
%%
% vis_type = nan(1,length(clean_units));
% vis_type(trans_cells) = 1;
% vis_type(sust_act_cells) = 2;
% vis_type(delay_act_cells) = 3;
% vis_type(delay_sup_cells) = 4;
% % barplot_by_layer(lightmod(visual_cells,end-1)',vis_type(visual_cells),ones(1,length(visual_cells)),'Light modulation index','Lightmodulation (evoked)')
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
plot(lightmod(:,lightcond),distfromlastch','.','color',area_color,'MarkerSize',24)
% hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
h = get(gca,'ytick');
set(gca,'yticklabel',h);
xlim([-1 1])
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--')
% legend('low','high')
xlabel('Light modulation index ','Fontsize',16)
ylabel(strcat('Depth (um)'),'Fontsize',16)
print(gcf, '-dpng','lightmodbydepth_bottom')
print(gcf,'-painters','-depsc','lightmodbydepth_bottom')

figure;
subplot(111)
plot(lightmod(:,lightcond),abs(distfromfirstch)','.','color',area_color,'MarkerSize',24)
% hold on; plot(lightmod_early(1:28,3),distfromlastch(1:28),'r.','MarkerSize',24)
h = get(gca,'ytick');
% set(gca,'yticklabel',h);
view(0,270)
xlim([-1 1])
ylim([0 max(abs(distfromfirstch))])
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
% legend('low','high')
xlabel('Light modulation index ','Fontsize',24)
ylabel(strcat('Depth (um)'),'Fontsize',24)
set(gca,'fontsize',18,'linewidth',2)
print(gcf, '-dpng','lightmodbydepth_top')
print(gcf,'-painters','-depsc','lightmodbydepth_top')

%by type
figure;
subplot(111)
scatter(lightmod(:,lightcond),abs(distfromfirstch)',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',area_color,'markerfacealpha',1)
hold on;
scatter(lightmod(enh_cells,lightcond),abs(distfromfirstch(enh_cells))',75,'filled','markerfacecolor',area_color,'markerfacealpha',.5)
scatter(lightmod(supp_cells,lightcond),abs(distfromfirstch(supp_cells))',75,'filled','markerfacecolor',area_color,'markerfacealpha',1)
h = get(gca,'ytick');
% set(gca,'yticklabel',h);
view(0,270)
xlim([-1 1])
ylim([0 max(abs(distfromfirstch))])
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
% legend('low','high')
xlabel('Light modulation index ','Fontsize',24)
ylabel(strcat('Depth (um)'),'Fontsize',24)
set(gca,'fontsize',18,'linewidth',2)
print(gcf, '-dpng','lightmodbydepth_top_bytype')
print(gcf,'-painters','-depsc','lightmodbydepth_top_bytype')

% 
% % trying by shank...
% shank = unique(shk);    % this is INCORRECT - shk is only from last experiment, not clean units
% dist = unique(distfromfirstch);
% for i = 1:length(dist)
%     for sh = 1:length(shank)
%         mean_lm(i,sh) = nanmean(lightmod((distfromfirstch==dist(i)&shk(clean_units)==shank(sh)),end));
%     end
% end
% mean_lm(isnan(mean_lm)) = 0;
% figure;
% for i = 1:size(mean_lm,2)
%     subplot(1,size(mean_lm,2),i)
%     bar(unique(distfromfirstch),mean_lm(:,i))
%     hold on
%     plot(distfromfirstch(ismember(clean_units,find(shk==i)))',lightmod(ismember(clean_units,find(shk==i)),end),'.','color',[0 .8 .7])
%     view(90,90)
%     ylim([-1 1])
%     xlim([0 max(distfromfirstch)])
%     h = get(gca,'xtick');
%     set(gca,'xticklabel',h*25);
%     title(sprintf('shank%d',i))
%     ylabel('Light modulation index ','Fontsize',12)
%     xlabel(strcat('Depth in LP (in um)'),'Fontsize',12)
% end
%     print(gcf, '-dpng','lightmodbyshank')


%%
if contains(exp_type,'halo','ignorecase',1)
color_mat = [0 0 0; .9 0 .3; 0.50, 0.0780, 0.10]; % for halo (red)
    nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];      % make red for halo
else
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4]; % % for lighter-shade dots
    nonvis_color_mat = [.5 .5 .5; .75 .8 1; .7 .8 .7];  % for lighter-shade dots
end


% % FR light vs no light, by shank
% plot_scatter(FRev(:,[1 2]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byshk', {'Most medial','','','Most lateral'}, 1)     % first lightcond pwr
% plot_scatter(FRev(:,[1 end]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byshk', {'Most medial','','','Most lateral'}, 1)   % last lightcond pwr
% % FR light vs no light - HIGH pwr, light onset, by shank
% plot_scatter(FRonset(:,[1 end]), shk(clean_units), {[.6 .6 .6],[0 0 1], [1 0 0], [0 1 0]}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FRonset_high_byshk', {'Most medial','','','Most lateral'}, 1)

% FR light vs no light - visual vs nonvisual units
% plot_scatter(FRev(:,[1 end-1]), ~isnan(sum(ppr_v,2)), {[.5 .5 .5],area_color}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis', {'other','activated'}, 0)     % first lightcond pwr
% plot_scatter(FRbl(:,[1 end-1]), ~isnan(sum(ppr,2)), {[.5 .5 .5],area_color}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_bl', {'other','activated'}, 0)     % last lightcond pwr
% plot_scatter(FRev(:,[1 end-1]), ones(1,size(FRev,1)), {area_color}, 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis', {}, 1)     % first lightcond pwr
% plot_scatter(FRbl(:,[1 end-1]), ones(1,size(FRbl,1), {,area_color}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_bl', {}, 0)     % last lightcond pwr


plot_scatter(FRev(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
% plot_scatter(FRev(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, 'Spks/s (light OFF)', 'Spks/s (light ON - high)', 'FR_high_byvis', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
plot_scatter(FRev(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_byvis', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% if strcmpi(exp_type,'trains')
%     plot_scatter(FRev(:,[1 end-1]), (vis_sig < .025)|(vis_sig_ons < .025), {[.7 .8 .7],[.7 0 1]}, 'Spks/s (light OFF)', 'Spks/s (light ON - med)', 'FR_med_byvis', {'Nonvisual','Visual'}, 1)     % last lightcond pwr
% end

% FR light vs no light - visual vs nonvisual units, blank trials
plot_scatter(FRbl(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_bl', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
plot_scatter(FRbl(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_byvis_bl', {'Nonvisual','Visual'}, 1)     % last lightcond pwr

% if trains experiment, also plot by "photo-activated" units
if contains(exp_type,'trains')    
    plot_scatter(FRev(:,[1 2]), nonans_v{lightcond-1}, {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_bytrains', {'Other','Hz-activated'}, 1)     % first lightcond pwr
    plot_scatter(FRev(:,[1 lightcond+1]), nonans_v{lightcond-1}, {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_bytrains', {'Other','Hz-activated'}, 1)     % last lightcond pwr
    plot_scatter(FRbl(:,[1 2]), nonans_bl{lightcond-1}, {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_bytrains_bl', {'Other','Hz-activated'}, 1)     % first lightcond pwr
    plot_scatter(FRbl(:,[1 lightcond+1]), nonans_bl{lightcond-1}, {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_bytrains_bl', {'Other','Hz-activated'}, 1)     % last lightcond pwr
end

plot_scatter(FRvison(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_onset', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
plot_scatter(FRvison(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_byvis_onset', {'Nonvisual','Visual'}, 1)     % last lightcond pwr

% FR light vs no light - visual vs nonvisual units, preferred trials
plot_scatter(FRpref(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON - low)', 'FR_low_byvis_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
plot_scatter(FRpref(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', 'FR_high_byvis_pref', {'Nonvisual','Visual'}, 1)     % last lightcond pwr

% change in preferred FR light vs no light - visual vs nonvisual units
plot_scatter(FRpref(:,[1 2])-FRbl(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-low)', 'FR_low_vischange_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
plot_scatter(FRpref(:,[1 lightcond+1])-FRbl(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked change in FR (Spks/s-light ON-high)', 'FR_high_vischange_pref', {'Nonvisual','Visual'}, 1)     % first lightcond pwr

% change in evoked FR light vs no light - visual vs nonvisual units
plot_scatter(FRev(:,[1 2])-FRbl(:,[1 2]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.25 .75], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-low)', 'FR_low_vischange', {'Nonvisual','Visual'}, 1)     % first lightcond pwr
plot_scatter(FRev(:,[1 lightcond+1])-FRbl(:,[1 lightcond+1]), (vis_sig < .025)|(vis_sig_ons < .025), {area_color,area_color}, [.5 1], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON-high)', 'FR_high_vischange', {'Nonvisual','Visual'}, 1)     % first lightcond pwr

% change in bursting
if contains(exp_type,'trains')
    diffvar = nonans_v{lightcond-1};
    diffvar20hz = nonans_v{lightcond};
    diffvarname = {'Other','Hz-activated'};
    burstrate_sig(end) = signrank(burstrate(ismember(exp_num(clean_units),cons_exps)),burstrate_light(ismember(exp_num(clean_units),cons_exps),i));     %re-do last sig values for trains exps (b/c inconsistent high-freq trains)
    burstrate_dir(end) = sign(nanmedian(burstrate_light(ismember(exp_num(clean_units),cons_exps),i)-burstrate(ismember(exp_num(clean_units),cons_exps))));
    burstresprate_sig(end) = signrank(burstresprate(ismember(exp_num(clean_units),cons_exps)),burstresprate_light(ismember(exp_num(clean_units),cons_exps),i));     %re-do last sig values for trains exps (b/c inconsistent high-freq trains)
    burstresprate_dir(end) = sign(nanmedian(burstresprate_light(ismember(exp_num(clean_units),cons_exps),i)-burstresprate(ismember(exp_num(clean_units),cons_exps))));
    plot_scatter([burstrate(ismember(exp_num(clean_units),cons_exps)) burstrate_light(ismember(exp_num(clean_units),cons_exps),end)], diffvar20hz, {area_color,area_color}, [.5 1], 'Bursting rate (light OFF)', 'Bursting rate (light ON-high)', 'Bursts_high', diffvarname, 1)     % first lightcond pwr
else
    diffvar = (vis_sig < .025)|(vis_sig_ons < .025);
    diffvarname = {'Nonvisual','Visual'};
    plot_scatter([burstrate burstrate_light(:,end)], diffvar, {area_color,area_color}, [.5 1], 'Bursting rate (light OFF)', 'Bursting rate (light ON-high)', 'Bursts_high', diffvarname, 1)     % first lightcond pwr
end
plot_scatter([burstrate burstrate_light(:,1)], diffvar, {area_color,area_color}, [.5 1], 'Bursting rate (light OFF)', 'Bursting rate (light ON-low)', 'Bursts_low', diffvarname, 1)     % first lightcond pwr
plot_scatter([burstrate burstrate_light(:,2)], diffvar, {area_color,area_color}, [.5 1], 'Bursting rate (light OFF)', 'Bursting rate (light ON-med)', 'Bursts_med', diffvarname, 1)     % first lightcond pwr

%% average PSTHs
% color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4]; % for graphing purposes (first is black, last is green)
psthV(:,:,:) = reshape(cell2mat(arrayfun(@(x) x.psthVisual([conds end],:), FRs(clean_units),'uniformoutput',0)),length(conds)+1,size(FRs(1).psthVisual,2),length(clean_units));
psthBl(:,:,:) = reshape(cell2mat(arrayfun(@(x) x.psthBlank([conds end],:), FRs(clean_units),'uniformoutput',0)),length(conds)+1,size(FRs(1).psthBlank,2),length(clean_units));
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

% rep_baseline = reshape(repmat(repmat(FRb,size(psthV,2),1),size(psthV,1),1),size(psthV,1),size(psthV,2),size(psthV,3));
% bs_psth = psthV-rep_baseline;       % baseline-subtracted psth
% % test_psth = bs_psth./repmat(repmat(max(max(bs_psth)),size(bs_psth,1),1),1,size(bs_psth,2));     % normalize
% % psthZ = zscore(bs_psth,0,2);        % zscore across timepoints
% mean_bs = repmat(mean(bs_psth(:,1:20,:),2),1,100,1);    % prestim currently hardcoded! 20 time bins x 25ms each = 500ms
% std_bs = repmat(std(bs_psth(:,1:20,:),[],2),1,100,1);
% psthZ = (bs_psth-mean_bs)./std_bs;      % setting 0 to the mean during prestim period only

if contains(exp_type,'halo')
    mean_bs = repmat(mean(psthV(:,1:10,:),2),1,100,1);    % prestim currently hardcoded! 10 time bins x 25ms each = 250ms (because halo experiments start during prestim period!)
else
    mean_bs = repmat(mean(psthV(:,1:20,:),2),1,100,1);    % prestim currently hardcoded! 20 time bins x 25ms each = 500ms
end
mean_bs_bl = repmat(mean(psthBl(:,1:20,:),2),1,100,1);
mean_bs(mean_bs==0) = nan;     % because otherwise could get infinity when normalizing by baseline
mean_bs_bl(mean_bs_bl==0) = nan;
% std_bs = repmat(std(psthV(:,1:20,:),[],2),1,100,1);
% psthZ = (psthV-mean_bs)./std_bs;      % setting 0 to the mean during prestim period only

norm_psth = psthV./mean_bs; % normalizes to prestim baseline, per condition (so that prestim=1)
norm_psth_bl = psthBl./mean_bs_bl;

% % define outliers
% max_visresps = squeeze(max(norm_psth(1,:,:),[],2)); % each unit's maximum (across time) normalized change in FR during NO LIGHT condition
% outlier_thresh = mean(max_visresps)+2*std(max_visresps);    % define outlier threshold as 2*standard dev in max vis responses amove the mean (across units) in max vis responses
% time_av = squeeze(mean(norm_psth,2));   % average across timepoints
% [~,outlier_units] = find(time_av>outlier_thresh);
% norm_psth(:,:,unique(outlier_units)) = nan;
% norm_psth_bl(:,:,outlier_units) = nan;
% % % alternative: use blank trials?
% % [a,b] = find(time_av_bl>max(prctile(norm_psth_bl(1,:,:),75,3)+iqr(norm_psth_bl(1,:,:),3)*10));
outlier_units = [];

% mean_visUp = mean(psthZ(:,:,clean_units(visUp)),3);
% test_mean = mean(psthZ(:,:,visual_cells),3);    % get average zscores across visually-responsive units
norm_mean = nanmean(norm_psth,3);
norm_se = nanstd(norm_psth,0,3)./sqrt(size(norm_psth,3));
norm_mean_bl = nanmean(norm_psth_bl,3);
norm_se_bl = nanstd(norm_psth_bl,0,3)./sqrt(size(norm_psth_bl,3));

% VISUAL trials
pop_fig3= figure;
xlim([-.475 2])   % ticks mark the END of 25ms bins
hold on;
for i = 1:size(norm_psth,1)
    shadedErrorBar([-.5:.025:1.975],norm_mean(i,:), norm_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
% ylim([-.5 .5])
yax = get(gca,'YLim');
% yax = [-min(abs(yax)) min(abs(yax))];
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
start_times = round([params(:).av_light_start],2)-unique([params(:).prestim]);  % round to nearest hundreth
stim_durs = round([params(:).light_dur],1);
stim_durs = stim_durs(stim_durs>0);
if length(unique(start_times))>1
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
    end
else
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_zscore';
print(pop_fig3,'-dpng',save_fig_name)
% print2eps(save_fig_name,pop_fig3)
print(pop_fig3,'-painters','-depsc',save_fig_name)

% BLANK trials
pop_fig4= figure;
xlim([-.475 2])   % ticks mark the END of 25ms bins
hold on;
for i = 1:size(norm_psth,1)
    shadedErrorBar([-.5:.025:1.975],norm_mean_bl(i,:), norm_se_bl(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
% ylim([-.5 .5])
yax = get(gca,'YLim');
% yax = [-min(abs(yax)) min(abs(yax))];
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
start_times = round([params(:).av_light_start],2)-unique([params(:).prestim]);  % round to nearest hundreth
stim_durs = round([params(:).light_dur],1);
stim_durs = stim_durs(stim_durs>0);
if length(unique(start_times))>1
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
    end
else
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_zscore_blanks';
print(pop_fig4,'-dpng',save_fig_name)
% print2eps(save_fig_name,pop_fig4)
print(pop_fig4,'-painters','-depsc',save_fig_name)

if ~contains(exp_type,'halo')       % temp
    % population psths: suppressed vs. activated cells (visual)
    pop_psth = figure;
    subplot(121)
    supp_cells_clean = supp_cells(~ismember(supp_cells,outlier_units));
    supp_mean = nanmean(norm_psth(:,:,supp_cells_clean),3);
    supp_se = nanstd(norm_psth(:,:,supp_cells_clean),0,3)./sqrt(length(supp_cells_clean));
    hold on;
    for i = 1:size(norm_psth,1)
        shadedErrorBar([-.475:.025:2],supp_mean(i,:), supp_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
    end
    title(sprintf('Suppressed cells (n=%d)',length(supp_cells_clean)),'fontsize',14);
    xlabel('Time from visual stim onset (s)','fontsize',14)
    ylabel('Normalized firing rate (spks/s)','fontsize',14)
    set(gca,'fontsize',16,'linewidth',2);
    xlim([-.475 2])   % ticks mark the END of 25ms bins
    yax = ylim;
    ylim(yax);
    patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)

    subplot(122)
    enh_cells_clean = enh_cells(~ismember(enh_cells,outlier_units));
    enh_mean = nanmean(norm_psth(:,:,enh_cells_clean),3);
    enh_se = nanstd(norm_psth(:,:,enh_cells_clean),0,3)./sqrt(length(enh_cells_clean));
    hold on;
    for i = 1:size(norm_psth,1)
        shadedErrorBar([-.475:.025:2],enh_mean(i,:), enh_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
    end
    title(sprintf('Activated cells (n=%d)',length(enh_cells_clean)),'fontsize',14);
    xlabel('Time from visual stim onset (s)','fontsize',14)
    ylabel('Normalized firing rate (spks/s)','fontsize',14)
    set(gca,'fontsize',16,'linewidth',2);
    xlim([-.475 2])   % ticks mark the END of 25ms bins
    yax = ylim;
    ylim(yax);
    patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    
    set(gcf, 'Position', [100, 100, 1000, 420])
    print(pop_psth, '-dpng','PopulationPSTH_subtypes_visual')
%     print2eps('PopulationPSTH_subtypes_visual',pop_psth)
    print(pop_psth,'-painters','-depsc','PopulationPSTH_subtypes_visual')


    % population psths: suppressed vs. activated cells (blanks)
    pop_psth = figure;
    subplot(121)
    supp_mean = nanmean(norm_psth_bl(:,:,supp_cells_clean),3);
    supp_se = nanstd(norm_psth_bl(:,:,supp_cells_clean),0,3)./sqrt(length(supp_cells_clean));
    hold on;
    for i = 1:size(norm_psth_bl,1)
        shadedErrorBar([-.475:.025:2],supp_mean(i,:), supp_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
    end
    title(sprintf('Suppressed cells (n=%d)',length(supp_cells_clean)),'fontsize',14);
    xlabel('Time from visual stim onset(s)','fontsize',14)
    ylabel('Normalized firing rate (spks/s)','fontsize',14)
    set(gca,'fontsize',16,'linewidth',2);
    xlim([-.475 2])   % ticks mark the END of 25ms bins
    yax = ylim;
    ylim(yax);
    patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)

    subplot(122)
    enh_mean = nanmean(norm_psth_bl(:,:,enh_cells_clean),3);
    enh_se = nanstd(norm_psth_bl(:,:,enh_cells_clean),0,3)./sqrt(length(enh_cells_clean));
    hold on;
    for i = 1:size(norm_psth,1)
        shadedErrorBar([-.475:.025:2],enh_mean(i,:), enh_se(i,:),{'Color',color_mat(i,:),'linewidth',2},1);
    end
    title(sprintf('Activated cells (n=%d)',length(enh_cells_clean)),'fontsize',14);
    xlabel('Time from visual stim onset (s)','fontsize',14)
    ylabel('Normalized firing rate (spks/s)','fontsize',14)
    set(gca,'fontsize',16,'linewidth',2);
    xlim([-.475 2])   % ticks mark the END of 25ms bins
    yax = ylim;
    ylim(yax);
    patch([.5 .5 1.5 1.5 .5],[yax(1) yax(2) yax(2) yax(1) yax(1)], [0 .1 1], 'LineStyle', 'none', 'FaceAlpha',.15 );
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    
    set(gcf, 'Position', [100, 100, 1000, 420])
    print(pop_psth, '-dpng','PopulationPSTH_subtypes_blanks')
%     print2eps('PopulationPSTH_subtypes_blanks',pop_psth)
    print(pop_psth,'-painters','-depsc','PopulationPSTH_subtypes_blanks')
end

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
%% pie graph of number of enhanced vs. suppressed vs. unaffected units
figure;
piegraph = pie([length(supp_cells) length(enh_cells) length(clean_units)-length(light_cells)]);
piegraph_labels = {'Suppressed','Activated','Other'};
piegraph_labels = piegraph_labels([~isempty(supp_cells) ~isempty(enh_cells) length(clean_units)-length(light_cells)>0]);
colormap([area_color; area_color; 1 1 1])
set(piegraph([1:2:end]),'edgecolor',area_color)
set(piegraph([2:2:end]),'fontsize',16)
for i=2:2:length(piegraph)
    set(piegraph(i),'string',sprintf('%s (%s)',piegraph_labels{i/2},piegraph(i).String))
    if strcmpi(piegraph_labels{i/2},'Activated')
        set(piegraph(i-1),'facealpha',.5) % enhanced units will be lighter colored
    end
end
print(gcf,'-dpng','Piegraph_lighteffects')
print(gcf,'-painters','-depsc','Piegraph_lighteffects')

%% OSI 
plot_scatter(OSI_CV(:,[1 lightcond+1]), ones(1,length(clean_units)), {area_color}, 1, 'OSI(CV) (light OFF)', 'OSI(CV) (light ON)', 'OSI(CV)', {'All cells'}, 1)     % first lightcond pwr
plot_scatter(OSI(:,[1 lightcond+1]), ones(1,length(clean_units)), {area_color}, 1,'OSI (light OFF)', 'OSI (light ON)', 'OSI', {'All cells'}, 1)     % first lightcond pwr
plot_scatter(DSI_CV(:,[1 lightcond+1]), ones(1,length(clean_units)), {area_color}, 1,'DSI(CV) (light OFF)', 'DSI(CV) (light ON)', 'DSI(CV)', {'All cells'}, 1)     % first lightcond pwr
plot_scatter(DSI(:,[1 lightcond+1]), ones(1,length(clean_units)), {area_color}, 1,'DSI (light OFF)', 'DSI (light ON)', 'DSI', {'All cells'}, 1)     % first lightcond pwr

plot_scatter(OSI_CV(tuned_cells,[1 lightcond+1]), lightmod(tuned_cells,lightcond)>0, {'k','b'}, [1 1], 'OSI(CV) (light OFF)', 'OSI(CV) (light ON - high)', 'OSI(CV)_bylightmod', {'Light-suppressed','Light-enhanced'}, 2)     % first lightcond pwr
plot_scatter(OSI(tuned_cells,[1 lightcond+1]), lightmod(tuned_cells,lightcond)>0, {'k','b'}, [1 1], 'OSI (light OFF)', 'OSI (light ON - high)', 'OSI_bylightmod', {'Light-suppressed','Light-enhanced'}, 2)     % first lightcond pwr
plot_scatter(DSI_CV(tuned_cells,[1 lightcond+1]), lightmod(tuned_cells,lightcond)>0, {'k','b'}, [1 1], 'DSI(CV) (light OFF)', 'DSI(CV) (light ON - high)', 'DSI(CV)_bylightmod', {'Light-suppressed','Light-enhanced'}, 2)     % first lightcond pwr
plot_scatter(DSI(tuned_cells,[1 lightcond+1]), lightmod(tuned_cells,lightcond)>0 , {'k','b'}, [1 1], 'DSI (light OFF)', 'DSI (light ON - high)', 'DSI_bylightmod', {'Light-suppressed','Light-enhanced'}, 2)     % first lightcond pwr

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
fileID = fopen('driver_results.txt','w');
fprintf(fileID,'Number of visually responsive cells: %d of %d\r\n',length(visual_cells),num_units);
fprintf(fileID,'Percent visually responsive: %.2f\r\n',100*length(visual_cells)/num_units);
fprintf(fileID,'Number of light-modulated cells: %d of %d\r\n',length(light_cells),num_units);
fprintf(fileID,'Percent light-modulated: %.2f\r\n', 100*length(light_cells)/num_units);
fprintf(fileID,'Number of outlier units: %d of %d\r\n',length(outlier_units),num_units);
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

fileID2 = fopen('driver_stats.txt','w');
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, ALL cells: %s (%s) \r\n',num2str(round(lightsig_all,3)),num2str(lightsig_all_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, ALL cells: %s (%s) \r\n',num2str(round(lightsig_bl,3)),num2str(lightsig_bl_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, VISUAL cells: %s (%s) \r\n',num2str(round(lightsig_vis,3)),num2str(lightsig_vis_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, VISUAL cells: %s (%s) \r\n',num2str(round(lightsig_blvis,3)),num2str(lightsig_blvis_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, NONVISUAL cells: %s (%s) \r\n',num2str(round(lightsig_nonvis,3)),num2str(lightsig_nonvis_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, NONVISUAL cells: %s (%s) \r\n',num2str(round(lightsig_blnonvis,3)),num2str(lightsig_blnonvis_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on pref FR (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_vispref,3)),num2str(lightsig_vispref_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on visual FR change (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_visFRdelt,3)),num2str(lightsig_visFRdelt_dir));
fprintf(fileID2,'Signed-rank test of sig light effect on pref visual FR change (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_prefFRdelt,3)),num2str(lightsig_prefFRdelt_dir));

fprintf(fileID2,'\r\n');
fprintf(fileID2,'Signed-rank test of sig OSI_CV change (low to high conditions) TUNED cells: %s (%s) \r\n',num2str(round(osiCVsig_tuned,3)),num2str(osiCVsig_tuned_dir));
fprintf(fileID2,'Signed-rank test of sig OSI change (low to high conditions), TUNED cells: %s (%s) \r\n',num2str(round(osisig_tuned,3)),num2str(osisig_tuned_dir));
fprintf(fileID2,'Signed-rank test of sig DSI_CV change (low to high conditions), TUNED cells: %s (%s) \r\n',num2str(round(dsiCVsig_tuned,3)),num2str(dsiCVsig_tuned_dir));
fprintf(fileID2,'Signed-rank test of sig DSI change (low to high conditions), TUNED cells: %s (%s) \r\n',num2str(round(dsisig_tuned,3)),num2str(dsisig_tuned_dir));
fprintf(fileID2,'Signed-rank test of sig OSI_CV change (low to high conditions), ALL cells:%s (%s) \r\n',num2str(round(osiCVsig_vis,3)),num2str(osiCVsig_vis_dir));
fprintf(fileID2,'Signed-rank test of sig OSI change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(osisig_vis,3)),num2str(osisig_vis_dir));
fprintf(fileID2,'Signed-rank test of sig DSI_CV change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(dsiCVsig_vis,3)),num2str(dsiCVsig_vis_dir));
fprintf(fileID2,'Signed-rank test of sig DSI change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(dsisig_vis,3)),num2str(dsisig_vis_dir));

fprintf(fileID2,'\r\n');
fprintf(fileID2,'Signed-rank test of burst rate change (low to high conditions) ALL cells: %s (%s) \r\n',num2str(round(burstrate_sig,3)),num2str(burstrate_dir));
fprintf(fileID2,'Signed-rank test of burst response rate change (low to high conditions) ALL cells: %s (%s) \r\n',num2str(round(burstresprate_sig,3)),num2str(burstresprate_dir));

fclose(fileID2);

fileID3 = fopen('driver_cleanunits.txt','w');
fprintf(fileID3,'%d\r\n',[unitinfo(clean_units).name]);
fclose(fileID3);

if contains(exp_type,'trains')
    hz_cells = find(nonans_v{lightcond-1});   % hz-activated at 10hz in visual trials
    fac_cells = find(ppr_v{1}(:,1)>1);
    dep_cells = find(ppr_v{1}(:,1)<=1);
    save('unit_types.mat','hz_cells','fac_cells','dep_cells')
else
    save('unit_types.mat','supp_cells','enh_cells')
end

end

function plot_scatter(data, color_var, colors, alpha, xlab, ylab, title, leg, lobf)
% if lobf = 1, make one lobf regardless of variables; if lobf = 2, make
% separate lobfs for each variable
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
    scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75, 'filled','markerfaceColor', colors{i},'markerfacealpha',alpha(i));
end
set(gca,'Fontsize',18,'linewidth',2)
xlabel(xlab,'Fontsize',18)
ylabel(ylab,'Fontsize',18)
xax = get(gca,'XLim');
yax = get(gca,'YLim');
x = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)),100);
y=x;
plot(x,y,'k--','color','k');
line([min(xax(1),yax(1)) max(xax(2),yax(2))],[0 0],'color','k')
line([0 0], [min(xax(1),yax(1)) max(xax(2),yax(2))],'color','k')
xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
if lobf == 1
    fittedX = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)), 200);
    coeffs = polyfitZero(data(:,1), data(:,2), 1);
    fittedY = polyval([0 coeffs], fittedX);
    plot(fittedX,fittedY,'color', colors{1},'linewidth',2)
elseif lobf == 2
    fittedX = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)), 200);
    for i = 1:length(vars)
        if sum(color_var==vars(i)) > 1    % comment out if you don't want to separately calculate lobfs
            coeffs(i,:) = polyfitZero(data(color_var==vars(i),1), data(color_var==vars(i),2), 1);     
            fittedY(i,:) = polyval([0 coeffs(i,:)], fittedX);
            plot(fittedX,fittedY(i,:),'color', colors{i},'linewidth',2)
        end
    end

else    % plot median
    for i = 1:length(vars)
        plot(median(data(color_var==vars(i),1)),median(data(color_var==vars(i),2)), 'marker','+', 'Color', colors{i},'markersize',28,'linewidth',4);
    end
end
if ~isempty(leg)
    l=legend(leg,'location','best');
end
set(l,'fontsize',18)
print(f, '-dpng',title)
% print2eps(title,f)        % doesn't seem to work with new matlab...
print(f,'-painters','-depsc',title)
end

% % calculating response latency for dLGN and TRN
% for b = 1:length(clean_units)
%     test(:,:,b) = unitinfo(clean_units(b)).rast(blank_trials{exp_num(clean_units(b))}(1:144),:);
% end
% 
% bs = mean(mean(test(:,1:500,:),2)); % mean across prestim timepoints and all blank trials
% bs = squeeze(bs);
% 
% for b = 1:length(clean_units)
%     bs_std(b,1) = std(reshape(test(:,1:500,b),1,numel(test(:,1:500,b))));
% end
% 
% test_mean(:,:) = mean(test);
% test_bs = test_mean./bs';
% bs_std = std(test_bs(1:500,:));
% bs_mean = mean(test_bs(1:500,:));
% 
% test_z = (mean(test_bs,2)-ones(size(test_bs,1),1))./bs_std;
% figure;plot(mean(test_z,2))
% test_meanz = mean(test_z,2);
% find(test_meanz(1000:end)>2,1,'first')