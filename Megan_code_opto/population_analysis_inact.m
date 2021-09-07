function population_analysis_inact(area, pop, exp_type)

% new version created 10/21/20 to better accomodate inactivation experiments with variable different opto conditions experiments (MAK)
% checked throuroughly 1/9/21
% inputs:
% area = e.g., 'LPlateral', 'LPmedial', 'LGN'
% pop = string indicating which population of interest (e.g., 'driver',
% 'modulator', 'ChR2', 'Halo', etc.)
% exp_type = 'ramp', 'trains', 'step', or 'size' 

% inputs to assign for different experiments:
% exp_paths - cell of all experimental directories you want to use
% shanks - cell of shanks (0:3 medial:lateral for bottom probes with
%   contacts pointing posterior) for each experiment (length(exp_paths) ==
%   length(shanks))
% probe (e.g., '128D_bottom' (just for z- and shank-locations, so 128DN
%   vs 128D doesn't matter). currently assumes all experiments used the
%   same probe type
% lightconds = cell array of light conditions (1 index) you want to analyze
%   from each experiment. order indicates order in which this script will
%   process them (e.g., [2 1 3])
%   -for SC/V1 combo experiments, V1inactivation should be 1st entry and SC
%   should be second. 
%   -for L5/L6 combo experiments, L5inactivation should be 1st entry and L6
%   should be second.
% lightcond = index for which light condition (from lightconds) to use to assess
% significance 

% If downloading data from Mendeley, start from line 286 ('Identify
% experiment parameters for multi-experiment analyses') (but will need to
% set your own main_dir and fig_dir

vissigdef = 'blanks';  % defined coparison for determining visual responsivity: 'prepost' for comparing pre- to post-stimulus in preferred visual trials, or 'blanks' to compare preferred visual to blank trials
%^use 'prepost' for dtA experiments because no blank trials were used

%% set up directories, identify experiments to analyze
main_dir = 'H:\LPproject\LPresults';
if strcmpi(pop,'driver')
    if strcmpi(exp_type,'step_halo') && strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
        exp_paths = {'J:\LPproject\NPRH18\LP\2019-11-12_18-06-34_LP_halo',...
            'J:\LPproject\NPRH19\LP\2019-11-13_17-32-45_LP_halo',...
            'J:\LPproject\NPRH21\LP\2019-12-13_17-37-53_Halo',...
            'J:\LPproject\NPRH23\LP\2019-12-19_12-09-54_Halo',...
            'J:\LPproject\NPRH24\LP\2019-12-18_12-16-55_Halo',...
            'J:\LPproject\NPRH5\LP\2019-08-16_15-31-01_LP_halo',...
            'J:\LPproject\NPRH27\LP\2020-01-28_13-01-16_Halo_REAL'};
        shanks = {[2:3],[2],[3],[2:3],[3],[2:3],[3]};
        probe = '128D_bottom';
        lightconds = {1,1,1,1,1,2,1};  % should be 2 b/c NPRH5 had two light conds (2nd had same start time as the rest of these experiments)
        lightcond = 1;

    elseif strcmpi(exp_type,'step_gtacr') && strcmpi(area, 'LPlateral')         % removed NPRG1 because I think expression is also in RSP...
        exp_paths = {'Z:\L5\stGtACR\NPRSC4\LPl\2020-07-21_15-45-08_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC8\LPl\2020-08-31_16-14-35_driftgrat_2LEDs'};
        shanks = {[0 1],[0:2],[1:2],[0:2]};
        probe = '128D_bottom';  
        lightconds = {1,1,2,2}; % for NPRSC6,8 SCinact was first cond)
        lightcond = 1;
        
        elseif strcmpi(exp_type,'step_gtacr') && strcmpi(area, 'LP')  %combining LPl & LPrm
        exp_paths = {'Z:\L5\stGtACR\NPRSC4\LPl\2020-07-21_15-45-08_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC5\LPrm\2020-07-22_17-55-27_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPrm\2020-09-02_19-08-32_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC8\LPl\2020-08-31_16-14-35_driftgrat_2LEDs'};
%             'Z:\L5\stGtACR\NPRSC8\LPrm\2020-08-31_20-19-41_driftgrat_2LEDs'};  % too few trials, very messy
        shanks = {[0 1],[0:2],[0:1],[1:2],[0:3],[0:2]}; % [0 1]};
        probe = '128D_bottom';  
        lightconds = {1,1,1,2,2,2};
        lightcond = 1;
        
        elseif strcmpi(exp_type,'step_gtacr_wRF') && strcmpi(area, 'LP')  %combining LPl & LPrm - excluding NPRSC8 because couldn't get RFs
        exp_paths = {'Z:\L5\stGtACR\NPRSC4\LPl\2020-07-21_15-45-08_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC5\LPrm\2020-07-22_17-55-27_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPrm\2020-09-02_19-08-32_driftgrat_2LEDs'};
        shanks = {[0 1],[0:2],[0:1],[1:2],[0:3]};
        probe = '128D_bottom';  
        lightconds = {1,1,1,2,2};
        lightcond = 1;
    
    elseif strcmpi(exp_type,'step_gtacr') && strcmpi(area, 'LPmedial')
        exp_paths = {'Z:\L5\stGtACR\NPRSC5\LPrm\2020-07-22_17-55-27_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\LPrm\2020-09-02_19-08-32_driftgrat_2LEDs'};
        shanks = {[0:1],[0:3]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {1,2};
        lightcond = 1;
        
    elseif strcmpi(exp_type,'step_L5SC') && strcmpi(area, 'LPlateral')          % previously wasn't including NPRSC4 because some halo in RSP...?
        exp_paths = {'Z:\L5\stGtACR\NPRSC4\LPl\2020-07-21_15-45-08_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC8\LPl\2020-08-31_16-14-35_driftgrat_2LEDs'};
        shanks = {[0 1],[0:2],[1:2],[0:2]};
        probe = '128D_bottom';  
        lightconds = {[1:3],[1:3],[2 1 3], [2 1 3]}; 
        lightcond = 3;      % previously 3 (changed 12/10/20)

    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LPmedial')
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
        exp_paths = {'J:\LPproject\NPRH18\LP\2019-11-12_18-06-34_LP_halo',...
            'J:\LPproject\NPRH21\LP\2019-12-13_17-37-53_Halo',...
            'J:\LPproject\NPRH23\LP\2019-12-19_12-09-54_Halo',...
            'J:\LPproject\NPRH24\LP\2019-12-18_12-16-55_Halo',...
            'J:\LPproject\NPRH5\LP\2019-08-16_15-31-01_LP_halo',...
            'J:\LPproject\NPRH27\LP\2020-01-28_13-01-16_Halo_REAL',...
            'J:\LPproject\NPRH37\LP\2020-06-19_17-08-03_Halo'};
        shanks = {[0],[0:2],[1],[0:2],[1],[2],[0:1]};
        probe = '128D_bottom';
        lightconds = {1,1,1,1,2,1,1};  % should be 2 b/c NPRH5 had two light
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LP')      %combining LPl & LPrm
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
        exp_paths = {'J:\LPproject\NPRH18\LP\2019-11-12_18-06-34_LP_halo',...
            'J:\LPproject\NPRH21\LP\2019-12-13_17-37-53_Halo',...
            'J:\LPproject\NPRH19\LP\2019-11-13_17-32-45_LP_halo',...
            'J:\LPproject\NPRH23\LP\2019-12-19_12-09-54_Halo',...
            'J:\LPproject\NPRH24\LP\2019-12-18_12-16-55_Halo',...
            'J:\LPproject\NPRH5\LP\2019-08-16_15-31-01_LP_halo',...
            'J:\LPproject\NPRH27\LP\2020-01-28_13-01-16_Halo_REAL',...
            'J:\LPproject\NPRH37\LP\2020-06-19_17-08-03_Halo'};
        shanks = {[0 2:3],[0:3],[2],[1:3],[0:3],[1:3],[2:3],[0:1]};
        probe = '128D_bottom';
        lightconds = {1,1,1,1,1,2,1,1};  % should be 2 b/c NPRH5 had two light
        lightcond = 1;

    elseif strcmpi(exp_type,'step_halo2LEDs') && strcmpi(area,'LP') % halo in V1 and VisL
        exp_paths = {'J:\LPproject\NPRH33\LP\2020-03-19_16-38-47_halo_2LEDs',...
            'J:\LPproject\NPRH34\LP\2020-03-18_18-07-51_Halo_2LEDs',...
            'J:\LPproject\NPRH35\LP\2020-03-19_21-08-49_halo_2LEDs'};
        shanks = {[0:3],[0:3],[2:3]};
        probe = '128D_bottom';
        lightconds = {[2 1 3],[2 1 3],[2 1 3]}; % in these, lightcond 1 = visL and 2 = V1
        lightcond = 3; 

    elseif strcmpi(exp_type,'step_DTA') && strcmpi(area,'LPmedial')
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
        exp_paths = {'Z:\dTA\DDTA3\LP\2020-01-31_15-27-25_driftgrat',...
            'Z:\dTa\DDTA5\LP\2020-03-09_15-24-23_driftgrat_DTA',...
            'Z:\dTa\DDTA6\LPrm\2020-08-27_12-22-45_driftgrat'};
        shanks = {[0:2],[0:2],[0:1]};
        lightconds = {[],[],[]};
        probe = '128D_bottom';
        lightcond = 1;  % should be 2 b/c NPRH5 had two light
    elseif strcmpi(exp_type,'step_DTA') && strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\DH3\2018-06-19_13-43-37_LP_Halo_real'};
%             shanks = {[2]};
        exp_paths = {'Z:\dTA\DDTA3\LP\2020-01-31_15-27-25_driftgrat',...
            'Z:\dTa\DDTA5\LP\2020-03-09_15-24-23_driftgrat_DTA',...
            'Z:\dTa\DDTA6\LPl\2020-08-27_15-04-45_driftgrat'};
        shanks = {[3],[3],[3]};
        probe = '128D_bottom';
        lightconds = {[],[],[]};  
        lightcond = 1;
    end

elseif strcmpi(pop,'modulator')
   if strcmpi(exp_type,'step_halo') && strcmpi(area, 'LPlateral')
       exp_paths = {...   
        'J:\LPproject\MH19\LP\2018-12-04_21-43-54_LP_Halo',...
        'J:\LPproject\MH21\LP\2018-12-03_17-56-28_LP_Halo',...
        'H:\LPproject\MH24\LP\2019-04-12_13-33-18_LP_halo_real',...
            'J:\LPproject\MH25\LP\2019-04-10_17-05-08_LP_Halo',...
            'Z:\L6\MH30\LP\2020-01-03_12-56-56_Halo',...
            'J:\LPproject\MH32\LP\2020-01-03_17-07-12_Halo',...
            'J:\LPproject\MH33\LP\2020-04-09_18-50-28_halo',...
            'J:\LPproject\MH34\LP\2020-04-10_12-35-28_halo',...
            'Z:\L6\MH35\LP\2020-04-10_16-53-34_halo',...
            'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {[0],[1],[3],[1],[3],[2:3],[0],[0],[1],[1 2]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,2,2,1,1,1,1,1,1};
        lightcond = 1;  
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area, 'LPmedial')
        exp_paths = {...
    'J:\LPproject\MH27\LP\2019-08-28_19-09-30_LP_halo',...   
        'J:\LPproject\MH28\LP\2019-08-28_16-02-12_LP_halo',...
        'J:\LPproject\MH31\LP\2020-01-02_16-27-43_Halo',...
            'J:\LPproject\MH32\LP\2020-01-03_17-07-12_Halo',...
            'Z:\L6\MH35\LP\2020-04-10_16-53-34_halo'};
        shanks = {[0:2],[0:1],[0:2],[0:1],[0]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,1,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LGN')
            % these are the ones with too high light power...
%             exp_paths = {'J:\LPproject\MH15\2018-05-23_16-36-18_LP_Halo',...
%                 'J:\LPproject\MH11\2018-05-09_12-08-51_Halo_LP',...
%                 'H:\LPproject\MH9\2018-03-27_19-39-26_Halo_LP'};
%                 shanks = {[3],[2 3],[2 3]};
        exp_paths = {...
            'J:\LPproject\MH19\LP\2018-12-04_21-43-54_LP_Halo',...
            'J:\LPproject\MH21\LP\2018-12-03_17-56-28_LP_Halo',...
            'J:\LPproject\MH25\LP\2019-04-10_17-05-08_LP_Halo',...
            'H:\LPproject\MH23\LP\2019-04-11_14-40-06_LP_Halo',...
            'J:\LPproject\MH33\LP\2020-04-09_18-50-28_halo',...
            'J:\LPproject\MH34\LP\2020-04-10_12-35-28_halo',...
            'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {[1,2],[2:3],[2],[0:3],[1:3],[1:3],[3]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,2,2,1,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo_new') && strcmpi(area,'LGN')
            % these are the ones with RFs in both LGN and V1
        exp_paths = {...
            'J:\LPproject\MH33\LP\2020-04-09_18-50-28_halo',...
            'J:\LPproject\MH34\LP\2020-04-10_12-35-28_halo',...
            'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {[1:3],[1:3],[3] };
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {1,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo_wRF') && strcmpi(area,'LP')  % these are the newer exps where I also did sparsenoise for RF mapping in V1
        exp_paths = {'J:\LPproject\MH27\LP\2019-08-28_19-09-30_LP_halo',...   
        'J:\LPproject\MH28\LP\2019-08-28_16-02-12_LP_halo',... 
            'Z:\L6\MH30\LP\2020-01-03_12-56-56_Halo',...
            'Z:\L6\MH31\LP\2020-01-02_16-27-43_Halo',...
            'J:\LPproject\MH33\LP\2020-04-09_18-50-28_halo',...
            'J:\LPproject\MH34\LP\2020-04-10_12-35-28_halo',...
            'Z:\L6\MH35\LP\2020-04-10_16-53-34_halo',...
            'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {[0:2],[0 1],[3],[0:2],[0],[0],[0 1],[1 2]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,1,1,1,1,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LP')  % combined LPl and LPrm
        exp_paths = {'J:\LPproject\MH19\LP\2018-12-04_21-43-54_LP_Halo',...
        'J:\LPproject\MH21\LP\2018-12-03_17-56-28_LP_Halo',...
        'H:\LPproject\MH24\LP\2019-04-12_13-33-18_LP_halo_real',...
            'J:\LPproject\MH25\LP\2019-04-10_17-05-08_LP_Halo',...
            'J:\LPproject\MH27\LP\2019-08-28_19-09-30_LP_halo',...   
        'J:\LPproject\MH28\LP\2019-08-28_16-02-12_LP_halo',... 
            'Z:\L6\MH30\LP\2020-01-03_12-56-56_Halo',...
            'Z:\L6\MH31\LP\2020-01-02_16-27-43_Halo',...
            'J:\LPproject\MH32\LP\2020-01-03_17-07-12_Halo',...
            'J:\LPproject\MH33\LP\2020-04-09_18-50-28_halo',...
            'J:\LPproject\MH34\LP\2020-04-10_12-35-28_halo',...
            'Z:\L6\MH35\LP\2020-04-10_16-53-34_halo',...
            'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {[0],[1],[3],[1],[0:2],[0 1],[3],[0:2],[0:3],[0],[0],[0 1],[1 2]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,2,2,2,2,1,1,1,1,1,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_gtacr') && strcmpi(area,'LGN')
        exp_paths = {'Y:\L6\stGtACR\MG2\LPl\2020-08-07_13-59-20_driftgrat_2lightpwrs',...
            'Z:\L6\stGtACR\MG4\LPl\2020-10-03_15-01-44_driftgrat_2LEDs',...
            'Z:\L6\stGtACR\MG6\LPl\2020-12-20_19-42-41_LP_haloREAL'};
        shanks = {[0:3],[2,3],[2 3] };
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {2,2,1};
        lightcond = 1;
     elseif strcmpi(exp_type,'step_gtacr') && strcmpi(area,'LPlateral')
        exp_paths = {'Z:\L6\stGtACR\MG6\LPl\2020-12-20_19-42-41_LP_haloREAL'};
        shanks = {[0 1] };
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_DTA') && strcmpi(area,'LPlateral')
        exp_paths = {'Z:\dTa\MDTA1\LP\2020-03-12_18-22-03_driftgrat_DTA',...
            'Z:\dTa\MDTA2\LP\2020-03-12_12-32-48_driftgrat_DTA',...
            'Z:\dTa\MDTA3\LP\2020-03-13_17-32-38_driftgrat_DTA'};
        shanks = {0,0,3};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {[],[],[]};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_DTA') && strcmpi(area,'LGN')
        exp_paths = {'Z:\dTa\MDTA1\LP\2020-03-12_18-22-03_driftgrat_DTA',...
            'Z:\dTa\MDTA2\LP\2020-03-12_12-32-48_driftgrat_DTA'};
        shanks = {[2 3],[1:3]};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {[],[]};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_DTA') && strcmpi(area,'LPmedial')
        exp_paths = {'Z:\dTa\MDTA3\LP\2020-03-13_17-32-38_driftgrat_DTA'};
        shanks = {0:2};
        probe = '128D_bottom';  % doesn't matter if it was DN or D
        lightconds = {[]};
        lightcond = 1;
   end

   elseif strcmpi(pop,'CTRL')
    if strcmpi(exp_type,'step_halo') && strcmpi(area,'LGN')
        exp_paths = {'Z:\Controls\MHCTRL8\LP\2020-05-08_13-40-38_Halo',...
            'Z:\Controls\MHCTRL9\LP\2020-07-02_17-39-17_Halo'};
        shanks = {[2:3],[2:3]};
        probe = '128DN_bottom';
        lightconds = {1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LPlateral')
        exp_paths = {'Z:\Controls\MHCTRL6\2019-07-15_14-58-36_LP_halo',...
            'Z:\Controls\MHCTRL8\LP\2020-05-08_13-40-38_Halo',...
            'Z:\Controls\MHCTRL9\LP\2020-07-02_17-39-17_Halo'};
        shanks = {[2:3],[0:1],[0:1]};
        probe = '128DN_bottom';
        lightconds = {2,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'all')
        exp_paths = {'Z:\Controls\MHCTRL6\2019-07-15_14-58-36_LP_halo',...
            'Z:\Controls\MHCTRL8\LP\2020-05-08_13-40-38_Halo',...
            'Z:\Controls\MHCTRL9\LP\2020-07-02_17-39-17_Halo'};
        shanks = {[2:3],[0:3],[0:3]};
        probe = '128DN_bottom';
        lightconds = {2,1,1};
        lightcond = 1;
    elseif strcmpi(exp_type,'L6_halo') && strcmpi(area,'LGN')
        exp_paths = {'Z:\L5&L6\MD9\LPl\2020-11-25_16-16-51_driftgrat_2LEDs'};
        shanks = {3};
        probe = '128DN_bottom';
        lightconds = {1};
        lightcond = 1;
    end

   elseif strcmpi(pop,'V1')

    if strcmpi(exp_type,'step_V1inact') && strcmpi(area,'LPlateral')
%             exp_paths = {'J:\LPproject\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
%             shanks = {[2]};
        exp_paths = {'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
                'Y:\V1Inactivation\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2'};
        shanks = {[1], [0 1],2};
        lightconds = {1,1,2};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_V1inact') && strcmpi(area,'LP')   % include LPrm
       exp_paths = {'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Z:\V1Inactivation\VSC4\LPrm\2020-06-04_20-16-50_HaloChR2'};
        shanks = {[1], [0 1],[0 1]};
        lightconds = {1,1,1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_cortexinact') && strcmpi(area,'LPlateral')
        exp_paths = {'Y:\V1Inactivation\VSC5\LPl\2020-07-28_14-55-27_driftgrat_2LEDs',...
            'Y:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs'};
        shanks = {[0:1], [0:1]};
        lightconds = {1,1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_SCcortexinact') && strcmpi(area,'LPlateral')
        exp_paths = {'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Y:\V1Inactivation\VSC5\LPl\2020-07-28_14-55-27_driftgrat_2LEDs',...
            'Y:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs'};
        shanks = { [1],[0 1],[0 1],[0 1]};
        lightconds = {1:3,1:3,1:3,1:3};
        lightcond = 3;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_cortexinact') && strcmpi(area,'LPmedial')
        exp_paths = {'Z:\V1Inactivation\VES3\LP\2020-05-07_15-47-43_ChR2_2LEDs'};
        shanks = {[0:3]};
        lightconds = {1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_SCV1inact') && strcmpi(area,'LPlateral')
        exp_paths = {'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Y:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs'};
        shanks = {[1], [0 1],[0 1]};
        lightconds = {1:3, 1:3,1:3};
        lightcond = 3;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_V1inact') && strcmpi(area,'LPmedial')
      exp_paths = {'Z:\V1Inactivation\VSC4\LPrm\2020-06-04_20-16-50_HaloChR2',...
          'Z:\V1Inactivation\VH2\2018-05-16_19-05-31_LP_PVChR2'};
        shanks = {[0:1],[1]};
        lightconds = {1,2};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_cortex2LEDs') && strcmpi(area,'LPmedial')
        exp_paths = {'Z:\V1Inactivation\VES3\LP\2020-05-07_15-47-43_ChR2_2LEDs'};
        shanks = {[0:3]};
        lightconds = {1:3};
        lightcond = 3;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'all') && strcmpi(area,'LP')
        exp_paths = {'Z:\V1Inactivation\VH3\LP\2018-05-17_15-16-48_VH3_LP_PVChR2',...
            'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Z:\V1Inactivation\VSC5\LPl\2020-07-28_14-55-27_driftgrat_2LEDs',...
            'Z:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs',...
            'Z:\V1Inactivation\VSC4\LPrm\2020-06-04_20-16-50_HaloChR2',...
          'Z:\V1Inactivation\VH2\2018-05-16_19-05-31_LP_PVChR2',...
          'Z:\V1Inactivation\VES3\LP\2020-05-07_15-47-43_ChR2_2LEDs'};
        shanks = {[2],[1],[0 1],[0 1],[0 1],[0:1],[1],[0:3]};
        lightconds = {2,1,1,1,1,1,2,1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(area,'LGN')   
       exp_paths = {'Z:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Z:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Z:\V1Inactivation\VSC5\LPl\2020-07-28_14-55-27_driftgrat_2LEDs',...
            'Z:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs'};
        shanks = {[2:3], [2:3],[3],[2 3]};
        lightconds = {1,1,1,1};
        lightcond = 1;
        probe = '128DN_bottom';
    end
elseif strcmpi(pop,'SC')
    if strcmpi(exp_type,'step_halo') && strcmpi(area,'LPlateral')
        exp_paths = {'Y:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Y:\V1Inactivation\VSC4\LPl\2020-06-04_17-08-19_ChR2Halo',...
            'Y:\V1Inactivation\VSC6\LPl\2020-07-30_12-39-22_drftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC4\LPl\2020-07-21_15-45-08_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC8\LPl\2020-08-31_16-14-35_driftgrat_2LEDs'};
        shanks = {[1], [0 1], [1],[0 1],[0:2],[1:2],[0:2]};
        lightconds = {2,2,2,2,2,1,1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'LPcaudal')
        exp_paths = {'Y:\V1Inactivation\VSC3\LPl\2020-06-04_12-22-21_HaloChR2',...
            'Y:\L5\stGtACR\NPRSC5\LPl\2020-07-22_13-49-52_driftgrat_2LEDs',...
            'Y:\L5\stGtACR\NPRSC6\LPl\2020-09-02_14-37-11_driftgrat_2LEDs'};
        shanks = {0,0:3,0};
        lightconds = {2,2,1};
        lightcond = 1;
        probe = '128DN_bottom';
    elseif strcmpi(exp_type,'step_halo') && strcmpi(area,'SC')
        exp_paths = {'Z:\L5\stGtACR\NPRSC5\SC\2020-07-22_17-55-27_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC6\SC\2020-09-02_19-08-32_driftgrat_2LEDs',...
            'Z:\L5\stGtACR\NPRSC8\SC\2020-08-31_20-19-41_driftgrat_2LEDs',...
            'Z:\V1Inactivation\VSC6\SC\2020-07-30_12-39-22_drftgrat_2LEDs',...
            'Z:\V1Inactivation\VSC7\SC\2020-07-29_13-24-00_driftgrat_2LEDs'};
        shanks = {0,0,0,0,0};
%         lightconds = {2,1,1,2,2};
        lightconds = {[2 1 3],1:3, 1:3, [2 1 3], [2 1 3]}; 
        lightcond = 1;
        probe = '64D_bottom';
    end
elseif strcmpi(pop,'drivermod')
    if strcmpi(exp_type,'gtacr_halo') && strcmpi(area,'LP')
        exp_paths = {'Z:\L5&L6\MD5\LPl\2020-11-07_15-25-10_driftgrat_2LEDS_REALREAL',...
            'Z:\L5&L6\MD5\LPrm\2020-11-07_19-45-43_driftgrat_2LEDs',...
            'Z:\L5&L6\MD7\LPl\2020-11-24_14-24-32_driftgrat_2LEDs',...
            'Z:\L5&L6\MD10\LPl\2020-12-21_15-56-41_driftgrat_2LEDsREAL'};
        shanks = {0:1,0:3,0:1,0:1};
        lightconds = {[2 1 3],[2 1 3],[2 1 3],[2 1 3]};  % reordered to be L5, L6, L5&6
        lightcond = 3;
        probe = '128DN_bottom';
    end
    
end

if ~exist(main_dir,'dir')
    mkdir(main_dir)
end
cd(main_dir)
fig_dir =  sprintf('%s\\%s\\%s_%s_%s_figs',main_dir,'InactFigures',area,pop,exp_type);
if strcmpi(vissigdef,'prepost');  fig_dir =  sprintf('%s\\%s\\%s_%s_%s_%sfigs',main_dir,'InactFigures',area,pop,exp_type,'dTAcomp'); end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

%% load experiment data
params = [];
unitinfo = [];
FRs = [];
tuning = [];
waveforms = [];
refr_idx = [];
cR = [];
uQ = [];
refV = [];
exp_num = [];

for i = 1:length(exp_paths)
    exp_path = exp_paths{i};
    % get experiment name
    out = regexp(exp_path,'\\','split');
    an_name{i} = out{end-1};       % animal name
    if contains(an_name{i},'LP','ignorecase',1) || contains(an_name{i},'V1','ignorecase',1) || contains(an_name{i},'TRN','ignorecase',1) || contains(an_name{i},'SC','ignorecase',1)
        an_name{i} = strcat(out{end-2},'_',an_name{i});
    elseif contains(an_name{i},'day','ignorecase',1) 
        an_name{i} = out{end-2};
    end
    if ~isempty(strfind(out{end},'shank'))      % if exp_path is a shank subdirectory w/in experimental folder
        inds = regexp(out{end-1},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
        exp_name = strcat(out{end-1}(1:inds(1)-1),sprintf('_%s',out{end}));     % experiment name includes shank
    else
        inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
        exp_name = out{end}(1:inds(1)-1);
        if sum(isstrprop(exp_name,'digit'))>4       % this means it's an OpenEphys file (because it names files starting with the date)
            inds = regexp(out{end},'_\D','start');
            if strfind(out{end}(inds(1):end),an_name{i})
                exp_name = out{end}(inds(1)+1:end);
            else
                exp_name = strcat(an_name{i},out{end}(inds(1):end));
            end
        end
    end
    exp_dir = strcat(main_dir,'\',exp_name);
    if ~exist(exp_dir,'dir') && exist(sprintf('%s\\%s\\%s',main_dir,an_name{i},exp_name))       % need better solution for this
        exp_dir = sprintf('%s\\%s\\%s',main_dir,an_name{i},exp_name);
    end
    fprintf(sprintf('Loading experiment %s\n',exp_name))
    
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
        cluster_file = sprintf('%s\\%s\\%s_cluster_quality.mat',main_dir,an_name{i},an_name{i});
    end
    exp = importdata(results_file,'-mat');     % load results mat
    clust = importdata(cluster_file,'-mat');
    clear cluster_file results_file

    exp_num = [exp_num i*ones(1,length(exp.unitinfo))];
    params = [params exp.params];
    unitinfo = [unitinfo exp.unitinfo];
    FRs = [FRs exp.FRs];
    if isfield(exp.tuning,'SF')
        exp.tuning = rmfield(exp.tuning,'SF');      % remove SF and TF fields, if present (for now...just so I can concatenate results from prior experiments)
        exp.tuning = rmfield(exp.tuning,'TF');
    end
    tuning = [tuning exp.tuning];
    waveforms = [waveforms exp.waveforms];
    uQ = [uQ clust.uQ'];
    cR = [cR clust.cR'];
    refV = [refV clust.refV];
    
end

%% Identify experiment parameters for multi-experiment analyses

% first, set up analysis window (esp. important for halo experiments where light
% starts before vis stim)
for n = 1:length(params)
    if ~isempty(lightconds{1})      % opto experiments
        if round(max(params(n).av_light_start(lightconds{n})),2) < params(n).prestim % if LED started BEFORE vis stim
            startT(n) = 1500*params(n).prestim+1;        % manually set start time to prestim amount of time after vis stim start (e.g.startT = 751, 250ms after vis stim start)
        else
            startT(n) = round(max(params(n).av_light_start(lightconds{n}))*1000);
        end
        endT(n) = round(min(params(n).av_light_start(lightconds{n})+params(n).light_dur(lightconds{n}+1))*1000);
    else            % non-opto experiments
        startT(n) = 1500*params(n).prestim+1;        % manually set start time to prestim amount of time after vis stim start (e.g.startT = 1001, 500ms after vis stim start)
        endT(n) = 4500*params(n).prestim+1;  % manually set end time to 2250
    end
end
window = [max(startT) min(endT)]; % use common time window across experiments (in ms)
ev_lighttime = (diff(window)+1)/1000; % in sec

% identify different trial types
for n = 1:length(params)
    if strcmp(params(n).IVs,'s_freq')
        vis_trials{n} = find((params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1)&(params(n).trial_type(:,strcmpi(params(n).IVs,'s_freq'))==.04)&(params(n).trial_type(:,strcmpi(params(n).IVs,'t_period'))==30));  % in case of multiple SFs and TFs, only compare regular 2Hz and .04cpd trials (for now)
    else
        vis_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==1);
    end
    blank_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'visual'))==0);
    nolight_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==0);
    conds{n} = unique(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit')));  % the actual values assigned to light conditions 
    for lc=1:length(lightconds{n})  
        light_trials{n}{lc} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1));
        vislight_trials{n}{lc} = intersect(vis_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1)));
        blanklight_trials{n}{lc} = intersect(blank_trials{n},find(params(n).trial_type(:,strcmpi(params(n).IVs,'light_bit'))==conds{n}(lightconds{n}(lc)+1)));
    end
    run_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==1);
    stat_trials{n} = find(params(n).trial_type(:,strcmpi(params(n).IVs,'running'))==0);
    fprintf('%s: %d total trials, %6.2f%% stationary \n',params(n).exp_name,size(params(n).trial_type,1),(length(stat_trials{n})/size(params(n).trial_type,1))*100)
end

%% identify which units should be included for analyses, pt. 1 (based on shank location and cluster quality)
incl_units = zeros(1,length(unitinfo)); 
if isfield(waveforms,'shank')
    shk = [waveforms.shank];
    if isempty(shk); shk = zeros(1,length(unitinfo)); end    % if single shank probe, shank field may be empty
else
    shk = zeros(1,length(unitinfo));
end

for n = 1:length(params)
    which_units = find(exp_num==n);
    incl_units(which_units(ismember(shk(exp_num==n),shanks{n}))) = 1;  % only include units from designated shank(s)
end
clean_units = find((uQ>16)&(refV<.5)&(incl_units));    % only include units with <.5% refractory period violations, unit quality > 16, and those from designated shanks 
incl_units = find(incl_units);      % necessary for later...use incl_units instead of clean_units to define LP boundaries

% also exclude low-FR cells (<.25Hz in blank and visual trials)
FRb = [FRs(incl_units).baseline];  %  baseline (from blank trials w/ no running, light or vis stim)
FRv = arrayfun(@(x) x.visual.ev(1,1),FRs(incl_units)); % visually-evoked FR from visual, stationary, no light trials
quiet_cells = max([FRv; FRb])<.25;
incl_units(quiet_cells) = [];
FRb(quiet_cells) = [];
FRv(quiet_cells) = [];

%% Calculate visual significance prior to getting distances from first
% and last ch and deciding which units are "clean". Use first and last
% visually significant channels to determine borders of LP
if contains(area,'LP','ignorecase',1)
    onsetT=.2;  % because LP can have more delayed onset vis responses
else
    onsetT=.1;
end

for i = 1:length(incl_units)      % for each unit
    nn = incl_units(i);
    tuning_curve{i} = tuning(nn).curve(:,:);
    oris = unique(params(exp_num(nn)).trial_type(:,strcmp(params(exp_num(nn)).IVs,'ori')));
        oris(oris>=999) = [];
    % wilcoxen rank-sum test to test for significant and visual- and light-modulation
    % for visual modulation, find preferred direction trials - only use THESE
    % trials to test for significant visual modulation (in case of extremely
    % tuned cells)
    [~,prefdir_deg(i)] = max(abs(tuning_curve{i}(1,:)-FRb(i)));   % direction w/ biggest change from baseline (defined from no-light condition)
    prefori_trials{i} = find(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori'))==oris(prefdir_deg(i)));
    if contains(exp_type,'dta','IgnoreCase',1) || strcmpi(vissigdef,'prepost') % if dta exp without blank trials, compare pre- and post-FRs (OR if i specifically indicate I wanna look at visual responsivity this way)
        vis_sig(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),window(1):window(1)+1000*params(exp_num(nn)).prestim-1),2),... % significance of visual response (evoked periods of 1000ms in blank vs. visual trials)
            sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1:1000*params(exp_num(nn)).prestim),2));
        vis_sig_ons(i) = signrank(sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2),...  % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials)
            sum(unitinfo(nn).rast(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),1:onsetT*1000),2));
    else % if comparing visual to blank trials, use UNPAIRED parametric text (ie wilcoxen ranksum)
        vis_sig(i) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),... % significance of visual response (evoked periods of 1000ms in blank vs. visual trials) % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2));
        vis_sig_ons(i) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2),...  % significance of visual response (100ms immediately following visual stimulus onset in blank vs. visual trials) % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),1000*params(exp_num(nn)).prestim+1:1000*(params(exp_num(nn)).prestim+onsetT)),2));   % changed from params(exp_num(nn)).onset to .2s 6/8/20 because 100ms may be too fast for vis-evoked response in LP
    end
    
%     if length(lightconds{exp_num(nn)})<2  % if only one light cond, can do ranksum test
    for lc = 1:length(lightconds{exp_num(nn)})  % changed 1/7/21 to specifically use preferred ori trials and stationary trials for assessing sig light effect on visual responses
        light_sig(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),...       % NEW 6/25/20 - separately assessing light significant for visual vs. blank trials  % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2));% significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
        light_sig_bl(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2),...       % NEW 6/25/20 - separately assessing light significant for visual vs. blank trials  % Stationary trials ONLY (new addition 1/7/21)
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2));% significance of light-evoked response (evoked periods of 500ms in nolight vs. light trials (each condition separately))
        if params(exp_num(nn)).av_light_start(lc) < params(exp_num(nn)).prestim    % if light started prestim, you can use visual+blank trials to test significance of light-onset response
            light_sig_ons(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(nolight_trials{exp_num(nn)},stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2),...
            sum(unitinfo(nn).rast(intersect(light_trials{exp_num(nn)}{lc},stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2)); % significance of light-evoked response at light onset (200ms after light onset in nolight vs light trials  - changed from params(exp_num(nn)).onset to .2s 12/28/20 % Stationary trials ONLY (new addition 1/7/21)
        else % otherwise, just use blank trials
            light_sig_ons(i,lc) = ranksum(sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2),...
            sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),round(1000*(params(exp_num(nn)).av_light_start(lc)))+1:round(1000*(params(exp_num(nn)).av_light_start(lc)+.2))),2)); % significance of light-evoked response at light onset (200ms after light onset in nolight vs light trials  - changed from params(exp_num(nn)).onset to .2s 12/28/20 % Stationary trials ONLY (new addition 1/7/21)
        end
    end
    if lightcond > 1 && (contains(exp_type,'SC','ignorecase',1) || contains(exp_type,'2LEDs','ignorecase',1))  % 
        rast{1} = sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
        tri_group{1} = zeros(length(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)})),1)';
        rast_bl{1} = sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
        tri_group_bl{1} = zeros(length(intersect(intersect(blank_trials{exp_num(nn)},nolight_trials{exp_num(nn)}),stat_trials{exp_num(nn)})),1)';
        for lc = 1:length(lightconds{exp_num(nn)})
            rast{lc+1} = sum(unitinfo(nn).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
            tri_group{lc+1} = lc*ones(length(intersect(intersect(prefori_trials{i},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)})),1)';
            rast_bl{lc+1} = sum(unitinfo(nn).rast(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)}),window(1):window(2)),2)';
            tri_group_bl{lc+1} = lc*ones(length(intersect(intersect(blank_trials{exp_num(nn)},light_trials{exp_num(nn)}{lc}),stat_trials{exp_num(nn)})),1)';
        end
        [light_sig(i,lightcond),~,unit_stats(i)] = kruskalwallis([rast{:}],[tri_group{:}],'off');
        mc(:,:,i) = multcompare(unit_stats(i),'ctype','dunn-sidak','display','off');
        [Q(:,:,i),vcrit(:,:,i)] = dunn([rast{:}],[tri_group{:}]+1,1);  % only comparisons to no light cond!
        [light_sig_bl(i,lightcond),~,unit_stats_bl(i)] = kruskalwallis([rast_bl{:}],[tri_group_bl{:}],'off');
        mc_bl(:,:,i) = multcompare(unit_stats(i),'ctype','dunn-sidak','display','off');
        [Q_bl(:,:,i),vcrit_bl(:,:,i)] = dunn([rast_bl{:}],[tri_group_bl{:}]+1,1); % only comparisons to no light cond!
%         light_sig_ons(i,1:length(lightconds{exp_num(nn)})) = nan(1,length(lightconds{exp_num(nn)})); % JUST A PLACEHOLDER  - not actually using for these experiments
    end
     
   % next, check tuning significance and get tuning curves
     if length(oris)>4        % DON'T calc tuning for experiments with only 4 (or fewer) orientations
        % evaluate significance or orientation tuning
        for lc = 1:length(lightconds{exp_num(nn)}) + 1  % FIXED 1/5/21 - was excluding no-light trials!
            ori_trials = cell(1, length(oris)/2);
            for o = 1:length(oris)/2
                if lc == 1
                    ori_trials{o} = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==conds{exp_num(nn)}(1)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in no-light trials
                else
                    ori_trials{o} = find((ismember(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'ori')),oris([o o+length(oris)/2])))&(params(exp_num(nn)).trial_type(:,strcmpi(params(exp_num(nn)).IVs,'light_bit'))==conds{exp_num(nn)}(lightconds{exp_num(nn)}(lc-1)+1)));   % trials at each orientation (**regardless of direction -seems to work better? captures more units) in particular light condition
                end
            end
            n_oritrials = cellfun(@(x) length(x),ori_trials);
            for o = 1:length(oris)/2
                if length(unique(n_oritrials)) > 1    % in cases of unequal # of trials per ori
                    ori_trials{o}(randperm(n_oritrials(o),n_oritrials(o)-min(n_oritrials))) = []; % if unequal # of trials per ori, randomly choose trial(s) to exclude
                end
                tuning_trials{exp_num(nn)}(:,o,lc) = sum(unitinfo(nn).rast(ori_trials{o},window(1):window(2)),2); % numbers of spikes across trials of given light condition for each orientation (by column)
            end
            tuned_sig(i,lc) = T2Hot1(tuning_trials{exp_num(nn)}(:,:,lc),0.05,zeros(1,length(oris)/2));     % if tuned_sig(lc) < .05, significantly tuned  (light trials, separately by condition)
        end
    end
end

%% Get depth information- use visual significance to approximate dorsal and ventral boundaries of LP
distfromlastch = nan(1,length(incl_units));
distfromfirstch = distfromlastch;
% [~, pthreshVis, ~, ~] = fdr_bh(vis_sig,.1); % correct for FDR (10%)
% [~, pthreshVisOns, ~, ~] = fdr_bh(vis_sig_ons,.1);
% vis_units = find(vis_sig<pthreshVis|vis_sig_ons<pthreshVisOns);      % benjamin-hochburg FDR correction
vis_units = find(vis_sig<.025|vis_sig_ons<.025); % alt: .025 (b/c bonferroni correction for 2 tests) instead of 

% get probe info 
p = eval(sprintf('probemap_%s_func',probe));
Zchan = flipud(sort(p.z(p.shaft==1))); % from top to bottom (bottom=0)
count=1;
for n = 1:length(params)        % for each exp
    for sh = 1:length(shanks{n})    % for each shank in exp
        shk_units{count} = find((exp_num(incl_units)==n) & (shk(incl_units)==shanks{n}(sh)));
        vischs = sort([waveforms(incl_units(intersect(shk_units{count},vis_units))).max_ch]);
        if strcmpi(area,'SC')
            vischs = sort([waveforms(incl_units(shk_units{count})).max_ch]);  % for SC recordings, include ALL units
        end
        firstch(count) = min(vischs);
        if ~strcmpi(area,'SC') && vischs(2)-min(vischs)>5      % added 6/9/20 in case a hippocampal unit is included (MAK)
            firstch(count) = vischs(2);
        end
        if firstch > 1
            if Zchan(firstch(count))==Zchan(firstch(count)-1) % for probes in hexagonal orientation, might leave out channel that is actually same height as "firstch"
                firstch(count) = firstch(count)-1;
            end
        end
        lastch(count) = max(vischs);
        if ~strcmpi(area,'SC') && lastch(count)-vischs(end-1)>5      % added 6/9/20 in too-deep unit is included (MAK)
            lastch(count) = vischs(end-1);
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

% adjust for new clean_units
distfromlastch = distfromlastch(ismember(incl_units,clean_units));          % only includes GOOD units
distfromfirstch = distfromfirstch(ismember(incl_units,clean_units));          % only includes GOOD units
tuning_curve=tuning_curve(ismember(incl_units,clean_units));
prefdir_deg = prefdir_deg(ismember(incl_units,clean_units));
prefori_trials = prefori_trials(ismember(incl_units,clean_units));
vis_sig = vis_sig(ismember(incl_units,clean_units));
vis_sig_ons = vis_sig_ons(ismember(incl_units,clean_units));
if exist('light_sig','var')  % if opto
    light_sig = light_sig(ismember(incl_units,clean_units),:);
    light_sig_bl = light_sig_bl(ismember(incl_units,clean_units),:);
    light_sig_ons = light_sig_ons(ismember(incl_units,clean_units),:);
    if exist('mc','var')
        mc = mc(:,:,ismember(incl_units,clean_units));
        mc_bl = mc_bl(:,:,ismember(incl_units,clean_units));
        Q = Q(:,:,ismember(incl_units,clean_units));
        Q_bl = Q_bl(:,:,ismember(incl_units,clean_units));
        vcrit = vcrit(:,:,ismember(incl_units,clean_units));
        vcrit_bl = vcrit_bl(:,:,ismember(incl_units,clean_units));
    end
end
tuned_sig = tuned_sig(ismember(incl_units,clean_units),:);
FRb = FRb(ismember(incl_units,clean_units));
exp_num = exp_num(clean_units);

%% test for light modulation

% first, need to load important FR, tuning and psth data
num_lcs = length(lightconds{1});  % currently, must choose to look at same number of lightconds from every experiment
FRev = nan(length(clean_units),num_lcs+1); 
FRearly = FRev;
FRlate = FRev;
FRonset = FRev;
FRbl = FRev;
FRvison = FRev;
OSI_CV = FRev;
OSI = FRev;
DSI_CV = FRev;
DSI = FRev;
psthV = nan(num_lcs+1,size(FRs(1).psthVisual,2),length(clean_units)); % #conds x #timebins x #units 
psthBl = psthV;
renormcurve = nan(num_lcs+1,size(tuning(1).curve,2),length(clean_units));
for n = 1:length(params)
    FRev(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.ev(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';   %vis-evoked
    FRearly(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evstart(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';   % early evoked period
    FRlate(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evlate(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % late evoked period
    FRonset(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.visual.evlightonset(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % vis-evoked light onset
    FRbl(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.blank.ev(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % blank, evoked period
    FRvison(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.onset(1,[1 lightconds{n}+1]),FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';  % vis stim onset
    % also get orientation values for later use
    OSI_CV(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.OSI_CV(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    OSI(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.OSI(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    DSI_CV(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.DSI_CV(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    DSI(exp_num==n,:) = reshape(cell2mat(arrayfun(@(x) x.DSI(1,[1 lightconds{n}+1]),tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,sum(exp_num==n))';
    renormcurve(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.curve([1 lightconds{n}+1],:), tuning(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(tuning(1).normcurve,2),sum(exp_num==n));
    % and get PSTHs for later use
    psthV(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.psthVisual([1 lightconds{n}+1],:), FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(FRs(1).psthVisual,2),sum(exp_num==n));
    psthBl(:,:,exp_num==n) = reshape(cell2mat(arrayfun(@(x) x.psthBlank([1 lightconds{n}+1],:), FRs(clean_units(exp_num==n)),'uniformoutput',0)),length(lightconds{n})+1,size(FRs(1).psthBlank,2),sum(exp_num==n));
end

% recalculate tuning using normalized tuning curves (additively or
% multiplicatively scaled to account for baseline FR change due to light)
renormcurve2 = renormcurve;     % multiplicatively rather than additively scaling
for i=2:size(renormcurve,1) % for each lightcond
    renormcurve(i,:,:) = squeeze(renormcurve(i,:,:)) - repmat(diff(FRev(:,[1 i]),[],2)',size(renormcurve,2),1); % using FRev instead of FRbl since halo is more effective in blank than visual trials?
    scalef = FRev(:,1)./FRev(:,i);
    renormcurve2(i,:,:) = squeeze(renormcurve2(i,:,:)).* repmat(scalef',size(renormcurve2,2),1);
end
for n=1:length(clean_units)
    for lc = 1:size(renormcurve,1)
        [newOSI(n,lc),newOSI_CV(n,lc),newDSI(n,lc),newDSI_CV(n,lc)] = calcOSIDSI(squeeze(renormcurve(lc,:,n)),oris');       % baseline-subtracted (baseline from each light condition separately)
%         if logical(exist('scalef','var')) && ~isinf(scalef) % in case of rare cases in which FRev==0 so scalef was infinite
            [newOSI2(n,lc),newOSI_CV2(n,lc),newDSI2(n,lc),newDSI_CV2(n,lc)] = calcOSIDSI(squeeze(renormcurve2(lc,:,n)),oris');       % normalized to pref deg (i.e. multiplicative scaling) -  **same as orig OSI vals - OSI calcs not affected by scaling!
%         end
    end
end

% calculate light modulation indexes
for ii = 2:size(FRev,2)
    lightmod(:,ii-1) = (diff(FRev(:,[1 ii]),[],2))./sum(FRev(:,[1 ii]),2);
    lightmod_early(:,ii-1) = (diff(FRearly(:,[1 ii]),[],2))./sum(FRearly(:,[1 ii]),2);
    lightmod_late(:,ii-1) = (diff(FRlate(:,[1 ii]),[],2))./sum(FRlate(:,[1 ii]),2);
    lightmod_onset(:,ii-1) = (diff(FRonset(:,[1 ii]),[],2))./sum(FRonset(:,[1 ii]),2);
    lightmod_bl(:,ii-1) = (diff(FRbl(:,[1 ii]),[],2))./sum(FRbl(:,[1 ii]),2);
end

%% get preferred stimulus FR for each condition (but preferred stimulus defined in no light condition) and prestim FR
FRpref = nan(size(FRev));
FRlightons = FRpref;
if exist('light_sig','var')
    lightstart = [1000*round(max(params(1).av_light_start(lightconds{1})),2)+1 1000*round(max(params(1).av_light_start(lightconds{1})),2)+200]; % first 200ms of light - CAUTION - currently assumes same light start and vis start for all exps
end
for i = 1:length(clean_units)
    for ii = 1:num_lcs+1
        if ii == 1
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},nolight_trials{exp_num(i)}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime; % only stationary trials (as for other FRs)
            if exist('light_sig','var')
                if lightstart(1) < 1000*params(1).prestim  % if light started before vis, can use blank+vis trials
                    FRlightons(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(nolight_trials{exp_num(i)},stat_trials{exp_num(i)}),lightstart(1):lightstart(2)),2))/((diff(lightstart)+1)/1000); % only stationary trials (as for other FRs), but including visual AND black trials
                else  % otherwise, just blanks
                    FRlightons(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},nolight_trials{exp_num(i)}),stat_trials{exp_num(i)}),lightstart(1):lightstart(2)),2))/((diff(lightstart)+1)/1000); % only stationary trials (as for other FRs), but including visual AND black trials
                end
            end
        else
            FRpref(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(prefori_trials{i},light_trials{exp_num(i)}{ii-1}),stat_trials{exp_num(i)}),window(1):window(2)),2))/ev_lighttime;  % only stationary trials (as for other FRs)
            if lightstart(1) < 1000*params(1).prestim  % if light started before vis, can use blank+vis trials
                FRlightons(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(light_trials{exp_num(i)}{ii-1},stat_trials{exp_num(i)}),lightstart(1):lightstart(2)),2))/((diff(lightstart)+1)/1000); % only stationary trials (as for other FRs), but including visual AND black trials
            else  % otherwise, just blanks
                FRlightons(i,ii) = mean(sum(unitinfo(clean_units(i)).rast(intersect(intersect(blank_trials{exp_num(i)},light_trials{exp_num(i)}{ii-1}),stat_trials{exp_num(i)}),lightstart(1):lightstart(2)),2))/((diff(lightstart)+1)/1000); % only stationary trials (as for other FRs), but including visual AND black trials
            end
            lightmod_lightons(:,ii-1) = (diff(FRlightons(:,[1 ii]),[],2))./sum(FRlightons(:,[1 ii]),2);
        end
    end
end

% and visual modulation (changed 1/28/21) - calculated from preferred vis
% trials (but onset is still from ALL vis trials)
% vismod = (FRpref(:,1) - FRb')./(FRpref(:,1) + FRb');        % using baseline (from blank trials w/ no running, light or vis stim)
vismod_on = (FRvison(:,1)-FRb')./(FRvison(:,1)+FRb');
vismod = reshape((FRpref(:)-FRbl(:))./(FRpref(:)+FRbl(:)),length(clean_units),num_lcs+1);

%% F1/F0 response analysis

binsize = .01;  %% currently hardcoded!!! 10ms bins (prev 25ms - changed 1/18/21)
% bs = squeeze(mean(psthV(:,1:200/(binsize*1000),:),2));  % baseline (mean across first 200ms (8*25ms binsize = 200ms), for each condition separately) **NOT using blanks in case of light effects on baseline FRs
Fratio = nan(length(clean_units),num_lcs+1);
zF1 = Fratio;
F1 = Fratio;
for i = 1:length(clean_units)
    vis_start = params(exp_num(i)).prestim;       % in sec
    vis_end = vis_start + params(exp_num(i)).stimtime;      % in sec
    if strcmpi(vissigdef,'prepost') % dta experiment (added 7/11/21)
        bs(:,i) = FRb(i);
    else
        bs(:,i) = mean(psthBl(:,vis_start/.025+11:round(window(2)/1000,2)/.025,i),2);      % get unit's baseline FR for SAME time window in blank trials of each light cond separately (so getting rid of baseline changes in FR due to light before calculating Fratios)
    end
    
     % only for preferred visual stim trials
    which_trials = ismember(1:size(unitinfo(clean_units(i)).rast,1),intersect(stat_trials{exp_num(i)},prefori_trials{i})); % preferred ori and stationary trials
    all_light = zeros(1,length(which_trials));
    if exist('light_sig','var')
        for ii = 1:num_lcs
            all_light(light_trials{exp_num(i)}{ii}) = ii; %light trials
        end
    end
    [~,pref_psthV] = make_psth_v2(binsize,0:binsize:size(unitinfo(clean_units(i)).rast,2)/1000,which_trials,unitinfo(clean_units(i)).rast,all_light);

    % **zF1 is unaffected by baseline subtraction, whereas baseline
    % subtraction is necessary for getting Fratios > 1
%     [Fratio(i,:), zF1(i,:)] = calc_F1F0(psthV(:,vis_start/binsize+21:vis_end/binsize,i)',binsize,2);        % **currently hardcoded for 2Hz tfreq - will need to change!! currently: +21 instead of +1 excludes first 20bins (20*.025binsize = .5sec excluded), i.e. excludes onset response)
%     [Fratio_bs(i,:), zF1_bs(i,:)] = calc_F1F0([psthV(:,vis_start/binsize+21:vis_end/binsize,i)-bs(:,i)]',binsize,2); % baseline-subtracted 

%     [Fratio(i,:), zF1(i,:), F1(i,:)] = calc_F1F0(psthV(:,vis_start/binsize+11:round(window(2)/1000,2)/binsize,i)',binsize,2);        % **currently hardcoded for 2Hz tfreq - will need to change!! currently: +11 instead of +1 excludes first 10bins (10*.025binsize = .25sec excluded), i.e. excludes onset response) and ends with the light offset period
%     [Fratio_bs(i,:), zF1_bs(i,:), F1(i,:)] = calc_F1F0([psthV(:,vis_start/binsize+11:round(window(2)/1000,2)/binsize,i)-bs(:,i)]',binsize,2); % baseline-subtracted 
    [Fratio_bs(i,:), zF1(i,:), F1(i,:)] = calc_F1F0([pref_psthV(:,vis_start/binsize+26:round(window(2)/1000,2)/binsize)-bs(:,i)]',binsize,2); % baseline-subtracted (currently using window starting 250ms after vis onset)
end

%% burst vs tonic firing stuff
burstrate = zeros(length(clean_units),1);
burstrate_light = zeros(length(clean_units),num_lcs);
burstnum = burstrate;
burstnumLight = burstrate_light;
burstresprate = burstrate;
burstresprate_light = burstrate_light;
if ~strcmpi(area,'SC')  % don't do for SC recordings
    for i = 1:length(clean_units)
        visnolight_trials = intersect(stat_trials{exp_num(i)},intersect(nolight_trials{exp_num(i)},vis_trials{exp_num(i)}));  % visual and STATIONARY trials
    %     visRast = unitinfo(clean_units(i)).rast(visnolight_trials,window(1):window(2));
        [burstrate(i),burstnum(i),burstresprate(i),avtonicrate(i),avburstrate(i)] = calc_bursting(unitinfo(clean_units(i)).rast(visnolight_trials,:),window(1):window(2),0);

        for ii = 1:num_lcs
            vislight_trials = intersect(stat_trials{exp_num(i)},intersect(light_trials{exp_num(i)}{ii},vis_trials{exp_num(i)})); % visual and STATIONARY trials
    %         vislightRast = unitinfo(clean_units(i)).rast(vislight_trials,window(1):window(2));
            [burstrate_light(i,ii),burstnumLight(i,ii),burstresprate_light(i,ii),avtonicrate_light(i,ii),avburstrate_light(i,ii)] = calc_bursting(unitinfo(clean_units(i)).rast(vislight_trials,:),window(1):window(2),0);
        end
    end
end

%% get waveform props
for i = 1:length(clean_units)
    [t2p_t(i),t2p_r(i),fwhm(i)] = get_waveform_props(waveforms(clean_units(i)).microV,params(exp_num(i)).amp_sr);
end

%% Figure time!
cd(fig_dir)

% define different cell types
% visual_cells = find((vis_sig < pthreshVis)|(vis_sig_ons < pthreshVisOns));  % use previously defined FDR-corrected pthresholds
visual_cells = find((vis_sig < .025)|(vis_sig_ons < .025)); 
nonvisual_cells = find(~ismember(1:length(vis_sig),visual_cells));

if exist('light_sig','var')
    % % use benjamini-hochburg correction for FDR to determine
    % % light-responsive cells
    % % defining light responsivity from BLANK trials to avoid more circular statistics
    [~, pthreshBl, ~, ~] = fdr_bh(light_sig_bl(:,lightcond),.1); % 10% FDR
    % pthreshBl = .01;        % alt: use set, strict cut-off of .01
    pthreshBl = max(pthreshBl,.01); % don't allow pthreshBl to be <.01
    [~, pthresh, ~, ~] = fdr_bh(light_sig(:,lightcond),.1); % 10% FDR
    pthresh = max(pthresh,.01); % don't allow pthreshBl to be <.01
    fprintf('new corrected pvalue for 10percent FDR - blank lightsig: %6.4f \n', pthreshBl)
    light_cells = find(light_sig_bl(:,lightcond)<=pthreshBl);  % currently looking across all included lightconds (might want to change??)
    supp_cells = find((light_sig_bl(:,lightcond)<=pthreshBl)&(lightmod_bl(:,lightcond)<0));
    enh_cells = find((light_sig_bl(:,lightcond)<=pthreshBl)&(lightmod_bl(:,lightcond)>0));
    other_cells = find(light_sig_bl(:,lightcond)>pthreshBl);
    
    % also look at significance of light at light onset (pre-vis stim) - esp.
    % for L6 experiments
    [~, pthreshOns, ~, ~] = fdr_bh(light_sig_ons(:,lightcond),.1); % 10% FDR
    pthreshOns = max(pthreshOns,.01); % don't allow pthresh to be <.01
    supp_ons_cells = find((light_sig_ons(:,lightcond)<=pthreshOns)&(lightmod_lightons(:,lightcond)<0));  
else
    supp_cells = [];
    supp_ons_cells = [];
    enh_cells = [];
    light_cells = [];
end

% look at units that were suppressed under different conditions (defined
% from BLANK trials)
if contains(exp_type,'SC','ignorecase',1) || contains(exp_type,'2LEDs','ignorecase',1)
    light_cells = find(light_sig_bl(:,lightcond)<.05);  % currently looking across all included lightconds (might want to change??)
    supp_cells = find(((light_sig_bl(:,lightcond)<.025)&(lightmod_bl(:,lightcond)<0)) | ((light_sig(:,lightcond)<.025)&(lightmod(:,lightcond)<0)));  
    enh_cells = find((light_sig_bl(:,lightcond)<.05)&(lightmod_bl(:,lightcond)>0));
    supp_pop = zeros(size(light_sig_bl,1),1);
%     supp_pop(supp_cells(light_sig_bl(supp_cells,2)<=.025 & lightmod_bl(supp_cells,2)<0)) = 2;  % SC modulated  % using .05 instead of pthreshbl because chances of false positives were already accounted for?
%     supp_pop(supp_cells(light_sig_bl(supp_cells,1)<=.025 & lightmod_bl(supp_cells,1)<0)) = 1;  % V1/L5 modulated
%     supp_pop(supp_cells(light_sig_bl(supp_cells,2)<=.025& lightmod_bl(supp_cells,2)<0 &light_sig_bl(supp_cells,1)<=.025& lightmod_bl(supp_cells,1)<0)) = 3;
%     supp_pop(supp_cells(mc_bl(2,end,supp_cells)<.05 & lightmod_bl(supp_cells,2)<0)) = 2;  % second row of mc output compares groups 1v3 (ie ctrl v SC)

    sig_mat = squeeze(Q(:,:,supp_cells)>vcrit(:,:,supp_cells))';
    sig_mat_bl = squeeze(Q_bl(:,:,supp_cells)>vcrit_bl(:,:,supp_cells))';
%     sig_mat = light_sig<.025;
%     sig_mat_bl = light_sig_bl<.025;
    mod_mat = lightmod(supp_cells,:)<0;
    mod_mat_bl = lightmod_bl(supp_cells,:)<0;
     supp_pop(supp_cells(sig_mat(:,1)+mod_mat(:,1)==2  | (sig_mat_bl(:,1)+mod_mat_bl(:,1)==2))) = supp_pop(supp_cells(sig_mat(:,1)+mod_mat(:,1)==2 | (sig_mat_bl(:,1)+mod_mat_bl(:,1)==2))) +1;
    supp_pop(supp_cells(sig_mat(:,2)+mod_mat(:,2)==2  | (sig_mat_bl(:,2)+mod_mat_bl(:,2)==2))) = supp_pop(supp_cells(sig_mat(:,2)+mod_mat(:,2)==2 | (sig_mat_bl(:,2)+mod_mat_bl(:,2)==2))) +2;
%     supp_pop(supp_cells((sig_mat(:,1)+mod_mat(:,1)==2 & sig_mat(:,2)+mod_mat(:,2)==2 & sig_mat(:,end)+mod_mat(:,end)==2) | (sig_mat_bl(:,1)+mod_mat_bl(:,1)==2 & sig_mat_bl(:,2)+mod_mat_bl(:,2)==2 & sig_mat_bl(:,end)+mod_mat_bl(:,end)==2))) = 3;
    supp_cells(supp_pop(supp_cells)==0) = []; % exclude cells that appear to have been false positives?
    save('supp_pop.mat','supp_pop')
end

tuned_cells = find(tuned_sig(:,1) < .05);

reg_cells = find(t2p_t>=.4);
FS_cells = find(t2p_t<.4);

%% define colors for plotting
% set up colors
% currently: colors according to population being inactivated (blue = L5,
% green = L6, SC = pink, V1 = teal
if strcmpi(pop,'drivermod')
    color_mat = [0 .2 .9; .166 .674 .188; .083 .567 .600]; % blue, green, teal (old blue: 0 .447 .741;)
elseif strcmpi(pop,'driver') 
    if contains(exp_type,'SC','ignorecase',1) || contains(exp_type,'2LEDs','ignorecase',1)  % if L5 and halo inactivation exps or experiments with multiple cortical areas
        color_mat = [0 .2 .9; 1 .6 1; .6 .52 .95]; % blue, pink, purple 
    else
        color_mat = [0 .2 .9];  % blue
    end
elseif strcmpi(pop,'modulator')
    color_mat = [.166 .674 .188]; % green
elseif strcmpi(pop,'V1')
    color_mat = [.083 .567 .6; 1 .6 1; .5 .52 .875]; % teal, pink, purple
elseif strcmpi(pop,'SC')
    color_mat = [1 .6 1;.083 .567 .6;.5 .52 .875]; % SC = pink
else % anything else (e.g., ctrl experiments)
    color_mat = [.5 .5 .5];  % grey
end
% color_mat = [0 0 0; .9 0 .3; 0.50, 0.0780, 0.10]; % for halo (red)
% nonvis_color_mat = [.5 .5 .5; 1 .75 .75; 0.7, 0.5, 0.5];      % make red for halo
% color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4]; % % for lighter-shade dots
% nonvis_color_mat = [.5 .5 .5; .75 .8 1; .7 .8 .7];  % for lighter-shade dots

% if contains(area,'lgn','ignorecase',1)
%     area_color = [.85 .325 .098];   % LGN=orange
% elseif contains(area,'LP','ignorecase',1)
%     area_color = [.494 .184 .556];   % LP = purple
% elseif contains(area,'TRN','ignorecase',1)
%     area_color = [.5 .5 .5];
% end

%% population-level paired tests of significance using Wilcoxen signed-rank (nonparametric)
for i=1:num_lcs
    if exist('light_sig','var')
        lightsig_all(i) = signrank(FRev(:,1),FRev(:,i+1));
        lightsig_all_dir(i) = sign(nanmedian(FRev(:,i+1)-(FRev(:,1)))); % 1 if light increased FR; -1 if it decreased FR (reflects sign of median of differences,
    % rather than sign of difference of medians
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
        lightsig_visFRdelt(i) = signrank(abs(FRev(:,1)-FRbl(:,1)),abs(FRev(:,i+1)-FRbl(:,i+1)));
        lightsig_visFRdelt_dir(i) = sign(nanmedian(abs(FRev(:,i+1)-FRbl(:,i+1))-abs(FRev(:,1)-FRbl(:,1))));
        if ~isempty(supp_cells) 
            lightsig_visFRdelt_supp(i) = signrank(abs(FRev(supp_cells,1)-FRbl(supp_cells,1)),abs(FRev(supp_cells,i+1)-FRbl(supp_cells,i+1)));
            lightsig_visFRdelt_supp_dir(i) = sign(nanmedian(abs(FRev(supp_cells,i+1)-FRbl(supp_cells,i+1))-abs(FRev(supp_cells,1)-FRbl(supp_cells,1))));
        end
%         lightsig_prefFRdelt(i) = signrank(FRpref(:,1)-FRbl(:,1),FRpref(:,i+1)-FRbl(:,i+1));
%         lightsig_prefFRdelt_dir(i) = sign(nanmedian(FRpref(:,i+1)-FRbl(:,i+1)-(FRpref(:,1)-FRbl(:,1))));
    end
    if ~isempty(supp_cells)     % changed from "new" to raw/regular 7/4/21
        osiCVsig_tuned(i) = signrank(OSI_CV(supp_cells,1),OSI_CV(supp_cells,i+1));
        osiCVsig_tuned_dir(i) = sign(nanmedian(OSI_CV(supp_cells,i+1)-(OSI_CV(supp_cells,1))));
        osisig_tuned(i) = signrank(OSI(supp_cells,1),OSI(supp_cells,i+1));
        osisig_tuned_dir(i) = sign(nanmedian(OSI(supp_cells,i+1)-(OSI(supp_cells,1))));
        dsiCVsig_tuned(i) = signrank(DSI_CV(supp_cells,1),DSI_CV(supp_cells,i+1));
        dsiCVsig_tuned_dir(i) = sign(nanmedian(DSI_CV(supp_cells,i+1)-(DSI_CV(supp_cells,1))));
        dsisig_tuned(i) = signrank(DSI(supp_cells,1),DSI(supp_cells,i+1));
        dsisig_tuned_dir(i) = sign(nanmedian(DSI(supp_cells,i+1)-(DSI(supp_cells,1))));
    end
    osiCVsig_vis(i) = signrank(OSI_CV(visual_cells,1),OSI_CV(visual_cells,i+1));
    osiCVsig_vis_dir(i) = sign(nanmedian(OSI_CV(visual_cells,i+1)-(OSI_CV(visual_cells,1))));
    osisig_vis(i) = signrank(OSI(visual_cells,1),OSI(visual_cells,i+1));
    osisig_vis_dir(i) = sign(nanmedian(OSI(visual_cells,i+1)-(OSI(visual_cells,1))));
    dsiCVsig_vis(i) = signrank(DSI_CV(visual_cells,1),DSI_CV(visual_cells,i+1));
    dsiCVsig_vis_dir(i) = sign(nanmedian(DSI_CV(visual_cells,i+1)-(DSI_CV(visual_cells,1))));
    dsisig_vis(i) = signrank(DSI(visual_cells,1),DSI(visual_cells,i+1));
    dsisig_vis_dir(i) = sign(nanmedian(DSI(visual_cells,i+1)-(DSI(visual_cells,1))));
    
    burstrate_sig(i) = signrank(burstrate(visual_cells),burstrate_light(visual_cells,i));
    burstrate_dir(i) = sign(nanmedian(burstrate_light(visual_cells,i)-burstrate(visual_cells)));
    Fratio_sig(i) = signrank(Fratio_bs(visual_cells,1),Fratio_bs(visual_cells,i+1));
    Fratio_dir(i) = sign(nanmedian(Fratio_bs(visual_cells,i+1)-(Fratio_bs(visual_cells,1))));
    zF1_sig(i) = signrank(zF1(visual_cells,1),zF1(visual_cells,i+1));
    zF1_dir(i) = sign(nanmedian(zF1(visual_cells,i+1)-(zF1(visual_cells,1))));
    F1_sig(i) = signrank(F1(visual_cells,1),F1(visual_cells,i+1));
    F1_dir(i) = sign(nanmedian(F1(visual_cells,i+1)-(F1(visual_cells,1))));
end

%% plot distance from top of LP by lightmod
% visual trials 
for i = 1:num_lcs
    figure;
    plot(lightmod(:,i),abs(distfromfirstch)','.','color',color_mat(i,:),'MarkerSize',24)
    view(0,270)
    xlim([-1 1])
    ylim([0 max(abs(distfromfirstch))])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    % legend('low','high')
    xlabel('LMI_v_i_s','Fontsize',24)
    ylabel(strcat('Depth (�m)'),'Fontsize',24)
    set(gca,'fontsize',18,'linewidth',2)
    print(gcf, '-dpng',sprintf('%s_%d','lightmodbydepth_top',i))
    print(gcf,'-painters','-depsc',sprintf('%s_%d','lightmodbydepth_top',i))

    % blank trials
    figure;
    subplot(111)
    plot(lightmod_bl(:,i),abs(distfromfirstch)','.','color',color_mat(i,:),'MarkerSize',24)
    view(0,270)
    xlim([-1 1])
    ylim([0 max(abs(distfromfirstch))])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    % legend('low','high')
    xlabel('LMI_s_p_o_n_t','Fontsize',24)
    ylabel(strcat('Depth (�m)'),'Fontsize',24)
    set(gca,'fontsize',18,'linewidth',2)
    print(gcf, '-dpng',sprintf('%s_%d','lightmodBlbydepth_top',i))
    print(gcf,'-painters','-depsc',sprintf('%s_%d','lightmodBlbydepth_top',i))

    %by type
    figure;
    subplot(111)
    scatter(lightmod(:,i),abs(distfromfirstch)',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',color_mat(i,:),'markerfacealpha',1)
%     scatter(lightmod(:,i),abs(distfromfirstch)',75,'filled','markerfacecolor',[.75 .75 .75],'markerfacealpha',1)
    hold on;
    scatter(lightmod(enh_cells,i),abs(distfromfirstch(enh_cells))',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',.5)
    scatter(lightmod(supp_cells,i),abs(distfromfirstch(supp_cells))',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',1)
    view(0,270)
    xlim([-1 1])
    ylim([0 max(abs(distfromfirstch))])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    xlabel('LMI_v_i_s','Fontsize',24)
    ylabel(strcat('Depth (�m)'),'Fontsize',24)
    set(gca,'fontsize',18,'linewidth',2)
    legend({'other','enhanced','suppressed'},'location','bestoutside')
    legend off
    print(gcf, '-dpng',sprintf('%s_%d','lightmodbydepth_top_bytype',i))
    print(gcf,'-painters','-depsc',sprintf('%s_%d','lightmodbydepth_top_bytype',i))

    % by type - blank trials
    figure;
    subplot(111)
    scatter(lightmod_bl(:,i),abs(distfromfirstch)',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',color_mat(i,:),'markerfacealpha',1)
%     scatter(lightmod_bl(:,i),abs(distfromfirstch)',75,'filled','markerfacecolor',[.75 .75 .75],'markerfacealpha',1)
    hold on;
    scatter(lightmod_bl(enh_cells,i),abs(distfromfirstch(enh_cells))',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',.5)
    scatter(lightmod_bl(supp_cells,i),abs(distfromfirstch(supp_cells))',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',1)
    view(0,270)
    xlim([-1 1])
    ylim([0 max(abs(distfromfirstch))])
    yax = get(gca,'YLim');
    line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
    xlabel('LMI_s_p_o_n_t','Fontsize',24)
    ylabel(strcat('Depth (�m)'),'Fontsize',24)
    set(gca,'fontsize',18,'linewidth',2)
    legend({'other','enhanced','suppressed'},'location','bestoutside')
    legend off
    print(gcf, '-dpng',sprintf('%s_%d','lightmodBlbydepth_top_bytype',i))
    print(gcf,'-painters','-depsc',sprintf('%s_%d','lightmodBlbydepth_top_bytype',i))
    
end

% if experiment included SC inactivation or 2LEDs, plot lightmod of V1/L5, SC/VisL, and
% both-activated units by depth
if contains(exp_type,'SC','ignorecase',1) || contains(exp_type,'2LEDs','ignorecase',1)
        figure;
        subplot(111)
        scatter(lightmod(:,end),abs(distfromfirstch)',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',[.5 .5 .5],'markerfacealpha',1)
        hold on;
        scatter(lightmod(supp_pop==1,end),abs(distfromfirstch(supp_pop==1))',75,'filled','markerfacecolor',color_mat(1,:),'markerfacealpha',1)
        scatter(lightmod(supp_pop==2,end),abs(distfromfirstch(supp_pop==2))',75,'filled','markerfacecolor',color_mat(2,:),'markerfacealpha',1)
        scatter(lightmod(supp_pop==3,end),abs(distfromfirstch(supp_pop==3))',75,'filled','markerfacecolor',color_mat(3,:),'markerfacealpha',1)
        
%         % trying to plot mean lightmods at each depth but hasn't really
%         % worked yet...
%         [dist_sort,inds_sort] = sort(abs(distfromfirstch(supp_pop>0)));
%         dist_range = min(dist_sort):min(diff(unique(dist_sort))):max(dist_sort);
%         dist_coarse = floor(dist_range/50);  % put into 50um bins
%         lightmod_sort = lightmod(supp_pop>0,end);
%         lightmod_sort = lightmod_sort(inds_sort);
%         supp_pop_sort = supp_pop(supp_pop>0);
%         supp_pop_sort = supp_pop_sort(inds_sort)';
%         mean_mod = zeros(3,length(unique(dist_coarse)));
%         for ii = 1:length(unique(dist_coarse))
%             mean_mod(1,ii) = mean(lightmod_sort((floor(dist_sort/50)==dist_coarse(ii)) & (supp_pop_sort==1)));
%             mean_mod(2,ii) = mean(lightmod_sort((floor(dist_sort/50)==dist_coarse(ii)) & (supp_pop_sort==2)));
%             mean_mod(3,ii) = mean(lightmod_sort((floor(dist_sort/50)==dist_coarse(ii)) & (supp_pop_sort==3)));
%         end
%         mean_mod(isnan(mean_mod)) = 0;
%         for i = 1:size(mean_mod,1)
%             plot(mean_mod(i,:),dist_range(1:2:end),'color',color_mat(i,:),'linewidth',2)
%         end

        view(0,270)
        xlim([-1 1])
        ylim([0 max(abs(distfromfirstch))])
        yax = get(gca,'YLim');
        line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
        xlabel('LMI_v_i_s','Fontsize',24)
        ylabel(strcat('Depth (�m)'),'Fontsize',24)
        set(gca,'fontsize',18,'linewidth',2)
        legend({'Other','V1-supp','SC-supp','V1+SC-supp'},'location','best','fontsize',12)
        print(gcf, '-dpng',sprintf('%s_%d','lightmodbydepth_top_bymodtype',i))
        print(gcf,'-painters','-depsc',sprintf('%s_%d','lightmodbydepth_top_bymodtype',i))

elseif contains(pop,'drivermod','ignorecase',1) 
    LED_titles = {'L5&L6CTs inactivated','L5CTs inactivated','L6CTs inactivated'};
    figure; 
    plotord = [3 1 2]; % plot all first, then L5CT, then L6CT
    for ii = 1:size(lightmod,2)
        i=plotord(ii);
        subplot(1,size(lightmod,2),ii)
        scatter(lightmod(:,i),distfromfirstch',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0],'markerfacealpha',1)
        hold on;
        scatter(lightmod(supp_cells,i),distfromfirstch(supp_cells)',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',1)
        scatter(lightmod(enh_cells,i),distfromfirstch(enh_cells)',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',.5)
        h = get(gca,'ytick');
        % set(gca,'yticklabel',h);
        view(0,270)
        xlim([-1 1])
        ylim([0 max(abs(distfromfirstch))])
        yax = get(gca,'YLim');
        line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
%         legend('other','suppressed','enhanced')
        xlabel('LMI_v_i_s','Fontsize',24)
        ylabel(strcat('Depth (�m)'),'Fontsize',24)
        title(LED_titles{ii})
        set(gca,'fontsize',18,'linewidth',2)
    end
    set(gcf,'Paperposition',[0 0 20 6])
    print(gcf, '-dpng','lightmodbydepth_fromtop_all')
    print(gcf,'-painters','-depsc','lightmodbydepth_fromtop_all')
elseif contains(area,'SC')
    LED_titles = {'SC inactivated','V1/L5 inactivated','SC+V1/L5 inactivated'};
    figure; 
    plotord = [1 2 3]; %SC, V1, SC+V1
    for ii = 1:size(lightmod,2)
        i=plotord(ii);
        subplot(1,size(lightmod,2),ii)
        scatter(lightmod(:,i),distfromfirstch',75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0],'markerfacealpha',1)
        hold on;
        scatter(lightmod(supp_cells,i),distfromfirstch(supp_cells)',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',1)
        scatter(lightmod(enh_cells,i),distfromfirstch(enh_cells)',75,'filled','markerfacecolor',color_mat(i,:),'markerfacealpha',.5)
        h = get(gca,'ytick');
        % set(gca,'yticklabel',h);
        view(0,270)
        xlim([-1 1])
        ylim([0 max(abs(distfromfirstch))])
        yax = get(gca,'YLim');
        line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
%         legend('other','suppressed','enhanced')
        xlabel('LMI_v_i_s','Fontsize',24)
        ylabel(strcat('Depth (�m)'),'Fontsize',24)
        title(LED_titles{ii})
        set(gca,'fontsize',18,'linewidth',2)
    end
    set(gcf,'Paperposition',[0 0 20 6])
    print(gcf, '-dpng','lightmodbydepth_fromtop_all')
    print(gcf,'-painters','-depsc','lightmodbydepth_fromtop_all')
end

% %% trying new thing...look at probability of being suppressed/activated as a function of distance from other suppressed cells
% if exist('light_sig','var')
%     if isempty(supp_cells) && ~isempty(supp_ons_cells)
%          supp_cells_tmp = supp_ons_cells;
%          enh_cells_tmp = enh_cells(~ismember(enh_cells,supp_ons_cells));        % necessary b/c some onset-suppressed cells are enhanced in blank trials!
%         fprintf('using prestim period to classify "suppressed" cells because otherwise there are none')
%     else
%         supp_cells_tmp = supp_cells;
%         enh_cells_tmp = enh_cells;
%     end
% 
%     max_ch = [waveforms(clean_units).max_ch];
%     max_ch_shuf = nan(1000,length(max_ch));
%     used_cells = [];
% 
%     for i=1:length(supp_cells_tmp)
%         this_cellch = max_ch(supp_cells_tmp(i));
%         shuf_inds = find(exp_num==exp_num(supp_cells_tmp(i)) & shk(clean_units)==shk(clean_units(supp_cells_tmp(i)))); % indexes to clean_units of other units recorded from same shank in same exp
%     %     shuf_inds(shuf_inds==supp_cells(i)) = [];  % exclude itself! (for shuffling below)
%         shuf_inds = shuf_inds(~ismember(shuf_inds,[supp_cells_tmp(i) used_cells])); % exclude itself AND already used cells! (for shuffling below)
%         these_chs = max_ch(shuf_inds); % corresponding channel positions of other units recorded from same shank in same exp
%         shuf_chs = nan(1000,length(these_chs));
%         for ii=1:1000  % shuffle the channel position assignments across all neighboring (same exp, same shank) units
%             shuf_chs(ii,:) = these_chs(randperm(length(these_chs)));
%         end
%         for n=1:max(distfromfirstch)/min(diff(unique(distfromfirstch)))+1   % bin by vertical distance between probe contacts
%             den(i,n) = sum(abs(max_ch-this_cellch)==n-1 & exp_num==exp_num(supp_cells_tmp(i)) & shk(clean_units)==shk(clean_units(supp_cells_tmp(i))) & ~ismember(1:length(clean_units),used_cells)); % total number of units n distance away from suppressed unit (denominator)
%             num_supp(i,n) = length(intersect(find(abs(max_ch-this_cellch)==n-1 & exp_num==exp_num(supp_cells_tmp(i)) & shk(clean_units)==shk(clean_units(supp_cells_tmp(i)))),supp_cells_tmp)); % how many of those were also suppressed (numerator)
%             num_enh(i,n) = length(intersect(find(abs(max_ch-this_cellch)==n-1 & exp_num==exp_num(supp_cells_tmp(i)) & shk(clean_units)==shk(clean_units(supp_cells_tmp(i)))),enh_cells_tmp));  % how many of those were enhanced (numerator)
%             if n==1 % make sure not to include itself!
%                 den(i,n) = den(i,n)-1;
%                 num_supp(i,n) = num_supp(i,n)-1;
%             end
% 
%             %shuffle distribution
%             for ii=1:1000
%                 den_shuf(i,n,ii) = sum(abs(shuf_chs(ii,:)-this_cellch)==n-1);  % shuffle denominator should be SAME as unshuffled denominator (above) b/c irrespective of unit index
%                 num_supp_shuf(i,n,ii) = length(intersect(shuf_inds(abs(shuf_chs(ii,:)-this_cellch)==n-1),supp_cells_tmp));
%                 num_enh_shuf(i,n,ii) = length(intersect(shuf_inds(abs(shuf_chs(ii,:)-this_cellch)==n-1),enh_cells_tmp));
%             end
%         end
%         % need to exclude this suppressed cell for further comparisons,
%         % otherwise pairs of two suppressed cells will be double-counted!
%         used_cells =[used_cells supp_cells_tmp(i)];
%         supp_cells_tmp(i) = nan;
%     end
%     dist_mat_supp = sum(num_supp)./sum(den);
%     dist_mat_enh = sum(num_enh)./sum(den);
%     shufdist_mat_supp = squeeze(sum(num_supp_shuf))./squeeze(sum(den_shuf));
%     shufdist_mat_enh = squeeze(sum(num_enh_shuf))./squeeze(sum(den_shuf));
%     shuf_mat_supp = mean(shufdist_mat_supp,2);
%     shuf_mat_enh = mean(shufdist_mat_enh,2);
%     shuf_supp_se = std(shufdist_mat_supp,[],2)./sqrt(1000);
%     shuf_enh_se = std(shufdist_mat_enh,[],2)./sqrt(1000);
%     distvec = 0:min(diff(unique(distfromfirstch))):max(distfromfirstch);
%     %     critval = .05/length(distvec);
%     critval = .025;     % don't think mult compare correction is necessary for shuffle distribution? but use .025 b/c two-tailed (12/28/20)
%     p_supp = nan(length(distvec),2);
%     p_enh = p_supp;
%     for x=1:length(distvec)
%         if ~isnan(dist_mat_supp(x))
%             p_supp(x,:) = [(sum(shufdist_mat_supp(x,:)>=dist_mat_supp(x)))/1000 sum(shufdist_mat_supp(x,:)<=dist_mat_supp(x))/1000];  % MUST be inclusive! (<= and >=)
%             p_enh(x,:) = [(sum(shufdist_mat_enh(x,:)>=dist_mat_enh(x)))/1000 sum(shufdist_mat_enh(x,:)<=dist_mat_enh(x))/1000];
%         end
%     end
%     sig_supp = min(p_supp,[],2)<critval;
%     sig_enh = min(p_enh,[],2)<critval;
% 
%     figure; hold on;
%     % shadedErrorBar(distvec,-shuf_mat_supp,shuf_supp_se)
%     % shadedErrorBar(distvec,shuf_mat_enh,shuf_enh_se)
%     % plot(distvec,-dist_mat_supp,'color',color_mat(lightcond,:),'linewidth',2)
%     % plot(distvec,dist_mat_enh,'color',color_mat(lightcond,:),'linewidth',2,'linestyle','--')
%     b1=bar(distvec,-shuf_mat_supp,1,'FaceColor',[1 1 1],'edgecolor',[.75 .75 .75],'facealpha',1,'linewidth',1.25);
%     b2 = bar(distvec,shuf_mat_enh,1,'FaceColor',[1 1 1],'edgecolor',[.75 .75 .75],'facealpha',1,'linewidth',1.25);
%     b3 = bar(distvec,-dist_mat_supp,1,'FaceColor',color_mat(lightcond,:),'edgecolor','k','facealpha',1,'linewidth',1.25);
%     b4 = bar(distvec,dist_mat_enh,1,'FaceColor',color_mat(lightcond,:),'edgecolor','k','facealpha',.5,'linewidth',1.25);
%     s1 = scatter(distvec(sig_supp),-sig_supp(sig_supp==1)+.1,'*','markeredgecolor',color_mat(lightcond,:),'markeredgealpha',1);
%     s2 = scatter(distvec(sig_enh)+5,sig_enh(sig_enh==1)-.1,'*','markeredgecolor',color_mat(lightcond,:),'markeredgealpha',.25);
%     ylim([-1 1])
%     yticklabels(num2str(abs(yticks'*100)))
%     xlim([0-min(diff(distvec))/2 distvec(max([find(~isnan(dist_mat_supp),1,'last') find(~isnan(dist_mat_enh),1,'last')]))+min(diff(distvec))/2])
%     xax=xlim;
%     line(xax,[0 0],'color','k')
%     l=legend([b1 b4 b3],{'shuffle','enhanced','suppressed'});
%     set(gca,'linewidth',2)
%     xlabel('Distance (in �m) from suppressed unit','fontsize',15)
%     ylabel('% suppressed/enhanced','fontsize',15)
%     print(gcf, '-dpng','propmodulated_bydist')
%     print(gcf,'-painters','-depsc','propmodulated_bydist')
% end

%% scatter plots

for i = 1:num_lcs
     types = zeros(1,size(FRev,1));
    types(supp_cells) = 2;
    types(enh_cells) = 1;
    
    % FR light vs no light - visual vs nonvisual units, visual trials
%     plot_scatter(FRev(:,[1 i+1]), ismember(1:length(clean_units),visual_cells), {color_mat(i,:),color_mat(i,:)}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis',i), {'Nonvisual','Visual'}, 2)     
    plot_scatter(FRev(visual_cells,[1 i+1]),  types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],  'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis',i), {'other','enhanced','suppressed'}, 1)     

    % FR light vs no light - visual vs nonvisual units, blank trials
%     plot_scatter(FRbl(:,[1 i+1]), ismember(1:length(clean_units),visual_cells), {color_mat(i,:),color_mat(i,:)}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis_bl',i), {'Nonvisual','Visual'}, 2)    
    plot_scatter(FRbl(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis_bl',i),{'other','enhanced','suppressed'}, 1)    

    % FR light vs no light - visual vs nonvisual units, preferred trials
    plot_scatter(FRpref(:,[1 i+1]), ismember(1:length(clean_units),visual_cells),  {color_mat(i,:),color_mat(i,:)}, [.25 .75], 'Spks/s (light OFF)', 'Spks/s (light ON)', sprintf('%s_%d','FR_byvis_pref',i), {'Nonvisual','Visual'}, 1)    

    % change in preferred FR light vs no light - visual vs nonvisual units
    plot_scatter(FRpref(:,[1 i+1])-FRbl(:,[1 i+1]), ismember(1:length(clean_units),visual_cells),  {color_mat(i,:),color_mat(i,:)}, [.25 .75], 'Vis-evoked FR\Delta (Spks/s-light OFF)', 'Vis-evoked FR\Delta (Spks/s-light ON)', sprintf('%s_%d','FR_vischange_pref',i), {'Nonvisual','Visual'}, 1)     

    % abs val change in evoked FR light vs no light - light modulated vs non-modulated
    % units
   
    plot_scatter(abs(FRev(:,[1 i+1])-FRbl(:,[1 i+1])), types, {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1], '|Vis-evoked FR\Delta| (Spks/s-light OFF)', '|Vis-evoked FR\Delta| (Spks/s-light ON)', sprintf('%s_%d','FR_vischange',i), {'Non-modulated','Light-activated','Light-suppressed'}, 2)    

    % %change in evoked FR light vs no light - light modulated vs non-modulated
    plot_scatter(abs(FRev(:,[1 i+1])-FRbl(:,[1 i+1]))./FRbl(:,[1 i+1]), types, {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1], '|Vis-evoked FR%\Delta| (Spks/s-light OFF)', '|Vis-evoked FR%\Delta| (Spks/s-light ON)', sprintf('%s_%d','FR_vis%change',i), {'Non-modulated','Light-activated','Light-suppressed'}, 2)     

    % change in bursting
%     plot_scatter([burstrate burstrate_light(:,i)], ismember(1:length(clean_units),visual_cells), {color_mat(i,:),color_mat(i,:)}, [.5 1], 'Bursting rate (light OFF)', 'Bursting rate (light ON)', sprintf('%s_%d','Bursts',i), {'Nonvisual','Visual'}, 1)     
%     plot_scatter([burstrate burstrate_light(:,i)],types, {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1],  'Bursting rate (light OFF)', 'Bursting rate (light ON)', sprintf('%s_%d','Bursts',i), {'Non-modulated','Light-activated','Light-suppressed'}, 1)     
    plot_scatter([burstrate(visual_cells) burstrate_light(visual_cells,i)],types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],  'Bursting rate (light OFF)', 'Bursting rate (light ON)', sprintf('%s_%d','Bursts',i), {'other','enhanced','suppressed'}, 1)     

    % change in linear responses (for visual cells only)
%     plot_scatter(log10(Fratio(visual_cells,[1 i+1])),types(visual_cells), {[.85 .85 .85],color_mat(i,:),color_mat(i,:)},[1 .5 1], 'Fratio (light OFF - log scale)', 'Fratio (light ON - log scale)', sprintf('%s_%d','Fratio',i), {'Non-modulated','Light-activated','Light-suppressed'}, 2) %
    plot_scatter(Fratio_bs(visual_cells,[1 i+1]), types(visual_cells), {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'Fratio (light OFF - log scale)', 'Fratio (light ON)', sprintf('%s_%d','Fratio_bs',i), {'Non-modulated','Light-activated','Light-suppressed'}, 1)    
    plot_scatter(zF1(visual_cells,[1 i+1]), types(visual_cells), {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'zF1 (light OFF)', 'zF1 (light ON)', sprintf('%s_%d','ZF1',i), {'Non-modulated','Light-activated','Light-suppressed'}, 1)
%     plot_scatter(F1(visual_cells,[1 i+1]), types(visual_cells), {[.75 .75 .75],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'F1 (light OFF - spks/s)', 'F1 (light ON - spks/s)', sprintf('%s_%d','F1',i), {'Non-modulated','Light-activated','Light-suppressed'}, 1)    
    plot_scatter(F1(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],  'F1 (light OFF - spks/s)', 'F1 (light ON - spks/s)', sprintf('%s_%d','F1',i),{'other','enhanced','suppressed'}, 1)       
% cumulative distribution plot of zF1 values for different cell types
    figure; plot(sort(zF1(intersect(other_cells,visual_cells),1)),100*(1:length(intersect(other_cells,visual_cells)))./length(intersect(other_cells,visual_cells)),'linewidth',2,'color',[.85 .85 .85])     % non-modulated visual cells
    hold on;
    plot(sort(zF1(supp_cells,1)),100*(1:length(zF1(supp_cells,1)))./length(zF1(supp_cells,1)),'linewidth',2,'color',color_mat(i,:))
    p=plot(sort(zF1(enh_cells,1)),100*(1:length(zF1(enh_cells,1)))./length(zF1(enh_cells,1)),'linewidth',2,'color',color_mat(i,:));
    if ~isempty(p); p.Color(4) = .5; end
    xlabel('zF1')
    ylabel('Proportion of units (%)')
    set(gca,'fontsize',16,'linewidth',2)
    l=legend;
    set(l,'String',{'Other visual cells','Suppressed','Enhanced'},'location','SouthEast')
    set(l,'fontsize',12)
    save_name = 'CumDist_zF1_bytype';
    print(gcf,'-dpng',save_name)
    print(gcf,'-painters','-depsc',save_name)
%     [p,anovatab,stats] = kruskalwallis([zF1(intersect(other_cells,visual_cells),1);zF1(supp_cells,1);zF1(enh_cells,1)],[ones(length(intersect(other_cells,visual_cells)),1);2*ones(length(supp_cells),1); 3*ones(length(enh_cells),1)],'off');
%     zF1_stats = multcompare(stats,'ctype','dunn-sidak','display','off');

    figure; plot(sort(Fratio_bs(visual_cells,1)),100*(1:length(visual_cells))./length(visual_cells),'linewidth',2,'color',[.85 .85 .85])
    hold on;
    plot(sort(Fratio_bs(supp_cells,1)),100*(1:length(Fratio_bs(supp_cells,1)))./length(Fratio_bs(supp_cells,1)),'linewidth',2,'color',color_mat(i,:))
    p=plot(sort(Fratio_bs(enh_cells,1)),100*(1:length(Fratio_bs(enh_cells,1)))./length(Fratio_bs(enh_cells,1)),'linewidth',2,'color',color_mat(i,:));
    if ~isempty(p); p.Color(4) = .5; end
    xlabel('MI')
    ylabel('Proportion of units (%)')
    set(gca,'fontsize',16,'linewidth',2)
    l=legend;
    set(l,'String',{'Other visual cells','Suppressed','Enhanced'},'location','SouthEast')
    set(l,'fontsize',12)
    save_name = 'CumDist_Fratio_bs_bytype';
    print(gcf,'-dpng',save_name)
    print(gcf,'-painters','-depsc',save_name)
    [p,anovatab,stats] = kruskalwallis([Fratio_bs(intersect(other_cells,visual_cells),1);Fratio_bs(supp_cells,1);Fratio_bs(enh_cells,1)],[ones(length(intersect(other_cells,visual_cells)),1);2*ones(length(supp_cells),1); 3*ones(length(enh_cells),1)],'off');
    F1_stats = multcompare(stats,'ctype','dunn-sidak','display','off');

    
    % change in tuning (for visual cells only)
    plot_scatter(OSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'OSI(CV) (light OFF)', 'OSI(CV) (light ON)', sprintf('%s_%d','OSI(CV)_bylightmod',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(OSI(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'OSI (light OFF)', 'OSI (light ON)', sprintf('%s_%d','OSI_bylightmod',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(DSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'DSI(CV) (light OFF)', 'DSI(CV) (light ON)', sprintf('%s_%d','DSI(CV)_bylightmod',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(DSI(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'DSI (light OFF)', 'DSI (light ON)', sprintf('%s_%d','DSI_bylightmod',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(newOSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1], 'OSI(CV) (light OFF)', 'OSI(CV) (light ON)', sprintf('%s_%d','OSI(CV)_bylightmod_baselinenorm',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(newOSI(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'OSI (light OFF)', 'OSI (light ON)', sprintf('%s_%d','OSI_bylightmod_baselinenorm',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(newDSI_CV(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'DSI(CV) (light OFF)', 'DSI(CV) (light ON)', sprintf('%s_%d','DSI(CV)_bylightmod_baselinenorm',i),{'other','enhanced','suppressed'}, 1)     % first lightcond pwr
    plot_scatter(newDSI(visual_cells,[1 i+1]), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, [1 .5 1],'DSI (light OFF)', 'DSI (light ON)', sprintf('%s_%d','DSI_bylightmod_baselinenorm',i), {'other','enhanced','suppressed'}, 1)     % first lightcond pwr

    plot_hist(OSI_CV(visual_cells,1), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, {color_mat(i,:),color_mat(i,:),color_mat(i,:)}, [1 .5 1], .05, 'count', 'OSI(CV)', '# of units', 'OSI(CV)_hist')
    plot_hist(DSI_CV(visual_cells,1), types(visual_cells), {[1 1 1],color_mat(i,:),color_mat(i,:)}, {color_mat(i,:),color_mat(i,:),color_mat(i,:)}, [1 .5 1], .05, 'count', 'DSI(CV)', '# of units', 'DSI(CV)_hist')

    close all
end

% if SC inactivation experiment, plot V1 (or L5 specifically) modulation by
% SC modulation - look at modulation in VISUAL trials
if contains(exp_type,'SC','ignorecase',1)
    plot_scatter(lightmod(:,1:2),  supp_pop, {[1 1 1],color_mat(1,:), color_mat(2,:), color_mat(3,:)} ,[1 1 1 1],'Lightmod (V1 inactivated)', 'Lightmod (SC inactivated)', 'Lightmod_V1vsVSC',{'Other','V1-suppressed','SC-suppressed','V1+SC-suppressed'},0)
    plot_scatter(lightmod_bl(:,1:2),  supp_pop, {[1 1 1],color_mat(1,:), color_mat(2,:), color_mat(3,:)} ,[1 1 1 1],'Lightmod (V1 inactivated)', 'Lightmod (SC inactivated)', 'Lightmod_bl_V1vsVSC',{'Other',''},0)

    % plot FRs (V1 inactivated) vs. FRs (V1&SC inactivated) to see if
    % there's an additive effect
    plot_scatter(FRev(:,[2 end]), supp_pop, {[1 1 1],color_mat(1,:), color_mat(2,:), color_mat(3,:)} ,[1 1 1 1],'FR_v_i_s (V1 inactivated-spks/s)', 'FR_v_i_s (V1&SC inactivated-spks/s)', 'FRs_V1vsSC',{'Other','V1-suppressed','SC-suppressed','V1+SC-suppressed'},2)

    plot_scatter([FRev(supp_pop==1,[2 end]); FRev(supp_pop==2,[3 end]); FRev(supp_pop==3,[2 end])],[ones(sum(supp_pop==1),1); 2*ones(sum(supp_pop==2),1); 3*ones(sum(supp_pop==3),1)],{color_mat(1,:), color_mat(2,:), color_mat(3,:)},[1 1 1],'FR_v_i_s - one LED (spks/s)', 'FR_v_i_s - both LEDs (spks/s)', 'FRs_OnevsTwo',{},2)
% for experiments with LEDs over different parts of cortex (currently only
% V1 & VisL)
elseif contains(exp_type,'2LEDs','ignorecase',1)
   plot_scatter(lightmod(:,1:2),  supp_pop, {[.9 .9 .9],color_mat(1,:), color_mat(2,:), color_mat(3,:)} ,[1 1 1 1],'Lightmod (V1 inactivated)', 'Lightmod (VisL inactivated)', 'Lightmod_V1vsVisL',{'Other','V1-suppressed','VisL-suppressed','V1+VisL-suppressed'},0)

    % plot FRs (V1 inactivated) vs. FRs (V1&SC inactivated) to see if
    % there's an additive effect
    plot_scatter(FRev(:,[2 end]), supp_pop, {[.9 .9 .9],color_mat(1,:), color_mat(2,:), color_mat(3,:)} ,[1 1 1 1],'Spks/s (V1 inactivated)', 'Spks/s (V1&VisL inactivated)', 'FRs_V1vsVisL',{'Other','V1-suppressed','VisL-suppressed','V1+VisL-suppressed'},1)

elseif contains(area,'SC')
        plot_scatter(lightmod(:,[2 1]),  types, {[1 1 1],color_mat(1,:), color_mat(1,:)} ,[1 .5 1],'Lightmod (V1 inactivated)', 'Lightmod (SC inactivated)', 'Lightmod_V1vsVSC',{'Other','SC-enhanced','SC-suppressed'},0)
    plot_scatter(FRev(:,[2 end]), types, {[1 1 1],color_mat(1,:), color_mat(1,:)} ,[1 .5 1],'FR_v_i_s (SC inactivated-spks/s)', 'FR_v_i_s (V1&SC inactivated-spks/s)', 'FRs_V1vsSC',{'Other','SC-enhanced','SC-suppressed'},1)

end

%%  make population-level average PSTHs (normalized to prestim baseline FR)
% thus, baseline = 1 
binsize = .025;

% determine baseline FRs (for visual and blank trials separately)
mean_bs = repmat(mean(psthV(:,1:200/(binsize*1000),:),2),1,size(psthV,2),1);    % define baseline from first 200ms
mean_bs_bl = repmat(mean(psthBl(:,1:200/(binsize*1000),:),2),1,size(psthBl,2),1);
mean_bs(mean_bs==0) = nan;     % because otherwise could get infinity when normalizing by baseline
mean_bs_bl(mean_bs_bl==0) = nan;

norm_psth = psthV./mean_bs; % normalizes to prestim baseline, per condition (so that prestim=1)
norm_psth_bl = psthBl./mean_bs_bl;

norm_mean = nanmean(norm_psth,3); % average across units
norm_se = nanstd(norm_psth,0,3)./sqrt(size(norm_psth,3));
norm_mean_bl = nanmean(norm_psth_bl,3);
norm_se_bl = nanstd(norm_psth_bl,0,3)./sqrt(size(norm_psth_bl,3));
% and separately for supp and enh cells
norm_supp_mean = nanmean(norm_psth(:,:,supp_cells),3); 
norm_supp_se = nanstd(norm_psth(:,:,supp_cells),0,3)./sqrt(length(supp_cells));
norm_enh_mean = nanmean(norm_psth(:,:,enh_cells),3); 
norm_enh_se = nanstd(norm_psth(:,:,enh_cells),0,3)./sqrt(length(enh_cells));

% set up parameters for plotting
for n = 1:length(params)
    tri_st(n) = params(n).prestim;
    tri_visdur(n) = params(n).stimtime;
    if exist('light_sig','var')
        start_times(n,:) = round(params(n).av_light_start(lightconds{n}),2)-params(n).prestim;  % find LED start times in diff exps, and round to nearest hundreth
        stim_durs(n,:) = round(params(n).light_dur(lightconds{n}+1),1);
    else
        start_times(n,:) = nan;
    end
end

if length(unique(tri_st))>1 || length(unique(tri_visdur))>1
    error('Disparate trial times across experiments')
else
    tri_st = tri_st(1);
    tri_visdur = tri_visdur(1);
end

if contains(exp_type,'halo','ignorecase',1)
    patch_color = [.9 .1 .1];  % red/orange patch for halo
else
    patch_color = [0 .1 1];  % default patch color is blue
end

% VISUAL trials
pop_fig= figure;
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean(1,:), norm_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean(i+1,:), norm_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_norm';
print(pop_fig,'-dpng',save_fig_name)
print(pop_fig,'-painters','-depsc',save_fig_name)

% BLANK trials
pop_fig2= figure;
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean_bl(1,:), norm_se_bl(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_mean_bl(i+1,:), norm_se_bl(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
set(gca,'fontsize',16,'linewidth',2);
save_fig_name = 'PopulationPSTH_blank_norm';
print(pop_fig2,'-dpng',save_fig_name)
print(pop_fig2,'-painters','-depsc',save_fig_name)

% VISUAL trials - by type
pop_fig3= figure;
subplot(121)
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_supp_mean(1,:), norm_supp_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_supp_mean(i+1,:), norm_supp_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
title(sprintf('Suppressed units (n=%d)',length(supp_cells)),'fontsize',14);
set(gca,'fontsize',16,'linewidth',2);

subplot(122)
xlim([-tri_st+binsize tri_visdur])   % ticks mark the END of 25ms bins
hold on;
shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_enh_mean(1,:), norm_enh_se(1,:),{'Color',[0 0 0],'linewidth',2},1);
for i = 1:num_lcs
    shadedErrorBar(-tri_st:binsize:tri_visdur-binsize,norm_enh_mean(i+1,:), norm_enh_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
end
yax = get(gca,'YLim');
line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
if length(unique(start_times))>1 && length(unique(start_times))*length(params)==length(start_times) % if multiple start times in EVERY experiment
    for ii = 1:length(unique(start_times))
        line([start_times(ii) start_times(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
        line([start_times(ii)+stim_durs(ii) start_times(ii)+stim_durs(ii)],yax,'Color',color_mat(ii,:),'LineStyle','--','linewidth',2)
    end
elseif exist('light_sig','var')
    patch([start_times(1) start_times(1) start_times(1)+stim_durs(1) start_times(1)+stim_durs(1) start_times(1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
end
ylim(yax)
xlabel('Time from visual stim onset (s)','fontsize',24)
ylabel('Normalized firing rate (spks/s)','fontsize',24)
title(sprintf('Activated units (n=%d)',length(enh_cells)),'fontsize',14);
set(gca,'fontsize',16,'linewidth',2);

set(gcf, 'Position', [100, 100, 1000, 420])
save_fig_name = 'PopulationPSTH_bytype_norm';
print(pop_fig3,'-dpng',save_fig_name)
print(pop_fig3,'-painters','-depsc',save_fig_name)

%% pie graph of number of enhanced vs. suppressed vs. unaffected units
% specifically looking among visually-responsive units
if exist('light_sig','var')
    figure;
    pie_data = [length(intersect(visual_cells,supp_cells)) length(intersect(visual_cells,enh_cells)) length(visual_cells)-length(intersect([supp_cells; enh_cells],visual_cells))];
    piegraph = pie(pie_data);
    piegraph_labels = {'Suppressed','Activated','Other'};
    piegraph_labels = piegraph_labels(pie_data>0);
    colormap([color_mat(lightcond,:); color_mat(lightcond,:); 1 1 1])
    % set(piegraph([1:2:end]),'edgecolor',color_mat(1,:))
    set(piegraph([2:2:end]),'fontsize',16)
    for i=2:2:length(piegraph)
        set(piegraph(i),'string',sprintf('%s (%s)',piegraph_labels{i/2},piegraph(i).String))
        if strcmpi(piegraph_labels{i/2},'Activated')
            set(piegraph(i-1),'facealpha',.5) % enhanced units will be lighter colored
        end
    end
    title(sprintf('%d Visual units',length(visual_cells)),'fontsize',16)
    print(gcf,'-dpng','Piegraph_lighteffects')
    print(gcf,'-painters','-depsc','Piegraph_lighteffects')

    if contains(exp_type,'SC','ignorecase',1) || contains(exp_type,'2LEDs','ignorecase',1)
        figure;
        pie2_data = [length(intersect(visual_cells,find(supp_pop==1))) length(intersect(visual_cells,find(supp_pop==2))) length(intersect(visual_cells,find(supp_pop==3))) length(intersect(visual_cells,find(supp_pop==0)))];
        piegraph2 = pie(pie2_data);
        colormap([color_mat(pie2_data(1:end-1)>0,:); 1 1 1])
        set(piegraph2([2:2:end]),'fontsize',16)
        for i=2:2:length(piegraph2)
            set(piegraph2(i),'string',sprintf('%s',piegraph2(i).String))
        end
        title(sprintf('%d Visual units',length(visual_cells)),'fontsize',16)
        print(gcf,'-dpng','Piegraph_bysupptype')
        print(gcf,'-painters','-depsc','Piegraph_bymodtype')
    end
end
    

%% save stuff

num_units = length(clean_units);
fileID = fopen('results.txt','w');
fprintf(fileID,'Number of visually responsive cells: %d of %d\r\n',length(visual_cells),num_units);
fprintf(fileID,'Percent visually responsive: %.2f\r\n',100*length(visual_cells)/num_units);
fprintf(fileID,'Number of light-modulated cells: %d of %d\r\n',length(light_cells),num_units);
fprintf(fileID,'Percent light-modulated: %.2f\r\n', 100*length(light_cells)/num_units);
fprintf(fileID,'Number of tuned cells: %d of %d\r\n',length(tuned_cells),num_units);
fprintf(fileID,'Percent significantly tuned: %.2f\r\n', 100*length(tuned_cells)/num_units);
fprintf(fileID,'Number of regular-spiking cells: %d of %d\r\n',length(reg_cells),num_units);
fprintf(fileID,'Percent regular-spiking: %.2f\r\n', 100*length(reg_cells)/num_units);
fprintf(fileID,'Number of fast-spiking cells: %d of %d\r\n',length(FS_cells),num_units);
fprintf(fileID,'Percent fast-spiking: %.2f\r\n', 100*length(FS_cells)/num_units);
fprintf(fileID,'Number of linear cells: %.2f\r\n', sum(Fratio_bs(:,1)>1));
fprintf(fileID,'Percent linear cells: %.2f\r\n', 100*sum(Fratio_bs(:,1)>1)/size(Fratio,1));
fclose(fileID);

good_units = [unitinfo(clean_units).name];
exps = {params.exp_name};
if exist('light_sig','var')
    fileID2 = fopen('stats.txt','w');
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, ALL cells: %s (%s) \r\n',num2str(round(lightsig_all,5)),num2str(lightsig_all_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, ALL cells: %s (%s) \r\n',num2str(round(lightsig_bl,5)),num2str(lightsig_bl_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, VISUAL cells: %s (%s) \r\n',num2str(round(lightsig_vis,5)),num2str(lightsig_vis_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, VISUAL cells: %s (%s) \r\n',num2str(round(lightsig_blvis,5)),num2str(lightsig_blvis_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), VISUAL trials, NONVISUAL cells: %s (%s) \r\n',num2str(round(lightsig_nonvis,5)),num2str(lightsig_nonvis_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on FR (low to high conditions), BLANK trials, NONVISUAL cells: %s (%s) \r\n',num2str(round(lightsig_blnonvis,5)),num2str(lightsig_blnonvis_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on pref FR (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_vispref,5)),num2str(lightsig_vispref_dir));
    fprintf(fileID2,'Signed-rank test of sig light effect on visual FR change (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_visFRdelt,5)),num2str(lightsig_visFRdelt_dir));
    if ~isempty(supp_cells)
        fprintf(fileID2,'Signed-rank test of sig light effect on visual FR change (low to high conditions), SUPPRESSED cells: %s (%s) \r\n',num2str(round(lightsig_visFRdelt_supp,5)),num2str(lightsig_visFRdelt_supp_dir));
    end
%     fprintf(fileID2,'Signed-rank test of sig light effect on pref visual FR change (low to high conditions): %s (%s) \r\n',num2str(round(lightsig_prefFRdelt,3)),num2str(lightsig_prefFRdelt_dir));

    fprintf(fileID2,'\r\n');
    if ~isempty(supp_cells)
        fprintf(fileID2,'Signed-rank test of sig OSI_CV change (low to high conditions) SUPPRESSED cells: %s (%s) \r\n',num2str(round(osiCVsig_tuned,5)),num2str(osiCVsig_tuned_dir));
        fprintf(fileID2,'Signed-rank test of sig OSI change (low to high conditions), SUPPRESSED cells: %s (%s) \r\n',num2str(round(osisig_tuned,5)),num2str(osisig_tuned_dir));
        fprintf(fileID2,'Signed-rank test of sig DSI_CV change (low to high conditions), SUPPRESSED cells: %s (%s) \r\n',num2str(round(dsiCVsig_tuned,5)),num2str(dsiCVsig_tuned_dir));
        fprintf(fileID2,'Signed-rank test of sig DSI change (low to high conditions), SUPPRESSED cells: %s (%s) \r\n',num2str(round(dsisig_tuned,5)),num2str(dsisig_tuned_dir));
    end
    fprintf(fileID2,'Signed-rank test of sig OSI_CV change (low to high conditions), ALL cells:%s (%s) \r\n',num2str(round(osiCVsig_vis,5)),num2str(osiCVsig_vis_dir));
    fprintf(fileID2,'Signed-rank test of sig OSI change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(osisig_vis,5)),num2str(osisig_vis_dir));
    fprintf(fileID2,'Signed-rank test of sig DSI_CV change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(dsiCVsig_vis,5)),num2str(dsiCVsig_vis_dir));
    fprintf(fileID2,'Signed-rank test of sig DSI change (low to high conditions), ALL cells: %s (%s) \r\n',num2str(round(dsisig_vis,5)),num2str(dsisig_vis_dir));

    fprintf(fileID2,'\r\n');
    fprintf(fileID2,'Signed-rank test of burst rate change (low to high conditions) ALL cells: %s (%s) \r\n',num2str(round(burstrate_sig,5)),num2str(burstrate_dir));
    fprintf(fileID2,'Signed-rank test of baseline-subtracted Fratio (low to high conditions) Visual cells: %s (%s) \r\n',num2str(round(Fratio_sig,5)),num2str(Fratio_dir));
    fprintf(fileID2,'Signed-rank test of zF1 (low to high conditions) Visual cells: %s (%s) \r\n',num2str(round(zF1_sig,5)),num2str(zF1_dir));
    fprintf(fileID2,'Signed-rank test of F1 (low to high conditions) Visual cells: %s (%s) \r\n',num2str(round(F1_sig,5)),num2str(F1_dir));

    fprintf(fileID2,'\r\n');
    fprintf(fileID2,'K-W test diff distributions - other vs supp, other vs enh, supp vs enh: %s  \r\n',num2str(round(F1_stats(:,end)',5)));

    fclose(fileID2);
    
    save('unit_info.mat','good_units','exps','exp_num','lightmod','lightmod_bl','supp_cells','supp_ons_cells','enh_cells','visual_cells','FRev','FRbl','newOSI_CV')  % note that good_units contains the actual cluster names, whereas supp_cells and enh_cells are indexes from 1:length(good_units)

else
    save('unit_info.mat','good_units','exps','exp_num','FRev','FRb','visual_cells')  % note that good_units contains the actual cluster names, whereas supp_cells and enh_cells are indexes from 1:length(good_units)
end

fileID3 = fopen('cleanunits.txt','w');
fprintf(fileID3,'%d\r\n',[unitinfo(clean_units).name]);
fclose(fileID3);


close all

end

function plot_scatter(data, color_var, colors, alpha, xlab, ylab, title, leg, lobf)
% if lobf = 1, make one lobf regardless of variables; if lobf = 2, make
% separate lobfs for each variable
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
    if sum(colors{i} == [1 1 1]) == 3  % if white circles, give them an outline
%         scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',[.5 .5 .5],'markerfacealpha',1);  %grey outline
        scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',colors{i+1},'markerfacealpha',1);
    else
        scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75, 'filled','markerfaceColor', colors{i},'markerfacealpha',alpha(i));
    end
end
set(gca,'Fontsize',18,'linewidth',2)
xlabel(xlab,'Fontsize',18)
ylabel(ylab,'Fontsize',18)
xax = get(gca,'XLim');
yax = get(gca,'YLim');

if prod(xax)<0 && abs(prod(xax))<1 && prod(yax)<0 && abs(prod(yax))<1
    xlim([-1 1])
    ylim([-1 1])
    xax = get(gca,'XLim');
    yax = get(gca,'YLim');
else
    xlim([min(xax(1),yax(1)), max(xax(2),yax(2))])
    ylim([min(xax(1),yax(1)), max(xax(2),yax(2))])
end
x = linspace(min(xax(1),yax(1)),max(xax(2),yax(2)),100);
y=x;
plot(x,y,'k--','color','k');
line([min(xax(1),yax(1)) max(xax(2),yax(2))],[0 0],'color','k')
line([0 0], [min(xax(1),yax(1)) max(xax(2),yax(2))],'color','k')
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
%             plot(fittedX,fittedY(i,:),'color', colors{i},'linewidth',2)
            scatter(fittedX,fittedY(i,:),20,'o','filled','markerfacecolor',colors{i},'markerfacealpha',alpha(i))
        end
    end

else    % plot median
    for i = 1:length(vars)
        plot(median(data(color_var==vars(i),1)),median(data(color_var==vars(i),2)), 'marker','+', 'Color', colors{i},'markersize',28,'linewidth',4);
    end
end
if ~isempty(leg)
    l=legend(leg,'location','best');
    set(l,'fontsize',12)
end
print(f, '-dpng',title)
print(f,'-painters','-depsc',title)
end

function plot_hist(data, color_var, colors, edgecolors, alpha, binW, norm, xlab, ylab, title)  % added 7/4/21
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
    histogram(data(color_var==vars(i)), 'binwidth', binW, 'normalization', norm, 'facecolor', colors{i}, 'facealpha', alpha(i), 'edgecolor', edgecolors{i})
end
set(gca,'Fontsize',18,'linewidth',2)
xlabel(xlab,'Fontsize',18)
ylabel(ylab,'Fontsize',18)
xax = get(gca,'XLim');
xlim([0 max(xax)+binW])
print(f, '-dpng',title)
print(f,'-painters','-depsc',title)
end

% function plot_avPSTH(timevec,psth,psth_se,patch_times,color_mat,patch_color,save_fig_name)
% % timevec = starttime:binsize:endtime-binsize
% % psth = average normalized psth across units - num_lcs x timebins (row 1 =
% % no light condition)
% % psth_se = same size as psth, but standard error across units
% % patch_times
% % color_mat = each row is a color (e.g., [1 1 1; 0 0 0]
% % patch_color = assign color of LED
% 
% figure;
% binsize = unique(diff(timevec));
% xlim([timevec(1)+binsize timevec(end)+binsize])   % add binsize so that ticks mark the END of 25ms bins
% hold on;
% shadedErrorBar(timevec,psth(1,:), psth_se(1,:),{'Color',color_mat(1,:),'linewidth',2},1);
% for i = 1:num_lcs
%     shadedErrorBar(timevec,norm_mean(i+1,:), norm_se(i+1,:),{'Color',color_mat(i,:),'linewidth',2},1);
% end
% yax = get(gca,'YLim');
% line([0 0],yax,'Color','k','LineStyle','--','linewidth',2)
% if size(patch_times,1) > 1  
%     for ii = 1:size(patch_times,1)
%         line([patch_times(ii,1) patch_times(ii,1)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
%         line([patch_times(ii,2) patch_times(ii,2)],yax,'Color',color_mat(ii+1,:),'LineStyle','--','linewidth',2)
%     end
% else
%     patch([patch_times(1,1) patch_times(1,1) patch_times(1,2) patch_times(1,2) patch_times(1,1)],[yax(1) yax(2) yax(2) yax(1) yax(1)], patch_color, 'LineStyle', 'none', 'FaceAlpha',.15 );
% end
% ylim(yax)
% xlabel('Time from visual stim onset (s)','fontsize',24)
% ylabel('Normalized firing rate (spks/s)','fontsize',24)
% set(gca,'fontsize',16,'linewidth',2);
% print(pop_fig,'-dpng',save_fig_name)
% print(pop_fig,'-painters','-depsc',save_fig_name)
% 
% end