function RFpopulation_bycond_analysis_V1(pop,mani,conds)

% pop = string indicating which population of interest (e.g., 'L5',
% 'L6', etc.)
% exp_type = e.g., "halo" or "gtacr"
% conds = 1xlength(exps) cell array of conds to use e.g., {1:4, 1:4, [1 2 5 6] }

fig_dir = sprintf('%s\\%s_%s_figs','H:\LPproject\LPresults\InactFigures_V1',pop,mani);
main_dir = 'H:\LPproject\LPresults';

% load important outputs of population_analysis_inact script
% 'good_units','exps','exp_num','lightmod','FRev','supp_cells','enh_cells'
load(fullfile(fig_dir,'results.mat'))

% load RFmap_analysis results
if length(conds)>1 && length(conds) ~= length(exps)
    warning('fix conds')
end

% load RF reference points from corresponding V1 experiments
s=dir(main_dir);
files = {s.name};
for i=1:length(exps)
    tmp = regexp(exps{i},'_','split');
    animal{i} = tmp{1};
    anV1_files = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,'V1','ignorecase',1),files));
    load(fullfile(main_dir,anV1_files{1},'RFresults.mat'),'rfMap');
    exp(i) = load(fullfile(main_dir,anV1_files{1},'RFresults.mat'),'good_units');
    load(fullfile(main_dir,anV1_files{1},'RFmap_analysis','RFs.mat'));
    all_rfs{i} = rfMap;
    RFunits{i} = good_RFunits;
%     RFs_off{i} = bw_off;
%     RFs_on{i} = bw_on;
    if exist('area','var')   % for experiments that  have RF info
        if length(conds) > 1
            all_area{i} = area(:,conds{i});
        %     all_overlap{i} = overlap;
            all_peakT{i} = peakT(:,conds{i});
        else
            all_area{i} = area;
            all_peakT{i} = peakT;
        end
       clear anV1_files rfMap good_RFunits area peakT
    else
       error('Need to rerun RFmap_analysis for %s',exps{i})
    end
        
end
clear area overlap peakT

count=1;
for i=1:length(exps)
    use_units{i} = find(ismember(RFunits{i},good_units(exp_num==i))); % indexes of RFunits to use
    area(count:count-1+length(use_units{i}),:) = all_area{i}(use_units{i},:);
%     overlap(count:count-1+length(use_units{i}),:) = all_overlap{i}(use_units{i},:);
    peakT(count:count-1+length(use_units{i}),:) = all_peakT{i}(use_units{i},:);
    supp(count:count-1+length(use_units{i})) = ismember(RFunits{i}(use_units{i}),good_units(intersect(find(exp_num==i),supp_cells)));
    count = count+length(use_units{i});
end

% peakT(area==0) = nan;
% area(area==0) = nan;
% badRFs = [sum(isnan(area(:,1:2:end)),2) sum(isnan(area(:,2:2:end)),2)];
use_OFF = sum(area(:,[1 3]),2)~=0;
use_ON = sum(area(:,[2 4]),2)~=0; % changed from [1 3] on 7/1/21

%% define colors for plotting (this section fixed 7/1/21)
finalfig_dir = 'H:\LPproject\LPresults\InactFigures\RFstuff';
cd(finalfig_dir)
% set up colors
% currently: colors according to population being inactivated (blue = L5,
% green = L6, SC = pink, V1 = teal
if strcmpi(pop,'drivermod')
    color_mat = [0 .2 .9; .166 .674 .188; .083 .567 .600]; % blue, green, teal (old blue: 0 .447 .741;)
elseif strcmpi(pop,'L5') 
    color_mat = [0 .2 .9];  % blue
elseif strcmpi(pop,'L6')
    color_mat = [.166 .674 .188]; % green
elseif strcmpi(pop,'V1')
    color_mat = [.083 .567 .6; 1 .6 1; .5 .52 .875]; % teal, pink, purple
elseif strcmpi(pop,'SC')
    color_mat = [1 .6 1;]; % SC = pink
else % anything else (e.g., ctrl experiments)
    color_mat = [.5 .5 .5];  % grey
end

%     figure; plot(area(use_ON,c),area(use_ON,c+2),'o')
plot_scatter([area(supp+use_ON'==2,2),area(supp+use_ON'==2,4)],ones(1,sum(supp+use_ON'==2)),{ color_mat},[1 1],'ON area - light OFF' , 'ON area - light ON','',{'inact'},1)
signrank(area(supp+use_ON'==2,2),area(supp+use_ON'==2,4))
plot_scatter([area(supp+use_OFF'==2,1),area(supp+use_OFF'==2,3)],ones(1,sum(supp+use_OFF'==2)),{ color_mat},[1 1],'OFF area - light OFF' , 'OFF area - light ON','',{'inact'},1)
signrank(area(supp+use_OFF'==2,1),area(supp+use_OFF'==2,3))
plot_scatter([peakT(supp+use_ON'==2,2),peakT(supp+use_ON'==2,4)],ones(1,sum(supp+use_ON'==2)),{ color_mat},[1 1],'ON peakT - light OFF' , 'ON peakT - light ON','',{'inact'},1)
signrank(peakT(supp+use_ON'==2,2),peakT(supp+use_ON'==2,4))
plot_scatter([peakT(supp+use_OFF'==2,1),peakT(supp+use_OFF'==2,3)],ones(1,sum(supp+use_OFF'==2)),{ color_mat},[1 1],'OFF peakT - light OFF' , 'OFF peakT - light ON','',{'inact'},1)
signrank(peakT(supp+use_OFF'==2,1),peakT(supp+use_OFF'==2,3))

save(sprintf('%s_%s_RFstuff.mat',pop,mani),'area','peakT', 'use_ON', 'use_OFF', 'supp')


end

function plot_scatter(data, color_var, colors, alpha, xlab, ylab, title, leg, lobf)
% if lobf = 1, make one lobf regardless of variables; if lobf = 2, make
% separate lobfs for each variable
vars = unique(color_var);
f=figure;
hold on
for i = 1:length(vars)
%     scatter(data(color_var==vars(i),1), data(color_var==vars(i),2), 75, 'filled','markerfaceColor', colors{i},'markerfacealpha',alpha(i));
    [uxy, ~, idx] = unique(data(color_var==vars(i),:),'rows');
    szscale = histc(idx,unique(idx));
    scatter(uxy(:,1),uxy(:,2),75, 'filled','markerfaceColor', colors{i},'markerfacealpha',alpha(i),'sizedata',szscale*25) ; % NEW - make markers bigger if multiple are on top of each other (i.e. identical xy vals)
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
% print(f, '-dpng',title)
% print(f,'-painters','-depsc',title)
end