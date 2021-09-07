function RFpopulation_analysis(region,pop,exp_type)

% region = e.g., 'LPlateral', 'LPmedial', 'LGN'
% pop = string indicating which population of interest (e.g., 'driver',
% 'modulator', etc.)
% exp_type = e.g., "step_halo" or "step_gtacr"

fig_dir = sprintf('%s\\%s_%s_%s_figs','H:\LPproject\LPresults\InactFigures',region,pop,exp_type);
main_dir = 'H:\LPproject\LPresults';

% load important outputs of population_analysis_inact script
% 'good_units','exps','exp_num','lightmod','FRev','supp_cells','enh_cells'
load(fullfile(fig_dir,'unit_info.mat'))

% load RFmap_analysis results


% load RF reference points from corresponding V1 experiments
s=dir(main_dir);
files = {s.name};
for i=1:length(exps)
    tmp = regexp(exps{i},'_','split');
    animal{i} = tmp{1};
    anV1_files = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,'V1','ignorecase',1),files));
    anLP_file = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,'LP','ignorecase',1),files));
    if length(anLP_file)>1
%         anLP_file = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,region,'ignorecase',1),files)); %in case of separate LPlateral and LPrm exps
        if contains(region,'LGN','ignorecase',1) || contains(region,'LPlateral','ignorecase',1)
            anLP_file = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,'LPl','ignorecase',1),files));
        elseif contains(region,'LPmedial','ignorecase',1)
            anLP_file = files(cellfun(@(x) contains(x,animal{i},'ignorecase',1)&contains(x,'LPrm','ignorecase',1),files));
        elseif sum(contains(anLP_file,tmp{2})  )    % if using all LP and >1 LP experiment
            anLP_file = {anLP_file{contains(anLP_file,tmp{2})}};
        end
    end
    fprintf('Loading RF results for %s \n',anLP_file{1})
    load(fullfile(main_dir,anLP_file{1},'RFresults.mat'),'rfMap');
    [nX_LP nY_LP] = size(rfMap{1}{1}); 
    clear rfMap
    vcount = 1;
    for ii = 1:length(anV1_files)
        load(fullfile(main_dir,anV1_files{ii},'RFmap_analysis','RFref.mat'));
        load(fullfile(main_dir,anV1_files{ii},'RFresults.mat'),'rfMap');
        [nX_V1 nY_V1] = size(rfMap{1}{1}); 
%         V1RFs{i}(ii,:) = mean([ref_off;ref_on]);      % get mean of off and on reference points
        for iii=1:size(ref_on,1)
            V1RFs{i}(vcount,:) = nanmean([ref_off(iii,:);ref_on(iii,:)]) + [.5*8*diff([nX_V1,nX_LP]) .5*8*diff([nY_V1,nY_LP])]; % NEW - this corrects for possible differences between V1 & LP sparsenoise exps. if nX_LP>nx_V1, add onto V1 ref coords (subtract if nX_V1>nX_LP)
            vcount = vcount+1;
        end
        clear rfMap nX_V1 nY_V1
    end
    load(fullfile(main_dir,anLP_file{1},'RFresults.mat'),'rfMap');
    exp(i) = load(fullfile(main_dir,anLP_file{1},'RFresults.mat'),'good_units');
    load(fullfile(main_dir,anLP_file{1},'RFmap_analysis','RFs.mat'));
    all_rfs{i} = rfMap;
    RFunits{i} = good_RFunits;
    RFs_off{i} = bw_off;
    RFs_on{i} = bw_on;
    all_area{i} = area(:,1:2);
        
end

for i=1:length(exps)
    use_units{i} = find(ismember(RFunits{i},good_units(exp_num==i))); % indexes of RFunits to use
    dis_shuf{i} = nan(1000,length(use_units{i}));
    max_area{i} = max(all_area{i}(use_units{i},:),[],2)'; % max b/w ON and OFF subfields
    for n =1:length(use_units{i})
        unit = use_units{i}(n); % index from 1:length(RFunits{i})
        offon = [~isempty(RFs_off{i}{unit}) ~isempty(RFs_on{i}{unit})];
        for ii=1:size(V1RFs{i},1)
            if offon(1)==1 && offon(2) ==0  % only sig off subfield
                diff_dis{i}{n}(ii) = norm(V1RFs{i}(ii,:)-[mean(RFs_off{i}{unit}(:,1)) mean(RFs_off{i}{unit}(:,2))]);
            elseif offon(1)==0 && offon(2)==1   % only sig on subfield
                diff_dis{i}{n}(ii) = norm(V1RFs{i}(ii,:)-[mean(RFs_on{i}{unit}(:,1)) mean(RFs_on{i}{unit}(:,2))]);
            elseif sum(offon)==2 % both sig
                diff_dis{i}{n}(ii) = min(norm(V1RFs{i}(ii,:)-[mean(RFs_off{i}{unit}(:,1)) mean(RFs_off{i}{unit}(:,2))]),norm(V1RFs{i}(ii,:)-[mean(RFs_on{i}{unit}(:,1)) mean(RFs_on{i}{unit}(:,2))]));
            else % neither sig
                diff_dis{i}{n}(ii) = nan;
            end
        end
        dis{i}(n) = min(diff_dis{i}{n});
        if isnan(dis{i}(n))
            test=1;
        end
        mod{i}(n) = lightmod(good_units==RFunits{i}(unit)&exp_num==i);
        rfs{i}{n} = all_rfs{i}{exp(i).good_units==RFunits{i}(unit)};
    end
    % make shuffled distribution
    if ~isempty(use_units{i})
        for ii=1:1000
            dis_shuf{i}(ii,:) = dis{i}(randperm(length(use_units{i})));
        end
    else
        dis{i} = [];
        mod{i} = [];
    end
end

figure; 
hold on
for i=1:length(exps)
    plot(dis{i},mod{i},'.','markersize',24)
end
xlim([0 100])
ylim([-1 1])
xax = get(gca,'XLim');
line(xax,[0 0],'color','k', 'linestyle','--','linewidth',2)
% legend(animal,'location','bestoutside')
xlabel('Distance from V1 RF center (in degrees)')
ylabel('Light modulation index')
% title(sprintf('%s\\_%s',region,pop))
cd(fig_dir)
print(gcf,'-dpng','RFdistbyLightmod')

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
    color_mat = [1 .6 1;]; % SC = pink
else % anything else (e.g., ctrl experiments)
    color_mat = [.5 .5 .5];  % grey
end

all_dis = [dis{:}];
all_mod = [mod{:}];
supp_units = good_units(supp_cells);
supp_exps = exp_num(supp_cells);
enh_units = good_units(enh_cells);
enh_exps = exp_num(enh_cells);
for i=1:length(exps)
    supp{i} = ismember(RFunits{i}(use_units{i}),supp_units(supp_exps==i))'; % index to use_units{i}
    enh{i} = ismember(RFunits{i}(use_units{i}),enh_units(enh_exps==i))';
end
all_supp = [supp{:}];
all_enh = [enh{:}];
edges = [0:10:round(max(all_dis))+10];
disbins = histcounts(all_dis,edges);
suppbins = histcounts(all_dis(all_supp),edges);
enhbins = histcounts(all_dis(all_enh),edges);
propsupp = suppbins./disbins;
propenh = enhbins./disbins;
%
for ii=1:1000
    for i=1:length(exps)
        suppshuf(i,:,ii) = histcounts(dis_shuf{i}(ii,supp{i}),edges);
        enhshuf(i,:,ii) = histcounts(dis_shuf{i}(ii,enh{i}),edges);
    end
end
distsupp_shuf = squeeze(sum(suppshuf))./repmat(disbins',1,1000);
distenh_shuf = squeeze(sum(enhshuf))./repmat(disbins',1,1000);
propsupp_shuf = nanmean(distsupp_shuf,2);
sesupp_shuf = std(distsupp_shuf,[],2)./sqrt(1000);
propenh_shuf = nanmean(distenh_shuf,2);
seenh_shuf = std(distenh_shuf,[],2)./sqrt(1000);
% critval = .05/(length(edges)-1); % or should it be .025 b/c 2-tailed?
critval = .025;
p_supp = nan(length(edges)-1,2);
p_enh = p_supp;
for x=1:length(edges)-1
    if disbins(x)
        if sum(suppbins) % only calc significance is there were suppressed units
            p_supp(x,:) = [(sum(distsupp_shuf(x,:)>=propsupp(x)))/1000 sum(distsupp_shuf(x,:)<=propsupp(x))/1000]; % has to be <= or else values of zero will always be significant...? 
        end
        if sum(enhbins) % only calc significance is there were enhanced units
            p_enh(x,:) = [(sum(distenh_shuf(x,:)>=propenh(x)))/1000 (sum(distenh_shuf(x,:)<=propenh(x)))/1000];
        end
    end
end
sig_supp = min(p_supp,[],2)<critval;
sig_enh = min(p_enh,[],2)<critval;

% also looks at how many VISUAL cells without significant RFs were
% suppressed/enhanced
for i=1:length(exps)
    exp_units{i} = good_units(intersect(find(exp_num==i),visual_cells));
    noRF_units{i} = exp_units{i}(~ismember(exp_units{i},RFunits{i}));
    supp_noRF{i} = noRF_units{i}(ismember(noRF_units{i},good_units(intersect(intersect(visual_cells,supp_cells),find(exp_num==i)))));
    enh_noRF{i} = noRF_units{i}(ismember(noRF_units{i},good_units(intersect(intersect(visual_cells,enh_cells),find(exp_num==i)))));
end
propsupp_noRF = length([supp_noRF{:}])/length([noRF_units{:}]);
propenh_noRF = length([enh_noRF{:}])/length([noRF_units{:}]);


figure; hold on;
% shadedErrorBar(edges(1:end-1)+5,-propsupp_shuf, sesupp_shuf)
% shadedErrorBar(edges(1:end-1)+5,propenh_shuf,seenh_shuf)
b1 = bar(edges(1:end-1)+5,-propsupp_shuf,1,'FaceColor',[1 1 1],'edgecolor',[.75 .75 .75],'facealpha',1,'linewidth',1.25);
bar(edges(1:end-1)+5,propenh_shuf,1,'FaceColor',[1 1 1],'edgecolor',[.75 .75 .75],'facealpha',1,'linewidth',1.25)
bar(edges(1:end-1)+5,-propsupp,1,'facecolor',color_mat(1,:),'edgecolor','k','facealpha',.5,'linewidth',1.25)
bar(edges(1:end-1)+5,propenh,1,'facecolor',color_mat(1,:),'edgecolor','k','facealpha',.1,'linewidth',1.25)
s1 = scatter(all_dis,all_mod,75,'filled','markerfacecolor',[1 1 1],'markeredgecolor',color_mat(1,:),'markerfacealpha',1);
bar(110,-propsupp_noRF,10,'facecolor',color_mat(1,:),'edgecolor','k','facealpha',.6,'linewidth',1.25)
bar(110,propenh_noRF,10,'facecolor',color_mat(1,:),'edgecolor','k','facealpha',.1,'linewidth',1.25)
hold on;
s2=scatter(all_dis(all_supp),all_mod(all_supp),75,'filled','markerfacecolor',color_mat(1,:),'markerfacealpha',1,'markeredgecolor',[.1 .1 .1]);
scatter(edges(sig_supp)+5,-sig_supp(sig_supp==1)+.05,'*','markeredgecolor',color_mat(1,:),'markeredgealpha',1)
s3=scatter(all_dis(all_enh),all_mod(all_enh),75,'filled','markerfacecolor',color_mat(1,:),'markerfacealpha',.4,'markeredgecolor',[.1 .1 .1]);
scatter(edges(sig_enh)+5,sig_enh(sig_enh==1)-.05,'*','markeredgecolor',color_mat(1,:),'markeredgealpha',.4)
xlim([0 125])
ylim([-1 1])
xax = xlim;
yax=ylim;
set(gca,'linewidth',2)
line([100 100],yax,'color','k','linewidth',1)
line(xax,[0 0],'color','k', 'linestyle','--','linewidth',2)
xticks(0:10:110)
ticklabs = xticklabels;
xticklabels({ticklabs{1:end-2} '','Other Vis'})
% plot(edges(1:end-1)+5,-propsupp,'color',color_mat(1,:),'linewidth',2)
% plot(edges(1:end-1)+5,propenh,'color',color_mat(1,:),'linestyle','--','linewidth',2)
xlabel('Distance from V1 RF center (in degrees)','fontsize',15)
ylabel('% suppressed/enhanced','fontsize',15)
yticklabels(num2str(abs(yticks')*100))
a2 = axes('YAxisLocation', 'Right');  % Create second Y axes on the right.
set(a2, 'color', 'none')  % Hide second plot.
set(a2, 'XTick', [])
set(a2,'XColor',[1 1 1])
set(a2, 'YLim', [-100 100])
set(a2,'YTickLabels',(num2str(yticks'/100)))
set(a2,'ycolor',[.15 .15 .15])
set(gca,'linewidth',2)
ylabel('LMI_v_i_s','fontsize',15)
l=legend([s1 s3 s2 b1],{'other','enhanced','suppressed','shuffle'},'location','northwest');
print(gcf,'-dpng','RFdistbyLightmod_bytype')
print(gcf,'-painters','-depsc','RFdistbyLightmod_bytype')

% % eventually, for looking at correlation coefficients between rfMaps...?:
% coeffs = cell(1,length(exps));
% for i=1:length(exps)
%     if ~isempty(use_units{i})
%         coeffs{i} = nan(length(rfs{i}{1})/2,length(use_units{i}));
%         good_units_black = find(cellfun(@(x) ~isempty(x),RFs_off{i}(use_units{i}))); % look at correlations for all units that passed the first pass (i.e. had peak in their timecourse) in no light and at least one lightcond (separately for white and black squares)
%         good_units_white = find(cellfun(@(x) ~isempty(x),RFs_on{i}(use_units{i}))); % look at correlations for all units that passed the first pass (i.e. had peak in their timecourse) in no light and at least one lightcond  (separately for white and black squares)
%         for n=1:length(good_units_black)
%             nn = good_units_black(n);
%             for ii = 3:2:length(rfs{i}{1})
%                 coeffs{i}(ii-2,nn) = corr2(rfs{i}{nn}{1},rfs{i}{nn}{ii});
%             end
%         end
%         for n=1:length(good_units_white)
%             nn = good_units_white(n);
%             for ii = 4:2:length(rfs{i}{1})
%                 coeffs{i}(ii-2,nn) = corr2(rfs{i}{nn}{2},rfs{i}{nn}{ii});
%             end
%         end
%     end
% end
% all_coeffs = cell2mat(coeffs)';
% figure;
% plot(all_coeffs(:,1)',all_mod,'o')
% ylim([-1 1])

diams = [967.56 1425 854.53 1322.64 995.12 1240.24 1130.24 1137.34]; % NPRH18 21 19 23 24 5 27 37
for i=1:length(exps)
    if sum(supp{i})
        max_dist(i) = max(dis{i}(supp{i}));
    else
        max_dist(i) = nan;
    end
end
figure; plot(diams,max_dist,'o')
end

