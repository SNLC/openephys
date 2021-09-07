function plot_unitRF(exp,unit,conds,condlabels,type)
% exp: eg.,'NPRH23_LP'
% unit: unit name (e.g., 198)
% conds: e.g., [2 4] would plot ON subfields in light OFF and ON conditions
% cond labels: e.g., {'LED OFF','LED ON'} (should be same length as conds)
% type: 'psth' or 'filt' (raw baseline-subtracted PSTH or spatial filter
% from SVD)
% alt: add 't' as timebin at which to plot (if using psth_norm instead of
% rfMap)

main_dir = 'H:\LPproject\LPresults';

load(fullfile(main_dir,exp,'RFresults.mat'))
load(fullfile(main_dir,exp,'RFmap_analysis','RFs.mat'))

i=find(good_units==unit);
ii=find(good_RFunits==unit);

numconds = length(stats(i).issig);
numstim = size(psth_norm{i},2)/numconds;
nX = size(rfMap{i}{conds(1)},1);
nY = size(rfMap{i}{conds(1)},2);

if rem(conds,2)
    borders = {all_bw_off{ii,:}};
else
    borders = {all_bw_on{ii,:}};
end

if strcmpi(type,'psth')
    ts = peakT(ii,conds);
    fprintf('%s%s\n','peak timepoints: ',num2str(ts))
    same = input('Do you want to use the same time bin for all conds? Y/N ','s');

    if strcmpi(same,'Y') && length(unique(ts))>1  % if you want to use same timepoint for both:
        t = input(strcat('Which time bin do you want to use? ',num2str(ts),'\n'));
        ts=repmat(t,1,numconds);
    end
end

f=figure('Position', [300, 0, 800*length(conds),600]); 
for n=1:length(conds)
    subplot(1,length(conds),n)
    if strcmpi(type,'psth')
        map = reshape(psth_norm{i}(ts(n),1+numstim*(conds(n)-1):numstim*(conds(n))),nX,nY)';
    elseif strcmpi(type,'zscore')
        map = reshape(zscore(rfMap{i}{conds(n)}(:)),nX,nY)';
    elseif strcmpi(type,'filt')
        map = rfMap{i}{conds(n)}';
    end
%     h(n) = imagesc(map(2:end-1,2:end-1));  % if V1 and LP sparsenoise exps were diff sizes...
    h(n) = imagesc(map);
%     if ~isempty(borders{ceil(conds(n)/2)})
    if ~isempty(borders) &&  ~isempty(borders{ceil(conds(n)/2)})
         hold on;
         plot(borders{ceil(conds(n)/2)}(:,1)/8+.5,borders{ceil(conds(n)/2)}(:,2)/8+.5,'linewidth',4,'color',[.25 .25 .25])
         plot(mean(borders{ceil(conds(n)/2)}(:,1)/8+.5),mean(borders{ceil(conds(n)/2)}(:,2)/8+.5),'.','color',[.25 .25 .25],'markersize',20)
%         plot(borders{n}(:,1)/8-.5,borders{n}(:,2)/8-.5,'linewidth',4,'color',[.25 .25 .25])   % if V1 and LP sparsenoise exps were diff sizes...
%          plot(mean(borders{n}(:,1)/8-.5),mean(borders{n}(:,2)/8-.5),'.','color',[.25 .25 .25],'markersize',20)
    end
    cax(n,:) = caxis();
    xticklabels ''
    yticklabels ''
    title(condlabels{n},'fontsize',16)
end

for n=1:length(conds)
    subplot(1,length(conds),n)
    xax=xlim;
    yax=ylim;
    caxis([min(cax(:,1)) max(cax(:,2))]);
    c(n) = colorbar('fontsize',12);
    c(n).Ticks = c(n).Ticks([1 end]);
%     c(n).TickLabels = {'min','max'};
    xticks(xax(1):2:xax(2)) 
    yticks(yax(1):2:yax(2)) 
    set(gca,'fontsize',18,'linewidth',3);
    set(c(n),'fontsize',18)    
end

cd('H:\LPproject\LPresults\UnitFigures')
print(f,'-dpng',strcat(exp,'_',num2str(unit),'_RF_bycond'))
print(f,'-painters','-depsc',strcat(exp,'_',num2str(unit),'_RF_bycond'))
    
end