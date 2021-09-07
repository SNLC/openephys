function barplot_by_layer(data2plot,layers,clusttype,ylabel,title,color_mat)

% data2plot should be an  m by n matrix of what you want to plot (e.g.
% OSI) where m is the number of conditions and n is the number of units
% layers should be n-length vector indicating which layer each unit belongs to
% clusttype is an n-length vector specifying whether unit is single unit
% (1) or multiunit (2)
% ylabel is a string name of the dependent variable
% color_mat - mx3 matrix of color specifications for each condition - changed
% from mani 6/13/20 so that the use can flexibly determine which colors to use

chs_23 = find(layers == 2.5);
chs_4 = find(layers == 4);
chs_5 = find(layers == 5);
chs_55 = find(layers == 5.5);
chs_6 = find(layers == 6);

data_23 = data2plot(:,chs_23);
data_4 = data2plot(:,chs_4);
data_5 = data2plot(:,chs_5);
data_55 = data2plot(:,chs_55);
data_6 = data2plot(:,chs_6);

% % means and SE:
% mean_data(1,:) = nanmean(data_23,2);
% mean_data(2,:) = nanmean(data_4,2);
% mean_data(3,:) = nanmean(data_5,2);
% mean_data(4,:) = nanmean(data_55,2);
% mean_data(5,:) = nanmean(data_6,2);
% se_data(1,:) = nanstd(data_23,0,2)/sqrt(size(data_23,2));
% se_data(2,:) = nanstd(data_4,0,2)/sqrt(size(data_4,2));
% se_data(3,:) = nanstd(data_5,0,2)/sqrt(size(data_5,2));
% se_data(4,:) = nanstd(data_55,0,2)/sqrt(size(data_55,2));
% se_data(5,:) = nanstd(data_6,0,2)/sqrt(size(data_6,2));

% medians and IQR:
mean_data(1,:) = nanmedian(data_23,2);
mean_data(2,:) = nanmedian(data_4,2);
mean_data(3,:) = nanmedian(data_5,2);
mean_data(4,:) = nanmedian(data_55,2);
mean_data(5,:) = nanmedian(data_6,2);
se_data(1,:) = iqr(data_23,2)/2;
se_data(2,:) = iqr(data_4,2)/2;
se_data(3,:) = iqr(data_5,2)/2;
se_data(4,:) = iqr(data_55,2)/2;
se_data(5,:) = iqr(data_6,2)/2;

fig = figure;
h = bar(mean_data);

% set(lh,'Location','BestOutside','Orientation','horizontal')

hold on;

numbars = size(mean_data, 1);
numconds = size(mean_data, 2);

groupwidth = min(0.8, numconds/(numconds+1.5));

for i = 1:numconds
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x(i,:) = (1:numbars) - groupwidth/2 + (2*i-1) * groupwidth / (2*numconds); % Aligning error bar with individual bar
errorbar(x(i,:), mean_data(:,i),se_data(:,i), se_data(:,i), 'color',[.2 .2 .2], 'linestyle', 'none','linewidth',5);
end

set(get(gca,'YLabel'),'String',ylabel,'Fontsize',16)
% set(gca,'XTicklabel','L2/3| L4| L5A| L5B| L6','Fontsize',24)
set(gca,'XTicklabel',{'L2/3', 'L4', 'L5A', 'L5B', 'L6'},'Fontsize',16)

% set(h(1),'FaceColor','k','EdgeColor','k');
% if length(h) > 1
%     set(h(2),'FaceColor','b','EdgeColor','b');
% end

% % color_mat = [0 0 0; 0 0 1; 0 .8 1; 0 0.5 .4; 0 .7 .2]; % for graphing purposes (first is black, last is green)
% if contains(mani,'halo','ignorecase',1)
% %     color_mat = [.9 0 .3; 0.50, 0.0780, 0.10]; % for halo (red)
%     color_mat = [ 0.6350, 0.0780, 0.1840];
% else
%     color_mat = [0 .8 1; 0 0 1; 0 0.5 .4];   % first is light blue, last is green
% end

for i = 1:length(h)
    set(h(i),'FaceColor',color_mat(i,:),'EdgeColor',color_mat(i,:));
end

hold on;
layer_ind = zeros(numconds,length(layers));
for i = 1:numconds
    layer_ind(i,chs_23) = x(i,1);
    layer_ind(i,chs_4) = x(i,2);
    layer_ind(i,chs_5) = x(i,3);
    layer_ind(i,chs_55) = x(i,4);
    layer_ind(i,chs_6) = x(i,5);
end
SUs = find(clusttype==1);
MUs = find(clusttype==2);

if length(h) > 1
    plot(layer_ind(:,SUs),data2plot(:,SUs),'.','color',[.5 .5 .5],'MarkerSize',6)
%     plot(layer_ind(:,SUs),data2plot(:,SUs),'o','color',[.5 .5 .5],'MarkerSize',3)
    plot(layer_ind(:,MUs),data2plot(:,MUs),'rs--','MarkerSize',24)
else
    plot(layer_ind(:,SUs),data2plot(:,SUs),'.','color',[.5 .5 .5],'MarkerSize',6)
%     plot(layer_ind(:,SUs),data2plot(:,SUs),'o','color',[.5 .5 .5],'MarkerSize',3)
    plot(layer_ind(:,MUs),data2plot(:,MUs),'rs','MarkerSize',24)
end

% xax = get(gca,'xlim');
xlim([min(h.XData)-.6 max(h.XData)+.6])
if max(data2plot(:))<=1 
    if min(data2plot(:))<0
        ylim([-1 1])
    else
        ylim([0 1])
    end
end
% line(xax,[-1/3 -1/3],'linestyle','--','color','k')
% line(xax,[1/3 1/3],'linestyle','--','color','k')

% paired ttest
% laynums = unique(layers);
% for l = 1:length(laynums)
%     [h(l),p(l)] = ttest(data2plot(1,find(layers==laynums(l))),data2plot(2,find(layers==laynums(l))));
% end

% NEW (6/14/19) - rotate plot so shows layers top-bottom on y-axis
view(90,90)
set(gca,'fontsize',18,'linewidth',2)

print(fig, '-dpng',sprintf('%s%s',title,'_bylayer'))
print2eps(sprintf('%s%s',title,'_bylayer'),fig)
% close all

return
