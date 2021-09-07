function boxplot_by_layer(data2plot,layers,clusttype,ylabel,title,color_mat)

% 12/12/20 - modified from barplot_by_layer to show median and 25th and 7th
% quartiles by layer, in addition to all data points (MAK)

% data2plot should be an  m by n matrix of what you want to plot (e.g.
% OSI) where m is the number of conditions and n is the number of units
% layers should be n-length vector indicating which layer each unit belongs to
% clusttype is an n-length vector specifying whether unit is of a particular class (e.g., single unit
% (1) or multiunit (2); inactivated or other)
% ylabel is a string name of the dependent variable
% color_mat - mx3 matrix of color specifications for each condition - changed
% from mani 6/13/20 so that the use can flexibly determine which colors to use

chs_23 = find(layers == 2.5);
chs_4 = find(layers == 4);
chs_5 = find(layers == 5);
chs_55 = find(layers == 5.5);
chs_6 = find(layers == 6);

data{1} = data2plot(:,chs_23);
data{2} = data2plot(:,chs_4);
data{3} = data2plot(:,chs_5);
data{4} = data2plot(:,chs_55);
data{5} = data2plot(:,chs_6);

fig = figure;
hold on;

numbars = length(data);
numconds = size(data2plot, 1);

groupwidth = min(0.8, numconds/(numconds+1.5));

for i = 1:numconds
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x(i,:) = (1:numbars) - groupwidth/2 + (2*i-1) * groupwidth / (2*numconds); % Aligning error bar with individual bar
    for ii=1:numbars
%         line([x(i,ii) x(i,ii)], quantile(data{ii}(i,:),[.25 .75]),  'color',[.2 .2 .2],'linewidth',5);  % MAK change 12/12/20 - plot actually 25th and 75th percentiles, rather than just iqr centered on median
        line([x(i,ii)-groupwidth x(i,ii)+groupwidth], repmat(quantile(data{ii}(i,:),[.25]),1,2),  'color',color_mat(i,:),'linewidth',2.5);
        line([x(i,ii)-groupwidth x(i,ii)-groupwidth], quantile(data{ii}(i,:),[.25 .75]),  'color',color_mat(i,:),'linewidth',2.5);
        line([x(i,ii)+groupwidth x(i,ii)+groupwidth], quantile(data{ii}(i,:),[.25 .75]),  'color',color_mat(i,:),'linewidth',2.5);
        line([x(i,ii)-groupwidth x(i,ii)+groupwidth], repmat(quantile(data{ii}(i,:),[.75]),1,2),  'color',color_mat(i,:),'linewidth',2.5);
        line([x(i,ii)-groupwidth x(i,ii)+groupwidth], repmat(nanmedian(data{ii}(i,:)),1,2),  'color',color_mat(i,:),'linewidth',5);
    end
end

set(get(gca,'YLabel'),'String',ylabel,'Fontsize',16)
set(gca,'XTicklabel',{'L2/3', 'L4', 'L5A', 'L5B', 'L6'},'Fontsize',16)

hold on;
layer_ind = zeros(numconds,length(layers));
for i = 1:numconds
    layer_ind(i,chs_23) = x(i,1);
    layer_ind(i,chs_4) = x(i,2);
    layer_ind(i,chs_5) = x(i,3);
    layer_ind(i,chs_55) = x(i,4);
    layer_ind(i,chs_6) = x(i,5);
end
all_layers = unique(layer_ind);

types = unique(clusttype);
newcolormat = [.5 .5 .5; color_mat];
for n = 1:length(all_layers)
    layer_units = find(layer_ind(1,:)==all_layers(n));
    plotvals = linspace(all_layers(n)-groupwidth/4,all_layers(n)+groupwidth/4,length(layer_units));
    plotvals = plotvals(randperm(length(layer_units)));
    for i=1:length(types)
        plot(plotvals(clusttype(layer_units)==types(i)),data2plot(:,layer_units(clusttype(layer_units)==types(i))),'.','color',newcolormat(i,:),'MarkerSize',8)
    end
%         plot(layer_ind(:,clusttype==types(i)),data2plot(:,clusttype==types(i)),'.','color',newcolormat(i,:),'MarkerSize',6)
    %     plot(layer_ind(:,SUs),data2plot(:,SUs),'o','color',[.5 .5 .5],'MarkerSize',3)
end

xlim([min(x(1,1))-.6 max(x(end,end))+.6])
xax = get(gca,'xlim');
if max(data2plot(:))<=1 
    if min(data2plot(:))<0
        ylim([-1 1])
    else
        ylim([0 1])
    end
end
line(xax,[0 0],'color','k','linestyle','--','linewidth',2)
% line(xax,[1/3 1/3],'linestyle','--','color','k')

% rotate plot so shows layers top-bottom on y-axis
view(90,90)
box on
set(gca,'fontsize',18,'linewidth',2)

print(fig, '-dpng',sprintf('%s%s',title,'_bylayer'))
print2eps(sprintf('%s%s',title,'_bylayer'),fig)
% close all

return