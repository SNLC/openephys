function RFmap_analysis_gaussian(exp_path)
% exp_path: eg., 'H:\LPproject\LPresults\NPRSC5_LPl'

cd(exp_path)
load('RFresults.mat')
% load('H:\LPproject\LPresults\MH18_LP\MH18_LP_cluster_quality.mat')
% good_units = intersect(find(SNR>=1.5&refr_idx<1),find(shank==0));
% allsig = reshape([stats.issig],4,length(stats))';
% sig_mat = [allsig(:,1)+allsig(:,3) allsig(:,2)+allsig(:,4)];
% good_units(sum(sig_mat(good_units,:),2)<1) = [];

if ~exist('RFmap_analysis','dir')
    mkdir('RFmap_analysis');
end
cd('RFmap_analysis')

s=dir(exp_path);
files = {s.name};
cq_file = files(cellfun(@(x) contains(x,'cluster_quality.mat','ignorecase',1),files));
load(fullfile(exp_path,cq_file{1}));
if length(uQ) ~= length(good_units)
    error('something wrong - units between driftgrat and sparsenoise exps dont match')
end
notclean_units = find((uQ<=16) | (refV'>=.5));    % only include units with <.5% refractory period violations, unit quality > 16, and those from designated shanks 

numconds = length(stats(1).issig);
numstim = size(psth_norm{1},2)/numconds;
FR = [cellfun(@(x) mean(x,2),cellfun(@(x) mean(abs(x(4:15,1:numstim))),psth_norm,'uniformoutput',0))' cellfun(@(x) mean(x,2),cellfun(@(x) mean(abs(x(4:15,1+numstim:2*numstim))),psth_norm,'uniformoutput',0))'];
notclean_units = union(notclean_units,find(max(FR,[],2)<.25));  % also exclude units with average FR change from baseline <.25spks/s in either ON or OFF stim - no light conds(more likely to turn out as false positive)

sig = reshape([stats.issig],length(stats(1).issig),length(good_units))';  % from timecourse
sig2 = reshape([stats.issig2],length(stats(1).issig2),length(good_units))'; % fron modified SNR
sig3 = reshape([stats.issig3],length(stats(1).issig3),length(good_units))'; % from statistical comparison with null distribution
Zs = reshape([stats.peakZscore],length(stats(1).peakZscore),length(good_units))'; % fron modified SNR
% noRFunits = find(sum(sig2(:,1:2),2)<1 | sum(sig(:,1:2),2)<1); % if neither no light stim (on- or off-subunits) passed both first and second passes
noRFunits = find(sum(sig2.*sig3,2)<size(sig3,2)/4);  % sig2*sig3 so that we're only counting ones that passed both second and third passes
good_units(union(noRFunits,notclean_units)) = [];
stats(union(noRFunits,notclean_units)) = [];
rfMap(union(noRFunits,notclean_units)) = [];
sig2(union(noRFunits,notclean_units),:) = [];
sig(union(noRFunits,notclean_units),:) = [];
sig3(union(noRFunits,notclean_units),:) = [];
Zs(union(noRFunits,notclean_units),:) = [];
shank(union(noRFunits,notclean_units)) = [];
max_ch(union(noRFunits,notclean_units)) = [];

shs = unique(shank);
for sh = 1:length(shs)
    sh_units{sh} = find(shank==shs(sh));
    [~,ord] = sort(max_ch(shank==shs(sh)));
    sh_units{sh} = sh_units{sh}(ord);
end

nX = size(rfMap{1}{1},1);
nY = size(rfMap{1}{1},2);
% if size(sig3,2)>2   % if opto conditions
%     sig_byfield = [sum(sig2(:,[1:2:end])+sig(:,[1:2:end]),2) sum(sig2(:,[2:2:end])+sig(:,[2:2:end]),2)];  % FIX
    sig_byfield = [sum(sig3(:,1:2:end),2) sum(sig3(:,2:2:end),2)]; % [off on]
% else
%     sig_byfield = [max(sig(:,1),sig3(:,1)) max(sig(:,2),sig3(:,2))];
% end

area = nan(length(good_units),size(sig2,2));
angle = area;
xsigma = area;
ysigma = area;
center = nan(length(good_units),size(sig2,2),2);
all_bw_off = cell(1,length(good_units));
all_bw_on = all_bw_off;

% for gaussian fitting:
[xi,yi] = meshgrid(1:8*nX,1:8*nY);
opts.iso = 0;
opts.tilted = 1;
t = linspace(0,2*pi,100);

for n = 1:length(good_units)
    unit = good_units(n);
%     reps = 10000;
%     test = zeros(reps,numel(psth_norm{unit}(:,1:120)));
%     test_psth = zeros(size(psth_norm{unit}(:,1:120),1),size(psth_norm{unit}(:,1:120),2),reps);
%     rfShuf = zeros(12,10,reps);
%     rfShuf_int = zeros(12*8,10*8,reps);
%     for n=1:reps
%         test(n,:) = randperm(numel(psth_norm{unit}(:,1:120)));
%         test_psth(:,:,n) = reshape(psth_norm{unit}(test(n,:)),size(psth_norm{unit},1),120);
%         [Ushuf,~,~] = svd(test_psth(:,:,n)','econ');
%         rfShuf(:,:,n) = reshape(Ushuf(:,1),12,10);
%         rfShuf_int(:,:,n) = imresize(rfShuf(:,:,n),[12*8 10*8],'bilinear');
%     end
%     pvals = zeros(size(RF));
%     for xx = 1:size(RF,1)
%         for yy = 1:size(RF,2)
%             if RF(xx,yy)<0
%                 pvals(xx,yy) = -sum(rfShuf(xx,yy,:)<=RF(xx,yy))/reps;
%             else
%                 pvals(xx,yy) = sum(rfShuf(xx,yy,:)>=RF(xx,yy))/reps;
%             end
%         end
%     end
%     pvals(abs(pvals)>=.05) = nan;   % do I need a multiple comparisons correction?
% 
%     % shuf_psth = mean(test_psth,3);
%     figure;imagesc(1:120,20:10:240,shuf_psth)
%     [U,S,V] = svd(shuf_psth','econ');
%     rfTest = U(:,1);
%     rfTest = reshape(rfTest,12,10);
%     figure;imagesc(1:12,1:10,rfTest')

    f1 = figure;
    for i = 1:length(rfMap{n})
%         if sig_byfield(n,rem(i+1,2)+1)
            RF = rfMap{n}{i}(:,:);
            ZRF = (RF-mean(RF(:)))./std(RF(:));
%             if (sig2(n,i)||sig(n,i)||Zs(n,i)>=6) && sig_byfield(n,rem(i+1,2)+1)>1     % only compare maps for which there was a significant RF in both light on and light off conds (was previously >=length(rfMap{n})/2)
            if sig_byfield(n,rem(i,2)==rem(1:2,2))
                zmap = abs(ZRF)>=2;
                pmap = stats(n).pvals{i};
            else
                zmap = zeros(size(RF));
                pmap = nan(size(RF));
            end
            subplot(length(rfMap{n}),5,1+(i-1)*5)
            imagesc(1:nX,1:nY,ZRF')
            RFnorm = RF./max(abs(RF(:)));
            RFnorm(isnan(pmap)) = 0;
%             RFnorm(~zmap) = 0;
            subplot(length(rfMap{n}),5,2+(i-1)*5)
            imagesc(1:nX,1:nY,RFnorm')
            testRF = imresize(RFnorm,[nX*8 nY*8],'bilinear'); % bilinear interpolation to 1deg resolution
            subplot(length(rfMap{n}),5,3+(i-1)*5)
            imagesc(1:nX,1:nY,testRF')
%             filtRF = imgaussfilt(testRF, 4); % smooth with 2d gaussian filter (2 stdev = 2*(1deg resolution?) = 2degs??)
            % NEW 12/7/20 - fit with 2D gaussian
            if sum(testRF(:))
                [results] = autoGaussianSurf(xi,yi,testRF',opts);
                fitRF = results.G';
                subplot(length(rfMap{n}),5,4+(i-1)*5)
                imagesc(1:nX,1:nY,fitRF')
                grid
        %         filtRF(abs(filtRF(:))<.5) = 0;
                subplot(length(rfMap{n}),5,5+(i-1)*5)
        %         imagesc(1:12,1:10,filtRF')
    %             colorbar
    %             I{i} = imbinarize(abs(filtRF)).*sign(filtRF);
    %             imagesc(1:nX,1:nY,I{i}')
    %             bw{i} = bwboundaries(I{i});
                bw{n}{i}(:,1) = results.x0 + results.sigmax*cos(t)*cos(results.angle) - results.sigmay*sin(t)*sin(results.angle); % x vals of ellipse 
                bw{n}{i}(:,2) = results.y0 + results.sigmay*sin(t)*cos(results.angle) + results.sigmax*cos(t)*sin(results.angle);  % y vals of ellipse
                plot(bw{n}{i}(:,1)./8,bw{n}{i}(:,2)./8,'r','linewidth',2)
                xlim([0 nX])
                ylim([0 nY])
                view(0,270)
                grid
            % get rid of "rfs" not part of the main one
%             if ~isempty(bw{i})
%                 if max(cellfun(@(x) numel(x),bw{i})) < min(cellfun(@(x) numel(x),bw{i}))*1.25 % if possible RFs are similarly sized, chose one with max pixel
%                     [a,b] = find(abs(filtRF)==max(abs(filtRF(:)))); % get max val pixel
%                     for nn=1:length(bw{i})
%                         [p,q] = meshgrid(round(min(bw{i}(:,1))):round(max(bw{i}(:,1))),round(min(bw{i}(:,2))):round(max(bw{i}(:,2))));
% %                         if find(p(:)==a(1)&q(:)==b(1))    % if peak pixel falls within this border
%                             I{i} = zeros(size(fitRF));
%                             I{i}(p(:),q(:)) = 1;
%                             which_bord(i) = nn;
%                         end
%                     end
%                 else  % choose the biggest one
%                     [~,which_bord(i)] = max(cellfun(@(x) numel(x),bw{i}));
%                     [p,q] = meshgrid(min(bw{i}{which_bord(i)}(:,1)):max(bw{i}{which_bord(i)}(:,1)),min(bw{i}{which_bord(i)}(:,2)):max(bw{i}{which_bord(i)}(:,2)));
%                     filt_grid = zeros(size(I{i}));
%                     filt_grid(p(:),q(:)) = 1;
%                 end
%                 [~,which_bord(i)] = max(cellfun(@(x) length(x),bw{i},'uniformoutput',1));      % if multiple borders, only use biggest one
%                 [p,q] = meshgrid(min(bw{i}{which_bord(i)}(:,1)):max(bw{i}{which_bord(i)}(:,1)),min(bw{i}{which_bord(i)}(:,2)):max(bw{i}{which_bord(i)}(:,2)));
%                 filt_grid = zeros(size(I{i}));
%                 filt_grid(p(:),q(:)) = 1;
%                 I{i} = I{i}.*filt_grid;
%                 perim(i) = size(bw{i}{which_bord(i)},1);
                area(n,i)= results.sigmax*results.sigmay*pi;
                angle(n,i) = rad2deg(results.angle);
                xsigma(n,i) = results.sigmax;
                ysigma(n,i) = results.sigmay;
                center(n,i,:) = [results.x0 results.y0];
            else
                bw{n}{i} = [nan nan];
            end
            clear results
            
%             if i>2
%                 overlap(n,i-2) = (sum(abs(I{i-2}(:)).*abs(I{i}(:))))/sum(abs(I{i-2}(:))); % as % of no-light condition
%             end
%         end
    end


    f2 = figure; 
    hold on
    for i=1:length(rfMap{n})
        subplot(2,length(rfMap{n})/2,ceil(i/2))
        if rem(i,2)
            c = 'b-';
            hold on
        else
            c = 'r-';
        end
        if ~isempty(bw{n}{i})
            plot(bw{n}{i}(:,1),bw{n}{i}(:,2),c,'linewidth',4)
        end
        view(0,270)
        xlim([0 size(fitRF,1)])
        ylim([0 size(fitRF,2)])
        title(strcat('Cond ',num2str(ceil(i/2))))
        if i==length(rfMap{n})
            legend({'OFF','ON'},'Location','Best')
        end
        
    end
    
   all_bw_off{n} = bw{n}{1}; 
   all_bw_on{n} = bw{n}{2}; 
    if length(rfMap{n})>2   % if opto exp
        for ii=1:2
            subplot(2,length(rfMap{n})/2,ii+length(rfMap{n})/2)
            for i=ii:2:length(rfMap{n})
                plot(bw{n}{i}(:,1),bw{n}{i}(:,2),'linewidth',4)
                hold on;
                view(0,270)
                xlim([0 size(fitRF,1)])
                ylim([0 size(fitRF,2)])
            end
            if rem(i,2)
                title('OFF subfields')
            else
                title('ON subfields')
            end
            if ii==2
                legend(num2str([1:length(rfMap{n})/2]'),'Location','Best')
            end
        end
    end
    save_name1 = sprintf('Fields_Cluster%d',good_units(n));
    save_name2 = sprintf('Optocomp_Cluster%d',good_units(n));
    print(f1,'-dpng',save_name1)
    print(f2,'-dpng',save_name2)
    close(f1,f2)
   
end
for sh=1:length(shs)
    fON = figure;
%     subplot(1,4,sh)
    chs = min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}));
    colorMap = zeros(length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))),3);
    colorMap(:,3) = linspace(0,1,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))));
    colorMap(:,2) = linspace(1,0,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))));
    colorMap(:,1) = abs(linspace(-.5,.25,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh})))));
    for nn = 1:length(sh_units{sh})
        if ~isempty(all_bw_on{sh_units{sh}(nn)})
            plot(all_bw_on{sh_units{sh}(nn)}(:,1),all_bw_on{sh_units{sh}(nn)}(:,2),'color',colorMap(chs==max_ch(sh_units{sh}(nn)),:),'linewidth',3)
        end
        hold on;
    end
    view(0,270)
    xlim([0 size(fitRF,1)])
    ylim([0 size(fitRF,2)])
    title(strcat('shank ',num2str(shs(sh))))
    set(gca,'fontsize',18,'linewidth',3);
    %     if sh == length(shs)
%         legend(num2str(max_ch(sh_units{sh})'),'Location','Southoutside','Orientation','horizontal')
%         legend(num2str(max_ch(sh_units{sh})'))
        c = colorbar('Ticks',[0 1],'TickLabels', {'ventral','dorsal'});
        colormap(flipud(colorMap))
        
%     else
%         legend('','Location','Southoutside','Orientation','horizontal')
%     end
    save_nameON{sh} = strcat('RFs_by_depth_ONfields_shank',num2str(shs(sh)));
%     set(fON,'Paperposition',[0 0 30 10])
    print(fON,'-dpng',save_nameON{sh})
    close(fON)
end

for sh=1:length(shs)
    fOFF = figure;
%     subplot(1,4,sh)
    chs = min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}));
    colorMap = zeros(length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))),3);
    colorMap(:,3) = linspace(0,1,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))));
    colorMap(:,2) = linspace(1,0,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh}))));
    colorMap(:,1) = abs(linspace(-.5,.25,length(min(max_ch(sh_units{sh})):max(max_ch(sh_units{sh})))));
    for nn = 1:length(sh_units{sh})
        if ~isempty(all_bw_off{sh_units{sh}(nn)})
            plot(all_bw_off{sh_units{sh}(nn)}(:,1),all_bw_off{sh_units{sh}(nn)}(:,2),'color',colorMap(chs==max_ch(sh_units{sh}(nn)),:),'linewidth',4)
        end
        hold on;
    end
    view(0,270)
    xlim([0 size(fitRF,1)])
    ylim([0 size(fitRF,2)])
    title(strcat('shank ',num2str(shs(sh))))
    set(gca,'fontsize',18,'linewidth',3);
%     if sh == length(shs)
%         legend(num2str(max_ch(sh_units{sh})'),'Location','Southoutside','Orientation','horizontal')
%         legend(num2str(max_ch(sh_units{sh})'))
        c = colorbar('Ticks',[0 1],'TickLabels', {'ventral','dorsal'});
        colormap(flipud(colorMap))
        
%     else
%         legend('','Location','Southoutside','Orientation','horizontal')
%     end
    
    save_nameOFF{sh} = strcat('RFs_by_depth_OFFields_shank',num2str(shs(sh)));
% set(fOFF,'Paperposition',[0 0 30 10])
    print(fOFF,'-dpng',save_nameOFF{sh})
%     close(fOFF)
end

if contains(exp_path,'V1','ignorecase',1)
    ref_on = [median(cellfun(@(x) median(x(:,1)),all_bw_on(cellfun(@(x) ~isempty(x),all_bw_on)))) median(cellfun(@(x) median(x(:,2)),all_bw_on(cellfun(@(x) ~isempty(x),all_bw_on))))]; 
    ref_off = [median(cellfun(@(x) median(x(:,1)),all_bw_off(cellfun(@(x) ~isempty(x),all_bw_off)))) median(cellfun(@(x) median(x(:,2)),all_bw_off(cellfun(@(x) ~isempty(x),all_bw_off))))];
    save('RFref.mat','ref_on','ref_off')
end

good_RFunits = good_units;
save('RFs.mat','good_RFunits','all_bw_off','all_bw_on')

end

% eventually, for looking at correlation coefficients between rfMaps:
% ok_units_black = find(sig(:,1)==1 & sum(sig(:,3:2:end),2)>=1); % look at correlations for all units that passed the first pass (i.e. had peak in their timecourse) in no light and at least one lightcond (separately for white and black squares)
% ok_units_white = find(sig(:,2)==1 & sum(sig(:,4:2:end),2)>=1); % look at correlations for all units that passed the first pass (i.e. had peak in their timecourse) in no light and at least one lightcond  (separately for white and black squares)
% good_units_black = find(ismember(1:size(sig,1),ok_units_black)' & sig2(:,1)==1);
% good_units_white = find(ismember(1:size(sig,1),ok_units_white)' & sig2(:,2)==1);
% coeffs = nan(size(sig));
% for n=1:length(good_units_black)
%     nn = good_units_black(n);
%     for i = 3:2:size(sig,2)
%         c = corrcoef(rfMap{nn}{1}(:),rfMap{nn}{i}(:));
%         coeffs(nn,i) = c(1,2);
%     end
% end
% for n=1:length(good_units_white)
%     nn = good_units_white(n);
%     for i = 4:2:size(sig,2)
%         c = corrcoef(rfMap{nn}{2}(:),rfMap{nn}{i}(:));
%         coeffs(nn,i) = c(1,2);
%     end
% end
% coeffs(:,1:2) = [];

% change in peak timepoints?
% [~,peakT_cell] = arrayfun(@(y) cellfun(@(x) max(x),y.timeCourse),stats,'uniformoutput',0);
% peakT = reshape([peakT_cell{:}],4,length(stats))';
