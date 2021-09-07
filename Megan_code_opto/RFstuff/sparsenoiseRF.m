function [psth_norm, rfMap, stats] = sparsenoiseRF(unit_times_ds,stim_times,binsize,T,stim_trials,run_trials,stim_key)
% description here:

% unit_times_ds = 
% stim_times = 
% binsize = 
% T = stimulus duration
% stim_key = # unique stimuli presentations x 4 (columns: x, y, bw, lcond)


%% extract spike times 
% spike_raster = make_raster_v2(unit_times_ds,field_trials,1000,4)';      % time x trials matrix of spike times (0s and 1s)
% spike_raster = make_raster_v2(unit_times_ds,stim_times(1:240:end),1000,60)';
spike_raster = make_raster_v2(unit_times_ds,stim_times,1000,T/1000)';
% [psth,~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,ones(1,size(spike_raster,2)),spike_raster',ones(1,size(spike_raster,2)));
% psth = psth';
% big_psth = psth(:);
% big_rast = spike_raster(:);     % turn into single time vector of spikes

% refresh = 4;        % vis stimuli per sec (hz)
% T = 60;        % time duration
% dt = 1/refresh;
% t = linspace(0,T,refresh*T+1);

N = sum(sum(spike_raster));
% psth = nan(nX*nY*nBW*length(diff_conds),(T/binsize)/1000);
% psth = nan(nX,nY,(T/binsize)/1000,nBW,length(diff_conds));
num_stim = size(stim_key,1);
X = unique(stim_trials(:,1));
nX = length(X);
Y = unique(stim_trials(:,2));
nY = length(Y);
BW = unique(stim_trials(:,3));
nBW = length(BW);
lconds = unique(stim_trials(:,4));
psth = nan(T/binsize,num_stim);
count = 1;
for lc = 1:length(lconds)
    for b = 1:nBW
        for y = 1:nY
            for x = 1:nX
                which_trials = zeros(1,size(spike_raster,2));
                which_stim = (stim_trials(:,1) == x)&(stim_trials(:,2) == y)&(stim_trials(:,3)==BW(b))&(stim_trials(:,4)==lconds(lc))&(run_trials'==0);    % added 1/6/20 - exclude running trials!
%                 which_stim = (stim_trials(:,1) == x)&(stim_trials(:,2) == y)&(stim_trials(:,3)==BW(b))&(stim_trials(:,4)==lconds(lc));    % include running trials
                which_trials(stim_trials(which_stim,5))=1;
                which_trials = logical(which_trials);
                [psth(:,count),~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,which_trials',spike_raster',ones(1,length(which_trials)));
                count = count+1;
            end
        end
    end
end

% % trying simple STA (11/30/18)
% test = reshape(psth(:,1:nX*nY)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% test_bl = test(:,:,1);  % first time bin = baseline
% % test_bs = test - repmat(test_bl,1,1,size(test,3)); % baseline-subtracted
% % test_norm = test_bs./max(abs(test_bs(:)));  % normalize to the maximum response (up or down) in any time bin
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
% %     imagesc(squeeze(test_norm(:,:,i)'))
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(1) = peakVar/noiseVar;
% [a,b] = max(abs(test_norm),[],3);
% [~,maxind] = max(abs(a(:)));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% 
% test = reshape(psth(:,nX*nY+1:2*nX*nY)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(2) = peakVar/noiseVar;
% 
% test = reshape(psth(:,2*nX*nY+1:3*nX*nY)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(3) = peakVar/noiseVar;
% 
% 
% test = reshape(psth(:,3*nX*nY+1:end)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(4) = peakVar/noiseVar;


% % find peak response for each subfield and each light condition
% % psth_bls = psth - repmat(psth(1,:),size(psth,1),1);
% for lc = 1:length(lconds)
%     for b = 1:nBW
%         K = psth(:,stim_key(:,3)==BW(b) & stim_key(:,4)==lconds(lc));
%         [max_bin{b,lc}, max_stim{b,lc}] = find(ismember(K, max(K(:))));
%         max_r(b,lc) = K(max_bin{b,lc}(1),max_stim{b,lc}(1));
%     end
% end
% max_r = max(max_r);        % get max for each light cond across subfields
% max_r_mat = zeros(size(psth));    
% for lc = 1:length(lconds)
%     max_r_mat(:,stim_key(:,4)==lconds(lc)) = max_r(lc);
% end
% psth_norm = psth./max_r_mat;
% % all_max = max(max_r(:));        % get cell's max firing rate
% % psth_norm = psth./all_max;
% [lowx,lowy] = find(psth_norm < .35);
% for i=1:length(lowx)
%     psth_norm(lowx(i),lowy(i)) = 0;
% end
% RFmap_on = zeros(nX,nY,T/binsize,length(lconds));
% RPmap_off = RFmap_on;
% for lc = 1:length(lconds)
%     for x=1:nX
%         for y = 1:nY
%             for t=1:T/binsize
%                 RFmap_on(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==1 & stim_key(:,4)==lconds(lc));
%                 RFmap_off(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==-1 & stim_key(:,4)==lconds(lc));
%             end
%         end
%     end
% end

psth_norm = zeros(size(psth));
% psth_norm = zeros(size(psth,1)-2,size(psth,2));
Model = zeros(size(psth_norm'));
baseline = zeros(1,size(psth_norm,2)/(nX*nY));
% for nn = 1:nBW
%     baseline([nn nn+2]) = mean(psth(1,(nn-1)*(nX*nY)+1:nn*nX*nY));     % baseline = start of first time bin for each bw condition in light OFF - use for exps with single pulses?? trials only
% end
for n = 1:size(psth_norm,2)/(nX*nY)
    baseline(n) = mean(mean(psth(1:2,(n-1)*(nX*nY)+1:n*nX*nY)));     % baseline = average of first 20ms
    psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(:,(n-1)*(nX*nY)+1:n*nX*nY)-baseline(n);
%     psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(3:end,(n-1)*(nX*nY)+1:n*nX*nY)-baseline(n);    % start from third time bin to avoid light artifacts
    [U,S,V] = svd(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY)','econ');    
    rfMapVec{n} = U(:,1);
    timeCourse{n} = V(:,1)';
    Scalar(n) = S(1,1);
    sings = S(find(S)); % singular values
    relVar = (sings.^2)./sum(sings.^2);
    if length(relVar)<size(S,1)
        relVar(length(relVar)+1:size(S,1)) = nan;  % in case there aren't enought singular values (I think only happens for bad units, but do this just so script doesn't error out...)
    end
    peakTimePoint(n) = find(abs(timeCourse{n})==max(abs(timeCourse{n})),1);    
    ts(n) = sign(timeCourse{n}(peakTimePoint(n)));
    timecoursePeak = max(ts(n)*timeCourse{n}(1:15)) > ts(n)*mean(timeCourse{n}(1:3))+5*std(timeCourse{n}(1:3));
    
    % check if first component is best one to use for RF
    whichmap(n) = 1;
    if n>2 && (relVar(2)>.1 || whichmap(rem(n,2)==rem(1:2,2))==2) % for light conds: if there's also a good second component, OR this was the case for the corresponding nolight condition
        corrs = nan(1,length(relVar));
        for i=1:sum(relVar>.1)
            corrmat = corrcoef(rfMapVec{rem(n,2)==rem(1:2,2)},U(:,i));
            corrs(i) = corrmat(1,2);
        end
        [~,whichmap(n)] = max(abs(corrs));  % find which good map is most similar to nolight cond map
    elseif relVar(2)>.1 && ~timecoursePeak % for nolight conds: if there's also a good second component but no peak timepoint in first component
        peak2TimePoint = find(abs(V(:,2))==max(abs(V(:,2))),1);
        ts2 = sign(V(peak2TimePoint,2));
        if max(ts2*V(1:15,2)) > ts2*mean(V(1:3,2))+5*std(V(1:3,2))  % first component didn't have a peak, but second does! 
            whichmap(n) = 2;
        end
    end
    if whichmap(n)>1
        rfMapVec{n} = U(:,whichmap(n));
        timeCourse{n} = V(:,whichmap(n))';
        Scalar(n) = S(whichmap(n),whichmap(n));
        peakTimePoint(n) = find(abs(V(:,whichmap(n)))==max(abs(V(:,whichmap(n)))),1);
        ts(n) = sign(V(peakTimePoint(n),whichmap(n)));
        timecoursePeak = max(ts(n)*timeCourse{n}(1:15)) > ts(n)*mean(timeCourse{n}(1:3))+5*std(timeCourse{n}(1:3));
    end

    if n>2 
%         if sign(timeCourse{n-2}(peakTimePoint(n-2))) ~= sign(ts(n)) % matches sign flip with no light condition
        if sign(timeCourse{rem(n,2)==rem(1:2,2)}(peakTimePoint(rem(n,2)==rem(1:2,2)))) ~= sign(timeCourse{n}(peakTimePoint(rem(n,2)==rem(1:2,2)))) 
            rfMapVec{n} = -rfMapVec{n};
            timeCourse{n} = -timeCourse{n};
        end            
    elseif ts(n)<0
        rfMapVec{n} = -rfMapVec{n};
        timeCourse{n} = -timeCourse{n};
    end
    rfMap{n} = reshape(rfMapVec{n},nX, nY);
    Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n) + baseline(n);
%     Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n);
    stats.timeCourse{n} = timeCourse{n};
    stats.relVar{n} = relVar(whichmap(n));
    stats.whichmap(n) = whichmap(n);
    maxZ(n) = (max(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    minZ(n) = (min(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    stats.peakZscore(n) = max(abs([minZ(n) maxZ(n)]));
%     if max(timeCourse{n}(:)) > mean(timeCourse{n}(1:3))+3*std(timeCourse{n}) % if peak of time course is >3STDs above baseline (baseline = mean of first 3 bins i.e. 30ms)
    if timecoursePeak && relVar(whichmap(n))>.1 % if peak of time course is >5STDs above baseline (baseline = mean of first 3 bins i.e. 30ms) AND has to explain at least 10% of variance
        stats.issig(n) = 1;
    else
        stats.issig(n) = 0;
    end
    
end
Residual = psth_norm' - Model;

% 
% 
% %% test if RF is significant (chi-squared test of independence across
% % % neighboring pixels of peak)
% for i=1:size(psth_norm,2)/(nX*nY)
%     tmp_psth = psth_norm(:,(i-1)*(nX*nY)+1:i*nX*nY);
%     small_stim_key = stim_key((i-1)*(nX*nY)+1:i*nX*nY,:);
% %     [~,b] = max(abs(tmp_psth(:)));
% %     peak = small_stim_key(ceil(b/size(tmp_psth,1)),1:2);
%     [~,max_which] =sort(abs(rfMap{i}(:)));
%     peak = small_stim_key(max_which(end),1:2);
% %     [~,max_inds] = sort(max(abs(tmp_psth)));
% %     peaks = small_stim_key(max_inds(end-2:end),1:2)
%     if peak(1)==1
%         first_xneighb = 1;
%         last_xneighb = peak(1)+1;
%     elseif peak(1) == nX
%         first_xneighb = peak(1)-1;
%         last_xneighb = nX;
%     else
%         first_xneighb = peak(1)-1;
%         last_xneighb = peak(1)+1;
%     end
%     if peak(2) == 1
%         first_yneighb = 1;
%         last_yneighb = peak(2) +1;
%     elseif peak(2) == nY
%         first_yneighb = peak(2)-1;
%         last_yneighb = nY;
%     else
%         first_yneighb = peak(2)-1;
%         last_yneighb = peak(2) +1;
%     end
%     [p,q] = meshgrid(first_xneighb:last_xneighb,first_yneighb:last_yneighb);
%     pairs = [p(:) q(:)];
%     which = zeros(1,length(pairs));
%     for ii = 1:length(pairs)
%         which(ii) = find((small_stim_key(:,1)==pairs(ii,1))&(small_stim_key(:,2)==pairs(ii,2)));
%     end
% %     [~,max_which] =sort(abs(rfMap{i}(:)));
% % [~,max_which] =sort(max(abs(tmp_psth(:,which))));
% %     neighb_psth = num2cell(tmp_psth(:,which(max_which(end-6:end))),1);      % int8 to save memory
%     try
%         neighb_psth = num2cell(int8(tmp_psth(:,which)),1);      % int8 to save memory
%         [~, ~,pval,~] = crosstab(neighb_psth{:});
%     catch
%         fprintf('9 is too many\n')
%         [~,max_which] =sort(max(abs(tmp_psth(:,which))));
%         neighb_psth = num2cell(int8(tmp_psth(:,which(max_which(end-7:end)))),1);      % int8 to save memory
%         try
%             [~, ~,pval,~] = crosstab(neighb_psth{:});
%         catch
%             fprintf('8 is still too many\n')
%             neighb_psth = num2cell(int8(tmp_psth(:,which(max_which(end-6:end)))),1); 
%             [~, ~,pval,~] = crosstab(neighb_psth{:});
%         end
%     end
% %     if pval < .05
%     if pval < .01
%         stats.issig(i) = 1;
%     else
%         stats.issig(i) = 0;
%         stats.pvals{i} = nan(size(rfMap{i}));
%     end
% end
% 

%% Trying new thing - look at variance across X and Y dimensions of calculated spatial filter and compare that against "shuffled" spatial filters
for i=1:length(rfMap) % if either peak was significant
%     if stats.issig(rem(i,2)==rem(1:2,2))
    if stats.issig(i)
        varY = var(rfMap{i},0,2); % get the variance across Y values for every X
        varX = var(rfMap{i},0,1); % get the variance across X values for every Y
        [peakVarY,whereVarY] = max(varY); % get X coordinate and value with maximum variance in Y
        [peakVarX,whereVarX] = max(varX); % get Y coordinate and value with maximum variance in X
        shufmat = zeros(nX,nY,1000);
        for n=1:1000
            tmp =randperm(numel(rfMap{i}));
            shufmat(:,:,n) = reshape(rfMap{i}(tmp),size(rfMap{i}));   % shuffle the spatial filter values across positions in nX by nY grid
        end
        shufvarY = squeeze(var(shufmat,0,2)); % get the variance across Y values for every X in every shuffle
        shufvarX = squeeze(var(shufmat,0,1)); % get the variance across X values for every Y in every shuffle
        shufSNR_Y = zeros(1,1000);
        if whereVarY>nX/2   % if peak in Y variances is in right half of grid
            spSNR_Y = peakVarY/mean(varY(1:2));   % get ratio of the peak Y variance to the mean of Y variances in first two columns (aka "spatial SNR")    
        else
            spSNR_Y = peakVarY/mean(varY(nX-1:nX));  % if the peak in Y variances is in the left half of the grid, get ratio of peak Y variance to the mean of Y variances in last two columns
        end
        [shufpeaksY,wherepeaksY] = max(shufvarY);   % get peak of Y variance for each shuffle
        shufSNR_Y(wherepeaksY>nX/2) = shufpeaksY(wherepeaksY>nX/2)./mean(shufvarY(1:2,wherepeaksY>nX/2)); % get same SNR val from every shuffle to get a null distribution of SNRs (using max variance from each shuffle)
        shufSNR_Y(wherepeaksY<=nX/2) = shufpeaksY(wherepeaksY<=nX/2)./mean(shufvarY(nX-1:nX,wherepeaksY<=nX/2));

        shufSNR_X = zeros(1,1000);
        if whereVarX>nY/2   % if peak in X variances is in lower half of grid
            spSNR_X = peakVarX/mean(varX(1:2));   % get ratio of the peak X variance to the mean of X variances in first two rows (aka "spatial SNR")    
        else
            spSNR_X = peakVarX/mean(varX(nY-1:nY));  % if the peak in X variances is in the upper half of the grid, get ratio of peak X variance to the mean of X variances in last two rows
        end
        [shufpeaksX,wherepeaksX] = max(shufvarX);   % get peak of X variance for each shuffle
        shufSNR_X(wherepeaksX>nY/2) = shufpeaksX(wherepeaksX>nY/2)./mean(shufvarX(1:2,wherepeaksX>nY/2)); % get same SNR val from every shuffle to get a null distribution of SNRs (using max variance from each shuffle)
        shufSNR_X(wherepeaksX<=nY/2) = shufpeaksX(wherepeaksX<=nY/2)./mean(shufvarX(nY-1:nY,wherepeaksX<=nY/2));

        probSigSNR_Y(i) = sum(shufSNR_Y>=spSNR_Y)/1000;
        probSigSNR_X(i) = sum(shufSNR_X>spSNR_X)/1000;

        if probSigSNR_Y(i) < .025 || probSigSNR_X(i) < .025 % consider it a significant receptive field if there's a significant (<.01 - strict cutoff) spatial SNR in x or y dimensions
    %     if probSigSNR_Y(i) < .05 && probSigSNR_X(i) < .05 % consider it a significant receptive field if spatial SNRs in both dimensions are significant (<.05)
            stats.issig2(i) = 1;
        else
            if max(stats.peakZscore(i)) > 5  % changed from 6 11/8/20
                stats.issig2(i) = 1;
            else
                stats.issig2(i) = 0;
            end
        end
    %     subplot(4,nCol,3*nCol+i);
    %     title(sprintf('%s map, peakZ = %.2f, sig = %d', key_leg{i}, stats.peakZscore(i),stats.issig(i)));
    else
        stats.issig2(i) = 0;
    end
end

%% make null distributions for statistical comparison
% add max z score contingency because really small RFs are unlikely to
% cross SNR significance threshold (new - 5/11/20)
tic
reps = 500;     % changed from 1000 5/11/20
which_conds = 1:size(psth_norm,2)/(nX*nY);
for i=1:size(psth_norm,2)/(nX*nY)       % new 7/9/19 - use trials specific to each condition to generate null distribution for that condition (e.g., black square, light off)
%     if sum(stats.issig(rem(i,2)==rem(which_conds,2))>=1)  || sum(stats.issig2(rem(i,2)==rem(which_conds,2))>=1) || max(stats.peakZscore(rem(i,2)==rem(which_conds,2))>6)    % include  chisquared test results
%     if sum(stats.issig2(rem(i,2)==rem(which_conds,2))>=1) ||  max(stats.peakZscore(rem(i,2)==rem(which_conds,2))>6)    % only do it if any statistically sig RFs were identified (because this takes time) - changed 5/13 to only do for sig subfields, but in opto exp, if one light cond came up sig, do for other one two
%     if stats.issig2(rem(i,2)==rem(1:2,2))
    if stats.issig2(i) || (sum(stats.issig2(rem(i,2)==rem(which_conds,2)))>=length(which_conds)/4 && stats.issig(i))  % do if sig spatial SNR, OR if it had a significant timecourse and at least half of other lightconds (same subunit, hence /4) had significant spatial SNR  - to enable comparison between conditions
% ^ 10/13 - added back in chi square test, and also do this if that test comes out significant
        rfMap_shuf = nan(nX,nY,reps);
        nshufstim = num_stim/(length(lconds)*nBW);
    %         trial_mat = repmat(1:nshufstim,1,length(which_trials)/nshufstim); % combining all bw and light conditions        
%         trial_mat = repmat(1+(i-1)*nshufstim:i*nshufstim,1,size(spike_raster,2)/num_stim);    % separately for each bw and light condition
        which_trials = zeros(1,size(spike_raster,2));
%         which_stim =
%         stim_trials(:,3)==stim_key(i*nshufstim,3)&(stim_trials(:,4)==stim_key(i*nshufstim,4)); % shuffling across trials of same LED AND luminance conditions
        which_stim = (stim_trials(:,4)==stim_key(i*nshufstim,4)); % shuffling across trials of same LED condition
        which_trials(stim_trials(which_stim,5))=1;
        which_trials = logical(which_trials);   
%         trial_mat = unique(stim_trials(which_stim,5));   % new 6/10/20
        find_stim = find(which_stim);   % get trial indices
        trial_mat = nan(1,nshufstim); % defines the x and ypos for each trial with nth bw and LED condition
        for x=1:length(find_stim)  % couldn't figure out how to do without for loop if trial numbers aren't quite equal across conditions...(MAK 9/14/20)
            trial_mat(x) = find((stim_key(1+(i-1)*nshufstim:i*nshufstim,1) == stim_trials(find_stim(x),1))&(stim_key(1+(i-1)*nshufstim:i*nshufstim,2) == stim_trials(find_stim(x),2))); % new MAK 9/14/20)
        end
        for ii = 1:reps
            shuf_trials = trial_mat(randperm(length(trial_mat)));    % randomly assign each stim presentation to one of nX*nY stimulus conditions AND shuffle time points (MAK 1/6/20)
            [psth_shuf,~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,1:sum(which_trials),spike_raster(randperm(size(spike_raster,1)),which_trials)',shuf_trials); %also shuffles in time **important when RFs are very large
%             [psth_shuf,~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,1:sum(which_trials),spike_raster(:,which_trials)',shuf_trials);
    %             shutf_t = nan(size(psth_shuf));
    %             for tt = 1:size(psth_shuf,1)    % NEW 7/20/19 - also shuffle in time
    %                 shuf_t(tt,:) = randperm(size(psth_shuf,2));     % shuffle time bin indices for each stim rep
    %             end
    %             psth_shuf = psth_shuf(shuf_t);
            bs = mean(mean(psth_shuf(:,1:2)));    % same as above - baseline is mean of first 20ms
    %         psth_shuf = psth_shuf(:,3:end) - bs;
            psth_shuf = psth_shuf - bs;
            [Ushuf,~,~] = svd(psth_shuf,'econ');        % get spatial filter for baseline-subtracted shuffled psth
            rfMapVec_shuf = Ushuf(:,1);
            rfMap_shuf(:,:,ii) = reshape(rfMapVec_shuf,nX, nY);
        end
    % make_psth_v2(binsize/1000,[0:binsize:T]./1000,which_trials',spike_raster',ones(1,length(which_trials)));
        pvals{i} = zeros(size(rfMap{i}));
        for xx = 1:nX
            for yy = 1:nY
                if rfMap{i}(xx,yy)<0
                    pvals{i}(xx,yy) = -sum(rfMap_shuf(xx,yy,:)<=rfMap{i}(xx,yy))/reps;
                else
                    pvals{i}(xx,yy) = sum(rfMap_shuf(xx,yy,:)>=rfMap{i}(xx,yy))/reps;
                end
            end
        end
    %     pvals{i}(abs(pvals{i})>=.05/(nX*nY)) = nan;   % because two-sided test. do I need a multiple comparisons correction?
        [sort_pvals,~] = sort(abs(pvals{i}(:)));      % treat all pvals as +
    %         testp = zeros(1,length(sort_pvals));

    % Benjamini-Hochberg FDR correction for multiple comparisons (false discovery rate set to 25%)
%         rnk = 1:numel(sort_pvals);
%         imq = rnk./(nX*nY)*.25;
%         max_p = sort_pvals(find(sort_pvals<imq',1,'last'));
%     %         max_p = max(abs(pvals_ord(find(testp))));
%         if isempty(max_p)
%             max_p = 0;
%         end
        [h, max_p] = fdr_bh(sort_pvals, .2); % 10% false discovery rate
        if sum(h)
            pvals{i}(abs(pvals{i})>max_p) = nan; 
            stats.issig3(i) = 1;
        else
            pvals{i} = nan(size(pvals{i}));
            stats.issig3(i) = 0;
        end
        stats.pvals{i} = pvals{i};
        
    else % if no sig RF detected for on/off responses
        stats.pvals{i} = nan(size(rfMap{1}));
        stats.issig3(i) = 0;
    end
toc

end

%%
figure;
timeBins = 0:binsize:T;
timeBins = timeBins(1:end-1)+binsize/2;
nCol = length(rfMap);
nRow = 4;
% timeBins=timeBins(3:end);       %TEMP

if length(rfMap)<=4 
    key_leg = {'OFF subunit, light OFF', 'ON subunit, light OFF', 'OFF subunit, light ON', 'ON subunit, light ON'};
else
    key_leg = {'OFF subunit, LED OFF', 'ON subunit, LED OFF', 'OFF subunit, LED1 ON', 'ON subunit, LED1 ON', 'OFF subunit, LED2 ON', 'ON subunit, LED2 ON', 'OFF subunit, LEDs1&2 ON', 'ON subunit, LEDs1&2 ON'};
end
leg_ax = [min(cellfun(@(x) min(x(:)),rfMap,'uniformoutput',1)) max(cellfun(@(x) max(x(:)),rfMap,'uniformoutput',1))];
for n = 1:size(psth_norm,2)/(nX*nY)
    
    subplot(nRow,nCol,n);
    plot(timeBins, timeCourse{n});
    xlim([timeBins(1) timeBins(end)]);
    title('response time course')
    xlabel('time')
    
    ax = subplot(nRow,nCol,nCol+n);
    imagesc(1:size(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY),2), timeBins,  psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY));
    %         cax = caxis();
    %         caxis(max(abs(cax))*[-1 1]);
    % colormap(colorcet('L3'));
    %         axis image
    colorbar
    title(sprintf('all PSTHs - %s',key_leg{n}));
    xlabel('space')
    ylabel('time')
    
%     subplot(nRow,nCol,2*nCol+n); 
%     imagesc(Residual((n-1)*(nX*nY)+1:n*nX*nY,:)');
%     %         cax = caxis();
%     %         caxis(max(abs(cax))*[-1 1]);     
%     %         axis image
%     colorbar
%     title(sprintf('residual - %s',key_leg{n}))
%     xlabel('space')
%     ylabel('time')
    
    subplot(nRow,nCol,2*nCol+n);
    imagesc(1:nX, 1:nY, rfMap{n}(:,:)');
    axis image
%     caxis(leg_ax)
    colorbar
    title(sprintf('%s map, peakZ = %.2f', key_leg{n}, stats.peakZscore(n)));
    
    subplot(nRow,nCol,3*nCol+n)
%     zmap = transpose(rfMap{n}-mean(rfMap{n}(:)))./std(rfMap{n}(:));
%     zmap(abs(zmap)<2)= nan;
%     imagesc(zmap)
%     title(sprintf('%s Zmap thresh (issig = %.2f)',key_leg{n},stats.issig2(n)));
    imagesc(1:nX,1:nY,stats.pvals{n}');
    axis image
    caxis([-.05 .05])
    colorbar
    title(sprintf('%s RF sig map (issig = %.2f)',key_leg{n},stats.issig2(n)));

end


    
% % 
% % psth_long = reshape(psth',1,size(psth,2)*size(psth,1));
% % STA = (1/N).*(X*psth_long);
% 
% 
% 
% big_rast = spike_raster(:);     % turn into single time vector of spikes
% light_mat = zeros(size(spike_raster));
% light_mat(:,stim_trials(:,4)==1) = 1;
% light_mat = light_mat(:);
% stim_mat = zeros(nX,nY,size(big_rast,1));  % initialize 3d matrix for stimulus (ms in all trials x xposition x yposition)
% for i=1:size(spike_raster,2)
%     stim_mat(stim_trials(i,1:2),(i-1)*T+1:i*T) = stim_trials(i,3);
% end
% N = sum(big_rast);          % total number of spikes
% spks = find(big_rast);
% spks_light = find(big_rast&light_mat);
% spks_nolight = find(big_rast&~light_mat);
% N_light = length(spks_light);
% N_nolight = length(spks_nolight);
% X_light = zeros(nX*nY*(T+1),N_light);
% for i=1:N_light
%     X_light(:,i) = reshape(permute(stim_mat(spks_light(i):spks_light(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% %     X_light(:,i) = reshape(permute(stim_mat(spks_light(i):spks_light(i)+T,:,:),[1 3 2]),nX*nY*(T+1),1);
% end
% %trying something different but not working yet...
% X_light = reshape(stim_mat(:,:,light_mat==1),nX*nY,sum(light_mat));
% STAlight = (1/N_light).*(X_light*big_rast(light_mat==1));
% STAlight = reshape(STAlight,nX,nY);
% 
% STAlight = (1/N_light).*(X_light*big_rast(spks_light));
% STAlight = reshape(STAlight,nX,nY,T+1);
% figure;imagesc(mean(STAlight(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STAlight(:,:,(n-1)*10+1:n*10),3))
% end
% 
% clear X_light
% X_nolight = zeros(nX*nY*(T+1),N_nolight);
% for i=1:N_nolight
%     X_nolight(:,i) = reshape(permute(stim_mat_ext(spks_nolight(i):spks_nolight(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% end
% STAnolight = (1/N_nolight).*(X_nolight*big_rast(spks_nolight));
% STAnolight = reshape(STAnolight,nX,nY,T+1);
% figure;imagesc(mean(STAnolight(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STAnolight(:,:,(n-1)*10+1:n*10),3))
% end
% 
% clear X_nolight
% X = zeros(nX*nY*(T+1),N);
% for i=1:N
%     X(:,i) = reshape(permute(stim_mat_ext(spks(i):spks(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% end
% STA = (1/N).*(X*big_rast(spks));
% STA = reshape(STA,nX,nY,T+1);
% figure;imagesc(mean(STA(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STA(:,:,(n-1)*10+1:n*10),3))
% end
% 
% % sta = ifft(fft(big_rast).*conj(fft(stim_mat,[],1)),[],1);
% % sta_nolight = ifft(fft(big_rast(light_mat==0)).*conj(fft(stim_mat(light_mat==0,:,:),[],1)),[],1);
% % N_nolight = sum(big_rast(light_mat==0));
% % sta_light = ifft(fft(big_rast(light_mat==1)).*conj(fft(stim_mat(light_mat==1,:,:),[],1)),[],1);
% % N_light = sum(big_rast(light_mat==1));
% 
% 
%     



return