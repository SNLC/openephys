function [Fratio,zF1,F1] = calc_F1F0(pref_psth,binsize,tfreq)

% inputs:
% psth - num_conds x N (timebins) PSTH of firing rates, calculated from the
% cells preferred visual stimulus 
% binsize - length of timebins of psth (in sec)
% tfreq - temporal frequency of visual stimulus (in Hz)

% modified 4/30/20 to also output standardized F1, which is independent of
% baseline firing rate (from Wypych et al., 2012; https://www.sciencedirect.com/science/article/pii/S0042698912002957#f0005)


% pref_psth_bs = pref_psth-repmat(FRs(clean_units(i)).psthBlank(1,21:end),3,1);
fs = 1/binsize;    % sampling frequency
N = size(pref_psth,1);  % number of samples
F_Nyq = fs/2;
t = [0:N]*binsize;   % time vec
T = max(t);
dw = 2*pi/(N*binsize);     % freq res
w = [0:N-1]*dw; % freq vector
F = fft(pref_psth);
Fnorm=abs(F./N); 

w2 = w;
w2(w>(N*dw/2)) = w2(w>(N*dw/2))-N*dw;
df = w2(w2>=0)/(2*pi); % equivalent to [0:1/T:fs-1/T]
F1 = Fnorm(df==tfreq,:);        % F1 is unaffected by baseline subtraction
Fratio = F1./Fnorm(1,:);
range = find(df==1/(T)):find(df<=F_Nyq,1,'last');
zF1 = (F1-mean(Fnorm(range,:)))./std(Fnorm(range,:));
% figure;semilogy(df,abs(Fnorm(1:length(df),1)),'r.-')

end