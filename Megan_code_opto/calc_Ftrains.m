function F_Z = calc_Ftrains(psth,binsize,tfreq)

% inputs:
% psth - num_conds x N (timebins) PSTH of firing rates, calculated from the
% cells preferred visual stimulus 
% binsize - length of timebins of psth (in sec)
% tfreq - temporal frequency of visual stimulus (in Hz)



% pref_psth_bs = pref_psth-repmat(FRs(clean_units(i)).psthBlank(1,21:end),3,1);
fs = 1/binsize;    % sampling frequency
N = size(psth,1);  % number of samples
F_Nyq = fs/2;
t = [0:N-1]*binsize;   % time vec
T = max(t);
dw = 2*pi/(N*binsize);     % freq res
w = [0:N-1]*dw; % freq vector
F = fft(psth);
Fnorm=abs(F./N); 

w2 = w;
w2(w>(N*dw/2)) = w2(w>(N*dw/2))-N*dw;
df = w2(w2>=0)/(2*pi);
% Fratio = Fnorm(df==tfreq)./max(Fnorm(2:find(df==tfreq)-1));
F_Z = (Fnorm(df==tfreq)-mean(Fnorm(2:find(df==tfreq)-1)))/std(Fnorm(2:find(df==tfreq)-1));
% figure;semilogy(df,abs(Fnorm(1:length(df),1)),'r.-')

end