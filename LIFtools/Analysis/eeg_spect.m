function [s,f,err] = eeg_spect(ctx, dt)


freq_range = [0.1 5];    		    % min and max freqs in periodogram
tapers = [3 5];                         % recommended in Chronuz; 5 tapers with bandwidth*time = 3
pad = 2;                                % Padding the frequency sampling by another 2 powers of 2;
Fs = 1 ./ dt;                           % best to use the underlying sampling;
err_bars = [1 0.01];                    % method 1 (theoretical) and p_sig = 0.01
fscorr = 0;                             %don't use finite size corrections

data = ctx - mean(ctx);
data = data';
[s,f,err]=mtspectrumc(data,tapers,pad,Fs,freq_range,err_bars,0);

plot(f,s)

