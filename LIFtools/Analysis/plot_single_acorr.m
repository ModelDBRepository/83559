function [ccf, x, n_pairs]  = plot_single_acorr(times,binsize, time_seconds, maxlag)

[ccf,x,f1,f2,n_pairs] = LIF_xcorr(times,times,binsize,[0 time_seconds],maxlag);

minccf = min(ccf);
zbin = find(x == 0);
ccf(zbin) = minccf; % take out the zero self-spike bin
maxccf = max(ccf);

plot(x,ccf)
axis([min(x) max(x) 0 maxccf])
% figure
% bar(x,ccf,1)
