function [y,bins,f1,f2,num_pairs] = LIF_xcorr(times1,times2,bin_size,time_window,varargin)

% LIF_XCORR cross-correlogram/covariogram histogram
%
%   LIF_XCORR(A,B,BINSIZE,T) where A and B are time-stamp (in seconds) arrays of spiking events,
%   BINSIZE is the size of the sample bin in seconds, T is a 2-element array specifying the start and end of the spike trains in seconds.
%   Computes the cross-correlation histogram between two spike trains, with A as the reference train. This is most useful when the spike trains
%   can be compared over their entirety. For computing correlations between sets of spike trains recorded over many trials, use either
%   JPSTH or COVARIOGRAM (the latter calls this function)
%
%   LIF_XCORR(...,MAXLAG) sets the size of the maximum 
%   interval returned to be +/- MAXLAG seconds; 
% 
%   LIF_XCORR(...,'cov') calculates the normalised covariogram. Put MAXLAG = [] if necessary. 
%   This allows a more accurate assessment of significant peaks and troughs (and can be compared to
%   covariogram bounds calculated using XCORR_BOUNDS)
%
%   [Y,BINS,F1,F2,NP] = LIF_XCORR(...) returns the binned interval counts Y and the bin centres B. 
%   Plot using BAR(BINS,Y,1). Can also optionally specify F1 and F2 to return the mean firing rates of trains A and B.
%   Specify NP to return the number of spike pairs used to calculate the correlogram
%   
%   NOTE#1: specify a maximum interval MAXLAG for best results. It ensures that the histogram is
%   approximately flat if the signals are not correlated. If the whole
%   signal is used for both reference and comparison, then edge effects
%   will result as the extreme intervals are necessarily low in frequency.
%
%   NOTE#2: try initial values of BINSIZE = 0.001 and MAXLAG = 0.1 and
%   adjust as necessary - BINSIZE should generally be the quantising step of the underlying spike-trains,
%   or a sufficiently small time-window to guarantee a single spike per bin in each train; so for simulated trains,
%   this could be the smallest absolute refractory period if known
%
%   NOTE#3: to get autocorrelogram, just enter the same spike-train for A and B
%
%   REFERENCE: Dayan, P & Abbott, L. F. (2001) Theoretical Neuroscience. Cambridge, MA: MIT Press
%
%   Mark Humphries 22/4/04

if time_window(2) < time_window(1)
    error('End of spike train time window must be after start')
end

time_seconds = time_window(2) - time_window(1);

% turn arrays right way round for further processing
[r1 c1] = size(times1);
[r2 c2] = size(times2);
if r1 > c1
    times1 = times1';
end
if r2 > c2
    times2 = times2';
end

if nargin >= 5 & isnumeric(varargin{1}) & ~isempty(varargin{1})
     max_lag = varargin{1};
else
    max_lag = 0;    % do not use max lag
end

bins = -time_seconds:bin_size:time_seconds;

num_spikes1 = length(times1);
num_spikes2 = length(times2);

f1 = num_spikes1/time_seconds;
f2 = num_spikes2/time_seconds;

if num_spikes1 < 2 | num_spikes2 < 2
    % then not sufficient to process
    y = zeros(1,length(bins));
    num_pairs = 0;
    return
end

% remove time-stamps that are within max_lag of each end from reference
% train - just use spikes within full window
temp = times1;
%if max_lag
%    temp(times1 < time_window(1)+max_lag | times1 > time_window(2)-max_lag) = [];
%end

mat1 = repmat(temp,length(times2),1);
mat2 = repmat(times2',1,length(temp));

diff = mat1 - mat2;             % matrix of all time differences between spikes

% remove all that are greater than max_lag
%if max_lag
%    diff(abs(diff)>max_lag) = [];
%end

%keyboard
y = hist(diff(:),bins);        % bin time-differences according to bin array

num_pairs = length(diff(:));

% normalise bins - Abeles...

if max_lag
     start_bin = find(bins == -max_lag);
     end_bin = find(bins == max_lag);
     bins = bins(start_bin:end_bin);
     y = y(start_bin:end_bin);
end


% subtract mean to create  normalised covariogram
if nargin == 6  & strcmp('cov',varargin{2})
    mean_bin_value = num_pairs/length(bins);
    std_bin = sqrt(mean_bin_value);
    
    norm_mean = mean_bin_value / length(bins);
	std_bin = std_bin / length(bins);
    y = y./length(bins) - norm_mean;
end







