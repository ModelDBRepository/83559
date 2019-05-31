function [y,x,variance,j_max,j_min] = covariogram(events1,times1,events2,times2,pairs,bin_size,T,varargin)

% COVARIOGRAM spike train covariogram with shuffle corrector
%
%   COVARIOGRAM(E1,T1,E2,T2,P,BINSIZE,T) where E1, E2 are arrays of all events to be correlated, T1, T2 are their 
%   corresponding time-stamp arrays, P is a 2 column list of the event indices which match up (e.g. those corresponding to the same trial), 
%   BINSIZE is the required bin-size in seconds, 
%   and T is a 2-element array specifying the start and end of the spike trains in seconds. Note that E1 and E2 must contain the same
%   set of event indices
%   
%   [Y,X,V,JMAX,JMIN] = COVARIOGRAM(...) Generates the cross-covariogram Y (with bin centres X) between the spike trains in E1,T1 and E2,T2. 
%   Thus, these arrays should contain individual trains which are either (a) responses from the same cell to repeated stimuli 
%   or (b) response from different cells to the same stimuli which can be considered closely related 
%   (as this process computes the mean PSTH for the whole set of trains in each array; for
%   example, if we are looking at the population response of a single channel in the GPR/GHS model)
%   
%   Values in the cross-covariogram are normalised, thus they are also known as correlation coefficients. 
%   The variance for each bin is returned in array V. The upper and lower bounds for
%   the coefficient values are dependent on the mean firing rates of the constituent spike trains. Thus, these bounds (JMAX and JMIN) are
%   also computed so that comparisons between covariograms can be made.
%
%   COVARIOGRAM(...,MAXLAG) sets the size of the maximum interval to be calculated to +/- MAXLAG seconds.
%
%   NOTE#1: specify a maximum interval MAXLAG for best results. It ensures that the histogram is
%   approximately flat if the signals are not correlated. If the whole
%   signal is used for both reference and comparison, then edge effects
%   will result as the extreme intervals are necessarily low in frequency.
%
%   NOTE#2: try initial values of BINSIZE = 0.001 and MAXLAG = 0.1 and
%   adjust as necessary - BINSIZE should be the quantising step of the underlying spike-trains to
%   ensure accurate results
%
%   Reference: Brody, C. D. (1999) "Correlations without synchrony" Neural Computation, 11, 1537-1551
%
%   Mark Humphries 22/12/04

time_seconds = T(2)-T(1);

if nargin >= 8
    max_lag = varargin{1};
    bins = -max_lag:bin_size:max_lag;
else
    max_lag = [];
    bins = -time_seconds:bin_size:time_seconds;    
end

[num_loops c] = size(pairs);
correlogram = zeros(num_loops,length(bins));
f1 = zeros(num_loops,1);
f2 = zeros(num_loops,1);
t1 = [];
t2 = [];
e1 = [];
e2 = [];

% loop and compute all individual correlograms and extract time/event arrays
for loop = 1:num_loops
    current_times1 = times1(events1 == pairs(loop,1));     % create numerical array of all times for this event in pair
    current_times2 = times2(events2 == pairs(loop,2));
    current_events1 = events1(events1 == pairs(loop,1));
    current_events2 = events2(events2 == pairs(loop,2));
    [correlogram(loop,:),x,f1(loop),f2(loop),n] = LIF_xcorr(current_times1,current_times2,bin_size,T,max_lag); 
    
    % store extracted time and event arrays
    t1 = [t1 current_times1];
    t2 = [t2 current_times2];
    e1 = [e1 current_events1];
    e2 = [e2 current_events2];
end

% compute shuffle corrector
[K, variance] = shuffle_corrector(e1,t1,e2,t2,pairs,bin_size,T,max_lag);

% covariogram
y = mean(correlogram) - K';

% upper and lower bounds
mean_f1 = mean(f1);
mean_f2 = mean(f2);

[j_max,j_min] = xcorr_bounds(mean_f1,mean_f2,bin_size);




