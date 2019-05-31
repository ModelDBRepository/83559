function [correlogram,variance] = shuffle_corrector(events1,times1,events2,times2,pairs,bin_size,T,varargin)

% SHUFFLE_CORRECTOR compute shuffle corrector for cross-correlogram (i.e. for covariogram generation) 
%
%   [S,V] = SHUFFLE_CORRECTOR(E1,T1,E2,T2,P,BINSIZE,T) where E1, E2 are arrays of all events to be correlated, T1, T2 are their 
%   corresponding time-stamp arrays, P is a 2 column list of the event indices which match up (e.g. those corresponding to the same trial), 
%   BINSIZE is the required bin-size in seconds, and T is a 2-element array specifying the start 
%   and end of the spike trains in seconds. (If the result of this function is to be used to
%   create the covariogram, then the BINSIZE MUST be the same as the bin size used for the correlogram), 
%   
%   Generates the shuffle corrector of the spike trains by generating their PSTHs and correlating them, returning S the array for the
%   shuffle corrector and V the corresponding array of PSTH bin variance (computed according to Brody's eq 2.4).
%
%   SHUFFLE_CORRECTOR(...,MAXLAG) sets the size of the maximum interval to be calculated to +/- MAXLAG seconds: 
%   if a value was specified for the parent correlogram, then this must be specified here too.
%
%   NOTE: it is assumed that all events/times passed are to be used for the PSTHs - P is just passed for conveniene
%   (i.e. that they've been extracted in the calling function such as COVARIOGRAM)
%
%   Reference: Brody, C. D. (1999) "Correlations without synchrony" Neural Computation, 11, 1537-1551
%
%   Mark Humphries 22/12/04

if nargin >= 8 & isnumeric(varargin{1}) & ~isempty(varargin{1})
    max_lag = varargin{1};
    bins = -max_lag:bin_size:max_lag;
else
    max_lag = 0;
    bins = -time_seconds:bin_size:time_seconds;    
end
num_bins = length(bins);

time_seconds = T(2) - T(1);
[num_sweeps c] = size(pairs);

% compute mean PSTHs
[P1,X,V1] = LIF_psth(events1,times1,bin_size,T,'m');
[P2,X,V2] = LIF_psth(events2,times2,bin_size,T,'m');

% initialise all necessary values
correlogram = zeros(num_bins,1);
var_corr = zeros(num_bins,1);
P1_v2 = zeros(num_bins,1);
P2_v1 = zeros(num_bins,1);

% window reference PSTH according to max_lag
tau = max_lag / bin_size;                   % lag in number of bins
t1 = tau + 1;                               % start time in number of bins (+1 so that can go back to index 1, rather than 0)
t2 = length(P1)-tau;                        % end time in number of bins

% compute shuffle corrector (eq 2.2)
% AND compute variance in the null hypothesis (eq 2.4)
for loop1 = t1:t2
    correlogram = correlogram + P1(loop1) .* P2(loop1-tau:loop1+tau);
    var_corr = var_corr + V1(loop1) .* V2(loop1-tau:loop1+tau);
    P1_v2 = P1_v2 + (P1(loop1)^2) .* V2(loop1-tau:loop1+tau);
    P2_v1 = P2_v1 + (P2(loop1)^2) .* V1(loop1-tau:loop1+tau);
end
    
variance = (var_corr + P1_v2 + P2_v1) ./ num_sweeps;