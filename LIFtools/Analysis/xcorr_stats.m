function [m,s,sp,st] = xcorr_stats(y,n_pairs,varargin)

% NOTE: this function is currently incomplete
% XCORR_STATS correlogram statistics
%
%   XCORR_STATS(Y,N) where Y is an array of correlogram bin values (raw counts, not normalised), 
%   N is the number of spike pairs used to generate the correlogram. 
%
%   [M,S,CP,CT] = XCORR_STATS(...) Computes the mean M, standard deviation S
%   and 99% confidence interval (p < 0.01) for peaks and troughs, CP and CT.
%
%   XCORR_STATS(...,P) specifies the confidence interval probability value P to be calculated
%   (e.g. P = 0.01 is the 99% confidence interval).
%
%   XCORR_STATS(...,'cov') will return the mean, standard deviation, and confidence intervals for the normalised covariogram. The
%   mean is always zero.
%
%   It is assumed that the two (normally mean) spike time-series T1 and T2 used to compute Y are independent. Thus,
%   as long as they have a constant firing rate, we have the null hypothesis that the time-series 
%   are both homogenous Poisson processes. Therefore, the correlogram of two completely independent time-series
%   should be flat (that is, every bin has the same value) as all intervals are equally likely (though this is not the
%   case when refractory periods are enforced, their effect should be minimal). In other words, the expected value E (mean) 
%   of every bin of Y is the (number of spikes)/(number of bins). We can then calculate the confidence intervals using the cumulative
%   probability distribution of a Poisson function with lambda = E (for whole number, e.g. total, spike counts, confidence intervals are set
%   at nearest integer spike count to the required P value).
%   When interpreting the results, it should be remembered that if (number bins in Y) > 100  then for
%   a 99% confidence interval we can reasonably expect at least one significant bin: meaningful peaks or troughs should therefore
%   consist of a number of consecutively significant bins.
%
%   References: Dayan, P. & Abbot, L. F. (2001) Theoretical Neuroscience. Cambridge, MA: MIT Press.
%               Abeles, M. (1982) "Quantification, smoothing, and confidence limits for single-units' histograms" Journal of
%                   Neuroscience Methods, 5, 317-325
%
%   Mark Humphries 21/5/2004

num_bins = length(y);

% mean count per bin
m = n_pairs / num_bins;

% if P specified then set
if nargin >= 3
    P = varargin{1};
else
    % default to 99% confidence limit
    P = 0.01;
end
P_inv = 1 - P;


% if nargin == 4 & strcmp('cov',varargin{2})  % does this work??
%     m = 0;
%     s = s / length(y);
% end

% find lower and upper confidence limits using cumulative propbability function of Poisson distribution
% requires only mean value (m) - calculated by P[x;m] = exp(-m)*m^x / x! where x is the integer bin count 
% value being tested
found = 0;
x = 0;
P_low = 0;

% find lower confidence limit spike count
while ~found
    P_low = cdf('poiss',x,m);  % cumulative probability
    if P_low > P
        found = 1;
        x = x - 1;  % confidence limit set to integer spike count before probability value
    else
        x = x + 1;
    end    
end

if x < 0
    x = 0;
end
st = x; % lower confidence limit is spike count x

% find upper confidence limit
P_high = 0; % start from where it left off
found = 0;
x = 0;
while ~found
    P_high = cdf('poiss',x,m);  % cumulative probability
    if P_high > P_inv
        found = 1;
    else
        x = x + 1;
    end    
end
sp = x;


% estimate Poisson cumulative probabilities using Gaussian distribution
% std of bin
s = sqrt(m);

spn = m + s * 2.58;
stn = m - s * 2.58;

    
% calculate probability of each spike count according to Poisson distribution
% P = zeros(num_bins,1);
% 
% for loop = 1:num_bins
%     P(loop) = (exp(-m)*m.^y(loop)) / factorial(y(loop));
% end





