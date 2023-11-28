function [corr,xaxis] = LIF_matlab_xcorr(varargin)

% LIF_MATLAB_XCORR computes the cross-correlation of firing rate time-series waveforms
%
%   LIF_MATLAB_XCORR(A,B,dt) where A,B are time-series of firing rate data (min length 2), and dt is the step-size of the sample bin in seconds. 
%   The time-series do not have to be converted to
%   firing rates [i.e. 1/t] but must be binned because (a) the cross-correlation should be computed on time-series 
%   sampled at identical rates (though not necessarily for the same period) and (b) smoothing the data 
%   through binning will assist in detecting correlated firing that would otherwise be masked by noise). 
%   
%   LIF_MATLAB_XCORR(A,B,dt,'cov') substracts the mean of the provided waveforms before 
%   computing the cross-correlation so that the cross-covariation is the resultant waveform (generally easier to interpret).
%
%   LIF_MATLAB_XCORR(A,dt) where A is a time-series of firing rate data and dt is the step-size of the sample bin in seconds, 
%   performs the auto-correlation function. If 'cov' is specified as the third argument then the resultant waveform is 
%   the auto-covariation function.
%
%   [C,D] = LIF_MATLAB_XCORR(...) returns C the cross- (or auto-) correlation (covariation) array, and D the time-points (i.e. the x-axis) 
%   of the time-series specified in C.
%
%   Mark Humphries 2/4/2004

if nargin < 2
    error('Not enough input arguments')
elseif nargin > 4
    error('Too many input arguments')
end

% assign time-series to interim variables
signal1 = varargin{1};
if length(varargin{2}) == 1;  
    % second input argument is time-step so auto-correlation (covariation)
    signal2 = varargin{1};
    dt = varargin{2};
else
    signal2 = varargin{2};
    dt = varargin{3};
end

% if doing covariation then subtract mean
if strcmp(varargin{end},'cov')
    signal1 = signal1-mean(signal1);
    signal2 = signal2-mean(signal2); 
end  
    
% do correlation
corr = xcorr(signal1,signal2,'unbiased');

% rotate array if in wrong direction
[r c] = size(corr);
if r > 1
    corr = corr';
end

% generate appropriate time-stamp array for correlation result
wvfrm_length = max(length(signal1),length(signal2));
end_point = wvfrm_length - 1;
xaxis = dt .* (-end_point:end_point);