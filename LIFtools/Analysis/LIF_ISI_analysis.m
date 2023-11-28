function [firing_rate,isihist,times,fshist,varargout] = LIF_ISI_analysis(times,num_hist_bins,varargin)

% LIF_ISI_ANALYSIS create inter-spike interval (ISI) data from spike trains
%
%   LIF_ISI_ANALYSIS(A,NUMBINS), where A is a time-stamp array of the spike train in seconds,
%   NUMBINS is the number of bins for the ISI histogram.
%
%   [C,D,E,F] = LIF_ISI_ANALYSIS..., returns: C, array ISI-based firing-rates - i.e. the raw ISI plot; D, array of ISI histogram 
%   including all ISIs across entire range (plot with bar); E, array of time-points for ISI plot (i.e. the spike times minus the last spike); and F the x-axis for the ISI histogram
%   
%   LIF_ISI_ANALYSIS(...,L) sets limits on the ISI histogram, where L is a 2-element array, specifying the min and max ISI in seconds. 
%   This option is of most use when batch-analysing spike-trains so that ISI histograms may be averaged.
%   If explicit limits are used, then the function also returns G, a binary variable indicating the existence of outliers to this range. 
%   
%   LIF_ISI_ANALYSIS(...,L,'m') returns the time-point array E as the mid-point between the pair of spikes which generated the corresponding ISI.
%   Set L = [] if not required.
%     
%   Note: if A has less than two elements, then zero arrays of suitable length will be returned
%
%   Mark Humphries 29/11/2004

% initialisation
[r c] = size(times);
if r > c
    times = times';
end

% check there's something to process
if length(times) < 2
    % nothing to process
    firing_rate = zeros(1,length(times));
    isihist = 0;
    fshist = 0;
    return
end


% calculate ISI data %
intervals = abs(diff(times));     % calculate intervals in seconds, use abs() in case times are pre- and post- stimulus

% either set range.....
if nargin >= 3 & ~isempty(varargin{1})
    L = varargin{1};
    if L(1) >= L(2)
        error('First ISI limit must be less than second ISI limit');
    end
    
    fshist = linspace(L(1),L(2),num_hist_bins); 
    % outliers?
    outliers = find(L(1) > intervals | L(2) < intervals);
    varargout{1} = ~isempty(outliers);
    
else
    % or dynamically size ISI hist
    max_diff = max(intervals);              
    fshist = linspace(0,max_diff,num_hist_bins);    % ISI histogram x-axis in seconds (lowest interval,biggest interval,steps between) 
end
% do histogram
isihist = hist(intervals,fshist)';              % create histogram of interval distribution

% firing rates
firing_rate = 1./intervals;         % take reciprocal of intervals to find instantaneous firing rate (IFR)  
times(end) = [];                    % can't plot for last spike....

if nargin >= 4 & findstr('m',varargin{2})
    times = times + intervals./2;       % plot at mid-points
end

%time_points = 0:dt:time_seconds;    % time-points for ISI-based firing rate interpolation
%IFRplot = kginterp1(times,firing_rate,time_points,'nearest');
