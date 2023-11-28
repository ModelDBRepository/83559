function [y,x,var_psth] = LIF_psth(events,times,bin_size,T,varargin)

% LIF_PSTH peri-stimulus time histogram generation
%
%   LIF_PSTH(A,B,BINSIZE,T) where A and B are arrays of spike train events-stamps and time-stamps (in seconds) 
%   for the simulation, respectively, BINSIZE is the required bin-size (in seconds), and T is the 2-element array of 
%   start and end times of the spike train(s) in seconds.
%   Generates the peristimulus histogram of the spike trains encoded in A and B.
%
%   [Y,X,V] = LIF_PSTH(...) Returns Y the array of the PSTH (spike counts or mean), X the bin centre times, and V the array of variance
%   in each bin (particularly useful for SHUFFLE_CORRECTOR).
%
%   LIF_PSTH(A,B,CELLS,BINSIZE,T,FLAG) where FLAG is (or any combination of):
%       'p' - will plot the PSTH with the correct formatting
%       'm' - will calculate mean number of spikes per bin (i.e. will normalise by number of sweeps)
%
%   LIF_PSTH(A,B,CELLS,BINSIZE,T,FLAG,STRING) adds the text in STRING to the end of the title   
%
%   NOTE#1: divide resulting raw output values by binsize to get firing rate estimates
%
%   NOTE#2: if the times in B encode a single continuous train (and there is, therefore, a single event number in A) then using
%   a large BINSIZE will result in a spike-count per time-interval (e.g. use BINSIZE = 1 to generate accurate firing rates)
%
%   NOTE#3: it is assumed that all events and times passed to this function will contribute to the PSTH 
%
%   Mark Humphries 21/12/04

if nargin >= 5 & findstr(varargin{1},'m')
    mean_bins = 1;
else
    mean_bins = 0;
end
     
if T(1) >= T(2)
    error('Start time must be less than end time')
end

time_seconds = T(2) - T(1);
x = T(1)+bin_size/2:bin_size:T(2)-bin_size/2;   % give x-axis values as centres for plotting
bin_edges = T(1):bin_size:T(2);                % specify bin-edges for histogram
num_bins = length(bin_edges)-1;

% calculate number of sweeps from set of event numbers passed
sweeps = unique(events);
num_sweeps = length(sweeps);

y = zeros(num_bins,1);
var_psth = zeros(num_bins,1);
event_count = zeros(num_sweeps,1);

% for every bin, sum all times from every event that fall within that bin
for loop = 1:num_bins
    b1 = T(1) + bin_size * (loop-1);
    b2 = T(1) + bin_size * loop;
    bin_times_idx = times > bin_edges(loop) & times <= bin_edges(loop+1);    
    y(loop) = sum(bin_times_idx);                   % total incidences in this bin across all sweeps
    these_events = events(bin_times_idx);           % retrieve all event indices
    mat_event = repmat([these_events inf]',1,num_sweeps);   % add Inf as a dummy event index - then event_count will always be array
    mat_sweeps = repmat(sweeps,y(loop)+1,1);
    event_count = sum(mat_event == mat_sweeps);
  
%     for loop2 = 1:num_sweeps
%         event_count(loop2) = sum(these_events == sweeps(loop2));    % array of event occurrence frequency in this bin
%     end
    var_psth(loop) = var(event_count);  
    
    if mean_bins
        % average by number of recording sweeps if required
        y(loop) = y(loop) ./ num_sweeps;
    end
end

% NOTE: could use Matlab's built-in histc function to do this, except that variance calculation would have to
% be done as above anyway.
%test = histc(times,bins) ./ num_sweeps;

if nargin >= 5 & findstr(varargin{1},'p')
    if nargin == 6
        text_add = varargin{2};
    else
        text_add = [];
    end
    
    figure
    bar(bins,y,1)
    title(['Peristimulus time histogram (bin size ' num2str(bin_size) ' seconds) ' text_add ])
    xlabel('time (seconds)');
    ylabel('spikes per bin');
end



