function [bursts,burst_isis] = sncburst(isis,times)

% SNCBURST SNc neuron burst criteria
%   [B,L] = SNCBURST(I,T) computes the start and end-points of SNc neuron bursts from the inter-spike interval 
%   array I (in seconds) and the time-stamp array T (in seconds). It returns a 2-column matrix with the start and end-point
%   of each burst B (times of the corresponding spikes) and a cell array L of the ISIs in each burst.
%   If no bursts are found the empty matrix is returned for both. If a final burst occurs with no end-point, it is omitted 
%   from the list (to avoid is effects on any further burst processing).
%
%   REFERENCE: Grace, A. A. & Bunney, B. S. (1984). The control of firing pattern in nigral dopamine neurons: burst firing. J Neurosci, 4, 2877-2890.
%
%   Mark Humphries 15/12/04

if max(isis) > 2    % then could be firing rate data by mistake
    warning('Check that the ISI data is in seconds')
end

    
start_thresh = 0.08;    %  80 ms
end_thresh = 0.16;      % 160 ms
counter = 0;
bursts = [];
this_burst = [];
burst_isis = {};
start_flag = 0;
end_flag = 0;

for loop = 1:length(isis)
    if isis(loop) < start_thresh & ~start_flag      % meets criterion and start not found yet
        t_start = times(loop);
        this_burst = isis(loop);     
        start_flag = 1; % start of next burst found
        counter = counter + 1;
    elseif isis(loop) < end_thresh & start_flag     % start found but end not yet found
        this_burst = [this_burst isis(loop)];       % accumulate intra-burst ISIs
    elseif isis(loop) > end_thresh & start_flag     % start found and end criterion met
        bursts = [bursts; t_start times(loop)];
        burst_isis{counter} = this_burst;
        start_flag = 0;
    end
end
