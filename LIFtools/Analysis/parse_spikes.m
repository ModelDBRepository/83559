function [new_times,idxs] = parse_spikes(times,varargin)

% PARSE_SPIKES remove anomalous spike events
%   [NT,I] = PARSE_SPIKES(T) will parse the received spike train time-stamp T, removing all anomalous spikes
%   here defined as any ISIs < 0.001 s (i.e. those spikes which imply a firing rate of > 1000 Hz). 
%   Returns the new time-stamps (NT) and index of removed spikes
%
%   PARSE_SPIKES(T,THETA) sets the anaomalous spike-interval limit to THETA seconds
%
%   Mark Humphries 3/12/2004

if nargin>=2 & isnumeric(varargin{1})
    threshold = varargin{1};
else
    threshold = 0.001;
end

ISIs = abs(diff(times));

idxs = find(ISIs < threshold);
idxs = idxs + 1;    % remove second spike

new_times = times;

new_times(idxs) = [];