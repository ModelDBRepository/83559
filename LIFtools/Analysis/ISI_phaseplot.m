function [h] = ISI_phaseplot(ISIs,varargin)

% ISI_PHASEPLOT 2D phase-plot
%   ISI_PHASEPLOT(I) plots the inter-spike intervals (ISIs) in array I on a 2D phase-plot for (t,t+1); if I is a cell array,
%   then each ISI set is plotted in turn.
%
%   ISI_PHASEPLOT(I,N) plot the phase-plot on Figure N - if it already exists, then it is assumed that the plot is to be added
%   thus allowing a phase-plot to be built up in batch processing (note: if N is a string or the empty matrix, then a new figure
%   window is opened).
%
%   ISI_PHASEPLOT(I,N,FLAG) where FLAG is
%       't': plots the trajectory of the phase-plot in addition to the scatter-points
%
%   Returns the handle to the figure - this can then be used to iteratively call the function with the same figure window
%
%   Mark Humphries 30/11/2004

if nargin >=2 & isnumeric(varargin{1}) & ~isempty(varargin{1})
    h = figure(varargin{1});
else
    h = figure;
end
 
if nargin >= 3 & findstr(varargin{2},'t')
    strPlot = '.-';
else
    strPlot = '.';
end
    
if iscell(ISIs)
    num_loops = length(ISIs);
else
    num_loops = 1;
    temp = {ISIs};
    ISIs = temp;
end

for loop = 1:num_loops
    isis = ISIs{loop};
    if length(isis > 1)
        isi_t1 = isis(1:end-1);
        isi_t2 = isis(2:end);
        plot(isi_t1,isi_t2,strPlot)
        hold on
    end
end
drawnow
title('ISI phase plot')
xlabel('ISI at t (s)')
ylabel('ISI at t+1 (s)')