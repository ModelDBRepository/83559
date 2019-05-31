function [bursts,burst_isis] = kbsta(times,varargin)

%KBSTA Knowledge Based Spike Train Analysis
%   [B,I] = KBSTA(T) where T is an array of spike time-stamps in seconds, finds bursts in the spike-train based on 
%   fuzzy-logic detection methods derived from human observer behaviour. Returns B, a 2-column array of start and end times (spikes) for each
%   burst, and a cell array I of the inter-spike intervals in each burst. If no bursts are found the empty matrix is returned for both. 
%   If a final burst occurs with no end-point, it is omitted from the list (to avoid is effects on any further burst processing).
%
%   KBSTA(T,R,ON,OFF) specifies the parameters of the burst detection method: R is the size of the firing-rate estimate window in seconds, ON and OFF are
%   the fuzzy set-membership detection thresholds for burst-on and burst-off, respectively. Set any of these to [] for default values
%   (which are R=0.2, ON = OFF = 0.3);
%
%   REFERENCE: Xu, Z. M., Ivanusic, J. J., Bourke, D. W., Butler, E. G. & Horne, M. K. (1999). Automatic detection 
%                of bursts in spike trains recorded from the thalamus of a monkey performing wrist movements. J Neurosci Methods, 91, 123-133.
%
%   Mark Humphries 15/12/04

% deal with empty times
if length(times) < 2
    bursts = [];
    burst_isis = {};
    return
end

% default parameters
rate_window = 0.2;              % size of forward/backward firing rate window in seconds
eta_on = 0.3;                   % increase to make burst onset detection more stringent
eta_off = 0.3;                  % increase to make burst offset detection more stringent

if nargin >= 2 & ~isempty(varargin{1}) rate_window = varargin{1}; end
if nargin >= 3 & ~isempty(varargin{2}) eta_on = varargin{2}; end
if nargin >= 4 & ~isempty(varargin{3}) eta_off = varargin{3}; end

% variables
n_spikes = length(times)-2;     % look at all spikes except first and last if possible

start_temp = find(times-rate_window > times(1));
end_temp = find(times + rate_window < times(end));
if isempty(start_temp) | isempty(end_temp)          % then is insufficient data to do burst detection
    bursts = [];
    burst_isis = {};
    return
end

n_start = start_temp(1);              % index of time array to start at - first point after rate-window width from beginning    
n_end = end_temp(end);

if n_start < 2 n_start = 2; end                             % can't start from first spike
if n_end > length(times)-1 n_end = length(times) - 1; end   % or end on last spike

% storage
FFR = zeros(n_spikes,1);
BFR = zeros(n_spikes,1);
Mon = zeros(n_spikes,1);
Moff = zeros(n_spikes,1);

%%%%%%%%%%%%%%%%%%%%% method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate all ISIs
num_hist_bins = 50;
[firing_rate,isihist,isi_times,fshist] = LIF_ISI_analysis(times,num_hist_bins);
isis = 1 ./ firing_rate;
mean_isi = mean(isis);

% calculate firing rate
firing_rate(n_start-1) = length(find(times < times(n_start) & times >= times(n_start)-rate_window)) / rate_window; % do initial backward firing rate
for loop = n_start:n_end
    FFR(loop) = length(find(times > times(loop) & times <= times(loop)+rate_window)) / rate_window;
    BFR(loop) = length(find(times < times(loop) & times >= times(loop)-rate_window)) / rate_window;
end

max_rate = max([max(FFR) max(BFR)]);

% calculate fuzzy set membership functions - and determine bursting as we go
start_flag = 0;
counter = 0;
bursts = [];
burst_isis = {};

for loop = n_start:n_end
    mu_FFR = mem_function(FFR(loop),max_rate);
    mu_BFR = mem_function(BFR(loop),max_rate);
    mu_FTI = mem_function(isis(loop),mean_isi);     % isi is for the current spike
    mu_BTI = mem_function(isis(loop-1),mean_isi);
    
    Mon(loop) = mu_FFR * (1-mu_BFR) * (1-mu_FTI) * mu_BTI;
    Moff(loop) = (1-mu_FFR) * mu_BFR * mu_FTI * (1-mu_BTI);
    if Mon(loop) > eta_on & ~start_flag
        % start of burst
        t_start = times(loop);
        this_burst = isis(loop);     
        start_flag = 1; % start of next burst found
        counter = counter + 1;
    elseif Moff(loop) < eta_off & start_flag
        % continuation of burst
        this_burst = [this_burst isis(loop)];
    elseif Moff(loop) > eta_off & start_flag     
        % end criterion met during burst
        bursts = [bursts; t_start times(loop)];
        burst_isis{counter} = this_burst;
        start_flag = 0;
    end
end

function mu = mem_function(xi,ci)
    % where XI is the detector value and CI the appropriate trial-dependent constant for that detector
    if xi > 0.9 * ci
        mu = 1;
    elseif xi < 0.1 * ci
        mu = 0;
    else
        mu = xi / ci;
    end