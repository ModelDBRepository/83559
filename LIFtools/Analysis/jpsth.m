function [N,x,D_uv,y] = jpsth(events1,times1,events2,times2,pairs,bin_size,T)

% JPSTH normalised joint peri-stimulus time histogram
%
%   [N,X,D,R] = JPSTH(E1,T1,E2,T2,P,B,T) where E1, E2 are arrays of all events to be correlated, T1, T2 are their corresponding time-stamp arrays,
%   P is a 2 column list of the event indices which match up (e.g. those corresponding to the same trial),  
%   B is the bin size (see note below), and T is the 2-element array of start and end times of the spike train(s) in seconds.
%
%   Returns the normalised JPSTH matrix N, the correspoinding axis time-stamps X, the predictor JPSTH matrix D, and the raw JPSTH matrix R.     
%
%   NOTE: for this to be most accurate, BINSIZE <= refractory period for the neuron or, if no refractory period, then BINSIZE = simulation
%   time-step because there must a maximum of one spike per bin from each individual spike train. Thus, when the bin average is taken
%   in the PSTH function, the result is a coefficient (i.e. in range 0,1)
%
%   REFERENCES: Aertsen, A. M. H. J., Gerstein, G. L., Habib, M. K. & Palm, G. (1989)
%               Brody, C. D. (1999) "Correlations without synchrony" Neural Computation, 11, 1537-1551 [for a clearer definition of normalised JPSTH]
%
%   Mark Humphries 22/12/2004

if T(1) >= T(2)
    error('Start time must be less than end time')
end

time_seconds = T(2) - T(1);
num_bins = floor(time_seconds / bin_size);
[num_pairs c] = size(pairs);

if c ~= 2
    error('Pairs matrix must have two columns');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate raw JPSTH; u constant over cols, v constant over rows (and extract times/events to be used)
t1 = [];
t2 = [];
e1 = [];
e2 = [];

y = zeros(num_bins,num_bins);

for loop = 1:num_pairs
    current_times1 = times1(events1 == pairs(loop,1));     % create numerical array of all times for this event in pair
    current_times2 = times2(events2 == pairs(loop,2));
    current_events1 = events1(events1 == pairs(loop,1));
    current_events2 = events2(events2 == pairs(loop,2));
    
    [s1 x] = spike_train_from_times(current_times1,bin_size,T);     % get time-array x here
	s2 = spike_train_from_times(current_times2,bin_size,T);
    
    s_mat1 = repmat(s1,num_bins,1);
    s_mat2 = repmat(s2',1,num_bins);
   
    y = y  + (s_mat1 .* s_mat2) / num_pairs;
    
    % store extracted time and event arrays
    t1 = [t1 current_times1];
    t2 = [t2 current_times2];
    e1 = [e1 current_events1];
    e2 = [e2 current_events2];
end

% compute individual PSTHs too
[psth1,x,v1] = LIF_psth(e1,t1,bin_size,T,'m');
[psth2,x,v2] = LIF_psth(e2,t2,bin_size,T,'m');

% expected JPSTH
mat1 = repmat(psth1',num_bins,1);
mat2 = repmat(psth2,1,num_bins);
E = mat1 .* mat2;

% cross-covariance (or predictor) JPSTH
D_uv = y - E;

% normalise
sd1 = sqrt(v1);
sd2 = sqrt(v2);

sd_mat1 = repmat(sd1',num_bins,1);
sd_mat2 = repmat(sd2,1,num_bins);
std_mat = sd_mat1 .* sd_mat2;

% normalised JPSTH
N = D_uv ./ std_mat;

% if any of std deviations were zero (i.e. if mean of bin was zero) then NaNs will result - explicitly
% set these to zero
N(isnan(N)) = 0;


