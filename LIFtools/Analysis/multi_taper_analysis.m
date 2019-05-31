function [results, spectrum, fs_array, ] = multi_taper_analysis(time_array, dt, period, varargin)

%   MULTI_TAPER_ANALYSIS Multi taper periodogram analysis of
%   unevenly-sampled data (e.g. spikes)
%
%   [R,P,F] = SCARGLE_ANALYSIS(I,T,DT,Q) where I is an array of physical variable data (if I is a cell array, each cell is treated
%   as a separate record), T is the sampling time-stamps for I in seconds (again if cell array etc), DT is sample bin  
%   in the underlying data in seconds (i.e. for spike-trains, the sampling
%   time-step), and Q is the time period of the sampled data in seconds.
%   
%   SCARGLE_ANALYSIS(I,T,DT,Q,NS) reports the output at NS evenly-spaced frequency samples, up to the Nyquist frequency
%   (default uses nearest integer power of 2 above the number of samples in the data)
%   
%   SCARGLE_ANALYSIS(I,T,DT,Q,NS,R) where R is a 2-element array, sets the minimum and maximum frequencies to analyse, which will be
%   spaced by NS. Default is [2/Q Nyquist] - if min or max values are outside these ranges, they will be adjusted so that they are in the default range, and a
%   warning given.
%   
%   SCARGLE_ANALYSIS(I,T,DT,Q,NS,R,C) returns the power value at the confidence interval C (default is p=0.01, 99% confidence interval) 
%   for each data-set   
%
%   The results matrix R and periodogram matrix P each have one row per row of data set I. The values returned in R consist of:
%   the frequency of maximum power; the probability of this occuring by chance (according to method of Horne & Baliunas, 1986); the
%   significance level of this result (less than value given, or NaN for non-significance); and the power-level at the specified confidence 
%   interval (default p=0.01, as above). The final three columns will report signficance
%   results for the other method of testing significance in periodograms [not yet implemented].
%   [In fact currently returns the power at the peak, mean and std dev of spectrum]
%
%   F is an array of the frequencies used to construct the periodogram: use as x-axis values when plotting periodogram
%
%   REFERENCES: (1) Scargle, J. D. (1982). Studies in astonomical time series analysis. II. Statistical aspects of spectral analysis 
%                   of unevenly spaced data. Astrophysical Journal, 263, 835-853.
%               (2) Horne, J. H. & Baliunas, S. L. (1986). A prescription for period analysis of unevenly sampled time series. 
%                   Astrophysical Journal,  302, 757-763.
%               (3) ESO online notes - Spectral analysis and unevenly spaced data
%
%   Mark Humphries 17/1/2005

% number of datasets
if iscell(time_array)
    N = length(time_array);
else
    N = 1;
end

if nargin>=8 & ~isempty(varargin{3})
    C = varargin{3};
else
    C = 0.01;
end

sampling_frequency = 1 / dt;                            % the sampling frequency
nyquist = sampling_frequency / 2;                       % maximum frequency (Nyquist)

% array of sample points
%time_array = linspace(0,duration,nbins)';

% arbitrary number of frequency samples
if nargin >= 5 & ~isempty(varargin{1})
    number_frequency_samples = varargin{1};
else
    % default is number of samples rounded to next highest power of 2 (like
    % DFT)
    last_t = period / dt;
    np = nextpow2(last_t);
    number_frequency_samples = 2 ^ np;
end

%%% storage
store_power = zeros(N,number_frequency_samples);
fs_array = zeros(N,number_frequency_samples);
results = zeros(N,7);

for j = 1:N
    fprintf(1, 'periodogram loop count %d\n', j);

    if N > 1
        % assign current bin plot
        y = binplots{j}';
        t = time_array{j}';
    else
        y = binplots;
        t = time_array;
    end
    
    
    % check that this has enough data to meaningfully do something
    if length(y) < 12   % arbitrary number
        warning('Insufficient data to compute periodogram');
        fzero = nan;
        prob = nan;
        sig_horne = nan;
        L = nan;
        sig_peak = nan;
        psdxbar = nan;
        psdstd = nan;
        store_power(j,:) = nan;
        fs_array(j,:) = nan; 
    else
        % calculate frequencies
        nbins = length(t);                                      % number of time-stamps  
        lowest_f = 2 / (period / dt);                                  % lowest frequency resolvable
        
        if nargin >= 6 & ~isempty(varargin{2})
            FR = varargin{2};
            if FR(1) < lowest_f
                FR(1) = lowest_f;
                warning('Lowest frequency set to correct value')
            end
            if FR(2) > nyquist
                FR(2) = nyquist;
                warning('Highest frequency set to Nyquist')
            end
        else
            FR = [lowest_f nyquist]; 
        end
            
        % array of frequencies to calculate power for, from lowest resolvable to nyquist....
        fs_array(j,:) = linspace(FR(1),FR(2),number_frequency_samples);
        
		%%%% constant for use in false-alarm probability calculation
		%%%% this is a numerical fit to the number of independent frequencies in the data
		%%%% as determined by Horne and Baliunas (1986)
		Ni = -6.362 + 1.193.* nbins + 0.00098.* nbins .* nbins;
	
        %%%%%%% statistics
		mean_y = mean(y);
		variance_y = std(y) ^ 2;
		
		%DC_mag = y - mean_y; % for Lomb periodogram
        DC_mag = y;
       
		for i = 1:number_frequency_samples
            angular_frequency = fs_array(j,i) * 2 * pi;
            
            %%%%%%%%% calculate tau %%%%%%%%%%%%%%%%%%
            sine_term = sum(sin(2 .* angular_frequency .* t));
            cosine_term =  sum(cos(2 .* angular_frequency .* t));
            tau = atan(sine_term / cosine_term) ./ (2 * angular_frequency);
            
            %%%%%%% calculate power %%%%%%%%%%%%%%%%%%
            C_term = cos(angular_frequency .* (t - tau));
            S_term = sin(angular_frequency .* (t - tau));
	
            cos_top = sum(DC_mag .*  C_term);
            cos_top = cos_top .^2;
            cos_bottom = sum(C_term .^ 2);
            
            sin_top = sum(DC_mag .*  S_term);
            sin_top = sin_top .^2;
            sin_bottom = sum(S_term .^ 2);
            
            %power = (1/ (2 * variance_y)) * (( cos_top / cos_bottom ) + (sin_top / sin_bottom));    % Lomb's version
            power = 0.5 * ((cos_top / cos_bottom ) + (sin_top / sin_bottom)); % Scargle's version
            if isinf(power)
                error(['Power at frequency ' num2str(fs_array(j,i)) ' is infinite. Check time-series.']);
            end
            
            % store stuff
            store_power(j,i) = power;      
        end
        
        %%%% find the false-alarm probability for the maximum frequency....
        current_power = store_power(j,:);
        
        % what is max frequency?
        fzero_index = find(current_power == max(current_power));
        if isempty(fzero_index)
            fzero = 0;
        else
            fzero = fs_array(j,fzero_index);
        end
        % what is probability?
        z = max(current_power);
        prob = 1 - (1 - exp(-z)) ^ Ni;
        
        % determine significance level
        temp = find(prob < P);  % find which probability levels it is less than
        if isempty(temp)
            sig_horne = nan;
        else
            sig_horne = P(temp(end));   % is less than last probability in list
        end
        
        %% determine confidence intervals for all P values
        L = -log(1-(1-C)^(1/Ni));       % re-arrangement of above probability expression
        
        %%%% find the significance level for the maximum frequency...
        psdxbar = mean(current_power);
        psdstd = std(current_power);
        
        %% significance level of peak....
        sig_peak = z > psdxbar + 3 * psdstd;
        above = (z - psdxbar) / psdstd;
    end
    
    % store this data
    results(j,1) = fzero;
    results(j,2) = prob;
    results(j,3) = sig_horne;
    results(j,4) = L;
    results(j,5) = sig_peak;
    results(j,6) = psdxbar;
    results(j,7) = psdstd;

end
