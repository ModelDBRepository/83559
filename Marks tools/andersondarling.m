function [A,H,L] = andersondarling(data,dist,varargin)

% ANDERSONDARLING Anderson-Darling goodness-of-fit for distributions
%   [A,H,L] = ANDERSONDARLING(D,T) given an array of data values D, and a
%   distribution type T (see below), computes the Anderson-Darling test
%   statistic A. Also returns the truth value H 
%   for rejection of the null hypothesis, for a given significance 
%   level (default is 0.05), and the approximate significance level L
%   attained
%
%   ANDERSONDARLING(...,ALPHA) sets the significance level to ALPHA
%
%   NOTES:
%   (1) Anderson-Darling critical values only defined for alpha value set
%   {0.1 0.05 0.025 0.01}. The critical values are for n->infinity. Minor
%   corrections required for n <= 100 (see Stephens, 1974).
%
%   (2) Currently, only distribution type 'normal' is defined
%
%   Reference: 
%   (1) will be D'Agostino & Stephens (1986) Goodness-of-fit
%   Techniques. M. Dekker: New York. (need this to get other distributions,
%   small n values, check correction, critical values etc)
%
%   (2) Stephens, M. A. (1974). EDF Statistics for Goodness of Fit and Some
%   Comparisons. Journal of the American Statistical Association, 69, 730-737.
%
%   (3) Nelson, L. S. (1998). The Anderson-Darling test for normality.
%   Journal of Quality Technology, 30, 298-299.
%
%   Mark Humphries 5/2/2006

alpha = 0.05;   % default significance level

if nargin >= 3 alpha = varargin{1}; end

% critical values for the normal distribution, with n -> infinity
% (Stephens, 1974)
%           alpha  critical value
norm_crit = [0.1   0.656; ...
             0.05  0.787; ...
             0.025 0.918; ...
             0.01  1.092];

% check that alpha-level is defined
sig_level = find(alpha == norm_crit(:,1));
if isempty(sig_level)
    error('Anderson-Darling test undefined for that alpha level');
end


% number of data points
n = length(data);

%% order data
[ranked_data idx] = sort(data);

% compute A, compare to critical values
switch dist
    case 'normal'
        mu = mean(data);
        sigma = std(data);
        
        % vectorised form of A-D statistic calculation
        cdf_1 = normcdf(ranked_data,mu,sigma);
        cdf_2 = flipud(cdf_1);
        idxs = 2 * (1:n) - 1;
        S = sum(idxs' .* (log(cdf_1) + log(1 - cdf_2)));

% Keep this for reader reference  
%         S = 0;
%         for loop = 1:n
%             cdf_1 = normcdf(ranked_data(loop),mu,sigma);
%             cdf_2 = normcdf(ranked_data(n-loop+1),mu,sigma);
%             S = S + (2*loop-1) * (log(cdf_1) + log(1-cdf_2)); 
%         end
        
        A = -n - S/n;
        
        % apply correction as both mean and std derived from sample (from
        % D'Agostino & Stephens, 1986)
        A = A * (1 + 0.75/n + 2.25/n^2);
        
        % compute approximate significance level (Nelson, 1998)
        L = 3.6789468 * exp(-A / 0.1749916);    
        
        % check for significance
        if A > norm_crit(2,sig_level)
            H = 1;
        else
            H = 0;
        end
        
    otherwise
    error('Distribution not supported');
end
