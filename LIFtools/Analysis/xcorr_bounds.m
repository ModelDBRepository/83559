function [j_max,j_min] = xcorr_bounds(F1,F2,bin_size)

% XCORR_BOUNDS cross-correlation upper and lower boundaries
%
%   XCORR_BOUNDS(F1,F2,BINSIZE) where F1 and F2 are the mean firing rates of the two correlated
%   spike trains, and BINSIZE is the bin-size of the cross-correlogram (really the cross-covariogram)
%   Returns max and min values for the normalised covariogram.
%
%   Reference: Dorn, J. D. & Ringach, D. L. (2003) "Estimating membrane voltage correlations from extracellular
%              spike trains", Journal of Neurophysiology, 89, 2271-2278.
%
%   Mark Humphries 21/4/2004


prob1 = F1 * bin_size;
prob2 = F2 * bin_size;

j_max = (min(prob1,prob2) - prob1*prob2) / sqrt(prob1*(1-prob1) * prob2*(1-prob2));
j_min = (-prob1*prob2) / sqrt(prob1*(1-prob1)*prob2*(1-prob2));
