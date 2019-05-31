function [jpsth_xcorr,corr_x,jpsth_coincidence] = jpsth_data(J,bin_size,corr_width)

% JPSTH_DATA extract cross-correlogram and coincidence data
%
%   JPSTH_DATA(J,B,W) where J is a normalised JPSTH matrix, B is the bin size of the JPSTH (in seconds), and W is the width of one tail
%   of the correlogram in bins (e.g.if W = 10, the resulting correlogram will be 21 bins wide: 2 * W + 1, for the centre bin).
%
%   [JX,X,JC] = JPSTH_DATA(...) returns the cross-correlogram JX, and the time-shift for each bin in the array X (plot with bar(X,JX)),
%   and the coincidence histogram JC (plot against x-axis for JPSTH)
%
%   REFERENCE: Aertsen, A. M. H. J., Gerstein, G. L., Habib, M. K. & Palm, G. (1989)
%   
%   Mark Humphries 23/12/2004

[num_bins, c] = size(J);

%%%%%%%%%%%%%% extract the cross-correlogram - sum over the diagonals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centre_bin = corr_width + 1;
num_corr_bins = 2 * corr_width + 1;

jpsth_xcorr = zeros(num_corr_bins,1);

% do centre first
diag = eye(num_bins); 
jpsth_xcorr(centre_bin) = (sum(sum(J .* diag))) / num_bins;

% coincidence histogram
jpsth_coincidence = J(diag==1); 

L_temp = J;
R_temp = J;
for loop = 1:corr_width
    eye_size = num_bins - loop;
    diag = eye(eye_size);   % create left-to-right identity matrix to retrieve centre diagonal
   
    % do right-of-centre
    R_temp = R_temp(1:end-1,2:end);        % create smaller matrix so left-of-centre becomes centre
    jpsth_xcorr(centre_bin+loop) = (sum(sum(R_temp .* diag))) / eye_size;      % normalise by number of bins
    % do left-of-centre
    L_temp = L_temp(2:end,1:end-1);        % create smaller matrix so right-of-centre becomes centre
    jpsth_xcorr(centre_bin-loop) = (sum(sum(L_temp .* diag))) / eye_size;      % normalise by number of bins
    
end
corr_x = -bin_size * corr_width:bin_size:bin_size * corr_width;


