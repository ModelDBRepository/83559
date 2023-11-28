function [coeffs_mapped,freqs,varargout] = LIF_wavelet(A,T,varargin)

%LIF_WAVELET wavelet analysis of LIF unit output
%
% LIF_WAVELET(A,T) where A is a single time-series of firing rates (either from ISIs or from moving-average) and T is the period of
%   simulation time (in seconds) covered by the time-series. Computes the wavelet function at scales 1:128 using the Morlet wavelet. 
%
%   LIF_WAVELET(A,T,FLAG) where FLAG is
%       'p': plots the wavelet analysis of the time-series, with calibrated colour bar and frequencies
%
%   [W,F] = LIF_WAVELET(...) returns the colour-mapped wavelet coefficient matrix W, and F the frequencies corresponding to the scales 
%   (values on Y axis when plotting W). W should be plotted using image(W). 
%
%   Mark Humphries 22/7/2004

scales = 1:256;
period = T / length(A);

coeffs = cwt(A,scales,'morl');

freqs = scal2frq(scales,'morl',period);

%% set colormap to jet
colormap('default')

NBC = 128;  %% number of colours

%% map values to full colour map
coeffs_mapped = wcodemat(coeffs,NBC,'mat',1);

if nargin >= 3 & findstr(varargin{1},'p')
    figure
    h = image(coeffs_mapped);
    p = get(h,'Parent');
    ticks = 16:16:256;
    set(p,'YTick',ticks)
    set(p,'YTickLabel',freqs(ticks));    
    varargout{1} = h;
end
    