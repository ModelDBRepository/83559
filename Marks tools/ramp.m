function value = ramp(v,varargin)

%RAMP piece-wise linear threshold / rectifcation
% RAMP(V,M,T) thresholds value V using a piecewise linear ramp function (or rectification function) with slope M and threshold T.
%   M and T are optional arguments: their default values are 1 and 0 respectively; if T only required, set M = [].
%
% Mark Humphries 19/8/2004

M = 1;
T = 0;

if nargin == 2 M = varargin{1};
elseif nargin == 3 & isempty(varargin{1}) T = varargin{2};
elseif nargin == 3
    M = varargin{1};
    T = varargin{2};
end

% limit input to (0,1)
if v < 0 value = T;
elseif v > 1/M+T value = 1;
else
    value = M*(v-T);
end