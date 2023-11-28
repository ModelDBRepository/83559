function varargout = plot_fI(tm,R,theta,max_A,abs_ref,varargin)

% PLOT_FI plots frequency-current relationship for basic LIF neuron
%
%   plot_fI(tm,R,theta,max_A,abs_ref) where
%       tm:         membrane time constant      (in seconds)
%       R:          resistance                  (in ohms)
%       theta:      threshold of firing         (in volts)          [value above resting potential (itself assumed to be 0)]
%       max_A:      maximum current input       (in amps)
%       abs_ref:    absolute refractory period  (in seconds)
%   
%   plot_fI(tm,R,theta,max_A,abs_ref,color) where
%       color:      standard color/shape specifying switch from PLOT command
%
%   Plots the frequency-current relationship for the LIF neuron specified by tm, R, theta, 
%   and the absolute refractory period. Displayed plot values given as "physiological" units.
%
%   [A,B] = plot_fi(...) returns the firing rate estimate A (in Hz) for each current injection step in array B (in amps)
%
%   Mark Humphries. Last rev: 17/12/204

if nargin == 6
    color = varargin{1};
else
    color = '';
end

% calculate f-I curve for 1000 current steps
step = max_A / 1000;
current = step:step:max_A;

out = 1 ./ (abs_ref + tm .* log(R .* current ./ (R.*current - theta))); 

% get rid of complex terms
complex_idx = imag(out);
out(complex_idx ~= 0) = 0;

% scale x-axis
nA_current = current .* 1e9;

% scale variables to physiological units for display
scale_tm = tm / 1e-3;                       % to milliseconds   
scale_R = R  / 1e6;                         % to mega Ohms
scale_theta = theta  / 1e-3;       % to millivolts
scale_abs_ref = abs_ref / 1e-3;   % to milliseconds


%%% plot result
figure
plot(nA_current,out,color)
title(['f-I plot for \tau_{m} = ' num2str(scale_tm) ' ms, R = ' num2str(scale_R) ' M \Omega, \theta = ' num2str(scale_theta) ' mV, \delta^{abs} = ' num2str(scale_abs_ref) 'ms.']);
xlabel('current (nA)')
ylabel('frequency (Hz)')
    
if nargout == 2
    varargout{1} = current;
    varargout{2} = out;
elseif nargout > 0
    disp('Incorrect number of output variables - 2 required')
end

    


