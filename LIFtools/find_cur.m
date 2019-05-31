function I = find_cur(tm,R,theta,abs_ref,fire,varargin)

% FIND_CUR determines the current injection needed for a given firing rate for basic LIF neuron
%
%   find_cur(tm,R,theta,abs_ref,fire)
%       tm:         membrane time constant      (in seconds)
%       R:          resistance                  (in ohms)
%       theta:      threshold of firing         (in volts)          [value above resting potential (itself assumed to be 0)]
%       abs_ref:    absolute refractory period  (in seconds)
%       fire:       target firing rate          (in Hz)
%
%   find_cur(tm,R,theta,abs_ref,fire,'d') 
%       the 'd' switch displays the required current in "physiological" units to screen
%
%   Returns the value of the current (in amps) required to drive the LIF neuron to achieve the specified
%   firing rate. If the specified firing rate cannot be achieved with the given parameters, then
%   returns NaN and a warning given.
%
%   Mark Humphries. Last rev: 17/12/2004

I = -theta / (R * exp(-(1/tm * (1/fire - abs_ref))) - R);

if I < 0
    I = nan;
end

if nargin == 6
    if isnan(I)
        disp(['Cannot achieve firing rate of ' num2str(fire) 'Hz with the specified LIF neuron']);
    else
        disp(['A current of ~' num2str(I * 1e9) 'nA is required for an output of ' num2str(fire) 'Hz (for specified LIF neuron)']); 
    end
end