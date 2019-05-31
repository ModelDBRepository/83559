function I = find_Vm_cur(R,inj,target,varargin)

% FIND_VM_CUR determines the current injection needed for a given membrane potential for a basic LIF neuron
%
%   find_Vm_cur(R,I,target)
%       R:          resistance                  (in ohms)
%       I:          constant injection current  (in volts)   [e.g. spontaneous currents]
%       target:     target membrane potential   (in volts)   [below resting potential, as that is zero]
%  
%   find_Vm_cur(R,I,target,'d') 
%       the 'd' switch displays the required current in "physiological" units to screen
%
%   Returns the value of the current (in amps) required to drive the LIF neuron to achieve the specified
%   membrane potential. 
%
%   Mark Humphries. Last rev: 17/12/2004

I = (target / R) - inj;

if nargin == 4
    disp(['A current of ~' num2str(I * 1e9) 'nA is required for a membrane potential of ' num2str(target) ' V (for specified LIF neuron)']); 
end