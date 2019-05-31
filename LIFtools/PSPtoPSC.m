function weight = PSPtoPSC(ts,tm,R,peak,type,varargin)

% PSPtoPSC calculation of PSC weight
%
%   PSPtoPSC(ts,tm,R,peak,type) where
%       ts:     synaptic time constant  (in seconds)
%       tm:     membrane time constant  (in seconds)
%       R:      resistance              (in ohms)
%       peak:   value of PSP peak       (in volts)
%       type:   form of PSC to fit      ('step','alpha','compact alpha')  
%
%   Returns the weight required to scale the specified PSC function 
%   to fit the PSP shape specified by ts and tm. Decay of the PSC is
%   dictated by ts, which is the rise time of the PSP; tm dictates the 
%   decay time of the PSP.
%
%   Type:  
%                'step'        1/ts * exp(-s/ts) H(s)
%               'alpha'        1/(ts-tr) * (exp(-s/ts) - exp(-s/tr)) H(s)   (not yet supported)
%       'compact alpha'        s/ts^2 * exp(-s/ts);                         (not yet supported) (suitable for tr -> ts, see Gerstner & Kistler (2002) "Spiking Neuron Models" Cambridge: CUP)

if nargin > 5
    ts2 = varargin{1};
end

% start with unit weight
c = 1;

switch type
case 'step'
	%% find peak of PSP for these time constants
	f = inline('- R*c/(ts-tm) * (exp(-s/ts) - exp(-s/tm))','s','c','tm','ts','R');
    
    peak_time = fminbnd(f,0,tm,[],c,tm,ts,R);
    %keyboard
    
	%% use to find suitable weight c for peak value
    weight = peak * (ts - tm) / (R * (exp(-peak_time/ts) - exp(-peak_time/tm)));
    %keyboard
case 'combined'
	%% find peak of PSP for these time constants
	f = inline('- R*c/(ts-tm) * (exp(-s/ts) - exp(-s/tm)) - R * c/(ts2-tm) * (exp(-s/ts2) - exp(-s/tm))','s','c','tm','ts','ts2','R');
    
    peak_time = fminbnd(f,0,tm,[],c,tm,ts,ts2,R);
    %keyboard
    
	%% use to find suitable weight c for peak value
    weight = peak / (R/(ts-tm) * (exp(-peak_time/ts) - exp(-peak_time/tm)) +...
                        R/(ts2-tm) * (exp(-peak_time/ts2) - exp(-peak_time/tm)));
    
case 'alpha'
    disp('Alpha function not yet supported')
    
case 'compact alpha'
	disp('Compact alpha function not yet supported')
end


% function f = get_PSP_peak(s,c,tm,ts,R)
%     %% minimising this function will return the value of s at which the PSP peak occurs
%     f = - (R*c/tm) * (1/(1-(ts/tm)) * (exp(-s/tm) - exp(-s/ts)));
%     
% function f = get_PSP_weight(c,s,tm,ts,R,peak)        
%     %% s is peak time; finding the zero-crossing of this function will give the PSC weight
%     f = peak - ((R*c/tm) * (1/(1-(ts/tm)) * (exp(-s/tm) - exp(-s/ts))));



    