function [in_n,in_t] = BG_GHS_inputs(dt,time_steps,n_inputs,switches,hz,n_neurons,neurons_per_channel,n_channels,ref_period,type)
% BG_GHS_INPUTS spike train inputs
%   BG_GHS_INPUTS(DT,T,I,S,HZ,N,C,NC,R,FLAG) where DT is the simulation time-step (in seconds), T is the number of time-steps in the simulationm
%   I is the number of inputs, S is the array of switching times, HZ is the matrix of firing frequencies (length(S) x NC), N is the number
%   of simulated input neurons, C is the number of neurons per channel, NC is the number of channels, and FLAG is one of:
%       'vivo'      -   Poisson process, irregular tonic firing
%       'organo'    -   organotypic-like correlated tonic firing, specify a single HZ value or just use first HZ value in matrix
%       'slow'      -   anaesthesia-like slow-wave firing
%
%   
%   Mark Humphries & Rob Stewart 23/12/2004

switch type
case 'vivo'
	old = 0;
	in_n = [];
	in_t = [];
	
	% in vivo normal tonic - Poisson processes
	for i = 1:length(switches)
        h = hz(i,:);            %channel frequencies for this switch point
        time_seconds = (switches(i)-old)*dt;
        num_isis = ceil(time_seconds * max(h) * 2);
        isi = [];
        for j = 1:n_channels
            isi = [isi,random('exp',1/h(j),num_isis,neurons_per_channel)];
        end    
        isi(isi<ref_period) = ref_period;
       
        intervals = round(isi/dt);
        if num_isis > 1
            spike_times = cumsum(intervals)+old; %only on a matrix
        else
            spike_times = intervals+old;
        end    
        spike_times = spike_times';
        [t1,t2] = find(spike_times<(switches(i)+2)); %shifted forward two (20/02/04)
        if ~isempty(t1)
            ind = sub2ind(size(spike_times),t1,t2);
            times = spike_times(ind);
            [in_t_temp,ord] = sort(times);   %sort into time order
            in_n_temp = t1(ord);
            in_t = [in_t;in_t_temp];
            in_n = [in_n;in_n_temp];
        end
        old = switches(i);
	end
	in_n = uint32(in_n+n_neurons-1);
	in_t = uint32(in_t-2); %not quite sure why this is two...

case 'organo'
% organotypic - correlated tonic 
% NOTE: in_t must be in time-order
isi = 1/hz(1)/dt;                                  % ISI in time-steps
regular_times = round(isi:isi:time_steps);      % template train spike-times

% create jittered versions for each input
reg_mat = repmat(regular_times,n_inputs,1);
jit_std = isi * 0.3;                            % 30% jitter - reduces with *increasing* Hz
jitter = round(randn(n_inputs,length(regular_times)) * jit_std);  % normally-distributed tonic firing - a la SNc neuron
trains = reg_mat + jitter;                      % jittered trains for each input
[in_t ord] = sort(trains(:));                   % sort in to time order

event_idxs = (1:n_inputs)+n_neurons-1;          % input indices are after all other indices (shift back by 1 for MEX indexing)
events = repmat(event_idxs',1,length(regular_times));
events = events(:);
in_n = events(ord);                             % put input indices in same order as time of spike

% remove all out-of-range events
out_of_range = find(in_t < 0 | in_t > time_steps);
in_t(out_of_range) = [];
in_n(out_of_range) = [];

in_n = uint32(in_n);
in_t = uint32(in_t);

case 'slow'
% in vivo anaesthetic - slow-wave, like slow-wave sleep, at ~1Hz
% up-period: correlated firing at ~12 Hz (Steriade et al 2001) [KG:
% shouldn't it be 24Hz because thisis the *mean* rate over both states?]
% down-period: complete silence
% desynchronised: during ECoG recording, wave spontaneously desynchonises, suggesting that underlying cortical firing
%   is no longer correlated; however, Kasanetz et al. (2002) data shows that striatal neurons remain in their up-state during
%   this period. Therefore, if up/down state transition does follow cortical slow-wave exactly, as seems to be the case, then
%   this suggests that desynchronised ECoG corresponds to decorrelated cortical firing at a frequency at least equivalent to the 
%   slow-wave up-period

% To model: alternate periods of silence and correlated firing (a la organotypic) at ~ 1Hz, randomly adding
% Poisson-generated section of desynchronisation if requested by user

time_seconds = time_steps * dt;
half_second = 0.5 * time_steps / time_seconds; 
in_t = [];
in_n = [];

slow_jitter = 2.5; % jitter - reduces with *increasing* Hz (try 2.6 with 24Hz input)


%%%%%% KG mod: vary intensity of each coherent up-state in cortex - No
%%%%%% spikes in each upstate is same for all inputs within each period
slow_std_hz = 0.15; % * 100 % std with normally distributed rate 

for loop = 1:time_seconds*2   % silent on odd, correlated on even (so starts on silence)
    if ~mod(loop,2)         % is even - do half second of correlated firing  
        
        %%%%%%  KG mod %%%%%%%%%%%%
        hz_eff = hz(1) .* (1 +  slow_std_hz .* randn);
        isi = 1/hz_eff/dt;                                       % ISI in time-steps
        %%%%%%% KG mod end %%%%%%%
        
        
        regular_times = round(isi:isi:half_second) + half_second.*(loop-1);      % template train spike-times for this second
        
       
        % create jittered versions for each input
		reg_mat = repmat(regular_times,n_inputs,1);
		jit_std = isi * slow_jitter;                            % 30% jitter - reduces with *increasing* Hz
		jitter = round(randn(n_inputs,length(regular_times)) * jit_std);  % normally-distributed tonic firing - a la SNc neuron
		trains = reg_mat + jitter;                      % jittered trains for each input
		[temp ord] = sort(trains(:));                   % sort in to time order
		in_t = [in_t; temp]; 
        
		event_idxs = (1:n_inputs)+n_neurons-1;          % input indices are after all other indices (shift back by 1 for MEX indexing)
		events = repmat(event_idxs',1,length(regular_times));
		events = events(:);
		temp = events(ord);                             % put input indices in same order as time of spike
        in_n = [in_n; temp];          
    end
end
% remove all out-of-range events
out_of_range = find(in_t < 0 | in_t > time_steps);
in_t(out_of_range) = [];
in_n(out_of_range) = [];

% sort ready for coincidence removal
[ts tn] = sort(in_t);
in_t = ts;
in_n = in_n(tn);
% remove duplicate events
No_ins = length(in_t);
if No_ins > 1
    No_coincidents = 0;
    coincident_indices = [];
    input_t1 = in_t(1);
    input_n1 = in_n(1);
    for j=2:No_ins
        input_t2 = in_t(j);
        input_n2 = in_n(j);
        if input_t2 == input_t1 & input_n1 == input_n2
            No_coincidents = No_coincidents + 1;
            coincident_indices = [coincident_indices j];
        end
        input_t1 = input_t2;
        input_n1 = input_n2;
    end
    in_t(coincident_indices) = [];
    in_n(coincident_indices) = [];
    fprintf(1, 'There were %d coincident events removed from inputs\n', No_coincidents);
end
in_n = uint32(in_n);
in_t = uint32(in_t);

end % switch
