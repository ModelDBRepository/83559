function r_fname = BATCH_BG_heterogenous_AMPA_NMDA(index,pars_file,seed1,pathroot,exp_name)

% BATCH_BG_hetereogenous_AMPA_NMDA(I,F,S,PATH,NAME), where I is the index into the input grid
% (give dummy value if not needed), 
% and F is a string specifying the filename of the parameters file necessary for 
% the current simulation, S is the random number seed, PATH is not used at the moment, and
% NAME is the prefix for the results file (essential for batch work on the cluster)  
%
  % Mark Humphries, Rob Stewart & Kevin Gurney, 29/1/2006


% warning off MATLAB:divideByZero;

%%%%% generate file name %%%%%%%%%%%%%%%%%
time_now = clock;
unique_name = datestr(time_now,30);  
r_fname = [exp_name '_' unique_name '_results'];
r_fname = [pathroot r_fname '.mat']; % not needed?
 
%% START PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
eval(pars_file);

hz = ones(3,n_channels)*tonic_rate;
switch input_method
    case 'tonic'
        hz = ones(3,n_channels)*tonic_rate;
        input_array = [];
    case 'switch'
        load('input_grid'); 
        hz(2:3,1) = rate_scaling .* input_array(index,1);    
        hz(3,2) = rate_scaling .* input_array(index,2); 
    case 'simultaneous'
        load('input_grid'); 
        hz(2,1) = rate_scaling .* input_array(index,1);  
        hz(2,2) = rate_scaling .* input_array(index,2);    
    otherwise
        error('unknown input method')
end


%% END PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys_time = 1e3;     % scale time values by this to get physiological units (milliseconds)
phys_R = 1e-6;             % scale R values by this to get physiological units (mega-ohms)

mean_tau_m = cell2mat(mean_tau_m);
mean_R = cell2mat(mean_R);

tau_AMPA = ones(n_neurons,1)*mean_tau_AMPA;
tau_NMDA = ones(n_neurons,1)*mean_tau_NMDA;
tau_GABAa = ones(n_neurons,1)*mean_tau_GABAa;
tau_m = zeros(n_neurons,1);

% do tau checks
check1 = find(mean_tau_m == mean_tau_AMPA);
check2 = find(mean_tau_m == mean_tau_GABAa);
if check1 error('AMPA time constant same as a membrane time constant'); end
if check2 error('GABAa time constant same as a membrane time constant'); end

% create tau m distributions: use Gamma distribution-based noise where quoted 'mean' values seem
% extreme

tau_m(SD1) = mean_tau_m(1) + std_tau_m(1) * randn(neurons_per_nucleus,1);
tau_m(SD2) = mean_tau_m(2) + std_tau_m(2) * randn(neurons_per_nucleus,1);

STN_distribution = random('Gamma',mean_tau_m(3)*phys_time,1,neurons_per_nucleus,1);   % generate physiological unit (integer) values 
tau_m(STN) =  STN_distribution / phys_time;                                           % then convert back
%tau_m(STN) = mean_tau_m(3) + std_tau_m(3) * randn(neurons_per_nucleus,1);

tau_m(GPe) = mean_tau_m(4) + std_tau_m(4) * randn(neurons_per_nucleus,1);
tau_m(GPi) = mean_tau_m(5) + std_tau_m(5) * randn(neurons_per_nucleus,1);

% create R distributions - again use Gamma for R (use same skew - assume capacitance fixed..)
% where appropriate
R = zeros(n_neurons,1); 
R(SD1) = mean_R(1) + std_R(1) * randn(neurons_per_nucleus,1);
R(SD2) = mean_R(2) + std_R(2) * randn(neurons_per_nucleus,1);

R(STN) = mean_R(3) + std_R(3) * randn(neurons_per_nucleus,1);

%% use same distribution - assuming relatively fixed C, increased Tm = increased R...
%temp = (mean_R(3) * phys_R) .* ones(neurons_per_nucleus,1) - (mean_tau_m(3)*phys_time);      % convert to phys values, and shift base 
%R(STN) = (temp + STN_distribution) / phys_R;                                                % apply tau_m skew and re-scale

R(GPe) = mean_R(4) + std_R(4) * randn(neurons_per_nucleus,1);
R(GPi) = mean_R(5) + std_R(5) * randn(neurons_per_nucleus,1);

% Discrete-time exact solution coefficients
e_AMPA = exp(-dt./tau_AMPA);
e_NMDA = exp(-dt./tau_NMDA);
e_GABAa = exp(-dt./tau_GABAa);
em = exp(-dt./tau_m);

% ae1 = R./tau_m .* (-tau_m ./(tau_m +tau_e1));
% ae2 = R./tau_m .* (-tau_m ./(tau_m +tau_e2));
% ai = R./tau_m .* (-tau_m ./(tau_m +tau_i));

a_AMPA = R./(tau_AMPA - tau_m);
a_NMDA = R./(tau_NMDA - tau_m);
a_GABAa = R./(tau_GABAa - tau_m);

a_AMPA = a_AMPA .* (e_AMPA-em); % pre-multiplied full solution!
a_NMDA = a_NMDA .* (e_NMDA-em); % pre-multiplied full solution!
%ae = ae1 + ae2;

a_GABAa = a_GABAa .* (e_GABAa-em); % pre-multiplied full solution!
as = R .* (1-em);

% structural constants
weights = [SD1_w SD2_w STN_GPew STN_GPiw GPe_STNw GPe_GPiw GPe_GPew GPi_GPiw EXT_w STN_ext_ratio]; 
delays = [SD12GPi_d SD22GPe_d STN2GPe_d STN2GPi_d GPe2STN_d GPe2GPi_d GPe2GPe_d GPi2GPi_d EXT2SD1_d EXT2SD2_d EXT2STN_d];
proportions = [SD12GPi_p; SD22GPe_p; GPe2STN_p; GPe2GPi_p; GPe2GPe_p; GPi2GPi_p];

n_inputs = neurons_per_nucleus;         %same input signal goes to SD1, SD2, STN
n_sources = n_neurons+n_inputs; %total spike sources

time_steps = round(time_seconds/dt);            

%%%%%% burst pseudo-current values %%%%%%%%
n_ca_cells = length(ca_cells);
burst_t1 = zeros(n_neurons,1);
burst_t2 = zeros(n_neurons,1);
alphaCA = zeros(n_neurons,1);
thetaCA = zeros(n_neurons,1);
c_cells = uint32(zeros(n_neurons,1));

c_cells(ca_cells) = 1;
burst_t1(ca_cells) = mean_t1 + std_t1 .* randn(n_ca_cells,1); % step-period of burst current
burst_t2(ca_cells) = mean_t2 + std_t2 .* randn(n_ca_cells,1); %  ramp-period of burst current

%keyboard

alphaCA(ca_cells) = mean_alphaCA + std_alphaCA .* randn(n_ca_cells,1);
thetaCA(ca_cells) = mean_thetaCA + std_thetaCA .* randn(n_ca_cells,1);

burst_t1 = uint32(round(burst_t1 ./ dt));
burst_t2 = uint32(round(burst_t2 ./ dt));
%keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% spontaneous values 
% hard-limit to above zero
limit = find(spon(STN) < 0);
spon(STN(limit)) = 0;
limit = find(spon(GPe) < 0);
spon(GPe(limit)) = 0;
limit = find(spon(GPi) < 0);
spon(GPi(limit)) = 0;

% scale STN spontaneous rate by dopamine impact
spon(STN) = spon(STN) .* (1 + stnda(3) * dop2);

% find striatal current values for down-state??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% shunting inhibition current
i_gate = zeros(n_neurons,1);
for loop = 1:n_neurons
    i_gate(loop) = find_Vm_cur(R(loop),spon(loop),shunt_to);
end

%%%%%%%% pre-multiply currents by exponential terms for efficiency
alphaCA = alphaCA .* as;      
spon = spon.*as;                
i_gate = i_gate .*as;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% spike input to model 
switches = switches ./ dt;
[in_n,in_t] = BG_GHS_inputs(dt,time_steps,n_inputs,switches,hz,n_neurons,neurons_per_channel,n_channels,ref_period,input_type);
n_in = length(in_n);

%%% current input to model
n_pulses = (time_seconds - t_offset) * pulse_Hz;
t_on = linspace(t_offset,time_seconds,n_pulses+1);   % on time
t_off = linspace(t_offset+t_width,time_seconds+t_width,n_pulses+1);    % off time
t_on = uint32(round(t_on ./ dt));    % convert to time-steps 
t_off = uint32(round(t_off ./ dt));
p_cells = uint32(zeros(n_neurons,1));
p_cells(pulse_cells) = 1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% network connectiviry 

% reset random seeds to ensure that, if run in batch mode, same network is constructed each
% iteration, given the same seed and network parameter combination
rand('state',seed1); randn('state',seed1);

scalar_constants = [p_connect, dt, dop1, dop2, scale, AMPA_PSP_size, NMDA_PSP_size, GABAa_PSP_size, n_channels,neurons_per_channel,neurons_per_nucleus,n_neurons,n_sources,max_scale];
[link_n,link_wAMPA,link_wNMDA,link_wGABAa,link_t,n_link,delays,bounds,n_paths,max_prox,max_soma,M_hat] = BG_net_shunt_AMPA_NMDA(scalar_constants,a_AMPA,a_NMDA,a_GABAa,as,weights,delays,proportions,mean_tau_AMPA,mean_tau_NMDA,mean_tau_GABAa,mean_tau_m,mean_R,stnda,gpeda);


%%%% queue parameters
ref_period_N = ref_period / dt;   % convert to time-steps
h = max(tau_m./dt)*3;                           %Horizon (scaled by time constants)
MaxEx = length(n_link);                         %Max external events per time step
MaxSelf = n_neurons;                            %Max self events per time step
MaxOut = (time_steps/ref_period_N * neurons_per_nucleus * n_nuclei)/2;        %Max number of network firing events

%%%% scalar parameters
dps = [sigma_bg,ref,dt,step_size,PSP_sigma,M_hat];                                   %Double parameters (25/02/04)
ips = uint32([h,MaxEx,MaxSelf,MaxOut,n_neurons,n_sources,n_in,time_steps,ref_period_N,trace_n-1]);   %Integer parameters (trace shifted to zero-base)

%%% initial values for state variables - must be arrays of size n_neurons
if load_state
   load('state');
   init_mem_pot = mem_pot;
   init_i_gabaa = i_gabaa;
   init_i_prox = i_prox;
   init_i_soma = i_soma;
   init_i_ampa = i_ampa;
   init_i_nmda = i_nmda;
else
   % initialise membrane potential to random voltage values (small deviations from resting potential
   % at zero are best to avoid triggering threshold-sensitive pseudo-currents e.g. bursting)
   % all others are best initialised to zero, otherwise inspired guessing at the appropriate current
   % (in Amperes) values is required
   init_mem_pot = randn(n_neurons,1) * 0.002;
   init_i_gabaa = zeros(n_neurons,1);
   init_i_prox = zeros(n_neurons,1);
   init_i_soma = zeros(n_neurons,1);
   init_i_ampa = zeros(n_neurons,1);
   init_i_nmda = zeros(n_neurons,1);
end

%%%% run model (this defines the neuron!)
[out_n,out_t,n_out,n_events,trace_vals,mem_pot,i_gabaa,i_soma,i_prox,i_ampa,i_nmda] = GHS_LIF_solver_shunt_AMPA_NMDA(in_n,in_t,link_n,link_wAMPA,link_wNMDA,link_wGABAa,link_t,...
                                                max_prox,max_soma,i_gate,delays,bounds,n_paths,spon,t_on,t_off,p_cells,...
                                                burst_t1,burst_t2,alphaCA,thetaCA,c_cells,...
                                                theta,mlimit,tau_AMPA,tau_NMDA,tau_GABAa,tau_m,R, weights, dps,ips,...
init_mem_pot, init_i_gabaa, init_i_prox, init_i_soma, init_i_ampa, init_i_nmda);

% get rid of useless space
out_t(out_t==0) = [];
out_n(out_n==0) = [];

% save all data
% file is 'results', many pars dumped - key vars are out_n, out_t which are
% the neuron index and time stamp of each spike event.
save(r_fname,'input_array','trace_vals','trace_n',...
        'R','theta','tau_m','tau_AMPA','tau_NMDA','tau_GABAa','spon','mlimit',...
        'dt','time_seconds','p_connect',...
        'n_channels','n_neurons','n_sources','n_inputs','neurons_per_nucleus','neurons_per_channel', 'n_nuclei',...
        'weights','delays','proportions','link_n','link_wAMPA','link_wNMDA','link_wGABAa','link_t',...
        't_on', 't_off', 'step_size', 'switches',...
        'out_n','out_t','GPi','STN','GPe','EXT','SD1','SD2','in_n','in_t');

if save_state
   % create a file to store state variables at end of simulation
   save('state','mem_pot','i_gabaa','i_soma','i_prox','i_ampa','i_nmda');
end
