%%%% BASIC PARAMETERS SCRIPT - Mark Humphries 18/8/2005

%%% Condition#1 of Magill et al. (2001) Ctx, DA 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   PARAMETERS THAT MOST OFTEN CHANGE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% random seed
rand('state',seed1); randn('state',seed1);

%%%% simulation time
time_seconds = 10;          %total simulation time in seconds (usually 10)  

do_urethane = 1;

%%%% dopamine
dop1 = 0.3; % tonic dopamine level (in [0 1] range) (0.3)
           % same weight calculation as in systems level model
dop2 = 0.3;
%%%% inputs
% In the following, for
% selection xpt             - 'vivo' (basic tonic with poisson stats)
% tonic rate calibration    - 'vivo' 
                                % (for tonic rates a rule of thumb is
                                % 30-30-12 for GPe-GPi-STN resp)
% slow wave xpt             - 'slow' (anaesthetic-like correlated slow-wave, uses the fixed train+jitter model)
% organo xpt                - 'organo' (correlated tonic, uses the fixed train+jitter model)

input_type = 'slow';        
switches = [0.2, 0.5, 1] .* time_seconds;  % in seconds

% In the following, for
% selection xpt             - 'switch' (switches salience at switch points using
                                % conventional 'salience grid' paradigm  
                                % alternative is 'simultaneous' (switch
                                % both channels at same time)
% tonic rate calibration    - 'tonic' (with 2-4Hz) - should probably set 'save-state' flag too
% slow wave xpt             - 'tonic' (with 24Hz after Steriade)
% organo xpt                - 'tonic' (rate?)


input_method = 'tonic'; 
tonic_rate = 32;             % in Hz 
rate_scaling = 0.2;          % scale rates in input grid by this amount

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify basic network parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_nuclei = 5; %SD1, SD2, STN, GPe, GPi
n_channels = 3;
neurons_per_channel = 64;
neurons_per_nucleus = neurons_per_channel*n_channels;
n_neurons = n_nuclei*neurons_per_nucleus;

% specify connection proportion
p_connect = 0.25; % cf 'rho'

%specify nuclei indices
SD1 = 1:neurons_per_nucleus;                            %Striatal D1 neurons
SD2 = neurons_per_nucleus+1:2*neurons_per_nucleus;      %Striatal D2 neurons
STN = 2*neurons_per_nucleus+1:3*neurons_per_nucleus;    %Subthalamic neurons
GPe = 3*neurons_per_nucleus+1:4*neurons_per_nucleus;    %Globus Pallidus internus
GPi = 4*neurons_per_nucleus+1:5*neurons_per_nucleus;    %Globus Pallidus externus
EXT = 5*neurons_per_nucleus+1:6*neurons_per_nucleus;    %Extrinsic input

trace_n = GPe(75); % index of neuron to record in detail membrane potential, epsp, J,...
                   % (number has to be less than No neurosn per nucleus)

%specify time parameters  
dt = 0.0001;                             %time steps in seconds 

% other neuron parameters
sigma_bg = 0.0003;                         %noise std dev - in V in solution of membrane eqn.
ref = 0;                                %refractory membrane reset potential (in V) - (V_reset)
ref_period = 0.002;                     %Absolute refractory period in seconds 
theta = ones(n_neurons,1)* 0.03;        %thresholds given in V above resting potential
theta(STN) = 0.02;                      % STN has lower firing threshold

mlimit = ones(n_neurons,1) * -0.02;       % limiting threshold below which membrane potential can't go: 
                                          % mimics limiting effect of GABA reversal potential

%time constants - one per neuron
mean_tau_AMPA = 0.002;   % AMPA component - 2ms excitatory current % 0.002
mean_tau_NMDA = 0.100;   % NMDA component - 120ms excitatory current 
mean_tau_GABAa = 0.003;   % GABAa - 2ms inhibitory current

mean_tau_m = {0.025, ...    % SD1
              0.025, ...    % SD2    
              0.006, ...    % STN
              0.014, ...    % GPe
              0.008};       % GPi/SNr 0.008

% adding noise to the time membrane constants
% all nuclei have gaussian dist. except STN which is gamma (outside
% parameter control block)
std_tau_m = cell2mat(mean_tau_m) .* 0.1;      % 10% std dev

% variable resistances
mean_R = {42e6, ...         % SD1 
          42e6, ...         % SD2
          18e6, ...         % STN
          88e6, ...         % GPe
          112e6};           % GPi/SNr
  
std_R = cell2mat(mean_R) .* 0.1;      % 10% std dev

% weights and associated parameters
% the factors w_ij / tau_s are calculated to give currents that yield PSPs
% of measured size. These are then *relative* weighting factors multiplied
% by w_ij / tau_s


SD1_w = -4; %-4;
SD2_w = -4; %-4;

STN_GPiw = 1.0; % 1 

%%%%%%%%% STN GP loop %%%%%%

STN_GPew = 1.0; %1.5

GPe_STNw = -1.0;  % (1) set to -2 to get consistent bursting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GPe_GPiw = -1.0; %-1;
 
% omit collaterals here
GPe_GPew = -1; % -1 collaterals
GPi_GPiw = -1; % -1 collaterals

EXT_w = 1.0; % ctx input 1.0
STN_ext_ratio = 1.0; % cortico-subthalamic weight is EXT_w * STN_ext_ratio (1)

%%%%%%%%% do urethane manipulation %%%%%%%%%%%
scale = 1.0; %for historical reasons in code

if do_urethane
    glut_scale = 0.65; % turn all glut wgts down (0.65? - KG)
    gaba_scale = 1.5; % turn all gaba  wgts up 
else
    glut_scale = 1.0; 
    gaba_scale = 1.0; 
end

SD1_w = SD1_w .* gaba_scale;
SD2_w = SD2_w .* gaba_scale;
GPe_STNw = GPe_STNw .* gaba_scale;
GPe_GPiw = GPe_GPiw .* gaba_scale;
GPe_GPew = GPe_GPew .* gaba_scale;
GPi_GPiw = GPi_GPiw .* gaba_scale;

STN_GPiw = STN_GPiw .* glut_scale;
STN_GPew = STN_GPew .* glut_scale;
EXT_w = EXT_w .* glut_scale;



%scale = (1/neurons_per_channel * p_connect);
AMPA_PSP_size = 0.003;                   % max size in V (for all populations 
NMDA_PSP_size = 0.0001;                   % (0.0001) 
GABAa_PSP_size = 0.003;

%PSP_sigma = PSP_size .* 0.1;        % noise std is 10% of PSP
PSP_sigma = 0.000; %0

% delays in seconds - axonal delays
SD12GPi_d = 0.004;  
SD22GPe_d = 0.005;
STN2GPe_d = 0.002;
STN2GPi_d = 0.0015;
GPe2STN_d = 0.004; %0.004
GPe2GPi_d = 0.003;
GPe2GPe_d = 0.001;
GPi2GPi_d = 0.001;
EXT2SD1_d = 0.01;
EXT2SD2_d = 0.01;
EXT2STN_d = 0.0025;

% dopamine coefficients
stnda = [0.25 0.5 0];  % STN coefficients for proportion of dopamine that affects [AMPA GABAa Spon] 
                         % currents (all < 1) (don't have to sum to 1)
                         % maybe spon = 0?
                         % nominally [0.25 0.5 0.1];
gpeda = [0.5 0.5];         % GPe coefficients for proportion of dopamine that affects [SD2 STN] input - NOT specific currents as above 
                         %(don't have to sum to 1) nominally [0.5 0.5]
%stnda = [0 0 0];   

%%%%% shunting inhibition (current size computed in main code)
% synaptic distributions [distal, prox, soma] - must sum to 1! (only
% required for inhibitory synapses)
SD12GPi_p = [1 0 0];    % distal only 
SD22GPe_p = [0.33 0.34 0.33];    % distal only 
GPe2STN_p = [0.3 0.4 0.3]; % [0.3 0.4 0.3];
GPe2GPi_p = [0 0.5 0.5];
GPe2GPe_p = [0 0.5 0.5];
GPi2GPi_p = [0 0.5 0.5];
max_scale = 0.5;      % scale of maximum proximal and somatic inputs (beta in maths)
shunt_to = -0.02;   % target membrane potential (in volts) of shunting current  
                    % the value to which V is driven by I_shunt with
                    % maximal inhibition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      INTRINSIC CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% spontaneous currents  (limited by Heaviside in main code)

spon = zeros(n_neurons,1); %  constant spontaneous input current for each neuron
%spon(STN) = randn(length(STN),1)*(1e-11) + (9e-11); 

if do_urethane
    spon(STN) = 5e-10; 
else
    spon(STN) = 11e-10; % 10.0e-10;       % 11e-10 with GP collaterals @ 1 
end

spon(GPe) = 3.8e-10; %2.7e-10;    %4.6e-10 with GP collaterals @ 1        
%spon(GPe) = randn(length(GPe),1)*(2e-11) + (4e-10);
%spon(GPi) = randn(length(GPi),1)*(2.5e-11) + (2.5e-10);          
spon(GPi) = 3.9e-10; %2.8e-10;    %3.2e-10 with GP and GPi collaterals @ 1 
spon(SD1) = -2.5e-10; % -2.5e-10
spon(SD2) = -2.5e-10; % -2.5e-10


%% burst-current parameters
mean_t1 = 0.2;          % seconds 0.1
std_t1 = 0.01;
mean_t2 = 1.0;             % 0.5
std_t2 = 0.22;             %~10% variance
mean_thetaCA = -0.01;      % volts (below rest) -0.01
std_thetaCA = -0.001;
mean_alphaCA = 9e-10;    % 9e-10. amperes - work out from membrane potential equation (e.g. set Hz=80)
                        % may need adjustmentin a circuit to set
                        % intra-burst freq = 80Hz (from Plenz and Kita)
std_alphaCA =  9e-11;   % 9e-11
ca_cells = STN;      % array of cells with burst pseudo-current


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     INPUTS          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% current pulse
pulse_Hz = 0.0;       % frequency of pulses 
t_offset = 1;        % start of pulse train (in seconds)
t_width = 0.5;      % width of pulse (in seconds)
step_size = -1e-8;          % in Amperes
pulse_cells = [];        % array of cells receiving pulse trains


% state variable flags
load_state = 1;      % if = 1, loads initial state from state.mat
save_state = 0;      % if = 1, saves final state to state.mat
