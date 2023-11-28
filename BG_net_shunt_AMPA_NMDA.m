function [link_n,link_wAMPA,link_wNMDA,link_wGABAa,link_t,n_link,delays,bounds,n_paths,max_prox,max_soma,M_hat] = BG_net_shunt_AMPA_NMDA(scalar_constants,ae_AMPA,ae_NMDA,ai_GABAa,as,weights,delay_times,proportions,tau_AMPA,tau_NMDA,tau_GABAa,tau_m,R,stnda,gpeda)
% BG_NET_SHUNT MetaSim function - for use with BG_HETEROGENOUS
%  BG_NET_SHUNT - creates structure of model

% constants
p_connect = scalar_constants(1);
dt = scalar_constants(2);
dop1 = scalar_constants(3);
dop2 = scalar_constants(4);
scale = scalar_constants(5);
AMPA_PSP_size = scalar_constants(6);
NMDA_PSP_size = scalar_constants(7);
GABAa_PSP_size = scalar_constants(8);
n_channels = scalar_constants(9);
neurons_per_channel = scalar_constants(10);
neurons_per_nucleus = scalar_constants(11);
n_neurons = scalar_constants(12);
n_sources = scalar_constants(13);
max_scale = scalar_constants(14);

%specify nuclei indices
SD1 = 1:neurons_per_nucleus;                            %Striatal D1 neurons
SD2 = neurons_per_nucleus+1:2*neurons_per_nucleus;      %Striatal D2 neurons
STN = 2*neurons_per_nucleus+1:3*neurons_per_nucleus;    %Subthalamic neurons
GPe = 3*neurons_per_nucleus+1:4*neurons_per_nucleus;    %Globus Pallidus internus
GPi = 4*neurons_per_nucleus+1:5*neurons_per_nucleus;    %Globus Pallidus externus
EXT = 5*neurons_per_nucleus+1:6*neurons_per_nucleus;    %Extrinsic input


%Memory efficient cell array connectivity structure
link_n = cell(n_sources,1);                 %link (target) neurons
link_wAMPA = cell(n_sources,1);                 %link weights
link_wNMDA = cell(n_sources,1);                 %link weights
link_wGABAa = cell(n_sources,1);                 %link weights
link_t = cell(n_sources,1);                 %link types

n_link = repmat(uint32(0),n_sources,1);     %number of links
n_link(SD1) = neurons_per_channel;
n_link(SD2) = neurons_per_channel;
n_link(STN) = 2*neurons_per_nucleus;
n_link(GPe) = 2*neurons_per_channel + neurons_per_nucleus;
n_link(GPi) = neurons_per_nucleus;
n_link(EXT) = 3*neurons_per_channel;

%%%%%%%%% SET THESE TO ZERO TO LESION STRUCTURE %%%%%%%%%%%%%%%%%%%
n_paths = repmat(uint32(0),n_sources,1);
n_paths(SD1) = 1; %GPi
n_paths(SD2) = 1; %GPe
n_paths(STN) = 2; %GPe and GPi
n_paths(GPe) = 3; %STN, GPi, and itself
n_paths(GPi) = 1;
n_paths(EXT) = 3; %SD1, SD2 and STN

%introduce delays and bounds as more efficient alternative to paths
delays = repmat(uint32(0),n_sources,double(max(n_paths)));
bounds = repmat(uint32(0),n_sources,double(max(n_paths))+1); %Slightly bigger for n_link add-on

delays(SD1,1) = delay_times(1) / dt;
delays(SD2,1) = delay_times(2) / dt; 
delays(STN,1) = delay_times(3) / dt;  
delays(STN,2) = delay_times(4) / dt; 
delays(GPe,1) = delay_times(5) / dt;  
delays(GPe,2) = delay_times(6) / dt; 
delays(GPe,3) = delay_times(7) / dt; 
delays(GPi,1) = delay_times(8) / dt; 
delays(EXT,1) = delay_times(9) / dt; 
delays(EXT,2) = delay_times(10) / dt; 
delays(EXT,3) = delay_times(11) / dt; 

bounds(:,1) = 0; %all start from zero in mex - all index first entry of a given pathway
bounds(SD1,2) = neurons_per_channel;
bounds(SD2,2) = neurons_per_channel;
bounds(STN,2) = neurons_per_nucleus;
bounds(STN,3) = 2*neurons_per_nucleus;
bounds(GPe,2) = neurons_per_channel;
bounds(GPe,3) = 2*neurons_per_channel;
bounds(GPe,4) = 2*neurons_per_channel+neurons_per_nucleus;
bounds(GPi,2) = neurons_per_nucleus;
bounds(EXT,2) = neurons_per_channel;
bounds(EXT,3) = 2*neurons_per_channel;
bounds(EXT,4) = 3*neurons_per_channel;

%%%% determine weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate weights -> [relative strength * PSP-generating current (given mean R, Tm, Ts) * scalar (for no. of synapses)]
SD1_w_GABAa = weights(1) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_GABAa,tau_m(5),R(5),GABAa_PSP_size,'step') *scale; 
SD2_w_GABAa = weights(2) * ones(1,neurons_per_nucleus) * PSPtoPSC(tau_GABAa,tau_m(4),R(4),GABAa_PSP_size,'step') *scale;
STN_GPe_w_AMPA = weights(3) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_AMPA,tau_m(4),R(4),AMPA_PSP_size,'step') * scale; 
STN_GPe_w_NMDA = weights(3) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_NMDA,tau_m(4),R(4),NMDA_PSP_size,'step') * scale; 
STN_GPi_w_AMPA = weights(4) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_AMPA,tau_m(5),R(5),AMPA_PSP_size,'step') * scale; 
STN_GPi_w_NMDA = weights(4) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_NMDA,tau_m(5),R(5),NMDA_PSP_size,'step') * scale; 

% discrete STN output
%STN_GPew = weights(3) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_e1,tau_m(4),R(4),PSP_size,'step')) * scale; 
%STN_GPiw = weights(4) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_e1,tau_m(5),R(5),PSP_size,'step')) * scale; 

GPe_STN_w_GABAa = weights(5) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_GABAa,tau_m(3),R(3),GABAa_PSP_size,'step') *scale; 
GPe_GPi_w_GABAa = weights(6) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_GABAa,tau_m(5),R(5),GABAa_PSP_size,'step')*scale;
GPe_GPe_w_GABAa = weights(7) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_GABAa,tau_m(4),R(4),GABAa_PSP_size,'step')*scale;
GPi_GPi_w_GABAa = weights(8) * ones(1,neurons_per_nucleus)* PSPtoPSC(tau_GABAa,tau_m(5),R(5),GABAa_PSP_size,'step')*scale;
EXT_w_AMPA = weights(9) * ones(1,3*neurons_per_nucleus)*scale;
EXT_w_AMPA(1:neurons_per_nucleus) = EXT_w_AMPA(1:neurons_per_nucleus)*PSPtoPSC(tau_AMPA,tau_m(1),R(1),AMPA_PSP_size,'step'); 
EXT_w_AMPA(1+neurons_per_nucleus:2*neurons_per_nucleus) = EXT_w_AMPA(1+neurons_per_nucleus:2*neurons_per_nucleus)*PSPtoPSC(tau_AMPA,tau_m(2),R(2),AMPA_PSP_size,'step'); 
EXT_w_AMPA(1+2*neurons_per_nucleus:3*neurons_per_nucleus) = EXT_w_AMPA(1+2*neurons_per_nucleus:3*neurons_per_nucleus)*PSPtoPSC(tau_AMPA,tau_m(3),R(3),AMPA_PSP_size,'step'); 


EXT_w_NMDA = weights(9) * ones(1,3*neurons_per_nucleus)*scale;
EXT_w_NMDA(1:neurons_per_nucleus) = EXT_w_NMDA(1:neurons_per_nucleus)*PSPtoPSC(tau_NMDA,tau_m(1),R(1),NMDA_PSP_size,'step'); 
EXT_w_NMDA(1+neurons_per_nucleus:2*neurons_per_nucleus) = EXT_w_NMDA(1+neurons_per_nucleus:2*neurons_per_nucleus)*PSPtoPSC(tau_NMDA,tau_m(2),R(2),NMDA_PSP_size,'step'); 
EXT_w_NMDA(1+2*neurons_per_nucleus:3*neurons_per_nucleus) = EXT_w_NMDA(1+2*neurons_per_nucleus:3*neurons_per_nucleus)*PSPtoPSC(tau_NMDA,tau_m(3),R(3),NMDA_PSP_size,'step'); 

%dopamine modulation via weight adjustment 
% 1. Striatum
EXT_w_AMPA(1:neurons_per_nucleus) = EXT_w_AMPA(1:neurons_per_nucleus)*(1+dop1); %D1
EXT_w_AMPA(1+neurons_per_nucleus:2*neurons_per_nucleus) = EXT_w_AMPA(1+neurons_per_nucleus:2*neurons_per_nucleus)*(1-dop2); %D2

EXT_w_NMDA(1:neurons_per_nucleus) = EXT_w_NMDA(1:neurons_per_nucleus)*(1+dop1); %D1
EXT_w_NMDA(1+neurons_per_nucleus:2*neurons_per_nucleus) = EXT_w_NMDA(1+neurons_per_nucleus:2*neurons_per_nucleus)*(1-dop2); %D2

% 2. STN (includes cortico-subthalamic weight ratio with cortico-striatal
% strength)
EXT_w_AMPA(1+2*neurons_per_nucleus:3*neurons_per_nucleus) = weights(10) * EXT_w_AMPA(1+2*neurons_per_nucleus:3*neurons_per_nucleus)*(1-stnda(1)*dop2);
EXT_w_NMDA(1+2*neurons_per_nucleus:3*neurons_per_nucleus) = weights(10) * EXT_w_NMDA(1+2*neurons_per_nucleus:3*neurons_per_nucleus)*(1-stnda(1)*dop2);
GPe_STN_w_GABAa = GPe_STN_w_GABAa .* (1-stnda(2)*dop2);

% 3. GPe
SD2_w_GABAa = SD2_w_GABAa .* (1 - gpeda(1)*dop2);
STN_GPe_w_AMPA = STN_GPe_w_AMPA .* (1 - gpeda(2)*dop2); 
STN_GPe_w_NMDA = STN_GPe_w_NMDA .* (1 - gpeda(2)*dop2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connections: 
%1. connect, on average, to proportion p_connect neurons in the corresponding channel
max_prox = zeros(n_neurons,1);
max_soma = zeros(n_neurons,1);

for c = 1:n_channels
    % input indexes
    in = (c-1)*neurons_per_channel+(1:neurons_per_channel);
   
    for n = 1:neurons_per_channel    
        
        out = (c-1)*neurons_per_channel + n;    % current out neuron

        % channel-wise input
        ext_probs = rand(3,neurons_per_channel);
        blnExt = ext_probs <= p_connect;
        EXT_STN_in = uint32((STN(in)-1) .* blnExt(1,:));
        EXT_SD1_in = uint32((SD1(in)-1) .* blnExt(2,:));
        EXT_SD2_in = uint32((SD2(in)-1) .* blnExt(3,:));
        EXT_out = [EXT_SD1_in,EXT_SD2_in,EXT_STN_in];

        % channel-wise output
		probs = rand(neurons_per_channel,4);        % 4 channel-wise connections to sparsify
		blnIn = probs <= p_connect;                 % logical array of selected unit
        D1_GPi_in = uint32((GPi(in)-1).*blnIn(:,1)'); %shifted back to zero-base for mex file
		D2_GPe_in = uint32((GPe(in)-1).*blnIn(:,2)');
        GPe_STN_in = uint32((STN(in)-1).*blnIn(:,3)');
        GPe_GPi_in = uint32((GPi(in)-1).*blnIn(:,4)');
        GPe_out = [GPe_STN_in GPe_GPi_in];

        % diffuse connections
        diffuse_probs = rand(neurons_per_nucleus,4);                % 3 diffuse connections to sparsify
        blnDiffuseIn = diffuse_probs <= (p_connect / n_channels);      % scale probability by number of channels
        % diffuse STN output
        STN_out = [uint32((GPe-1).*blnDiffuseIn(:,1)') uint32((GPi-1).*blnDiffuseIn(:,2)')];         
        
        % diffuse GPe-GPe
        blnDiffuseIn(out,3) = 0;     % do not connect to self            
        GPe_out = [GPe_out uint32((GPe-1).*blnDiffuseIn(:,3)')];
        
        % diffuse GPi-GPi
        blnDiffuseIn(out,4) = 0;     % do not connect to self            
        GPi_out = uint32((GPi-1).*blnDiffuseIn(:,4)');
        
        % discrete STN output
%         stn_probs = rand(neurons_per_channel,2);                % 2 discrete connections
%         blnSTN_in = stn_probs <= p_connect;
%         STN_out = [uint32((GPe(in)-1).*blnSTN_in(:,1)') uint32((GPi(in)-1).*blnSTN_in(:,2)')];         
        
        % assign connections        
        syn = rand(neurons_per_channel,4) .* blnIn;     % all of blnIn as it defines all inhibitory connections
        syn_diffuse = rand(neurons_per_nucleus,2) .* blnDiffuseIn(:,3:4);
        
        link_n{SD1(out)} = D1_GPi_in;
        link_wGABAa{SD1(out)} = SD1_w_GABAa(in) .* ai_GABAa(GPi(in))' .* blnIn(:,1)';
        link_wAMPA{SD1(out)} = 0;
        link_wNMDA{SD1(out)} = 0;
        type1 = zeros(neurons_per_channel,1);
        type1(syn(:,1) > 0 &  syn(:,1) <= proportions(1,1)) = 1;                                         % distal
        type1(syn(:,1) > proportions(1,1) &  syn(:,1) <= (proportions(1,1) + proportions(1,2))) = 2;     % proximal
        type1(syn(:,1) > (proportions(1,1) + proportions(1,2))) = 3;                                     % somatic    
        link_t{SD1(out)} = uint32(type1); 
        
        link_n{SD2(out)} = D2_GPe_in;
        link_wGABAa{SD2(out)} = SD2_w_GABAa(in) .*ai_GABAa(GPe(in))' .* blnIn(:,2)';
        link_wAMPA{SD2(out)} = 0;
        link_wNMDA{SD2(out)} = 0;
        type2 = zeros(neurons_per_channel,1);
        type2(syn(:,2) > 0 &  syn(:,2) <= proportions(2,1)) = 1;                                         % distal
        type2(syn(:,2) > proportions(2,1) &  syn(:,2) <= (proportions(2,1) + proportions(2,2))) = 2;     % proximal
        type2(syn(:,2) > (proportions(2,1) + proportions(2,2))) = 3;                                     % somatic    
        link_t{SD2(out)} = uint32(type2); 

        
        link_n{GPe(out)} = GPe_out;
        link_wGABAa{GPe(out)} = [GPe_STN_w_GABAa(in) .* ai_GABAa(STN(in))' .* blnIn(:,3)',...
                GPe_GPi_w_GABAa(in) .* ai_GABAa(GPi(in))' .* blnIn(:,4)',...
                GPe_GPe_w_GABAa .* ai_GABAa(GPe)' .* blnDiffuseIn(:,3)'];
        link_wAMPA{GPe(out)} = 0;
        link_wNMDA{GPe(out)} = 0;
        type3 = zeros(neurons_per_channel,1); % TO STN
        type3(syn(:,3) > 0 &  syn(:,3) <= proportions(3,1)) = 1;                                         % distal
        type3(syn(:,3) > proportions(3,1) &  syn(:,3) <= (proportions(3,1) + proportions(3,2))) = 2;     % proximal
        type3(syn(:,3) > (proportions(3,1) + proportions(3,2))) = 3;                                     % somatic    
        type4 = zeros(neurons_per_channel,1);   % TO GPI
        type4(syn(:,4) > 0 &  syn(:,4) <= proportions(4,1)) = 1;                                         % distal
        type4(syn(:,4) > proportions(4,1) &  syn(:,4) <= (proportions(4,1) + proportions(4,2))) = 2;     % proximal
        type4(syn(:,4) > (proportions(4,1) + proportions(4,2))) = 3;                                     % somatic    
        type5 = zeros(neurons_per_nucleus,1);   % TO ITSELF
        type5(syn_diffuse(:,1) > 0 &  syn_diffuse(:,1) <= proportions(5,1)) = 1;                                         % distal
        type5(syn_diffuse(:,1) > proportions(5,1) &  syn_diffuse(:,1) <= (proportions(5,1) + proportions(5,2))) = 2;     % proximal
        type5(syn_diffuse(:,1) > (proportions(5,1) + proportions(5,2))) = 3;                                     % somatic    
        
        link_t{GPe(out)} = [uint32(type3)' uint32(type4)' uint32(type5)']; 
        
        link_n{GPi(out)} = GPi_out;
        link_wGABAa{GPi(out)} = GPi_GPi_w_GABAa .* ai_GABAa(GPi)' .* blnDiffuseIn(:,4)';
        link_wAMPA{GPi(out)} = 0;
        link_wNMDA{GPi(out)} = 0;        
        type6 = zeros(neurons_per_nucleus,1);
        type6(syn_diffuse(:,2) > 0 &  syn_diffuse(:,2) <= proportions(6,1)) = 1;                                         % distal
        type6(syn_diffuse(:,2) > proportions(6,1) &  syn_diffuse(:,2) <= (proportions(6,1) + proportions(6,2))) = 2;     % proximal
        type6(syn_diffuse(:,2) > (proportions(6,1) + proportions(6,2))) = 3;                                     % somatic    
        link_t{GPi(out)} = uint32(type6)';
        
        % keyboard
        % diffuse STN
        link_n{STN(out)} = STN_out;
        link_wAMPA{STN(out)} = [STN_GPe_w_AMPA .*ae_AMPA(GPe)'.*blnDiffuseIn(:,1)',...
                                STN_GPi_w_AMPA .*ae_AMPA(GPi)'.*blnDiffuseIn(:,2)']; 
        link_wNMDA{STN(out)} = [STN_GPe_w_NMDA .*ae_NMDA(GPe)'.*blnDiffuseIn(:,1)',...
                                STN_GPi_w_NMDA .*ae_NMDA(GPi)'.*blnDiffuseIn(:,2)']; 
        link_wGABAa{STN(out)} = 0;
        link_t{STN(out)} = zeros(neurons_per_nucleus*2,1);      % STN output excitatory
        
        % keyboard
        % discrete STN
%         link_n{STN(out)} = STN_out;
%         link_w{STN(out)} = [STN_GPew.*ae(GPe(in))'.*blnSTN_in(:,1)', STN_GPiw.*ae(GPi(in))'.*blnSTN_in(:,2)']; 
%         link_t{STN(out)} = zeros(neurons_per_nucleus*2,1);      % STN output excitatory
        
        link_n{EXT(out)} = EXT_out;
        link_wAMPA{EXT(out)} = [EXT_w_AMPA(in) .* ae_AMPA(SD1(in))' .* blnExt(2,:),...
                                EXT_w_AMPA(in+neurons_per_nucleus) .* ae_AMPA(SD2(in))' .* blnExt(3,:),...
                                EXT_w_AMPA(in+2*neurons_per_nucleus) .* ae_AMPA(STN(in))' .* blnExt(1,:)];
        link_wNMDA{EXT(out)} = [EXT_w_NMDA(in) .* ae_NMDA(SD1(in))' .* blnExt(2,:),...
                                EXT_w_NMDA(in+neurons_per_nucleus) .* ae_NMDA(SD2(in))' .* blnExt(3,:),...
                                EXT_w_NMDA(in+2*neurons_per_nucleus) .* ae_NMDA(STN(in))' .* blnExt(1,:)];
        link_wGABAa{EXT(out)} = 0;
        link_t{EXT(out)} = zeros(neurons_per_channel*3,1);      % no shunting for external drive 
        
        % collate synapse types - divide by weight to leave only on
        % numerator so that they have effect!
%           weights = [SD1_w 
%                      SD2_w 
%                      STN_GPew 
%                      STN_GPiw 
%                      GPe_STNw 
%                      GPe_GPiw 
%                      GPe_GPew 
%                      GPi_GPiw 
%                      EXT_w 
%                      STN_ext_ratio]; 
        if weights(5) ~= 0
            max_prox(STN(in)) = max_prox(STN(in)) + (type3 == 2) .* (GPe_STN_w_GABAa(in)' ./ ( abs(weights(5)) .* (1-stnda(2)*dop2) ) ) .* ai_GABAa(STN(in));    % maximum weight value possible
            max_soma(STN(in)) = max_soma(STN(in)) + (type3 == 3) .* (GPe_STN_w_GABAa(in)' ./ ( abs(weights(5)) .* (1-stnda(2)*dop2) ) ) .* ai_GABAa(STN(in));    % maximum weight value possible
        end
        
        if weights(2) ~= 0 
            max_prox(GPe(in)) = max_prox(GPe(in)) + (type2 == 2) .* (SD2_w_GABAa(in)' ./ (abs(weights(2)).*  (1-gpeda(1)*dop2) ) ) .* ai_GABAa(GPe(in));
            max_soma(GPe(in)) = max_soma(GPe(in)) + (type2 == 3) .* (SD2_w_GABAa(in)' ./ (abs(weights(2)) .*  (1-gpeda(1)*dop2) ) )  .* ai_GABAa(GPe(in));
        end
        
        if weights(1) ~= 0
            max_prox(GPi(in)) = max_prox(GPi(in)) + ((type1 == 2) .* (SD1_w_GABAa(in)' ./ abs(weights(1))) .* ai_GABAa(GPi(in)));
            max_soma(GPi(in)) = max_soma(GPi(in)) + ((type1 == 3) .* (SD1_w_GABAa(in)' ./ abs(weights(1))) .* ai_GABAa(GPi(in)));
        end
        if weights(6) ~=0
            max_prox(GPi(in)) = max_prox(GPi(in)) + ((type4 == 2) .* (GPe_GPi_w_GABAa(in)'./ abs(weights(6))) .* ai_GABAa(GPi(in)));
            max_soma(GPi(in)) = max_soma(GPi(in)) + ((type4 == 3) .* (GPe_GPi_w_GABAa(in)'./ abs(weights(6))) .* ai_GABAa(GPi(in)));        
        end
        
    %keyboard
        
        % collate synapse types
%         pad = zeros(neurons_per_channel*(n_channels-1),1); % use to pad all those for which type is specified on a per channel basis
%         max_prox(STN) = max_prox(STN) + [(type3 == 2); pad]   .* GPe_STNw' .* ai(STN);    % maximum weight value possible
%         max_prox(GPe) = max_prox(GPe) + ([(type2 == 2); pad] .* SD2_w' .* ai(GPe)) + ((type5 == 2) .* GPe_GPew' .* ai(GPe));    
%         max_prox(GPi) = max_prox(GPi) + ([(type1 == 2); pad] .* SD1_w' .* ai(GPi)) + ([(type4 == 2); pad] .* GPe_GPiw' .* ai(GPi)) + ((type6 == 2) .* GPi_GPiw' .* ai(GPi));
% 	
%         max_soma(STN) = max_soma(STN) + [(type3 == 3); pad] .* GPe_STNw' .* ai(STN);    % maximum weight value possible
%         max_soma(GPe) = max_soma(GPe) + ([(type2 == 3); pad] .* SD2_w' .* ai(GPe)) + ((type5 == 3) .* GPe_GPew' .* ai(GPe));    
%         max_soma(GPi) = max_soma(GPi) + ([(type1 == 3); pad] .* SD1_w' .* ai(GPi)) + ([(type4 == 3); pad] .* GPe_GPiw' .* ai(GPi)) + ((type6 == 3) .* GPi_GPiw' .* ai(GPi));

    end    
end

% add diffuse inhibitory output to max prox/soma calculations
GG_sum_type = zeros(3,neurons_per_nucleus);
SS_sum_type = zeros(3,neurons_per_nucleus);     % use S prefix fo SNr-SNr collaterals (named GPi-GPi above) 

GPe_idx = double(bounds(GPe(1),3))+1:double(bounds(GPe(1),4));
GPi_idx = double(bounds(GPi(1),1))+1:double(bounds(GPi(1),2));

for loop = 1:neurons_per_nucleus
    out_type = double(link_t{GPe(loop)}(GPe_idx));    % get output type for non-summed connection        
    GG_sum_type(1,:) = GG_sum_type(1,:) + (out_type == 1);    % distal
    GG_sum_type(2,:) = GG_sum_type(2,:) + (out_type == 2);    % proximal
    GG_sum_type(3,:) = GG_sum_type(3,:) + (out_type == 3);    % somatic
    
    out_type = double(link_t{GPi(loop)}(GPi_idx)); 
    SS_sum_type(1,:) = SS_sum_type(1,:) + (out_type == 1);    % distal
    SS_sum_type(2,:) = SS_sum_type(2,:) + (out_type == 2);    % proximal
    SS_sum_type(3,:) = SS_sum_type(3,:) + (out_type == 3);    % somatic
end
%           weights = [SD1_w 
%                      SD2_w 
%                      STN_GPew 
%                      STN_GPiw 
%                      GPe_STNw 
%                      GPe_GPiw 
%                      GPe_GPew 
%                      GPi_GPiw 
%                      EXT_w 
%                      STN_ext_ratio]; 

if weights(7) ~= 0
    max_prox(GPe) = max_prox(GPe) + GG_sum_type(2,:)' .* ai_GABAa(GPe) .* (GPe_GPe_w_GABAa' ./ abs(weights(7)));   % NOTE: this will have to change should weights ever be output specific
    max_soma(GPe) = max_soma(GPe) + GG_sum_type(3,:)' .* ai_GABAa(GPe) .* (GPe_GPe_w_GABAa' ./ abs(weights(7)));
end

if weights(8) ~= 0
    max_prox(GPi) = max_prox(GPi) + SS_sum_type(2,:)' .* ai_GABAa(GPi) .* (GPi_GPi_w_GABAa' ./ abs(weights(8)));   % NOTE: this will have to change should weights ever be output specific
    max_soma(GPi) = max_soma(GPi) + SS_sum_type(3,:)' .* ai_GABAa(GPi) .* (GPi_GPi_w_GABAa' ./ abs(weights(8)));
end



% scale max_prox and max_soma - to be the sum of simultaneous IPSPs or
% less..
max_prox = max_prox .* max_scale;
max_soma = max_soma .* max_scale;

% remove all zeros - necessary to prevent divide by zero in C code
% is safe because all corresponding link_t{} values are also zero and therefore inhibtory currents
% for these cells never change!
max_prox(max_prox==0) = 1;
max_soma(max_soma==0) = 1;  

% find the median value
M_hat = median([max_prox(max_prox<1)' max_soma(max_soma<1)']);

% keyboard
% hist([max_prox(max_prox<1)' max_soma(max_soma<1)'],25)



