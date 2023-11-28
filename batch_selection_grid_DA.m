function [sel_results_list,a_fname,txt_fname] = batch_selection_grid(prefix,pars_file,r_seed,pathroot,exp_name)

% handles a batch of selection experiments on a single model
%
%
% Mark Humphries 27/1/2006

thresh = 5;

load input_grid
[No_sims c] = size(input_array);


batch_gpis = struct('means', {}, 'stds', {}, 'units', {}, 'switch_sum', {});
batch_sum = [];
sel_results_list = {};

for i = 1:No_sims
    fprintf(1, 'simulation %d\n', i);
    
    % changes this
    rfname = BATCH_BG_heterogenous_AMPA_NMDA(i,pars_file,r_seed,pathroot,exp_name);
    
    
    
    % keep this
    fprintf(1, 'post processing %d\n', i);
    [gpi_means, gpi_stds, GPi_ch, summary] = mean_outputs(i,rfname,thresh);
    batch_gpis(i).means = gpi_means;
    batch_gpis(i).stds = gpi_stds;
    batch_gpis(i).units = GPi_ch;
    batch_gpis(i).switch_sum = summary;
    
    batch_sum = [batch_sum; summary];
    
    sel_results_list{i} = rfname;
end

% save all classifications to text file...
time_now = clock;
unique_name = datestr(time_now,30);  

txt_fname = [prefix '_' unique_name '_sel_sum.txt'];
save ([pathroot txt_fname],'batch_sum','-ascii')


% save analysis to path
a_fname = [prefix '_' unique_name '_analysis.mat'];
save([pathroot a_fname],'batch_gpis','batch_sum');
