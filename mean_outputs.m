function [gpi_means, gpi_stds, GPi_ch, summary] = mean_outputs(index,fname,thresh)

% MEAN_OUTPUTS compute mean GPi outputs and categorises for selection....
% now incorporates same categorisation as FIND_SWITCHED_BATCH

load(fname)

bubble_scale = 4;
t_out_t = double(out_t) .* dt;

% get mean firing rates for all channels within each time segement
No_t_segs = length(switches);
GPi_mean_FRs = zeros(neurons_per_nucleus, No_t_segs);

for loop = 1:neurons_per_nucleus
    gpi_times = t_out_t(out_n == GPi(loop));
    t_start = 0;
    for j = 1:No_t_segs
        t_end = switches(j) .* dt;
        duration_seg = t_end - t_start;
        N_spikes_in_seg = sum(gpi_times > t_start & gpi_times < t_end);
        GPi_mean_FRs(loop, j) = N_spikes_in_seg ./ duration_seg;
        t_start = t_end;
    end
end

% break up into channels
GPi_ch = cell(1,n_channels);
for i = 1:n_channels
    start_i = neurons_per_channel .* (i-1) +1;
    end_i = start_i + neurons_per_channel - 1;
    GPi_ch{i} = GPi_mean_FRs(start_i:end_i, :);
end

% compute  stats (3 channels)
gpi_means = zeros(n_channels, No_t_segs);
gpi_stds = zeros(n_channels, No_t_segs);
for i = 1:n_channels
    gpi = GPi_ch{i};
    gpi_means(i,:) = mean(gpi);
    gpi_stds(i,:) = std(gpi);
end

% find summary - complete analysis
    % encode outcome using a number which, in binary, is given by the truth values in a 
    % truth table with columns c1_t1, c1_t2, c2_t2
    % which are channel 1 - time segment 1 and 2, and channel 2 time
    % segment 2, selected respectively. Thus:
    % 0 is no selection
    % 1 is single channel slection (ch 2)
    % 2, 3 are not valid
    % 4 is interference
    % 5 is switching
    % 6 is single channel selection (ch 1)
    % 7 is dual channel selection

    gpi_eqms1 = gpi_means(1:2, No_t_segs - 1);  % first 2 channels over penultimate segment
    gpi_eqms2 = gpi_means(1:2, No_t_segs);      % first 2 channels over last segment
    if gpi_eqms1(1) < thresh
        c1_t1 = 1;
    else
        c1_t1 = 0;
    end
    sel = c1_t1;
    if gpi_eqms2(1) < thresh
        c1_t2 = 1;
    else
        c1_t2 = 0;
    end
    sel = bitshift(sel, 1);
    sel = bitor(sel, c1_t2);
    if gpi_eqms2(2) < thresh
        c2_t2 = 1;
    else
        c2_t2 = 0;
    end
    sel = bitshift(sel, 1);
    sel = bitor(sel, c2_t2);
    if (sel == 2) | (sel == 3)
        ok = 0;
    else
        ok = 1;
    end
    
   
    % additional boolean coding of special cases for bubble plot
    if (sel == 1) | (sel == 6)
        single = bubble_scale;
    else
        single = 0;
    end
    if sel == 5
        switched = bubble_scale;
    else
        switched = 0;
    end
    if sel == 0
        none = bubble_scale;
    else 
        none = 0;
    end
    if sel == 7
        dual = bubble_scale;
    else
        dual = 0;
    end
    if sel == 4
        interf = bubble_scale;
    else
        interf = 0;
    end
    
    in_c1 = input_array(index, 1);
    in_c2 = input_array(index, 2);
    summary = [in_c1 in_c2 sel ok none single switched interf dual];
end



