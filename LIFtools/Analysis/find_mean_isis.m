function rs = find_mean_isis(cum_spikes)

f_length = 50;
filter = gausswin(f_length);
events = find(cum_spikes);
intervals = diff(events);
rates = 1 ./ intervals;
rates = [0 rates];
rates =  cum_spikes(events) .* rates;
rs = zeros(1, length(cum_spikes));
rs(events) = rates;
rs = conv(rs, filter);
