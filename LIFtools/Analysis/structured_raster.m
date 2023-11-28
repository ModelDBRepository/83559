function structured_raster(n_in, t_in, out_n, t_out, dt, neurons_per_nucleus,SD1, SD2, STN, GPe, GPi, duration)

t_out_t = double(t_out) .* dt; % timestamps for network events
t_in_t = double(t_in) .* dt;    % timestamps for input streams

figure
subplot(3,1,1)
plot(t_in_t, n_in, '.');
axis([0 duration min(n_in) max(n_in)]);
title('inputs');

ts = [];
ns = [];
for i=1:neurons_per_nucleus
    ts = [ts t_out_t(out_n == SD1(i))];
    ns = [ns out_n(out_n == SD1(i))];
end
subplot(3,1,2)
plot(ts, ns, '.')
axis([0 duration min(SD1) max(SD1)]);
title('Striatum D1');

ts = [];
ns = [];
for i=1:neurons_per_nucleus
    ts = [ts t_out_t(out_n == SD2(i))];
    ns = [ns out_n(out_n == SD2(i))];
end
subplot(3,1,3)
plot(ts, ns, '.')
axis([0 duration min(SD2) max(SD2)]);
title('Striatum D2');

ts = [];
ns = [];
for i=1:neurons_per_nucleus
    ts = [ts t_out_t(out_n == STN(i))];
    ns = [ns out_n(out_n == STN(i))];
end
figure
subplot(3,1,1)
plot(ts, ns, '.')
axis([0 duration min(STN) max(STN)]);
title('STN');

ts = [];
ns = [];
for i=1:neurons_per_nucleus
    ts = [ts t_out_t(out_n == GPe(i))];
    ns = [ns out_n(out_n == GPe(i))];
end
subplot(3,1,2)
plot(ts, ns, '.')
axis([0 duration min(GPe) max(GPe)]);
title('GPe');

ts = [];
ns = [];
for i=1:neurons_per_nucleus
    ts = [ts t_out_t(out_n == GPi(i))];
    ns = [ns out_n(out_n == GPi(i))];
end
subplot(3,1,3)
plot(ts, ns, '.')
axis([0 duration min(GPi) max(GPi)]);
title('GPi');

