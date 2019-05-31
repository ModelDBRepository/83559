function h = do_single_raster(times)

times = times';
t = [times times];
ysp = [0 1];
timeline = ysp(2)/2;

figure;
hold on
plot(t, ysp, 'k');
mint = min(times) - 1;
maxt = max(times) + 1;
axis([mint maxt  -1 2]);
plot([mint maxt], [timeline timeline], 'k');
hold off
