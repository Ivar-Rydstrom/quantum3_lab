x = table2array(readtable('pulse_duration_raw.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('sin_squared_decay', x(:,1), x(:,2), [2 900 0.02 pi/2 0]');

f = 89;
T = 1/f;
t = x(:,1) .* T;
x_plot = linspace(10, 500, 1000);
y_plot = p(1) .* exp(-x_plot./p(2)) .* (sin(p(3).*x_plot + p(4))).^2 + p(5);
x_plot = x_plot.*T;
plot(t, x(:,2), 'ok', x_plot, y_plot, '-r')

l_ = legend("Data", "Decaying Sine Squared Fit");
t_ = title("Intervention Pulse Duration");
x_ = xlabel("Pulse Duration (ms)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
