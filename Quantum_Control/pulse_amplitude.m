clf;
clear;

x = table2array(readtable('pulse_amp_raw.csv', 'NumHeaderLines',1));
[p,fitdata] = lsqfit113('sin_squared_rise', x(:,1), x(:,2), [1.6 10 0.82 4.21 0.12 0 0.01]', [1 0 1 1 1 0 0]');
x_plot = linspace(0, 10, 1000);
y_plot = p(1) .* exp(-x_plot./p(2)) .* (sin(p(3).*x_plot + p(4))).^2 + p(5);

plot(x(:,1), x(:,2), 'ok');
hold on;
%plot(x_plot, p(1).*exp(-x_plot./p(2)));

l_ = legend("Data", "Sine Squared Exponential Rise Fit");
t_ = title("Intervention Pulse Amplitude");
x_ = xlabel("Pulse Amplitude (V)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');

hold off;