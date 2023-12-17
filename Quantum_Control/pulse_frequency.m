x = table2array(readtable('pulse_freq_raw.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('abs_sinc', x(:,1), x(:,2), [1.8 0.0015 -0.3 0.1]');

xplot = linspace(-1600, 2200, 1000);
x_ = sym(xplot);
yplot = double(p(1) .* abs(sinc(p(2) * x_ + p(3))) + p(4));

hold on;
plot((x(:,1)+89E3)*1E-3, x(:,2), 'ok', (xplot).*1E-3 + 89, yplot, '-r')

hold off;


l_ = legend("Data", "| Sinc | Fit");
t_ = title("Intervention Pulse Frequency");
x_ = xlabel("Pulse Frequency (kHz)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
