x = table2array(readtable('pulse_amp_raw.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('sin_squared_decay', x(:,1), x(:,2), [2 -1 2*pi/8 3*pi/2 0]');

plot(x(:,1), x(:,2), 'ok', x(:,1), fitdata, '-r')

l_ = legend("Data", "Sine Squared Fit");
t_ = title("Intervention Pulse Duration");
x_ = xlabel("Pulse Amplitude (V)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
