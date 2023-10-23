x = table2array(readtable('pulse_freq_raw.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('sinc_squared', x(:,1), x(:,2), [1.5 0.01 -0.5 0.4]');

plot((x(:,1)+89E3)*1E-3, x(:,2), 'ok', (x(:,1)+89E3)*1E-3, fitdata, '-r')

l_ = legend("Data", "Sinc Fit");
t_ = title("Intervention Pulse Frequency");
x_ = xlabel("Pulse Frequency (kHz)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
