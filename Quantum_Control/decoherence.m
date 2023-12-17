x = table2array(readtable('decoherence.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('exp_decay', x(:,1), x(:,2), [1 50 5 0]');

plot(x(:,1), x(:,2), 'ok', x(:,1), fitdata, '-r')

l_ = legend("Data", "Exponential Fit");
t_ = title("Decoherence (T2)");
x_ = xlabel("Tau (ms)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
