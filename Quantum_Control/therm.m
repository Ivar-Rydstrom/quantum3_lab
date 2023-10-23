x = table2array(readtable('therm.csv', 'NumHeaderLines',1));
clf;
[p,fitdata] = lsqfit113('exponential_rise', x(:,1), x(:,2), [200 100 2 0]');

plot(x(:,1), x(:,2), 'ok', x(:,1), fitdata, '-r')

l_ = legend("Data", "Exponential Fit");
t_ = title("Thermalization (T1)");
x_ = xlabel("Hold Time (s)");
y_ = ylabel("Echo Amplitude (V)");

fontsize(l_,15,'points');
fontsize(t_,15,'points');
fontsize(x_,15,'points');
fontsize(y_,15,'points');
