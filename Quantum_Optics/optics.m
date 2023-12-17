clear;
clf;

x = table2array(readtable('data.csv', 'NumHeaderLines',1));

polarplot(x(:,2).*pi./180,x(:,1).*1E4*10, '*', 'LineStyle','None');

title("Downconverted Counts Per Second")

ax = gca;
ax.FontSize = 20;