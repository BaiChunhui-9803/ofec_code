filepath = "//172.24.24.151/e/DiaoYiya/papre_con_data/review_data/running_time2/dim-time_v3.txt";
filename = "dim-time_v3";
data =readmatrix(filepath, "NumHeaderLines",1);
x= data(:,1);
y= data(:,2);

f = figure('visible','on');
plot(x,y);

fondsize = 15;

set(gca,'FontSize',fondsize);
%title('Dimensions - Time(s)')
%xlabel('Number of Samples','FontSize',fondsize)
%xlabel('Number of Peaks','FontSize',fondsize)
xlabel('Number of Dimensions','FontSize',fondsize)
ylabel('Time (s)', 'FontSize',fondsize);
setExportFigureType(filename,'origin',0.15);