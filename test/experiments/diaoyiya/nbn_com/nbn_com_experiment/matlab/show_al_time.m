readdir = "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\alg_time\";
filename= "2202"
formatSpec="%f";
bridgeIdsName= readdir + filename + "_tsp_algTimeMatlab.txt";
sizeA= [1 2]; 
fileID = fopen(bridgeIdsName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);

time= fscanf(fileID,formatSpec,sizeM);


x = time(:,1)';
y1 = time(:,2)';    % 第一组y数据
y2 = time(:,3)'; % 第二组y数据
% 绘制图形
figure; % 创建新图形窗口

% 绘制图形
figure; % 创建新图形窗口

% 绘制第一组数据，设置线条格式
plot(x, y1, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 6); % 红色实线，带圆圈标记
%'b--s'
hold on; % 保持当前图形，以便在同一坐标轴上绘制第二组数据

% 绘制第二组数据，设置不同的线条格式
plot(x, y2, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6); % 蓝色虚线，带正方形标记

hold off; % 释放图形，不再绘制新的数据

ax = gca;

% 设置x轴刻度标签的大小
set(ax, 'FontSize', 12); % 设置刻度标签的字体大小为12，你可以根据需要调整这个值


% 设置坐标轴标签
ylabel('Time (s)', 'FontSize', 12); % x 轴标签
xlabel('Number of Samples', 'FontSize', 12); % y 轴标签



% 添加图例
legend("CNBSI", "CNBSRP");

% 添加标题
title('Running time of the two NBN algorithms', 'FontSize', 13);