dir = '\\172.24.6.26\e\DiaoYiya\paper_cswidn\cswidn_sampling_phase_ananlysis2\';

alg_data = importdata([dir 'alg_data.txt']);
com_data = importdata([dir 'com_data.txt']);
con_data = importdata([dir 'con_data.txt']);
con_sub = importdata([dir 'con_sub.txt']);
mix_data = importdata([dir 'mix_data.txt']);

width = 800; % 假设宽度为 800 像素
height = width * 11 / 25; % 根据宽高比计算高度
figure('Position',[0,0,width,height]);
hold on;

yLenName=  [""];
showId=2;

if showId==5

alg_data(1,:) = []; % 删除第一行
com_data(1,:) = []; % 删除第一行
con_data(1,:) = []; % 删除第一行
con_sub(1,:) = []; % 删除第一行
mix_data(1,:) = []; % 删除第一行

end

alg_data(:,5)=1.0- alg_data(:,5);
com_data(:,5)=1.0- com_data(:,5);
con_data(:,5)=1.0- con_data(:,5);
con_sub(:,5)=1.0- con_sub(:,5);
mix_data(:,5)=1.0- mix_data(:,5);


% time = data(:, 1);
% Ir = data(:, 2);
% In = data(:, 3);
% opt = data(:, 4);
% diff = data(:, 5);


% Plot second column of com_data with dotted blue line
h1= plot(com_data(:,1), com_data(:,showId),'m');

% Plot second column of con_data with solid green line
h2= plot(con_data(:,1), con_data(:,showId),'b');

% Plot second column of con_sub with dash-dot cyan line
h3= plot(con_sub(:,1), con_sub(:,showId),'g');

% Plot second column of mix_data with a thicker magenta solid line
h4= plot(mix_data(:,1), mix_data(:,showId),'k');
% Plot second column of alg_data with dashed red line
h5 = plot(alg_data(:,1), alg_data(:,showId),'r');

h6= xline(1,'r--'); % Red vertical line for first interval start
xline(29,'r--'); % Red vertical line for first interval end
h7=xline(37,'b--'); % Blue vertical line for second interval end
xline(6,'b--'); % Blue vertical line for second interval start
h8=xline(46,'g--'); % Green vertical line for third interval end
xline(14,'g--'); % Green vertical line for third interval start

legend([h1, h2, h3, h4, h5, h6,h7,h8], {'$  \mathbf S_{com}(\mathbf o) $', '$ \mathbf S_{con}(\mathbf o) $', '$\mathbf S_{con}(\mathbf o, 1e^{-6})$', '$ \mathbf S_{o}$', '$ \mathbf S_{alg}$', 'Injection Time of $o_{l_1}$','Injection Time of $o_{l_2}$','Injection Time of $o_{l_3}$'},'Interpreter','latex', 'FontSize', 12,'Location','eastoutside');
%legend([h1, h2, h3, h4, h5, h6,h7,h8], {'com data', 'con data', 'con sub data', 'mix data', 'alg data', 'p1','p2','p3'});
legend show;
xlabel('Time $t$','Interpreter','latex', 'FontSize', 12);
ylabel('Ruggedness $I_r$','Interpreter','latex','FontSize', 12);
%ylabel('Neutrality $I_n$','Interpreter','latex','FontSize', 12);
%ylabel('Number of optima','Interpreter','latex','FontSize', 12);
%ylabel('$\| \mathbf  S(t) , \mathbf S(t-1)  \|$','Interpreter','latex','FontSize', 12);
hold on;
%title('Comparison of Second Columns from Different Datasets', 'FontSize', 13);
% hFigure = figure;
% % 进行绘图操作
% set(hFigure,'PaperPositionMode','auto');
% set(hFigure,'PaperUnits','inches');
% paperWidth = 10; % 设置纸张宽度（单位：英寸）
% paperHeight = 4; % 设置纸张高度（单位：英寸）
% set(hFigure,'PaperSize',[paperWidth,paperHeight]); % 设置纸张大小


%set(hFigure,'PaperPositionMode','auto');
%set(hFigure,'PaperUnits','inches');
%set(hFigure,'PaperSize',[width height]); % 设置纸张大小，这里的 width 和 height 是你想要的图像大小，单位为英寸
%print('5.png','-dpng','-r600');


