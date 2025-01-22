
filedir="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_cswidn\cswidn_pro_info\";
formatSpec= "%f";
sizeA= [1 2]; 
fileID = fopen(filedir+"edges.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
edges= fscanf(fileID,formatSpec,sizeM);

fileID = fopen(filedir+"positions.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
points= fscanf(fileID,formatSpec,sizeM);
points=points';

fileID = fopen(filedir+"lakeIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
lakes= fscanf(fileID,formatSpec,sizeM);
lakes= lakes+1;


fileID = fopen(filedir+"riverIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
riverIds= fscanf(fileID,formatSpec,sizeM);
riverIds= riverIds+1;

fileID = fopen(filedir+"tankIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
tankIds= fscanf(fileID,formatSpec,sizeM);
tankIds= tankIds+1;


fileID = fopen(filedir+"junctionIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
junctionIds= fscanf(fileID,formatSpec,sizeM);
junctionIds= junctionIds+1;

fileID = fopen(filedir+"sensorIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
sensorIds= fscanf(fileID,formatSpec,sizeM);
sensorIds= sensorIds+1;


fileID = fopen(filedir+"optIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
optIds= fscanf(fileID,formatSpec,sizeM);
optIds= optIds+1;


fileID = fopen(filedir+"path.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
pathIds= fscanf(fileID,formatSpec,sizeM);
pathIds= pathIds+1;


edge_a = edges(1,:)+1;
edge_b= edges(2,:)+1;
cswidnNet = digraph(edge_a,edge_b);

gray_color = [0 0 0];
cswidn_plot = plot(cswidnNet,'XData',points(:,1),'YData',points(:,2),'NodeColor', gray_color , 'LineWidth',1 ,'EdgeColor',gray_color,'DisplayName', 'Network');
cswidn_plot.NodeLabel = [];
cswidn_plot.MarkerSize=3;

markerSize = 100;
lenthwidth=3;

hold on;


lakePos = points(lakes,:); % 根据 optIds 选择对应的点
lake_points=  scatter(lakePos(:,1),lakePos(:,2),'filled','DisplayName','Lake');
lake_points.MarkerFaceColor= 'none';
lake_points.MarkerEdgeColor= 'blue';
lake_points.Marker= 'o';
lake_points.SizeData= markerSize;

hold on;
riverPos = points(riverIds,:); % 根据 optIds 选择对应的点
river_points=  scatter(riverPos(:,1),riverPos(:,2),'filled','DisplayName','River');
river_points.MarkerFaceColor= 'none';
river_points.MarkerEdgeColor= 'green';
river_points.Marker= 'diamond';
river_points.SizeData= markerSize;

hold on;





% junctionPos = points(junctionIds,:); % 根据 optIds 选择对应的点
% junction_points=  scatter(junctionPos(:,1),junctionPos(:,2),'filled','DisplayName','Junction');
% junction_points.MarkerFaceColor= 'none';
% junction_points.MarkerEdgeColor= 	"#EDB120";
% junction_points.Marker= '^';
% junction_points.SizeData= 30;

hold on;
tankPos = points(tankIds,:); % 根据 optIds 选择对应的点
tank_points=  scatter(tankPos(:,1),tankPos(:,2),'filled','DisplayName','Tank');
tank_points.MarkerFaceColor= 'none';
tank_points.MarkerEdgeColor= "#D95319";
tank_points.Marker= '^';
tank_points.SizeData= markerSize;

hold on;

sensorPos = points(sensorIds,:); % 根据 optIds 选择对应的点
sensor_points=  scatter(sensorPos(:,1),sensorPos(:,2),'filled','DisplayName','Sensor');
sensor_points.MarkerFaceColor= 'none';

sensor_points.MarkerEdgeColor= '#FF00FF';
sensor_points.Marker= 'hexagram';
sensor_points.SizeData= markerSize;



hold on;



optPos = points(optIds,:); % 根据 optIds 选择对应的点
opt_points=  scatter(optPos(:,1),optPos(:,2),'filled','DisplayName','Source');
opt_points.MarkerFaceColor= 'none';
opt_points.MarkerEdgeColor= 'r';
opt_points.Marker= 'square';
opt_points.SizeData= markerSize;
opt_points.LineWidth=1.5;
opt_points.UserData.index = [1 2 3];
hold on;


numPoints = size(optPos,1);


% 设置文本对象数组用于显示下标
textObjects = gobjects(numPoints,1);
for i = 1:numPoints
    % 在点的位置显示 LaTeX 编码的下标
    textObjects(i) = text(optPos(i,1)-2.5,optPos(i,2),['$o_{l_{',num2str(i),'}}$'],'FontSize',15,'Color','k','Interpreter','latex');
end
%x_{l_i}
% pathPos = points(pathIds,:); % 根据 optIds 选择对应的点
% path_points=  scatter(pathPos(:,1),pathPos(:,2),'filled','DisplayName','Source');
% path_points.MarkerFaceColor= 'none';
% path_points.MarkerEdgeColor= 'g';
% path_points.Marker= 'o'; 
% path_points.SizeData= markerSize;
% path_points.LineWidth=1.5;


           lgd = legend('show');
           lgd.FontSize = 15;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);


