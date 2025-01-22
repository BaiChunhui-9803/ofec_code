
%filedir = "//172.24.24.151/f/DiaoYiya/2018/DiaoYiya/tspLonGenerator/";
filedir= "//172.24.24.151/e/DiaoYiya/code/data/ofec-data/paper_com_experiment_data/tsp_typical_comparison_exp3/";
filedir = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\lon_data";
outputdir= "//172.24.24.151/f/DiaoYiya/paper_com_experiment_data/total_tsp/lon_figure_final/";
outputdir= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/lon_figure_final/"
tspname = "u574";
curoutdir= outputdir + tspname+"/";
mkdir(curoutdir);
    
    force_layout_iteration= 100;
%subfix= "_lon_001.txt";

str = ["1281" "2202" "6702" "6717" "u574" ];
str= ["5955"];


str = ["2202"  "u574" "5955"];
strSize = size(str);
strSize= strSize(1,2);

for taskId= 1:strSize
clf;
tspname = str(1,taskId);


subfix= "_lon_filter001.txt";
filepath = filedir+ tspname+"/"+ tspname +subfix;

curoutfilename= curoutdir + tspname + "_filter001";



fileID = fopen(filepath);
task_name = fgetl(fileID);
info_s = fgetl(fileID);
graph_dim = str2double(info_s);
formatSpec="%d";
sizeA = [2 graph_dim];
graph_edge= fscanf(fileID,formatSpec,sizeA);
info_s = fgetl(fileID);
info_s = fgetl(fileID);
node_size = str2double(info_s);
sizeA = [1 node_size];
node_funnel = fscanf(fileID,formatSpec,sizeA);
info_s = fgetl(fileID);
info_s = fgetl(fileID);
node_size = str2double(info_s);
sizeA = [1 node_size];
node_sizes = fscanf(fileID,formatSpec,sizeA);
info_s = fgetl(fileID);
info_s = fgetl(fileID);
node_size = str2double(info_s);
sizeA = [1 node_size];
node_fit = fscanf(fileID,formatSpec,sizeA);
graph_s= graph_edge(1,:);
graph_t=graph_edge(2,:);
%s = [1 1 1 1 1 6 6 6 6 6];
%t = [2 3 4 5 6 7 8 9 10 11];
G = graph(graph_s,graph_t);

h = plot(G, 'Layout','force', 'Iterations',force_layout_iteration,'visible','off');
h.NodeLabel = {};
hold on;
maxNodeFunnel = max(node_funnel);
mapColor = jet(maxNodeFunnel+1);

mapColor(2,:) = [1 0 0 ];
node_color = zeros(node_size,3);
gray_color = [0.5 0.5 0.5];
for idx=1:node_size
    node_sizes(idx)  = node_sizes(idx)+1;
    if node_funnel(idx)==-1
        node_color(idx,:)= gray_color;
    else
        node_color(idx,:) = mapColor(node_funnel(idx)+1,:);
    end
end



h3= scatter3( h.XData, h.YData, node_fit,node_sizes, node_color,'filled');
hold on;
h2 = plot(G, 'XData', h.XData, 'YData', h.YData, 'ZData', node_fit, 'LineWidth',0.1, 'NodeColor', node_color,'EdgeColor',[0.5 0.5 0.5]);
h2.NodeLabel = {};

hold on;

% 假设已经有了 h、G、h2、h3 等变量
% 获取 node_fit 的值
node_fit_values = node_fit;

% 找到 node_fit 最小的一些点的索引
[min_fit, min_indices] = min(node_fit_values); % num_points_to_select 是你想要选择的最小点的数量
min_value = min(node_fit);
indices = find(node_fit == min_value);


% 获取这些点的坐标
selected_XData = h3.XData(indices);
selected_YData = h3.YData(indices);
selected_ZData = h3.ZData(indices);


% 用黑色正方形画出来
h4 = scatter3(selected_XData, selected_YData,selected_ZData, 'filled','DisplayName','Optima $\mathbf o$');
h4.MarkerFaceColor= 'none';
h4.MarkerEdgeColor= 'k';
h4.Marker= 'square';
h4.LineWidth = 2.0;
h4.SizeData= 150;


h1.DisplayName = '';
h2.DisplayName = '';
h3.DisplayName = '';
%h4.DisplayName = 'YourDesiredDisplayNameForH4';



           hold on;


           %legend('show');
%           hl = legend('show');
%set(hl, 'Interpreter','latex')
%set(gca,'XTick',[],'YTick',[],'ZTick',[]);

 
view(3);
setExportFigureType(curoutfilename,'origin',0.15);
view(2);
setExportFigureType(curoutfilename,'to',0.15);
view([180 0]);
setExportFigureType(curoutfilename,'view180',0.15);
view([270 0]);
setExportFigureType(curoutfilename,'view270',0.15);





end
