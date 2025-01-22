function [] = ShowLON(figureId, filepath, savepath)
force_layout_iteration=100;
fileID = fopen(filepath);
task_name = fgetl(fileID);
info_s = fgetl(fileID);
graph_dim = str2double(info_s);
graph_edge = [];
formatSpec="%d";
if graph_dim~=0
sizeA = [2 graph_dim];
graph_edge= fscanf(fileID,formatSpec,sizeA);
info_s = fgetl(fileID);
end
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
maxNodeFunnel = max(node_funnel);
mapColor = autumn(maxNodeFunnel+1);
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

if graph_dim~=0 
    cur_figure = figure(figureId);
clf(cur_figure);
graph_s= graph_edge(1,:);
graph_t=graph_edge(2,:);
%s = [1 1 1 1 1 6 6 6 6 6];
%t = [2 3 4 5 6 7 8 9 10 11];
G = graph(graph_s,graph_t);
h1 = plot(G, 'Layout','force', 'Iterations',force_layout_iteration);
%hold on;
x= h1.XData;
y=h1.YData;    
cur_figure = figure(figureId);
clf(cur_figure);
h3= scatter3( x,y, node_fit,node_sizes, node_color,'filled');
hold on;
h2 = plot(G, 'XData',x, 'YData', y, 'ZData', node_fit, 'LineWidth',0.1, 'NodeColor', node_color,'EdgeColor',[0.5 0.5 0.5]);
ExportThreeFigures(savepath);
else
    h3= scatter3( [0], [0], node_fit,node_sizes, node_color,'filled');
    ExportThreeFigures(savepath);
end
end