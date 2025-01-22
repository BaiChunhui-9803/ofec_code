topDir ="D:";
readdir = topDir+ "/student/2018/DiaoYiya/tspLonGenerator/";
%savedir = dir+"";
filename= "att532";
lon_filename = ["_lon_filter001" "_lon_filter005"];
filedirpah = filename+"/"+ filename;
saveDir = "E:/DiaoYiya/experiment_figure/TSP_figure/LON/";
mkdir(saveDir);
force_layout_iteration= 100;
tsp_filenames= dir(readdir);
dir_size=  size(tsp_filenames,1);

for dirId = 1:dir_size
filename= tsp_filenames(dirId).name;
if(tsp_filenames(dirId).isdir)
display(filename);

    for fileId=1:2

%fileId =2 ;
figureId=1;
filepath =  readdir+ filedirpah +lon_filename(fileId) +".txt";
%filepath = "D:/student/2018/DiaoYiya/tspLonGenerator/att532/att532_lon_filter001.txt";
saveFilename =filename+ lon_filename(fileId);
%savefilepath = saveDir+ saveFilename;
%display(savefilepath);

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
f1 = figure(figureId);
clf(f1);
h = plot(G, 'Layout','force', 'Iterations',force_layout_iteration,'visible','off');
hold on;
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



h3= scatter3( h.XData, h.YData, node_fit,node_sizes, node_color,'filled');
hold on;
h2 = plot(G, 'XData', h.XData, 'YData', h.YData, 'ZData', node_fit, 'LineWidth',0.1, 'NodeColor', node_color,'EdgeColor',[0.5 0.5 0.5]);


%ExportThreeFigures(saveDir +saveFilename);

    end
end
end
%,node_size,node_color);

%G = graph(graph_t, graph_s);

%formatSpec="%d %d";
%sizeA= [2 graph_dim];

