readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_2/';
readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/nbn_data_lonSample3/';
%readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_lon_nbndata2/';

outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/totalNBN_lon_figure_final/';
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/lonNBN_lon_figure_final/';
outputdir= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/nbn_total_figure_final/"

mkdir(outputdir);
str = [ "2202" "u574","5955" ];
str = ["u574"];
strSize = size(str);
strSize= strSize(1,2);

basicNodeSize = 4;
numColor = 1000;
sizeA= [1 2]; 

gray_color= [0.6 0.6 0.6];
%gray_color2= [0.5 0.5 0.5];

for taskId= 1:strSize
clf;
tspname = str(1,taskId);
nbnfilename= tspname;
nbnfilename= tspname +"_lonSamples";
lonFunInfoName= tspname+ "_lonFunnelIdsMatlab.txt";

fileID = fopen(readdir+lonFunInfoName,'r');
lonFunInfoSize= fscanf(fileID,formatSpec,sizeA);
lonFunInfo= fscanf(fileID,formatSpec,lonFunInfoSize);

    nbn_filenpath=readdir+ nbnfilename+"_network.txt";    
    network_mat =readmatrix(nbn_filenpath, "NumHeaderLines",1);
maxNodeFunnel = max(lonFunInfo);
mapColor = jet(maxNodeFunnel+1);
    
mapColor(2,:) = [1 0 0 ];
edge_a = network_mat(:,1)';
edge_b = network_mat(:,5)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

%mapColor = jet(numColor);
node_fit = network_mat(:,4);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
        
for idx=1:num_node
 %  node_color(idx,:) = mapColor(numColor-node_ColorIdx(idx)+1,:);
 %  node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
    node_color(idx,:) = gray_color;
end

node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end

cursize= lonFunInfoSize(1,2)-1;
for idx=1:cursize
    if lonFunInfo(1, idx)>= 0

            node_color(idx+1,:) = mapColor(lonFunInfo(1, idx)+1,:);
                node_size(idx+1)=15;
                if lonFunInfo(1, idx)== 1
              %          disp(idx)
              %  node_size(idx)=50;
                end
    elseif lonFunInfo(1,idx)==-1
            node_color(idx+1,:) = gray_color;
    elseif lonFunInfo(1,idx)==-2
            node_color(idx+1,:) = gray_color;
    end
end

%gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];

%f = figure('visible','off');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',1.0 ,'EdgeColor',gray_color);
nbn_plot.EdgeAlpha=0.2;

hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
    

hold on;    

min_value = max(node_fit);
indices = find(node_fit == min_value);

if tspname=="u574"
% 删除最后一个元素
indices(end) = [];
end

% 获取这些点的坐标
selected_XData = nbn_points.XData(indices);
selected_YData = nbn_points.YData(indices);
selected_ZData = nbn_points.ZData(indices);

% 用黑色正方形画出来
h4 = scatter3(selected_XData, selected_YData,selected_ZData, 'filled','DisplayName','Optima $\mathbf o$');
h4.MarkerFaceColor= 'none';
h4.MarkerEdgeColor= 'k';
h4.Marker= 'square';
h4.LineWidth = 2.0;
h4.SizeData= 100;


curoutfilename= outputdir + tspname ;



view(3);

set (gca,'position',[0.1,0.1,0.9,0.9] );
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
   ax = gca;
   ax.OuterPosition = [0 0 1 1];
   outerpos = ax.OuterPosition; % 获取外部框位置

%setExportFigureType(curoutfilename,'origin',0.15);
%view(2);
%setExportFigureType(curoutfilename,'to',0.15);
%view([180 0]);
%setExportFigureType(curoutfilename,'view180',0.15);




 
    
end