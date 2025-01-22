readdir= "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore5\";
readDriTrait = "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore5_network_analysis/";
readdir = "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_more_hnsw_rnd_2\";
readdir= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_more_hnsw_rnd_3_linux\";
readdir = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\eax_lkh_more_hnsw_rnd_3_linux\";
optIdDir= readdir;

readDriTrait="\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore5_network_analysis\";
readDriTrait= readdir;
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/eax_lkh_lon_nbndata_rnd_optIds/';
outputdir= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/nbn_eax_trait_figures2\";
mkdir(outputdir);
str = ["1281" "2202" "6702" "6717" "u574","5955" ];
str = ["u574" "2202"  "5955" ];
str = ["5955"  ];
strSize = size(str);
strSize= strSize(1,2);
for taskId= 1:strSize
clf;
tspname = str(1,taskId);
traitname = ["_lkh_trait", "_eax_trait" ];
traitname= "_lkh_trait";
matlabTraitDataFilenName= readDriTrait+ tspname+ traitname+ "_0";
% whether found opt
found_opt=0;

traitFileName= matlabTraitDataFilenName+"_matlab.txt";
shortestPathName= matlabTraitDataFilenName +"_shortestPath_matlab.txt";
offset= 0.01;
sizeA= [1 2]; 

basicNodeSize = 8;
numColor = 1000;
mapColor = jet(numColor);
solMapColor = jet(numColor);


    
nbn_filenpath =strcat(readdir,tspname);
nbn_filenpath=strcat(nbn_filenpath,"_nbn.txt");
      %  nbn_mat  = readmatrix(nbn_filenpath);
nbn_filenpath =strcat(readdir,tspname);
nbn_filenpath=strcat(nbn_filenpath,"_network.txt");    
        
network_mat =readmatrix(nbn_filenpath, "NumHeaderLines",1);
curoutfilename= outputdir + tspname + "_eax_lkh";

            
edge_a = network_mat(:,1)';
edge_b = network_mat(:,5)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);


node_fit = network_mat(:,4);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
        
for idx=1:num_node
%    node_color(idx,:) = mapColor(numColor-node_ColorIdx(idx)+1,:);
       node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];

f = figure('visible','on');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.5 ,'EdgeColor',gray_color,'DisplayName', 'NBN');
nbn_plot.EdgeAlpha=0.2;
hold on;


%nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%nbn_points.Marker= 'o';
%nbn_points.MarkerFaceAlpha=0.5;

%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;

formatSpec="%f";









displayName = ["LKH's failture", "LKH's success", "EAX's failture", "EAX's success"  ];





%algTraitName = traitname(1,tId);
%traitname = ["eax_success" ];

disp(traitFileName);
fileID = fopen(traitFileName,'r');
formatSpec = '%f';
%    for tId=1:4
        
 %
     %  curTraitName= traitname(1, tId);
%       curTraitName= algTraitName;
      % curDisName = displayName(tId);
      % disp(curTraitName);
       sizeM= fscanf(fileID,formatSpec,sizeA);
       if (sizeM(1,2)>0)
           data= fscanf(fileID,formatSpec,sizeM);
           sizeM= sizeM(1,2);
           data_idx = data(1,:);
           data_idx = data_idx+1;
           data_value= 1-data(2,:);
           %data_value=numColor- int64(data_value.*numColor);
           data_value=numColor- int64(data_value.*numColor)+1;
           %data_value=data_value +1;
           solPos= zeros(sizeM,5); 
           for idx=1:sizeM
               solPos(idx,:)=  network_mat(data_idx(idx),:);
           end
           solPos(:,4)=solPos(:,4)+offset;
           
           solColor = zeros(sizeM,3);
           for idx=1:sizeM
                solColor(idx,:) = solMapColor(data_value(idx),:);
%                 solColor(idx,:) = black_color;
           end
           
           solSize = zeros(sizeM,1);
           for idx=1:sizeM
                solSize(idx)=20;
           end
           
           hold on;
           
           solPoints=  scatter3(solPos(:,2),solPos(:,3),solPos(:,4),10,solColor,'filled','hexagram','DisplayName','Trajectory $\it T$');
           %solPoints=  scatter3(solPos(:,2),solPos(:,3),solPos(:,4),solSize,solColor, "filled");
           solPoints.LineWidth= 1;
           solPoints.Marker= "hexagram";
           
       %    curoutfilename =  outputdir + tspname + "_trait_"+ curTraitName;
           hold on;
           



fileID = fopen(shortestPathName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
shorestPathIds= fscanf(fileID,formatSpec,sizeM);
sizeM= sizeM(1,2);
shorestPathIds=shorestPathIds+1;
shorestPathPos= zeros(sizeM,5);
for idx=1:sizeM
    shorestPathPos(idx,:)=  network_mat(shorestPathIds(idx),:);
end
shorestPathPos(:,4)=shorestPathPos(:,4)+offset; 
path=shorestPathIds;
%path =[];
% 创建路径子图
subgraph_path = subgraph(nbn_graph, path);



           if tspname=="5955"
 nearestId=[5815];
sizeM=1;
nearPos= zeros(sizeM,5);

for idx=1:sizeM
    nearPos(idx,:)=  network_mat(nearestId(idx),:);
end

nearPos(:,4)=nearPos(:,4)+offset; 
hold on;
near_points=  scatter3(nearPos(:,2),nearPos(:,3),nearPos(:,4),'filled','DisplayName','Decep.');
near_points.MarkerFaceColor= 'none';
near_points.MarkerEdgeColor= 'b';
near_points.Marker= 'diamond';
near_points.SizeData= 60;

            hold on;
end




bestIdName= readDriTrait + tspname + "_localOptima_matlab.txt";
%bestIdName = "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore5_network_analysis\"
disp(bestIdName);
%bestIdName=optIdName;
sizeA= [1 2]; 
fileID = fopen(bestIdName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
bestIds= fscanf(fileID,formatSpec,sizeM);
bestIds=bestIds+1;
sizeM= sizeM(1,2);
bestPos= zeros(sizeM,5);
for idx=1:sizeM
    bestPos(idx,:)=  network_mat(bestIds(idx),:);
end
bestPos(:,4)=bestPos(:,4)+offset; 
hold on;
best_points=  scatter3(bestPos(:,2),bestPos(:,3),bestPos(:,4),'filled','DisplayName','Local optima');
best_points.MarkerFaceColor= 'none';
best_points.MarkerEdgeColor= 'k';
best_points.Marker= 'o';

best_points.SizeData= 50;

optIdName= optIdDir + tspname + "_bestIds.txt";
sizeA= [1 2]; 
fileID = fopen(optIdName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
optIds= fscanf(fileID,formatSpec,sizeM);
optIds=optIds+1;
sizeM= sizeM(1,2);
optPos= zeros(sizeM,5);
for idx=1:sizeM
    optPos(idx,:)=  network_mat(optIds(idx),:);
end
optPos(:,4)=optPos(:,4)+offset; 
hold on;
opt_points=  scatter3(optPos(:,2),optPos(:,3),optPos(:,4),'filled','DisplayName','Optima $\mathbf o$');
opt_points.MarkerFaceColor= 'none';
opt_points.MarkerEdgeColor= 'r';
opt_points.Marker= 'o';
opt_points.SizeData= 60;
opt_points.LineWidth=1;

hold on;


if found_opt==0
% 绘制路径子图
nbn_plot_path = plot(subgraph_path, 'XData', network_mat(path, 2), 'YData', network_mat(path, 3), 'ZData', network_mat(path, 4)+offset, 'NodeColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', '$    \tilde{P}(\it T, \mathbf o)  $');
           nbn_plot_path.LineWidth= 1;
           nbn_plot_path.Marker= "square";
           nbn_plot_path.NodeLabel = [];
 %          nbn_plot_path.set

end
% 绘制路径子图
nbn_plot_path = scatter3(network_mat(shorestPathIds, 2), network_mat(shorestPathIds, 3),  network_mat(shorestPathIds, 4)+offset,30,'k', 'DisplayName', ' $\it{P}(\it T, \mathbf o)$ ');
           nbn_plot_path.LineWidth= 1;
           nbn_plot_path.Marker= "square";
           

           hold on;


           %legend('show');
           lgd = legend('show');
           lgd.FontSize = 15;
set(lgd, 'Interpreter','latex')
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
           view(3);
%setExportFigureType(curoutfilename,'origin',0.15);
%view(2);
%setExportFigureType(curoutfilename,'to',0.15);
%view([180 0]);
%setExportFigureType(curoutfilename,'view180',0.15);
         %  solPoints.visible='false';
  
           
       end
%    end



end