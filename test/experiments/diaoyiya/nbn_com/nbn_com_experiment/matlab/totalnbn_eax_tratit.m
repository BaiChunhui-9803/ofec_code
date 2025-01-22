readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_2/';
readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_lon_nbndata2/';
readdir= '//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_compareResults4_network/';
readDriTrait= "\\172.24.242.8\share\Student\2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_4_traitdata/";
readdir="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_hnswEqualBetter_network/";
readdir="\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\eax_lkh_nbn_hnswEqualBetter_network\";
readdir ="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_com_experiment_data\totalTsp3\eax_lkh_nbn_hnswRnd_network\";
readdir="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_com_experiment_data\totalTsp3\eax_lkh_lon_nbndata_rnd\";
readdir ="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_com_experiment_data\totalTsp3\eax_lkh_lon_nbndata_rnd\";
readdir= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore\";
readdir= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore_network\";
readdir ="\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\eax_lkh_nbn_rnd_lkhmore5\";
readDriTrait=readdir;
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/totalNBN_eax_figure_test/';
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/eax_lkh_nbn_hnswEqualBetter_network_figure/';
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/eax_lkh_lon_rnd_figure/';
mkdir(outputdir);
str = ["1281" "2202" "6702" "6717" "u574","5955" ];
str = ["u574" "2202"  "5955" ];
str = ["5955"  ];

strSize = size(str);
strSize= strSize(1,2);

basicNodeSize = 8;
numColor = 1000;
mapColor = jet(numColor);
solMapColor = jet(numColor);

for taskId= 1:strSize
clf;
tspname = str(1,taskId);
    
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
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.5 ,'EdgeColor',gray_color);
nbn_plot.EdgeAlpha=0.2;
hold on;


%nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%nbn_points.Marker= 'o';
%nbn_points.MarkerFaceAlpha=0.5;

%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;

formatSpec="%f";

bestIdName= readDriTrait + tspname + "_bestIds.txt";
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
bestPos(:,4)=bestPos(:,4)+0.05; 
hold on;
best_points=  scatter3(bestPos(:,2),bestPos(:,3),bestPos(:,4),'filled');
best_points.MarkerFaceColor= 'none';
best_points.MarkerEdgeColor= 'k';
best_points.Marker= 'square';
best_points.SizeData= 50;











traitname = ["lkh_fail", "lkh_success", "eax_fail", "eax_success"  ];
displayName = ["LKH's failture", "LKH's success", "EAX's failture", "EAX's success"  ];

%traitname = ["eax_success" ];
traitFileName= readDriTrait+ tspname+ "_traits.txt";
fileID = fopen(traitFileName,'r');
formatSpec = '%f';
    for tId=1:1
        
 %      curTraitName= traitname(1, tId);
       curTraitName= traitname(tId);
       curDisName = displayName(tId);
      % disp(curTraitName);
       sizeM= fscanf(fileID,formatSpec,sizeA);
       if (sizeM(1,2)>0)
           data= fscanf(fileID,formatSpec,sizeM);
           sizeM= sizeM(1,2);
           data_idx = data(1,:);
           data_idx = data_idx+1;
           data_value= 1-data(2,:);
           %data_value=numColor- int64(data_value.*numColor);
           data_value=numColor- int64(data_value.*numColor);
           data_value=max(1,data_value);
           %data_value=data_value +1;
           solPos= zeros(sizeM,5); 
           for idx=1:sizeM
               solPos(idx,:)=  network_mat(data_idx(idx),:);
           end
           solPos(:,4)=solPos(:,4)+0.02; 
           
           solColor = zeros(sizeM,3);
           for idx=1:sizeM
               
                solColor(idx,:) = solMapColor(data_value(idx),:);
%                 solColor(idx,:) = black_color;
           end
           
           solSize = zeros(sizeM,1);
           for idx=1:sizeM
                solSize(idx)=10;
           end
           
           hold on;
           
           solPoints=  scatter3(solPos(:,2),solPos(:,3),solPos(:,4),10,solColor,'filled','hexagram','DisplayName',curDisName);
           %solPoints=  scatter3(solPos(:,2),solPos(:,3),solPos(:,4),solSize,solColor, "filled");
           solPoints.LineWidth= 1;
           solPoints.Marker= "hexagram";
           
           curoutfilename =  outputdir + tspname + "_trait_"+ curTraitName;
          
           

           view(3);
%setExportFigureType(curoutfilename,'origin',0.15);
%view(2);
%setExportFigureType(curoutfilename,'to',0.15);
%view([180 0]);
%setExportFigureType(curoutfilename,'view180',0.15);
  %         solPoints.visible='false';
 %set(solPoints, 'Visible', 'off');
           
       end
    end



end