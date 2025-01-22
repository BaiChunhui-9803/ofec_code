nbndir ="\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_cswidn\cswidn_algorithm_data_problemPhase_network/";
nbndir= "\\172.24.6.26\e\DiaoYiya\paper_cswidn\cswidn_algorithm_data_problemPhase_network\";
algdir = "\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_cswidn\cswidn_algorithm_runProcess_algorithmData/";
algdir = "\\172.24.6.26\e\DiaoYiya\paper_cswidn\cswidn_algorithm_runProcess_algorithmData\";
phasetime= "48";
runtime= "20";
nbnfilename = "CSIWDN_phase_"+ phasetime;
algfilename= "algData_phase_"+ phasetime+"_numRun_"+ runtime;
sizeA= [1 2]; 
formatSpec ="%f";
fileID = fopen(algdir+algfilename +"_solIds.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
solId= fscanf(fileID,formatSpec,sizeM);
solId=solId+1;

fileID = fopen(algdir+algfilename +"_solIter.txt",'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
solIter= fscanf(fileID,formatSpec,sizeM);




    basicNodeSize = 40;
    numColor = 1000;
      
    nbn_filenpath =strcat(nbndir,nbnfilename);
    nbn_filenpath=strcat(nbn_filenpath,"_network.txt");    
    network_mat =readmatrix(nbn_filenpath, "NumHeaderLines",1);
     
 solPos = network_mat(solId,:); % 根据 optIds 选择对应的点
 


edge_a = network_mat(:,1)';
edge_b = network_mat(:,5)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = network_mat(:,4);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
        
for idx=1:num_node
    %node_color(idx,:) = mapColor(numColor-node_ColorIdx(idx)+1,:);
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];





sol_ColorIdx = uint32(rescale(solIter,1,numColor));
num_node =  size(solIter,2);
sol_color=  zeros(num_node,3);
        
for idx=1:num_node
    %node_color(idx,:) = mapColor(numColor-node_ColorIdx(idx)+1,:);
    sol_color(idx,:) = mapColor(sol_ColorIdx(idx),:);
end
sol_size = zeros(num_node,1);
for idx=1:num_node
    sol_size(idx)=basicNodeSize;
end

%f = figure('visible','off');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.5 ,'EdgeColor',gray_color);
%nbn_plot.EdgeAlpha=0.2;
hold on;
nbn_points=  scatter3(solPos(:,2),solPos(:,3),solPos(:,4),sol_size,sol_color,'filled');





