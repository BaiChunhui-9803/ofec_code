numColor = 1000;
basicNodeSize = 20;
conT = 1000;      
figureId=4;
filename="";
outputdir = "E:\Diao_Yiya\paper\nbn_experiment\nbn_testErrData\";
filepath2="";
filepath="E:\Diao_Yiya\paper\nbn_experiment\nbn_testErrData\";
display(filepath+filename+"_nbn.txt");
pronamefile= filepath+filename+"_nbn.txt";

filetasks = readlines(pronamefile);
        nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
        network_mat =readmatrix(filepath2+ filename + "_network.txt", "NumHeaderLines",1);
        nbn_solInfo = readmatrix(filepath + filename + "_nbnSolInfo.txt", "NumHeaderLines",1);


edge_a = nbn_mat(:,1)';
edge_b = nbn_mat(:,3)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = nbn_solInfo(:,2);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
numConstrait=0;
numDynamic=0;
for idx=1:num_node
    if nbn_solInfo(idx,4)~=1
        numConstrait= numConstrait+1;
    end
    if nbn_solInfo(idx,5)~=0
        numDynamic= numDynamic +1;
    end
end
for idx=1:num_node
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];


if numConstrait~=0
    nodeContrait_color= node_color;
    for idx=1:num_node
    if nbn_solInfo(idx,4)~=1
    nodeContrait_color(idx,:) =black_color;
    end
    end
end

if numDynamic~=0
    node_dynamic_color = node_color;
    node_ColorIdx = uint32(rescale(nbn_solInfo(:,5),1,numColor));
    for idx=1:num_node
    node_dynamic_color(idx,:) =mapColor(node_ColorIdx(idx),:);
    end
end


f= figure(figureId);
clf(f);


nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
ExportThreeFigures(outputdir+ filename +"_fitness_nbn" );

