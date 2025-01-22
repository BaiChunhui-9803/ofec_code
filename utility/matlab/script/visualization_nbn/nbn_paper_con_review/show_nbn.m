 
nbn_filenpath = "//172.24.24.151/e/DiaoYiya/papre_con_data/review_data/shubert/Proname_Classic_Shubert_dim_2_origin_nbn.txt";
nbn_filenpath2 = "//172.24.24.151/e/DiaoYiya/papre_con_data/review_data/shubert/Proname_Classic_Shubert_dim_2_network.txt";
nbn_mat  = readmatrix(nbn_filenpath);
network_mat =readmatrix(nbn_filenpath2, "NumHeaderLines",1);
numColor = 1000;
basicNodeSize = 20;        
edge_a = nbn_mat(:,1)';
edge_b = nbn_mat(:,3)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = nbn_mat(:,5);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
        
for idx=1:num_node
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.9 0.9 0.9 ];
black_color = [0 0 0];

f = figure('visible','on');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', node_color , 'LineWidth',0.5 ,'EdgeColor',black_color);
nbn_plot.EdgeAlpha= 1;
nbn_plot.LineStyle = '-';
nbn_plot.Marker=".";
nbn_plot.MarkerSize=0.01;
%nbn_plot.Visible='off';
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;
