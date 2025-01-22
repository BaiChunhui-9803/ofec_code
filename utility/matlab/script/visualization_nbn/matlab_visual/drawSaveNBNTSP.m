function [] = drawSaveNBNTSP(filepath, filepath2,outputdir, filename,figureId,numColor,basicNodeSize)
       % nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
        network_mat =readmatrix(filepath2+ filename + "_network.txt", "NumHeaderLines",1);

num_edge =size(network_mat,1);
edge_a = zeros(num_edge,1);
edge_b = network_mat(:,5)';

for idx= 1:num_edge
    edge_a(idx) = idx;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = network_mat(:,4);
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
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];

f= figure(figureId);
clf(f);


nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
ExportThreeFigures(outputdir+ filename +"_fitness_nbn" );

set(nbn_points, 'visible','off');
delete(nbn_points);
delete(nbn_plot);

            clear("black_color");
            clear("gray_color");
            clear("node_size");
            clear("node_color");
            clear("num_node");
            clear("node_ColorIdx");


                     clear("node_fit");
         clear("mapColor");
         clear("edge_a");
         clear("num_edge");
         clear("edge_b");

      clear("nbn_graph");
         clear("network_mat");
                  clear("nbn_mat");

                  clf(f);

end