function[] = drawNBNfun(figureId,edge_a,edge_b,x,y,z,node_fit,numColor,basicNodeSize)
     %   nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
     %   network_mat =readmatrix(filepath2+ filename + "_network.txt", "NumHeaderLines",1);
     %   nbn_solInfo = readmatrix(filepath + filename + "_nbnSolInfo.txt", "NumHeaderLines",1);


%edge_a = nbn_mat(:,1)';
%edge_b = nbn_mat(:,3)';
%num_edge =size(edge_a,2);
%for idx= 1:num_edge
%    edge_a(idx) = edge_a(idx)+1;
%    edge_b(idx) = edge_b (idx) + 1;
%end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
%node_fit = nbn_solInfo(:,2);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,2);
node_color=  zeros(num_node,3);


for idx=1:num_node
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
%black_color = [0 0 0];

f= figure(figureId);
clf(f);


nbn_plot = plot(nbn_graph,'XData',x,'YData',y,'ZData',z,'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
hold on;
nbn_points=  scatter3(x,y,z,node_size,node_color,'filled');

end