function []= showGraph3DHighligh(n,EdgeB,EdgeA, X,Y,Z,ColorVal, HnodeFrom, HnodeTo)
g = graph(EdgeB,EdgeA);
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
figure(n);
num_color=4;
num_node=5;
colorArr = zeros(num_node,num_color);
for n = 1:num_node
    for c = 1:num_color
        colorArr(n,c) = rand(1, 1, 'double');
    end
end
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', colorArr);
h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', [1 0 0 0.5; 0 1 0 1; 0 1 0 0.4; 0 1 0 1; 0 1 0 0.3]);
%highlight(h,popNode,'NodeColor','red');
colorbar;
end 