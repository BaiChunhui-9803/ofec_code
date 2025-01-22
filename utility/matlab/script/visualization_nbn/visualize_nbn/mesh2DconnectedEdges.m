function []= mesh2DconnectedEdges(n,EdgeB,EdgeA, xi,yi,zi,colorR,colorG,colorB, markerPtrs,Psize)
cut = 3000;
max_xi= max(xi);
min_xi = min(xi);
len_xi = (max_xi - min_xi)/cut;
max_yi = max(yi);
min_yi = min(yi);
len_yi = (max_yi- min_yi)/cut;
% Make the scattered interpolant.
F = scatteredInterpolant(xi, yi, zi);
% Get a grid of points at every pixel location in the RGB image.
[xGrid, yGrid] = meshgrid(min_xi: len_xi: max_xi,  min_yi: len_yi: max_yi);
Z = F(xGrid,yGrid);
colorArr3(:,1) = colorR(:);
colorArr3(:,2) = colorG(:);
colorArr3(:,3) = colorB(:);
g = graph(EdgeB,EdgeA);
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
figure(n);
surf(xGrid,yGrid,Z, 'EdgeColor', 'none');
colorbar;
hold on;
h = plot(g,'XData',xi,'YData',yi,'ZData',zi,'NodeColor', colorArr3);
highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize);
end 