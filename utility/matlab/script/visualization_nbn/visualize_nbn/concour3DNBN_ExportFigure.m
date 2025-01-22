function []= concour3DNBN_ExportFigure(n,EdgeB,EdgeA, xi,yi,zi,colorR,colorG,colorB, cut, markerPtrs,Psize,contour3Level,filename, filetype, fontSize)
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
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
g = graph(EdgeB,EdgeA);
fh = figure(n);
clf(n);
fh.WindowState = 'maximized';
%figure ('DefaultAxesFontSize',30); 
%clear current figure;
%surf(xGrid,yGrid,Z, 'EdgeColor', 'none');
contour3(xGrid,yGrid,Z,contour3Level);
%colormap jet;
%colorbar('westoutside');
hold on;
h = plot(g,'XData',xi,'YData',yi,'ZData',zi,'NodeColor', colorArr3);
highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize);
hold on;
ax = gca; % current axes
ax.FontSize = fontSize;
ExportFigure(n,filename, filetype);
end 