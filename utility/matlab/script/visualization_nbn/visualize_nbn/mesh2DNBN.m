function []= mesh2DNBN(n,EdgeB,EdgeA, xi,yi,zi,colorR,colorG,colorB, cut, markerPtrs,Psize)
n=1;
cut=1000;
Psize=10;
writematrix(EdgeB,'EdgeB.xls');
writematrix(EdgeA,'EdgeA.xls');
writematrix(xi,'xi.xls');
writematrix(yi,'yi.xls');
writematrix(zi,'zi.xls');
writematrix(colorR,'colorR.xls');
writematrix(colorG,'colorG.xls');
writematrix(colorB,'colorB.xls');
writematrix(markerPtrs,'markerPtrs.xls');



readmatrix(EdgeB,'EdgeB.xls');
readmatrix(EdgeA,'EdgeA.xls');
readmatrix(xi,'xi.xls');
readmatrix(yi,'yi.xls');
readmatrix(zi,'zi.xls');
readmatrix(colorR,'colorR.xls');
readmatrix(colorG,'colorG.xls');
readmatrix(colorB,'colorB.xls');
readmatrix(markerPtrs,'markerPtrs.xls');

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
clf(n);
%clear current figure;
surf(xGrid,yGrid,Z, 'EdgeColor', 'none');
colormap hsv;
colorbar;
hold on;
h = plot(g,'XData',xi,'YData',yi,'ZData',zi,'NodeColor', colorArr3);
highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize);
end 