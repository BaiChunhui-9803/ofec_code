function []= concour2D_ExportFigure(n, xi,yi,zi,cut,contour3Level,filename, filetype, fontSize)
% n=1;
% cut=1000;
% Psize=10;
% EdgeB = readmatrix('EdgeB.txt');
% EdgeA = readmatrix('EdgeA.txt');
% xi = readmatrix('xi.txt');
% yi = readmatrix('yi.txt');
% zi = readmatrix('zi.txt');
% colorR = readmatrix('colorR.txt');
% colorG = readmatrix('colorG.txt');
% colorB = readmatrix('colorB.txt');
% markerPtrs = readmatrix('markerPtrs.txt');
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
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
fh = figure(n);
close(fh);
fh = figure(n);
fh.WindowState = 'maximized';
clf(n);
%clear current figure;
%surf(xGrid,yGrid,Z, 'EdgeColor', 'none');
contour(xGrid,yGrid,Z,contour3Level);
%colormap jet;
%colorbar;
hold on;
ax = gca; % current axes
ax.FontSize = fontSize;
ExportFigure(n,filename, filetype);
end 