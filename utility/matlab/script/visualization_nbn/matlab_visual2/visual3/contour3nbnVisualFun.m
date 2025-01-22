function []= contour3nbnVisualFun(NBN_visual_data, outputDir, filename,Psize,cut,contour3Level)

outputfilepath = append(outputDir,'/',filename);

%xi = NBN_visual_data.Fit2d_pos(:,1);
%yi = NBN_visual_data.Fit2d_pos(:,2);
%zi = NBN_visual_data.Fit2d_pos(:,3);
max_xi= max(NBN_visual_data.Fit2d_pos(:,1));
min_xi = min(NBN_visual_data.Fit2d_pos(:,1));
len_xi = (max_xi - min_xi)/cut;
max_yi = max(NBN_visual_data.Fit2d_pos(:,2));
min_yi = min(NBN_visual_data.Fit2d_pos(:,2));
len_yi = (max_yi- min_yi)/cut;
% Make the scattered interpolant.
F = scatteredInterpolant(NBN_visual_data.Fit2d_pos(:,1), NBN_visual_data.Fit2d_pos(:,2), NBN_visual_data.Fit2d_pos(:,3));
% Get a grid of points at every pixel location in the RGB image.
[xGrid, yGrid] = meshgrid(min_xi: len_xi: max_xi,  min_yi: len_yi: max_yi);
Z = F(xGrid,yGrid);
contour3(xGrid,yGrid,Z,contour3Level);
colormap('jet');
ax = gca;
ax.DataAspectRatio = [1 1 1];
set(gca,'visible','on');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
figurefilepath = append(outputfilepath,'_','contour3');
   fh= gcf;
   fh.WindowState = 'maximized';
   view(2);
   setExportFigureTypeJPG(figurefilepath,'top','eps','eps',0.1);
hold on;
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
nbn_plot= plot(NBN_visual_data.Fit2Dgraph,'XData',NBN_visual_data.Fit2d_pos(:,1),'YData',NBN_visual_data.Fit2d_pos(:,2),'ZData',NBN_visual_data.Fit2d_pos_z, 'NodeColor',  NBN_visual_data.idConColorFit2dPos);
highlight(nbn_plot,NBN_visual_data.NBN_opt_idxsFit2D,'Marker','+','MarkerSize',Psize,'NodeColor','k');
%colorbar;
figurefilepath = append(outputfilepath,'_','contour3NBN');
ExportThreeFiguresJPG(figurefilepath,'eps','eps');
end