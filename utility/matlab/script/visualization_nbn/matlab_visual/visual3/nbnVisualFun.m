function []= nbnVisualFun(NBN_visual_data, outputDir, filename,Psize)
fh = figure(NBN_visual_data.figureId);
%figure('visible','off');
fh.WindowState = 'maximized';
nbn_plot = plot(NBN_visual_data.NBNgraph,'XData',NBN_visual_data.NBNpos(:,1),'YData',NBN_visual_data.NBNpos(:,2),'ZData',NBN_visual_data.NBNpos(:,3),'NodeColor', NBN_visual_data.idConColor);
nbn_plot.Marker = 'o';
nbn_plot.LineStyle = ':';
dimO  = size(NBN_visual_data.NBN_opt_idxs);
if(dimO(1)==1&&dimO(1)==1&&NBN_visual_data.NBN_opt_idxs==-1)
else 
    highlight(nbn_plot,NBN_visual_data.NBN_opt_idxs,'Marker','+','MarkerSize',Psize);
end
ax = gca;
ax.DataAspectRatio = [1 1 1];
set(gca,'visible','off');
%colorbar;
outputfilepath = append(outputDir,'/',filename);
figurefilepath = append(outputfilepath,'_','nbn');
ExportThreeFiguresJPG(figurefilepath,'eps','eps');
clf(fh);
nbn_plot = plot(NBN_visual_data.NBNgraph,'XData',NBN_visual_data.NBNpos(:,1),'YData',NBN_visual_data.NBNpos(:,2),'ZData',NBN_visual_data.NBNpos(:,3),'NodeColor', NBN_visual_data.conFitColor);
nbn_plot.Marker = 'o';
nbn_plot.LineStyle = ':';
if(dimO(1)==1&&dimO(1)==1&&NBN_visual_data.NBN_opt_idxs==-1)
else 
    highlight(nbn_plot,NBN_visual_data.NBN_opt_idxs,'Marker','+','MarkerSize',Psize);
end

ax = gca;
ax.DataAspectRatio = [1 1 1];
set(gca,'visible','off');
%colorbar;
outputfilepath = append(outputDir,'/',filename);
figurefilepath = append(outputfilepath,'_','nbn');
ExportThreeFiguresJPG(figurefilepath,'eps','eps');


%   fh= gcf;
%   fh.WindowState = 'maximized';
%   view(2);
%   ax = gca;
%ax.DataAspectRatio = [1 1 1];
%set(gca,'visible','off');
%figurefilepath = append(outputfilepath,'_','nbnFitColor');
%   setExportFigureTypeJPG(figurefilepath,'top','eps','eps',0.1);
end