function []= meshNbnVisualFun(NBN_visual_data,NBN_originData, outputDir, filename,numColor,Psize)
outputfilepath = append(outputDir,'/',filename);
[k,numSample]  = size(NBN_originData.mesh2D_colors);
numDiv =sqrt(numSample);
divSpace= 1.0/(numDiv-1);
[X,Y] = meshgrid(0:divSpace:1);
Z=zeros(numDiv,numDiv);
Zcolorval =zeros(numDiv,numDiv);
for idx=1:numDiv
    for idy=1:numDiv
        id = (idx-1)*numDiv + idy;
        x= idy;
        y=idx;
        Z(x,y) = NBN_originData.mesh2D_zval(1,id);
        Zcolorval(x,y,:) =NBN_originData.mesh2D_colors(1,id)*NBN_originData.mesh2D_colors(2,id); 
    end
end
maxZ =  max(Z,[],'all') ;
minZ =  min(Z,[],'all');
rangeZ = maxZ- minZ;
for idx=1:numDiv
    for idy=1:numDiv
        Z(idx,idy)=(Z(idx,idy)-minZ)/rangeZ;
    end
end

Zcolor = ValueToColor(Zcolorval,numColor);
f = figure(NBN_visual_data.figureId);
clf(f);
mesh(X,Y,Z,Zcolor);   
hold on;
ax = gca;
%ax.DataAspectRatioMode ='auto';
ax.DataAspectRatio = [1 1 1];
%axis(ax,'tight')
set(gca,'visible','on');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
view(2);
figurefilepath = append(outputfilepath,'_','mesh');
fh= gcf;
fh.WindowState = 'maximized';
view(2);
setExportFigureTypeJPG(figurefilepath,'top','eps','eps',0.1);
%xi = NBN_visual_data.Fit2d_pos(:,1);
%yi = NBN_visual_data.Fit2d_pos(:,2);
%zi = NBN_visual_data.Fit2d_pos(:,3);
%zi= NBN_visual_data.Fit2d_pos_oz;
[numPos,tmp]= size(NBN_visual_data.Fit2d_pos_oz);
for idx=1:numPos
    NBN_visual_data.Fit2d_pos_oz(idx)  = (NBN_visual_data.Fit2d_pos_oz(idx)-minZ)/rangeZ +0.005;
end
nbn_plot= plot(NBN_visual_data.Fit2Dgraph,'XData',NBN_visual_data.Fit2d_pos(:,1),'YData',NBN_visual_data.Fit2d_pos(:,2),'ZData', NBN_visual_data.Fit2d_pos_oz, 'NodeColor',  NBN_visual_data.idConColorFit2dPos);
highlight(nbn_plot,NBN_visual_data.NBN_opt_idxsFit2D,'Marker','+','MarkerSize',Psize,'NodeColor','k');
figurefilepath = append(outputfilepath,'_','meshNBN');
filename = figurefilepath;
filetype = 'eps';
fileSuffix= 'eps';
  fh= gcf;
   fh.WindowState = 'maximized';
   view([-37.5000 30]);
   setExportFigureTypeJPG(filename,'origin',filetype,fileSuffix,0.1);
   view(2);
   setExportFigureTypeJPG(filename,'top',filetype,fileSuffix,0.1);
      view([180 0]);
   setExportFigureTypeJPG(filename,'front',filetype,fileSuffix,0.1);
   
%ExportThreeFigures(figurefilepath,'eps','eps');
%hold on;
%view([-37.5000 30]);
end