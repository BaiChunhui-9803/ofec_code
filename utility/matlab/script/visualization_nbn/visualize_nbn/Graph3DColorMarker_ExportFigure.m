function []= Graph3DColorMarker_ExportFigure(n,EdgeB,EdgeA, X,Y,Z,colorR,colorG,colorB, markerPtrs,Psize,filename, filetype, fontSize,colorbar_label, axis_label)
colorArr3(:,1) = colorR(:);
colorArr3(:,2) = colorG(:);
colorArr3(:,3) = colorB(:);
g = graph(EdgeB,EdgeA);
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
fh = figure(n);
%fh.WindowState = 'maximized';
fh.Position = [0,0,1080,1080];
clf(n);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', colorArr3);
highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize)

%h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', [1 0 0 0.5; 0 1 0 1; 0 1 0 0.4; 0 1 0 1; 0 1 0 0.3]);
%highlight(h,popNode,'NodeColor','red');
if colorbar_label==1
    colorbar;
end 
hold on;
if axis_label == 1
    ax=gca
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    ax.ZAxis.Visible='off';
    axis off;
else
    ax = gca; % current axes
    ax.FontSize = fontSize;
end
ExportFigure(n,filename, filetype);
end 