function []= Graph3DColorMarker(n,EdgeB,EdgeA, X,Y,Z,colorR,colorG,colorB, markerPtrs,Psize)
colorArr3(:,1) = colorR(:);
colorArr3(:,2) = colorG(:);
colorArr3(:,3) = colorB(:);
g = graph(EdgeB,EdgeA);
%popNode =uint32(HnodeFrom):uint32(HnodeTo);
fh = figure(n);
fh.WindowState = 'maximized';
clf(n);
%h = plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', colorArr3);
highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize)
%h = plot(g,'XData',X,'YData',Y,'ZData',Z,'NodeColor', [1 0 0 0.5; 0 1 0 1; 0 1 0 0.4; 0 1 0 1; 0 1 0 0.3]);
%highlight(h,popNode,'NodeColor','red');
colorbar;
end 