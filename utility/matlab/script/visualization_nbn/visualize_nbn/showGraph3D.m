function []= showGraph3D(n,EdgeB,EdgeA, X,Y,Z,ColorVal)
g = graph(EdgeB,EdgeA);
figure(n);
plot(g,'XData',X,'YData',Y,'ZData',Z, 'NodeCData', ColorVal);
colorbar;
end 