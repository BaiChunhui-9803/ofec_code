[X,Y] = meshgrid(-8:0.017:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;
C = X.*Y;


cmap = colormap(jet);
%s.FaceAlpha = 0.5;
max_c = max(C,[],'all');
min_c = min(C,[],'all');
len_c= max_c  - min_c;
total_color = 255;
[rows,column]=size(C);
norC = C;
for i=1:rows
for j=1:column
    norC(i,j) = (C(i,j)-min_c)/len_c*total_color;
end
end
norC = round(norC);
for i=1:rows
for j=1:column
    norC(i,j) = norC(i,j)+1;
end
end
C_rgb =  ind2rgb(norC, cmap);
for i=1:100
for j=1:100
    C_rgb(i,j,1) = 0;
    C_rgb(i,j,2) = 0;
    C_rgb(i,j,3) = 0;
end
end
s = mesh(X,Y,Z,C_rgb);
%colormap jet;
%colorbar;
cc = s.CData;
%C_rgb =  ind2rgb(cc, cmap);
s.CDataMapping = 'scaled';
s.EdgeColor = 'interp';

xx = reshape(X.',1,[]);
yy = reshape(Y.',1,[]);
zz = reshape(Z.',1,[]);

EdgeA =uint32(2000):uint32(400000);
EdgeB = randi([1 1000000],1,400000-2000+1);
totalNum =  400000-2000+1;
colorR =  rand(totalNum, 1, 'double');
colorG =  rand(totalNum, 1, 'double');
colorB =  rand(totalNum, 1, 'double');
colorArr3(:,1) = colorR(:);
colorArr3(:,2) = colorG(:);
colorArr3(:,3) = colorB(:);
g = graph(EdgeB,EdgeA);
h = plot(g,'XData',xx,'YData',yy,'ZData',zz,'NodeColor', colorArr3);
%markerPtrs=[1 3 9 10 22 33 55 ];
%Psize=10;
%highlight(h,markerPtrs,'Marker','h','MarkerSize',Psize)