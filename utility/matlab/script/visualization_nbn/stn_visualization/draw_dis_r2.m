[X,Y] = meshgrid(0.1:0.01:1.0);
A = (2*3.1415926*X);
B= Y.^2;
C= B./X;
D= Y.^2./(X.*2);
D = exp(-D);
A = A.^(-0.5);
Z= A.*D;
%Z = (2*3.1415926*X).^(-0.5)*exp(-(Y.^2./(X.*2))) ;
contour3(X,Y,Z,100);
          fh= gcf;
   view([-37.5000 30]);
      ax = gca;
   ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.ZAxis.FontSize = 12;
   xlabel('$\| \mathbf a,\mathbf b\|$','Interpreter','latex',"FontSize",15) ;
ylabel('$r$ ','Interpreter','latex',"FontSize",15) ;
zlabel(' $p_m(\mathbf a \leftarrow \mathbf b)$','Interpreter','latex',"FontSize",15)