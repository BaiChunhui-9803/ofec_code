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
