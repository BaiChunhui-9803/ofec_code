[X,Y] = meshgrid(-8:.5:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;
dim = size(X);
dimxx= dim(1,1)*dim(1,1);
           % div = sqrt();
x= reshape(X,1,dimxx);
y = reshape(Y,1,dimxx);
z = reshape(Z,1,dimxx);
conT=1000;
drawMeshFun(x, y, z, conT);
