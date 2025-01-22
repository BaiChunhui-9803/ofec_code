function [] = drawMesh2(meshData)
x = meshData(:,3)';
y=meshData(:,4)';
z= meshData(:,5)';
sampleSize = size(x);
div = sqrt(sampleSize(1,2));
xx= reshape(x,div,div);
yy = reshape(y,div,div);
zz = reshape(z,div,div);
mesh(xx,yy,zz);
end