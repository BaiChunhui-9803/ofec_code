function [] = drawNBN(nbnData)
df_x= nbnData(:,4);
df_y = nbnData(:,5);
numPoints=0;
sz_nbn= size(nbnData);
for idx=1:sz_nbn(1,1)
    if nbnData(idx,1)~=nbnData(idx,3)
        numPoints = numPoints+1;
    end
end
x = zeros(numPoints,1);
y = zeros(numPoints,1);
numPoints=0;
for idx=1:sz_nbn(1,1)
    if nbnData(idx,1)~=nbnData(idx,3)
        numPoints= numPoints+1;
        x(numPoints,1)= df_x (idx,1);
        y(numPoints,1) = df_y (idx,1);
    end
end
scatter(x,y,10);
end