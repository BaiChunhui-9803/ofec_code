%colorVal=rand(30,1); %data to be plotted
colorVal = [1,2,3,4,5; -1,2,3,-4,5];
sizeCV= size(colorVal);

colorVal(1,3) = 9;
%colorVal(1,30)
numColor= 200;
%colorIdx = colorVal;

[dimx,dimy] = size(colorVal);
minVal = 1e9;
for x=1:dimx
    for y=1:dimy
        if colorVal(x,y)>=0
            minVal = min(colorVal(x,y),minVal);
        end
        %colorIdx(x,y) = floor(((colorVal(x,y)-minValue)/range)*(numColor-1))+1;
    end
end
colorIdx = max(colorVal,minVal);
maxValue = max(colorIdx,[],'all') ;
range = maxValue - minValue;
for x=1:dimx
    for y=1:dimy
        colorIdx(x,y) = floor(((colorVal(x,y)-minValue)/range)*(numColor-1))+1;
    end
end
dim = size(colorVal);
map = jet(numColor);
color =  0;
if dimy==1
    color = zeros(dimx,3);
    for x= 1:dimx
        if colorVal(x,1)>=0
            color(x,:) = map(colorIdx(x,1),:);
        else
            color(x,:)=	[0 0 0];
        end
    end
elseif dimx==1  
    color = zeros(dimy,3);
    for y= 1:dimy
        if colorVal(y,1)>=0
            color(y,:) = map(colorIdx(1,y),:);
        else
            color(y,:)=	[0 0 0];
        end
    end    
else 
    color = zeros(dimx,dimy,3);
    for x=1:dimx
        for y=1:dimy
            if colorVal(x,y)>=0
                 color(x,y,:) = map(colorIdx(x,y),:);
            else
                 color(x,y,:)=	[0 0 0];
            end
        end
    end
end
