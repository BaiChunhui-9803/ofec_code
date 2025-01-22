function color = ValueToColor(colorVal,numColor)

[dimx,dimy] = size(colorVal);
minValue = 1e9;
for x=1:dimx
    for y=1:dimy
        if colorVal(x,y)>=0
            minValue = min(colorVal(x,y),minValue);
        end
        %colorIdx(x,y) = floor(((colorVal(x,y)-minValue)/range)*(numColor-1))+1;
    end
end
colorIdx = max(colorVal,minValue);
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
        if colorVal(y)>=0
            color(y,:) = map(colorIdx(y),:);
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

end