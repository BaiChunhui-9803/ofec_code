function[] = drawSaveMeshFun2(x,y,z,conT,figureId,filepath)
            sampleSize = size(x);
            div = sqrt(sampleSize(1,2));
            x= reshape(x,div,div);
            y = reshape(y,div,div);
            z = reshape(z,div,div);
f = figure(figureId);
clf(f);
         contour(x,y,z,conT);
            colormap("jet");
            exportgraphics(gcf,filepath);
end