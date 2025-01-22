function[] = drawMeshFun(x,y,z,conT)
            sampleSize = size(x);
            div = sqrt(sampleSize(1,2));
            x= reshape(x,div,div);
            y = reshape(y,div,div);
            z = reshape(z,div,div);
            f= figure(1);
            clf(f);
            pro_contour= contour3(x,y,z,conT);
            colormap("jet");
end