function[] = drawSaveMesh(filepath, outputdir, filename,figureId, conT)
            mesh_name =  replace(filename,"_Dim_2","_mesh");
            meshData = readmatrix(filepath+mesh_name+".txt");
            x = meshData(:,3)';
            y=meshData(:,4)';
            z= meshData(:,5)';
            
            sampleSize = size(x);
            div = sqrt(sampleSize(1,2));
            xx= reshape(x,div,div);
            yy = reshape(y,div,div);
            zz = reshape(z,div,div);
         
f= figure(figureId);
clf(f);
            pro_contour= contour3(xx,yy,zz,conT);
            ExportThreeFigures(outputdir+ filename +"_pro_contour" );
        %    hold on;
        %    fit_nbn_plot = plot(nbn_graph,'XData',nbn_nbnfit(:,2),'YData',nbn_nbnfit(:,3),'ZData',nbn_nbnfit(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
        %    hold on;
        %    fit_node_scatter = scatter3(nbn_nbnfit(:,2),nbn_nbnfit(:,3),nbn_nbnfit(:,4),node_size,node_color,'filled');
%clf(f);
%delete(pro_contour);
end