function[]= calUADfunciton(x,y, map, numColor, picOutputDir,dataOutputDir, filename, method_name)
        %histogram(scores);
        %outputfilepath = outputDir + a +'iforest.png';
        %exportgraphics(gcf,outputfilepath);
        sz= size(x);
        data_mat = zeros(sz(1,1),2);
        data_mat(:,1)= x;
        data_mat(:,2)= y;

        [Mdl,tf,colorval] = ocsvm(data_mat);
        norScore = rescale(colorval, 0,1);
        colorval3= zeros(sz(1,1),3);
        for id = 1:sz(1,1)
            coloridx= int32(norScore(id)*numColor)+1;
            if coloridx>numColor
                coloridx= numColor;
            end
            colorval3(id,:)= map(coloridx,:);
        end
        scatter(y,x,10,colorval3);
        outputfilepath = picOutputDir + filename +'_df_'+method_name+'.png';
        exportgraphics(gcf,outputfilepath);

        outputfilepath = dataOutputDir + filename +'_df_'+method_name+'.txt';
        writematrix(colorval,outputfilepath);
end