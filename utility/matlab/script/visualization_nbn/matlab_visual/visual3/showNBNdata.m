function info= showNBNdata(figureId,outputDir,fileDir,filename,Psize)
    filepath = append(fileDir,'/',filename);
    info.ErrFlag = 0;
    info.ErrInfo = 'ok';
           NBN_originData = read_NBN_DataFun(filepath);
        NBN_originData.figureId= figureId;
        if NBN_originData.ErrFlag==1
            info.ErrFlag= 1;
            info.ErrInfo ='input Error';
            return;
        end 
        NBN_visual_data = transfer_NBN_visual_data(NBN_originData);
        NBN_visual_data.figureId= figureId;
        fh= figure(figureId);
        figure('visible','off');
        clf(fh);
        nbnVisualFun(NBN_visual_data, outputDir, filename,Psize);
        if NBN_originData.flag_continousFit2D~=0
            fh= figure(figureId);
            clf(fh);
            cut =   1000;
            contour3Level=200;
            contour3VisualFun(NBN_visual_data, outputDir, filename,Psize,cut,contour3Level);
            fh= figure(figureId);
            %figure('visible','off');
            clf(fh);
            numColor=1000;
            meshVisualFun(NBN_visual_data, NBN_originData,outputDir, filename,numColor,Psize)
        end
    
    
    try
        NBN_originData = read_NBN_DataFun(filepath);
        NBN_originData.figureId= figureId;
        if NBN_originData.ErrFlag==1
            info.ErrFlag= 1;
            info.ErrInfo ='input Error';
            return;
        end 
        NBN_visual_data = transfer_NBN_visual_data(NBN_originData);
        NBN_visual_data.figureId= figureId;
        fh= figure(figureId);
        figure('visible','off');
        clf(fh);
        nbnVisualFun(NBN_visual_data, outputDir, filename,Psize);
        if NBN_originData.flag_continousFit2D~=0
            fh= figure(figureId);
            clf(fh);
            cut =   1000;
            contour3Level=200;
            contour3VisualFun(NBN_visual_data, outputDir, filename,Psize,cut,contour3Level);
            fh= figure(figureId);
            %figure('visible','off');
            clf(fh);
            numColor=1000;
            meshVisualFun(NBN_visual_data, NBN_originData,outputDir, filename,numColor,Psize)
        end
    catch ME
        info.ErrFlag = 1;
        info.ErrInfo = ME.message;
    end    
end