
%dir = "/mnt/share151/";
dirName = "D:/";
dirName = "/mnt/share151/";
filepath =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8/";
filepath2= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8_df/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/lastConProPicTask2.txt";
error_file_dir = dirName + "student/2018/DiaoYiya/nbnConPro_experiment_data/error_info3/";
mkdir(error_file_dir);

filetasks = readlines(pronamefile);
%filetasks= dir(filepath);
numFile= size(filetasks,1);
outputdir = "E:/DiaoYiya/experiment_figure/conProFigure2/";
outputdir =  "/mnt/share151e/DiaoYiya/experiment_figure/conProFigure3/";
mkdir(outputdir);
numThread =13; %设置处理器数目
delete(gcp('nocreate'));
par = parpool('local', numThread); %开启并行计算
error_number= zeros( numFile,2);
error_info = string(error_number);
numColor = 1000;
basicNodeSize = 20;
conT = 1000;


parfor k=1:numFile
%for k=1:numFile
%k=309;
%k=1;
try
    figureId= mod(k,numThread)+200;
    %f = figure(figureId);
    fileinfo=filetasks(k);
    fileinfo2 = split(fileinfo," ");
    filename = fileinfo2(1);
    display(fileinfo2);
    if fileinfo2(2)=="df"
        drawSavedf(filepath2,outputdir,filename,figureId,numColor,basicNodeSize);
    elseif fileinfo2(2)=="nbn"
        drawSaveNBN(filepath,filepath2,outputdir,filename,figureId, numColor,basicNodeSize);
    elseif fileinfo2(2)=="contour"
        drawSaveMesh(filepath,outputdir,filename,figureId,conT);

    end


catch ME
   fileID = fopen(error_file_dir+ fileinfo(1)+"_errorInfo.txt",'w');
    fprintf(fileID,'%s\n',ME.message);
    fclose(fileID);
end
end

%delete(par);
