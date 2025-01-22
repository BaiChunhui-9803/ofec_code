fileDir = '//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_figure2';
fileDir = 'Z:/student/2018/diaoyiya/nbn_experiment/corProNBN_df_cp2/';
dirOfDirFile= dir(fileDir);
numOfDirFiles = length(dirOfDirFile);
outputDir= '//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_DisFit_figure/';
mkdir(outputDir);
coreNum =1; %设置处理器数目
delete(par) %关闭并行计算
par = parpool('local', coreNum); %开启并行计算
parfor iter_file= 1:numOfDirFiles
    filename = dirOfDirFile(iter_file).name;
    if (strcmp(filename,'.'))|| (strcmp(filename,'..'))
        disp(filename);
    else
        filepath = append(fileDir,'/',filename);
        subfilename = split(filename,".");
        disp(subfilename(1));
        disp(filepath);            
        filepath=[outputDir,filename,'nbn2d.png'];
        readfilepath = [fileDir, filename];
        figureId = mod(iter_file,coreNum)+1;
        f= figure(figureId);
        f = figure('visible','off');
        data = readmatrix(readfilepath);
        colorval = data(:,3);
x = data(:,1)';
y=data(:,2)';
colorval3 = ValueToColor(colorval,1e3);
scatter(y,x,10,colorval3);
   a = string(subfilename(1,1)) ;
        outputfilepath = outputDir + a +'.png';
        exportgraphics(gcf,outputfilepath);
    end
end

%delete(par) %关闭并行计算