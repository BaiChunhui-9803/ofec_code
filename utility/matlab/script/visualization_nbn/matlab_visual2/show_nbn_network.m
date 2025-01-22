
%dir = "/mnt/share151/";
dirName = "D:/";
dirName = "/mnt/share151/";
filepath =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8/";
filepath2= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8_df/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/lastConProPicTask.txt";
error_file_dir = dirName + "student/2018/DiaoYiya/nbnConPro_experiment_data/error_info/";
mkdir(error_file_dir);
filepath = '//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/';
filepath2 = 'E:/nbn_data/visualization_nbn_image/';
%filetasks = readlines(pronamefile);
filetasks= dir(filepath);
numFile= size(filetasks,1);

outputdir = "E:/DiaoYiya/experiment_figure/nbnFigure/";
%outputdir = "/mnt/share151e/DiaoYiya/experiment_figure/conProFigure2/";
mkdir(outputdir);
numThread =5; %设置处理器数目
delete(gcp('nocreate'));
%par = parpool('local', numThread); %开启并行计算
error_number= zeros( numFile,2);
error_info = string(error_number);
numColor = 1000;
basicNodeSize = 20;
conT = 1000;


%parfor k=1:numFile
for k=1:numFile
%k=309;
%k=1;
%try
    figureId= mod(k,numThread)+200;
    if filetasks(k).isdir==0
    %f = figure(figureId);
    fileinfo=filetasks(k).name;
%    display(fileinfo);
       if endsWith(fileinfo,"_nbn.txt")
               fileinfo2 = split(fileinfo,"_nbn");
    filename = fileinfo2(1);
    filename =filename{1,1};
    display(fileinfo);
    display(filename);
    
    
    nbn_filenpath =[];
    
        nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
        network_mat =readmatrix(filepath2+ filename + "_network.txt", "NumHeaderLines",1);
       end

    end
    %    drawSaveNBN(filepath,filepath2,outputdir,filename,figureId, numColor,basicNodeSize);
%catch ME
 %  fileID = fopen(error_file_dir+ fileinfo(1)+"_errorInfo.txt",'w');
  %  fprintf(fileID,'%s\n',ME.message);
   % fclose(fileID);
%end
end

%delete(par);
