
%dir = "/mnt/share151/";
dirName1 = "//172.24.242.8/share/";
dirName1 = "/mnt/Data/";
dirName2 = "//172.24.24.151/e/";
dirName2 = "/mnt/share151e/";
filepath =dirName1+"Student/2018/YiyaDiao/NBN_data/tsp_nbn/tsp_HyperSampling_nbn/";
filepath2= dirName2+"DiaoYiya/experiment_data/tsp_nbn_data/tsp_network_data/";
filepath2= dirName2+"tsp_nbn_data/tsp_network_data/";
pronamefile= dirName2 +"DiaoYiya/experiment_data/tsp_nbn_data/tsp_names.txt";
error_file_dir = dirName2 + "DiaoYiya/experiment_figure/TSP_figure/tsp_nbn/error_info2/";
mkdir(error_file_dir);

filetasks = readlines(pronamefile);
%filetasks= dir(filepath);
numFile= size(filetasks,1);
%\TSP_figure\tsp_nbn\figure
outputdir = dirName2+"DiaoYiya/experiment_figure/TSP_figure/tsp_nbn/figure3/";
mkdir(outputdir);
numThread =1; %设置处理器数目
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
try
    figureId= mod(k,numThread)+200;
    %f = figure(figureId);
    fileinfo=filetasks(k);
    fileinfo2 = split(fileinfo," ");
    if(fileinfo2=="") 
        continue;
    end
        display(fileinfo2);
    subDirName = fileinfo2(1);
    filename = fileinfo2(2);

    
    
    
    %if fileinfo2(3)=="df"
        drawSavedf(filepath2,outputdir,filename,figureId,numColor,basicNodeSize);
   % elseif fileinfo2(3)=="nbn"
        %display(filepath+filename+"/");
        
    %    drawSaveNBNTSP(filepath+subDirName+"/",filepath2,outputdir,filename,figureId, numColor,basicNodeSize);
  %  elseif fileinfo2(3)=="contour"
  %      drawSaveMesh(filepath,outputdir,filename,figureId,conT);

  %  end


catch ME
   fileID = fopen(error_file_dir+ fileinfo(1)+"_errorInfo.txt",'w');
    fprintf(fileID,'%s\n',ME.message);
    fclose(fileID);
end
end

%delete(par);
