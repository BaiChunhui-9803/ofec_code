
%dir = "/mnt/share151/";
dirName = "D:/";
filedir =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_2018linux/";
meshFileDir= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total4/";
filenames_mat = readlines(dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/nbnConFilenames.txt");
nbnNetworkFile = dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network/";
numOfDirFiles = size(filenames_mat);
numOfDirFiles = numOfDirFiles(1,1);


nbn_name = "_nbn";
mesh_name = "_mesh";
txt_name =".txt";
picture_name ="_pic";
df_name = "_df";
outputdir =  dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_totalfigure/";
mkdir(outputdir);
numThread =40; %设置处理器数目
delete(par) %关闭并行计算
delete(gcp('nocreate'));
par = parpool('local', numThread); %开启并行计算

%delete(par) %关闭并行计算
error_number= zeros( length(files),2);
error_info = string(error_number);
parfor k=1:numOfDirFiles
    figureId= mod(numThread,k)+1;
    f = figure(figureId);
    filename=filenames_mat(k,1);
    display(filename);
    try
%    if endsWith(filename,"nbn.txt")

%draw mesh 
      %  if  endsWith(filename, "_2")
      %      clf(f);
      %      meshfilename =filename+"_mesh.txt";
      %      mesh_mat= readmatrix(meshFileDir+meshfilename);
      %      drawMesh(mesh_mat,100);
      %      meshpicturename=  replace(meshfilename,txt_name,picture_name);
      %      ExportThreeFigures(outputdir+meshpicturename);
      %  end

        
        nbn_mat  = readmatrix(filedir+filename+"_nbn.txt");
        network_mat =readmatrix(nbnNetworkFile+ filename + ".txt");
        

     %   clf(f);
     %   drawNBN(nbn_mat);
     %   nbnpicturename=filename + "_distanefitness";
     %  setExportFigureType(outputdir+nbnpicturename,'origin',0.1);




        %        ExportThreeFigures(outputdir+nbnpicturename);
%    end
    catch ME
        display(ME.message);
%        error_info(k,1)= filename;
%        error_info(k,2)= ME.message;
    end
end
delete(par) %关闭并行计算
writematrix(error_info,'error.txt');