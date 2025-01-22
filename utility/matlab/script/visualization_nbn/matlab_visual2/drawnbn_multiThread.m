
%dir = "/mnt/share151/";
dirName = "D:/";
dirName = "/mnt/share151/";
filepath =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8/";
filepath2= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8_df/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/lastConProPicTask.txt";
error_file_dir = dirName + "student/2018/DiaoYiya/nbnConPro_experiment_data/error_info/";
mkdir(error_file_dir);
filepath = '//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/';
filepath = 'E:/DiaoYiya/experiment_data/con_visualization_nbn_stn/';
filepath2 = 'E:/nbn_data/visualization_nbn_image/';
filepath2 = 'E:/DiaoYiya/experiment_data/con_visualization_nbn_stn/';
%filetasks = readlines(pronamefile);
filetasks= dir(filepath);
numFile= size(filetasks,1);

outputdir = "E:/DiaoYiya/experiment_figure/nbnFigure/";
outputdir = 'E:/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn/';
error_file_dir = "E:/DiaoYiya/experiment_figure/nbnFigureErr/";
error_file_dir = 'E:/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn_error/';
mkdir(error_file_dir);
%outputdir = "/mnt/share151e/DiaoYiya/experiment_figure/conProFigure2/";
mkdir(outputdir);
numThread =40; %设置处理器数目
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
    
    
    nbn_filenpath =strcat(filepath,filename);
    nbn_filenpath=strcat(nbn_filenpath,"_nbn.txt");
    
        nbn_mat  = readmatrix(nbn_filenpath);
        
       
    nbn_filenpath =strcat(filepath2,filename);
    nbn_filenpath=strcat(nbn_filenpath,"_network.txt");     
        
        
        network_mat =readmatrix(nbn_filenpath, "NumHeaderLines",1);
        
        
edge_a = nbn_mat(:,1)';
edge_b = nbn_mat(:,3)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = nbn_mat(:,5);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
        
for idx=1:num_node
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];


clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;


    outputfilename =strcat(outputdir,filename);
    outputfilename=strcat(outputfilename,"_fitness_nbn");   

           fh= gcf;
   fh.WindowState = 'maximized';
   view([-37.5000 30]);
   setExportFigureType(outputfilename,'origin',0.1);
   view(2);
   setExportFigureType(outputfilename,'top',0.1);
      view([180 0]);
   setExportFigureType(outputfilename,'front',0.1);
        
       end

    end
    %    drawSaveNBN(filepath,filepath2,outputdir,filename,figureId, numColor,basicNodeSize);
catch ME
   fileID = fopen(error_file_dir+ filename+"_errorInfo.txt",'w');
    fprintf(fileID,'%s\n',ME.message);
    fclose(fileID);
end
end

delete(par);
