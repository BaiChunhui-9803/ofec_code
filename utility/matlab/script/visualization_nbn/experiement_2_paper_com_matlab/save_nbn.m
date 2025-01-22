
%dir = "/mnt/share151/";
dirName = "D:/";
dirName = "/mnt/share151/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/lastConProPicTask.txt";
%filepath ='//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stn_data/visualization_nbn_stns_filternbn/';
%filepath = '//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborData_3_opt_alg/2/';
filepath = '//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/oneMax/onemaxNeighborData_4_opt_alg/2/';
filepath3 = filepath;
filepath2 = filepath;
%filepath3= '//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stns1/';
%mkdir(filepath3);
%filetasks = readlines(pronamefile);
filetasks= dir(filepath);
numFile= size(filetasks,1);

outputdir= '//172.24.24.151/e/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn100/';
outputdir = '//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/oneMax/onemax_Neighbor/';
error_file_dir = "E:/DiaoYiya/experiment_figure/nbnFigureErr/";
error_file_dir = 'E:/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn_error2/';
error_file_dir = '//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/oneMax/onemax_Neighbor/';
mkdir(error_file_dir);
%outputdir = "/mnt/share151e/DiaoYiya/experiment_figure/conProFigure2/";
mkdir(outputdir);
numThread =2; %设置处理器数目
%delete(gcp('nocreate'));
%par = parpool('local', numThread); %开启并行计算
error_number= zeros( numFile,2);
error_info = string(error_number);
numColor = 1000;
basicNodeSize = 20;
conT = 1000;


%parfor k=1:numFile
%for k=1:numFile
for k=1:numFile
%k=309;
%k=1;
%try
    %figureId= mod(k,numThread)+200;
    if filetasks(k).isdir==0
    %f = figure(figureId);
    fileinfo=filetasks(k).name;
%    display(fileinfo);
       if endsWith(fileinfo,"_algSols_nbn.txt")
    fileinfo2 = split(fileinfo,"_nbn");
    filename = fileinfo2(1);
    filename =filename{1,1};
    display(fileinfo);
    display(filename);
    
    
    nbn_filenpath =strcat(filepath,filename);
    nbn_filenpath=strcat(nbn_filenpath,"_nbn.txt");
    
        nbn_mat  = readmatrix(nbn_filenpath);
        
           nbn_filenpath =strcat(filepath3,filename);
   nbn_filenpath=strcat(nbn_filenpath,"_solAlgInfo.txt");
        
       nbn_info = readmatrix(nbn_filenpath);
        
       
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
gray_color = [0.9 0.9 0.9 ];
black_color = [0 0 0];

num_algNode= size(nbn_info,1);



f = figure('visible','off');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', node_color , 'LineWidth',0.5 ,'EdgeColor',gray_color);
nbn_plot.EdgeAlpha= 0.01;
nbn_plot.LineStyle = '-';
nbn_plot.Marker=".";
nbn_plot.MarkerSize=0.01;
%nbn_plot.Visible='off';
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;






set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    outputfilename =strcat(outputdir,filename);
    outputfilename=strcat(outputfilename,"_fitness_nbn");   
   filenameStr=sprintf('%s',outputfilename);
 %  filetypeStr=sprintf('%s',filetype);
   veiwNameStr=sprintf('%s','origin');
   outputfilename=[filenameStr,'-',veiwNameStr,'.jpg' ];
%   exportgraphics(gcf,outputfilename,'Resolution',600);
exportgraphics(gcf,outputfilename);


%   setExportFigureType(outputfilename,'origin',0.15);
%   view(2);
%   setExportFigureType(outputfilename,'top',0.15);
%      view([180 0]);
%   setExportFigureType(outputfilename,'front',0.15);
 
outputfilename=[];
show_node_color=[];



MshowIterId=[];

    fileinfo2 = [];
    filename = [];

nbn_mat = [];
nbn_info = [];
nbn_filenpath= [];
network_mat = [];
edge_a = [];
edge_b = [];
num_edge = [];
nbn_graph = [];

mapColor = [];
node_fit = [];
node_ColorIdx = [];
num_node = [];
node_color= [];
node_size = [];
gray_color = [];
black_color = [];
nbn_plot= [];
node_info_num =  [];
numNode=  [];
showNodeIdxs= [];
showIter=[];
showNBNinfoIdx=  [];

show_network =  [];
show_node_size=  [];
show_node_color=   [];

maxColorIdx=  [];
showMapColor =  [];

gray_color2 =  [];
showNBNnodeSize=  [];
startColor= [];
endcolor = [];
%	out << it.m_id << "\t" << it.m_visited_times << "\t" << it.m_start << "\t" << it.m_end << "\t" << it.fitness << "\t" << it.m_belongAlg << std::endl;
color1 = [];
color2 =  [];
color3= [];

       end

    end
    %    drawSaveNBN(filepath,filepath2,outputdir,filename,figureId, numColor,basicNodeSize);
%catch ME
%   fileID = fopen(error_file_dir+ filename+"_errorInfo.txt",'w');
%    fprintf(fileID,'%s\n',ME.message);
%    fclose(fileID);
%end
end

%delete(par);
