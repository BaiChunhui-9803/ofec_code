
%dir = "/mnt/share151/";
dirName = "D:/";
dirName = "/mnt/share151/";
filepath =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8/";
filepath2= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8_df/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/lastConProPicTask.txt";
filepath = '//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/';
filepath = 'E:/DiaoYiya/experiment_data/con_visualization_nbn_stn/';
filepath = '//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_stns_filternbn/';
filepath ='//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stn_data/visualization_nbn_stns_filternbn/';

filepath2 = 'E:/nbn_data/visualization_nbn_image/';
filepath2 = 'E:/DiaoYiya/experiment_data/con_visualization_nbn_stn/';
filepath2 = '//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_stns_filternbn/';
filepath2 ='//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stn_data/visualization_nbn_stns_filternbn';
filepath2 ='//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stn_data/visualization_nbn_stns_filternbn/';
filepath3 ='//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stn_data/visualization_nbn_stns/';
%filepath3= '//172.24.24.151/e/DiaoYiya/experiment_data/visualization_nbn_stns1/';
%mkdir(filepath3);
%filetasks = readlines(pronamefile);
filetasks= dir(filepath);
numFile= size(filetasks,1);

outputdir = "E:/DiaoYiya/experiment_figure/nbnFigure/";
outputdir = 'E:/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn2/';
outputdir= '//172.24.24.151/e/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn100/';
error_file_dir = "E:/DiaoYiya/experiment_figure/nbnFigureErr/";
error_file_dir = 'E:/DiaoYiya/experiment_data/conVisualizationFigure_nbn_stn_error2/';
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
for k=1:3
%k=309;
%k=1;
%try
    %figureId= mod(k,numThread)+200;
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
        
            nbn_filenpath =strcat(filepath3,filename);
    nbn_filenpath=strcat(nbn_filenpath,"_nbnCurInfo.txt");
        
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

a
f = figure('visible','on');
clf;
nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
nbn_plot.EdgeAlpha= 0.01;
nbn_plot.LineStyle = '--';
nbn_plot.Marker=".";
nbn_plot.MarkerSize=0.01;
%nbn_plot.Visible='off';
hold on;
%nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%ExportThreeFigures();
%outputfilename = outputdir+ filename +"_fitness_nbn" ;

node_info_num = size(nbn_info,1);
numNode= 0;
showNodeIdxs= [];
showIter=[];
showNBNinfoIdx= 3;

for idx=1:node_info_num
    if nbn_info(idx,showNBNinfoIdx)~= -1
        numNode= numNode+1;
        showNodeIdxs(end+1) = idx;
        showIter(end+1)=nbn_info(idx,showNBNinfoIdx);
    end
    

end

showIter = sort(showIter);
showIterId= zeros(1,size(showIter,2));
for idx=1:size(showIter,2)
    showIterId(idx)=idx;
end
MshowIterId = containers.Map(showIter,showIterId);
mapColor = jet(size(showIter,2));
showNBNinfoIdx= 2;
show_network = zeros(numNode,size(network_mat,2));
show_node_size= zeros(numNode,size(node_size,2));
show_node_color= zeros(numNode,size(node_color,2));

maxColorIdx= max(nbn_info(:,showNBNinfoIdx));
showMapColor = jet(maxColorIdx+1);

gray_color2 = [0.8 0.8 0.8];
showNBNnodeSize= 20;


startColor= [1 1 0];
endcolor = [0 0 0 ];
%	out << it.m_id << "\t" << it.m_visited_times << "\t" << it.m_start << "\t" << it.m_end << "\t" << it.fitness << "\t" << it.m_belongAlg << std::endl;
color1 = [0.8500 0.3250 0.0980];
color2 = [0.4660 0.6740 0.1880];
color3= [0.3010 0.7450 0.9330];

showNBNinfoIdx=3;

for idx=1:numNode
    show_network(idx,:)= network_mat(showNodeIdxs(idx),:);
    show_node_size(idx,:)= showNBNnodeSize;
    if nbn_info(showNodeIdxs(idx),showNBNinfoIdx)==-1
        show_node_color(idx,:)= gray_color2;
    else 
        colorId= MshowIterId(nbn_info(showNodeIdxs(idx),showNBNinfoIdx));
        show_node_color(idx,:)= mapColor(colorId,:);
    end    
    
    
%         
%    if nbn_info(showNodeIdxs(idx),showNBNinfoIdx)<0 
%        show_node_color(idx,:) = 	gray_color2;
%    elseif nbn_info(showNodeIdxs(idx),showNBNinfoIdx) == 0
%        show_node_color(idx,:) = 	color1;
%    elseif nbn_info(showNodeIdxs(idx),showNBNinfoIdx) == 1
%        show_node_color(idx,:) = 	color2;
%    elseif nbn_info(showNodeIdxs(idx),showNBNinfoIdx) == 2    
%        show_node_color(idx,:) = 	color3;
%    end
end





scatter3(show_network(:,2),show_network(:,3),show_network(:,4),show_node_size,show_node_color,'filled');
hold on;

colormap(mapColor);
c = colorbar;
c.Label.String = 'Iteration';
c.Label.FontSize = fontSize;
c.Ticks = [];
c.Position =[0.89 0.1 0.1 0.9];
%c.Layout.Tile = 'east';
hold on;
    outputfilename =strcat(outputdir,filename);
    outputfilename=strcat(outputfilename,"_fitness_nbn");   

           fh= gcf;
   view([-37.5000 30]);
   
   
algId=nbn_info(1,4);

algname ='';
    if algId==-0
        algname='ILS';
    elseif algId == 1
        algname='DE';
    else
        algname='PSO';
    end   
   
   
x = [0.75 0.77];
y = [0.82 0.82];
annotation('textarrow',x,y,'String',algname, 'HeadStyle',"vback3",'HeadWidth',0.1,'Headlength',0.1,'Color', [0 0 0], 'FontSize', fontSize);
% x = [0.9 0.94];
% y = [0.88 0.88];
% annotation('textarrow',x,y,'String','DE', 'HeadStyle',"vback3",'HeadWidth',0.1,'Headlength',0.1,'Color', color2, 'FontSize', fontSize);
% 
% x = [0.9 0.94];
% y = [0.84 0.84];
% annotation('textarrow',x,y,'String','PSO', 'HeadStyle',"vback3",'HeadWidth',1.0,'Headlength',0.1,'Color', color3, 'FontSize', fontSize);
% 
%    

set(gca,'XTick',[],'YTick',[],'ZTick',[]);

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
