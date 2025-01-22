%fileID = fopen('//172.29.71.217/e/DiaoYiya/experiment_data/visualization_merge_stn/proName_ProName_BBOB_F09_FileName__Dim_3_evalSize_1500000.txt','rb');

filedir = '//172.29.71.217/e/DiaoYiya/experiment_data/visualization_lon_filtered_data/';
filedir = 'E:/nbn_data/visualization_lon_data9/';
%filedir = 'E:/nbn_data/visualization_save_lon2/';
filedir ='E:/nbn_data/visualization_lon_data10/';
savedir=  'E:/nbn_data/visualization_lon_save_image_228/';
mkdir(savedir);

Files=dir(filedir);
for k=1:length(Files)
%for k=1:3
   if(Files(k).isdir~=1)
     FileNames=Files(k).name;
     newStr = split(FileNames,'.');
     filename= newStr{1,1};
     display(filename);
     

filepath=[filedir,filename,'.txt'];
fileID = fopen(filepath,'r');
sizeNum=1;
format('longEng');
formatSpec = '%e';
matSize = fscanf(fileID,formatSpec,[2,1]);
matSize= matSize';
nodeInfo = fscanf(fileID,formatSpec,matSize);
matSize = fscanf(fileID,formatSpec,[2,1]);
matSize= matSize';
edgeInfo = fscanf(fileID,formatSpec,matSize);
nodeNum = size(nodeInfo);
nodeNum= nodeNum(1,2);
s = edgeInfo(1,:)';
t= edgeInfo(2,:)';
edgeNum = size(edgeInfo);
edgeNum = edgeNum(1,2);
for idx=1:edgeNum
    s(idx,1)= s(idx,1)+1;
    t(idx,1)= t(idx,1)+1;
end
weight =  edgeInfo(3,:)';
basicNodeSize = 25;
node_funnel= nodeInfo(4,:)';

maxNodeFunnel = max(node_funnel);
mapColor = autumn(maxNodeFunnel+1);
node_color = zeros(nodeNum,3);
gray_color = [0.5 0.5 0.5];
black_color = [0 0 0];
node_sizes = zeros(nodeNum,1);
for idx=1:nodeNum
    node_sizes(idx)  = basicNodeSize;
    if node_funnel(idx)==-1
        node_color(idx,:)= gray_color;
    else
        node_color(idx,:) = mapColor(node_funnel(idx)+1,:);
    end
end



%f = figure('visible','off');
clf;
G = digraph(s,t,weight,nodeNum);
h = plot(G,'Layout','force','NodeLabel',{},'NodeColor', gray_color , 'LineWidth',1.0 ,'EdgeColor',black_color);
h.ZData = nodeInfo(3,:)';
   view([-37.5000 30]);
h.LineWidth= 1.0000;
hold on;
h3= scatter3( h.XData,h.YData, h.ZData,node_sizes, node_color,'filled');

%h.EdgeColor= colorArr3;

%h.LineWidth= LWidths;
%h.ArrowSize=4;
%pictureName = filename + '.jpg';
%saveFilePath =[savedir,filename, '.jpg'];
%exportgraphics(gcf,saveFilePath);


saveFilePath =[savedir,filename];
   setExportFigureType(saveFilePath,'origin',0.10);


%exportgraphics(gcf,saveFilePath,'Resolution',600);
   end
end
