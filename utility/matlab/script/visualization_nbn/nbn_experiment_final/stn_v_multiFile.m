%fileID = fopen('//172.29.71.217/e/DiaoYiya/experiment_data/visualization_merge_stn/proName_ProName_BBOB_F09_FileName__Dim_3_evalSize_1500000.txt','rb');

filedir = 'E:/nbn_data/visualization_stn_merge2/';
savedir=  'E:/nbn_data/visualization_save_stn_merge_2/';
mkdir(savedir);
Files=dir(filedir);
%for k=1:3
for k=1:length(Files)
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
A = fscanf(fileID,formatSpec);
algName= fscanf(fileID,'%s', 3);
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
weight =  edgeInfo(3,:)';
G = digraph(s,t,weight,nodeNum);
edgeNum = size(edgeInfo);
edgeNum = edgeNum(1,2);
colorArr3 = zeros(edgeNum,3);
%lineWidthArr = zeros(edgeNum,1);
LWidths = max(1,5*weight/max(weight));

color1 = [0.8500 0.3250 0.0980];
color2 = [0.4660 0.6740 0.1880];
color3= [0.3010 0.7450 0.9330];
for idx = 1:edgeNum
   if edgeInfo(4,idx)<0 
       colorArr3(idx,:) = 	[0.5 0.5 0.5];
   elseif edgeInfo(4,idx) == 0
       colorArr3(idx,:) = 	color1;
   elseif edgeInfo(4,idx) == 1
       colorArr3(idx,:) = 	color2;
   elseif edgeInfo(4,idx) == 2    
       colorArr3(idx,:) = 	color3;
   end
   % lineWidthArr(idx)=1.5000;
   %     H(r,c) = 1/(r+c-1);
end

clf;
f = figure('visible','off');
clf;
h = plot(G,'Layout','force','NodeLabel',{});
h.EdgeColor= colorArr3;
h.LineWidth= 1.5000;
h.LineWidth= LWidths;
h.ArrowSize=4;
h.NodeColor = [0 0 0 ];
h.MarkerSize= 0.001;

nodeWeight = nodeInfo(2,:)';
nodeColorArr3 = zeros(nodeNum,3);
nodeWeight = max(10,25*nodeWeight/max(nodeWeight));
startColor= [1 1 0];
endcolor = [0 0 0 ];
%	out << it.m_id << "\t" << it.m_visited_times << "\t" << it.m_start << "\t" << it.m_end << "\t" << it.fitness << "\t" << it.m_belongAlg << std::endl;
	
for idx = 1:nodeNum
   if nodeInfo(4,idx)==1 
       nodeColorArr3(idx,:) = 	endcolor;
   elseif nodeInfo(3,idx) == 1
       nodeColorArr3(idx,:) = 	startColor;
   elseif nodeInfo(6,idx) == 0
       nodeColorArr3(idx,:) = 	color1;
   elseif nodeInfo(6,idx) == 1
       nodeColorArr3(idx,:) = 	color2;
   elseif nodeInfo(6,idx) == 2    
       nodeColorArr3(idx,:) = 	color3;
   else
      nodeColorArr3(idx,:) = 	[0.5 0.5 0.5];
   end
   % lineWidthArr(idx)=1.5000;
   %     H(r,c) = 1/(r+c-1);
end

hold on;
h3= scatter3( h.XData,h.YData, h.ZData,nodeWeight, nodeColorArr3,'filled');
fontSize=12;

x = [0.86 0.9];
y = [0.9 0.9];
annotation('textarrow',x,y,'String','ILS', 'HeadStyle',"vback3",'HeadWidth',0.1,'Headlength',0.1,'Color', color1, 'FontSize', fontSize);

x = [0.86 0.9];
y = [0.86 0.86];
annotation('textarrow',x,y,'String','DE', 'HeadStyle',"vback3",'HeadWidth',0.1,'Headlength',0.1,'Color', color2, 'FontSize', fontSize);


x = [0.86 0.9];
y = [0.82 0.82];
annotation('textarrow',x,y,'String','PSO', 'HeadStyle',"vback3",'HeadWidth',1.0,'Headlength',0.1,'Color', color3, 'FontSize', fontSize);
%pictureName = filename + '.jpg';
saveFilePath =[savedir,filename, '.png'];
exportgraphics(gcf,saveFilePath,'Resolution',600);
%   exportgraphics(gcf,filepath,'Resolution',600);
   end

end
