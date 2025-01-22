%fileID = fopen('//172.29.71.217/e/DiaoYiya/experiment_data/visualization_merge_stn/proName_ProName_BBOB_F09_FileName__Dim_3_evalSize_1500000.txt','rb');

filedir = 'E:/nbn_data/visualization_stn_merge/';
savedir=  'E:/nbn_data/visualization_save_stn_merge/';



filename = 'proName_ProName_BBOB_F09_FileName__Dim_3_evalSize_1500000';
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
LWidths = max(0.01,5*weight/max(weight));
for idx = 1:edgeNum
   if edgeInfo(4,idx)<0 
       colorArr3(idx,:) = 	[0.5 0.5 0.5];
   elseif edgeInfo(4,idx) == 0
       colorArr3(idx,:) = 	[1 0 0.5];
   elseif edgeInfo(4,idx) == 1
       colorArr3(idx,:) = 	[0 0 1];
   elseif edgeInfo(4,idx) == 2    
       colorArr3(idx,:) = 	[0 1 0];
   end
   % lineWidthArr(idx)=1.5000;
   %     H(r,c) = 1/(r+c-1);
end

h = plot(G,'Layout','force');
h.EdgeColor= colorArr3;
h.LineWidth= 1.5000;
h.LineWidth= LWidths;
h.ArrowSize=4;
%pictureName = filename + '.jpg';
saveFilePath =[savedir,filename, '.jpg'];
exportgraphics(gcf,saveFilePath);