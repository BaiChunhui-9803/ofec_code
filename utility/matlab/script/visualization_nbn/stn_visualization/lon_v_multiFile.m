%fileID = fopen('//172.29.71.217/e/DiaoYiya/experiment_data/visualization_merge_stn/proName_ProName_BBOB_F09_FileName__Dim_3_evalSize_1500000.txt','rb');

filedir = '//172.29.71.217/e/DiaoYiya/experiment_data/visualization_lon_filtered_data/';
savedir=  'E:/nbn_data/visualization_lon_save_image4/';
mkdir(savedir);

Files=dir(filedir);
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
G = digraph(s,t,weight,nodeNum);


h = plot(G,'Layout','force');

%h.EdgeColor= colorArr3;
h.LineWidth= 1.5000;
%h.LineWidth= LWidths;
%h.ArrowSize=4;
%pictureName = filename + '.jpg';
saveFilePath =[savedir,filename, '.jpg'];
exportgraphics(gcf,saveFilePath);

   end

end
