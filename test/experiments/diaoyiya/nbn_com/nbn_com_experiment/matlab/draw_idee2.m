
readdir= 'F:\code\ofec_data\idee_data_filter/';
readdir = "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\idee_image_data/";
outputdir ='\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/idee_figure_eax_lkh_trait/';
readDriTrait= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\total_data_idee/";
readdir = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\IDEE_total_data\idee_image_data\";
readDriTrait = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\IDEE_total_data\total_data_idee\";

mkdir(outputdir);
str = [ "u574","5955" ,"2202" "6702" "6717" "1281" ];
str = ["u574"];
%str = [ "u574","5955" ,"2202"  ];
strSize = size(str);
strSize= strSize(1,2);


traitname = ["_lkh_trait", "_eax_trait" ];
traitname= "_eax_trait";
taskId = 1;
tspname = str(1,taskId);
curname= tspname+ traitname+ "_0";
matlabTraitDataFilenName= readDriTrait+curname;
traitFileName= matlabTraitDataFilenName+"_matlab.txt";
curoutfilename = outputdir+ curname;
display(curoutfilename);

formatSpec="%f";
basicNodeSize = 8;
numColor = 1000;
mapColor = bone(numColor);
solMapColor = jet(numColor);
sizeA= [1 2]; 

%for taskId= 1:strSize
clf;
f = figure;  % 创建一个图形窗口

hFig = figure('Visible', 'off'); 

nodeFitFileName= readDriTrait + tspname + "_fitness_matlab.txt";
fileID = fopen(nodeFitFileName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
%        sizeM= sizeM(1,2);
nodeFit= fscanf(fileID,formatSpec,sizeM);
%nodeFit= -nodeFit;
black_color= [0 0 0];
%sizeM= sizeM(1,2);
%bestPos= zeros(sizeM,5);

ve_mat = readmatrix(readdir+ tspname + "IDEE_embedded_ve.txt");

maxDim = max(ve_mat(:,2));
[X,Y] = meshgrid(0:1:maxDim);
maxDim= maxDim+1;
sizeVe= size(ve_mat,1);
Z= zeros(maxDim, maxDim);
id2pos= zeros(sizeVe,2);
for idx=1:sizeVe
    Z(ve_mat(idx,3)+1,ve_mat(idx,2)+1)= nodeFit(idx);
    id2pos(idx,1)= ve_mat(idx,2);
    id2pos(idx,2)= ve_mat(idx,3);
end
clf;
hold on;
colormap jet;


x = linspace(1, 423, 1); % x轴的值
y = linspace(1, 423, 1); % y轴的值
ideemap = surf(X,Y,Z,'EdgeColor','none');
%imagesc(x, y, Z);

bestIdName= readDriTrait + tspname + "_bestSolIds.txt";
sizeA= [1 2]; 
fileID = fopen(bestIdName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
bestIds= fscanf(fileID,formatSpec,sizeM);
bestIds=bestIds+1;
bestIds(bestIds >= 190969) = [];
bestIdsSize = length(bestIds); % 获取数组 bestIds 的元素数量

sizeM= bestIdsSize;
bestPos= zeros(sizeM,5);


           bestPos= zeros(sizeM,2);
           bestZ= zeros(sizeM,1);
           
           for idx=1:sizeM
               bestPos(idx,:)=  id2pos(bestIds(idx),:);
               bestZ(idx,1)  = nodeFit(bestIds(idx));
           end
           


hold on;

best_points=  scatter3(bestPos(:,1),bestPos(:,2),bestZ,'filled');
best_points.MarkerEdgeColor= 'w';
best_points.MarkerFaceColor="w";
%best_points.Marker='o';
best_points.Marker='square';
best_points.SizeData=150;
best_points.LineWidth=6;
%opt_points.MarkerFaceColor= 'none';
%opt_points.MarkerEdgeColor= 'r';
%opt_points.Marker= 'o';
%opt_points.SizeData= 60;

%best_points.MarkerFaceAlpha= 0.5;




disp(traitFileName);
fileID = fopen(traitFileName,'r');
       sizeM= fscanf(fileID,formatSpec,sizeA);
           data= fscanf(fileID,formatSpec,sizeM);
           sizeM= sizeM(1,2);

    
           traitIds = data(1,:);
           if tspname== '5955'
           traitIds(traitIds >= 190969) = [];
           end
traitIdsSize = length(traitIds); % 获取数组 bestIds 的元素数量

sizeM= traitIdsSize;
           traitIds = traitIds+1;
           traitPos= zeros(sizeM,2);
           traitZ= zeros(sizeM,1);
           
           for idx=1:sizeM
               traitPos(idx,:)=  id2pos(traitIds(idx),:);
               traitZ(idx,1)  = nodeFit(traitIds(idx));
           end
           

trait_points=  scatter3(traitPos(:,1),traitPos(:,2),traitZ,'filled');
trait_points.MarkerEdgeColor= 'k';
trait_points.MarkerFaceColor="none";
trait_points.Marker='o';
trait_points.SizeData=20;
trait_points.LineWidth=1;

hold on;

min_value = max(bestZ);
indices = find(traitZ == min_value);



foundIdsSize = length(indices); % 获取数组 bestIds 的元素数量

sizeM= foundIdsSize


           foundPos= zeros(sizeM,2);
           foundZ= zeros(sizeM,1);
           
           for idx=1:sizeM
               foundPos(idx,:)=  traitPos(indices(idx),:);
               foundZ(idx,1)  = nodeFit(indices(idx));
           end
 if sizeM>0         

% 用黑色正方形画出来
h4 = scatter3(bestPos(:,1),bestPos(:,2),bestZ,'filled');
h4.MarkerEdgeColor= 'k';
h4.MarkerFaceColor='k';
h4.Marker='o';
h4.SizeData=2;
h4.LineWidth=5;

 end


hold on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%           legend('show');
% 设置图形的位置和大小
set(f, 'Position', [0 0 1 1]);  % 这将图形充满整个窗口

 %         view(3);
%setExportFigureType(curoutfilename,'origin',0);
view(2);
setExportFigureType(curoutfilename,'to',0);
%view([180 0]);
%setExportFigureType(curoutfilename,'view180',0);
           %solPoints.visible='false';

%end