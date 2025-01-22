
readdir= 'H:\code\ofec_data\idee_data_filter/';
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/idee_figure_eax_lkh_trait_4/';
outputdir ='\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/idee_figure_eax_lkh_trait_withoutLKH/';
readDriTrait= "\\172.24.242.8\share\Student\2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_4_traitdata/";

readdir = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\IDEE_lon_eax_data\idee_data_filter\";
readDriTrait = "\\172.24.24.151\f\DiaoYiya\paper_com\experiments_data\IDEE_lon_eax_data\eax_lkh_nbnData_4_traitdata\";


mkdir(outputdir);
str = [ "u574","5955" ,"2202" "6702" "6717" "1281" ];
%str = ["5955"];
str = [ "u574","5955" ,"2202"  ];
strSize = size(str);
strSize= strSize(1,2);
formatSpec="%f";

basicNodeSize = 8;
numColor = 1000;
mapColor = bone(numColor);
solMapColor = jet(numColor);
sizeA= [1 2]; 

for taskId= 1:strSize
clf;
hFig = figure('Visible', 'off'); 
tspname = str(1,taskId);
nodeFitFileName= readDriTrait + tspname + "_nodeFit.txt";
fileID = fopen(nodeFitFileName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
%        sizeM= sizeM(1,2);
nodeFit= fscanf(fileID,formatSpec,sizeM);
black_color= [0 0 0];
%sizeM= sizeM(1,2);
%bestPos= zeros(sizeM,5);

ve_mat = readmatrix(readdir+ tspname + "_IDEE_embedded_ve.txt");

maxDim = max(ve_mat(:,2));
[X,Y] = meshgrid(0:1:maxDim);
maxDim= maxDim+1;
sizeVe= size(ve_mat,1);
Z= zeros(maxDim, maxDim);
id2pos= zeros(sizeVe,2);
for idx=1:sizeVe
    Z(ve_mat(idx,2)+1,ve_mat(idx,3)+1)= nodeFit(idx);
    id2pos(idx,1)= ve_mat(idx,2);
    id2pos(idx,2)= ve_mat(idx,3);
end
clf;
hold on;
colormap jet;
ideemap = surf(X,Y,Z,'EdgeColor','interp');





bestIdName= readDriTrait + tspname + "_bestIds.txt";
sizeA= [1 2]; 
fileID = fopen(bestIdName,'r');
sizeM= fscanf(fileID,formatSpec,sizeA);
bestIds= fscanf(fileID,formatSpec,sizeM);
bestIds=bestIds+1;
sizeM= sizeM(1,2);
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
best_points.MarkerFaceColor="none";
best_points.Marker='square';
best_points.SizeData=150;
best_points.LineWidth=6;
% 
% best_points=  scatter3(bestPos(:,1),bestPos(:,2),bestZ,'filled');
% best_points.MarkerFaceColor= 'w';
% best_points.Marker= 'square';
% best_points.SizeData=500;
% %best_points.MarkerFaceAlpha= 0.5;

hold on;

traitname = ["lkh_fail", "lkh_success", "eax_fail", "eax_success"  ];
traitFileName= readDriTrait+ tspname+ "_traits.txt";
fileID = fopen(traitFileName,'r');





  for tId=1:4
       curTraitName= traitname(1, tId);
      % disp(curTraitName);
       sizeM= fscanf(fileID,formatSpec,sizeA);
       if (sizeM(1,2)>0)
           data= fscanf(fileID,formatSpec,sizeM);
           sizeM= sizeM(1,2);
           data_idx = data(1,:);
           data_idx = data_idx+1;
           data_value= 1-data(2,:);
           %data_value=numColor- int64(data_value.*numColor);
           data_value=numColor- int64(data_value.*numColor)+1;
           %data_value=data_value +1;
           
           solPos= zeros(sizeM,2);
           solZ= zeros(sizeM,1);
           
           for idx=1:sizeM
               solPos(idx,:)=  id2pos(data_idx(idx),:);
               solZ(idx,1)  = nodeFit(data_idx(idx));
           end
           
           
         %  solPos(:,4)=solPos(:,4)+0.02; 
           
           solColor = zeros(sizeM,3);
           for idx=1:sizeM
                solColor(idx,:) = solMapColor(data_value(idx),:);
                 solColor(idx,:) = black_color;
           end
           
           solSize = zeros(sizeM,1);
           for idx=1:sizeM
                solSize(idx)=100;
           end
           
           hold on;



           trait_points=  scatter3(solPos(:,1),solPos(:,2),solZ, "filled");
trait_points.MarkerEdgeColor= 'k';
trait_points.MarkerFaceColor="none";
trait_points.Marker='o';
trait_points.SizeData=20;
trait_points.LineWidth=1;
           
           curoutfilename =  outputdir + tspname + "_trait_"+ curTraitName;
           display(curoutfilename);

 %         view(3);
%setExportFigureType(curoutfilename,'origin',0);
view(2);
setExportFigureType(curoutfilename,'to',0);
%view([180 0]);
%setExportFigureType(curoutfilename,'view180',0);
           %solPoints.visible='false';
  
           
       end
    end



end