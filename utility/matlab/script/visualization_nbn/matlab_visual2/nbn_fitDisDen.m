
filename = "ProName_MMOP_CEC2013_F12_Dim_2";
filename = "ProName_Classic_Ackley_Dim_2";
filedir= "Z:/student/2018/diaoyiya/nbn_experiment/";
filedir="//172.24.242.8/share/Student/2018/YiyaDiao/nbn_experiment/";
filedir ="//172.24.24.151/d/student/2018/DiaoYiya/";
filename_opt= filename+"_opt";

nbn_mat  = readmatrix(filedir+"conProNBN_total4/"+filename_opt+"_FitDisRadiusFilterByThreahold.txt");


df_x= nbn_mat(:,3);
df_y = nbn_mat(:,4);
df_color = nbn_mat(:,5);
filter_flag = nbn_mat(:,6);
df_color = df_color.*2;
numColor = 1e3;
map = jet(numColor);
sz_xy  = size(df_x);
numPoints= sz_xy(1,1);
colorIdx = zeros(numPoints,1);
numFilterPoitnts = 0;
for idx= 1:sz_xy(1,1)
    colorIdx(idx) = uint32(df_color(idx)*numColor);
    if(colorIdx(idx)==0)
        colorIdx(idx)=1;
    end
   if(filter_flag(idx)==1)
        numFilterPoitnts=numFilterPoitnts+1;
    end
end

colorVal = zeros(sz_xy(1,1),3);
for idx=1:sz_xy(1,1)
    colorVal(idx,:) = map(colorIdx(idx),:);
end
f1= figure(1);
clf(f1);
scatter(df_y,df_x,10,colorVal);


filter_x= zeros(numFilterPoitnts,1);
filter_y = zeros(numFilterPoitnts,1);
filter_colorVal= zeros(numFilterPoitnts,3);
pidx=0;

fitPos= readmatrix(filedir+"conProNBN_total4/"+filename_opt+"_FitDisRadiusFitPos.txt");
filter_fitpos= zeros(numFilterPoitnts,4);

dd_x= zeros(numFilterPoitnts,1);
dd_y = zeros(numFilterPoitnts,1);
for idx= 1:sz_xy(1,1)
    if(filter_flag(idx)==1)
        pidx=pidx+1;
        filter_x(pidx) = df_x(idx);
        filter_y(pidx) = df_y(idx);
        colorVal(idx,:) =[1,0,0]; 
       % dd_x(pidx,)
        filter_colorVal(pidx,:) = colorVal(idx,:);
        filter_fitpos(pidx,:)= fitPos(idx,:);

    end
end
f2= figure(2);
clf(f2);

scatter(df_y,df_x,1,colorVal);
hold on;
scatter(filter_y,filter_x,10,filter_colorVal);
%scatter(df_y,df_x,1,colorVal);


f3 = figure(3);
clf(f3);
%meshData=  readmatrix(savepath+ filename+"_mesh.txt");
meshData = readmatrix(filedir+"conProNBN_total4/"+filename+"_mesh.txt");
x = meshData(:,3)';
y=meshData(:,4)';
z= meshData(:,5)';
sampleSize = size(x);
div = sqrt(sampleSize(1,2));
xx= reshape(x,div,div);
yy = reshape(y,div,div);
zz = reshape(z,div,div);
f4 = figure(4);
contour3(xx,yy,zz,50);
hold on;
scatter3(filter_fitpos(:,2)',filter_fitpos(:,3)',filter_fitpos(:,4)',20,filter_colorVal);

