% test nbn matrix
%readfilepath="Z:/student/2018/diaoyiya/nbn_experiment/ProName_MMOP_CEC2013_F11_Dim_5_SampleSize_2000006_Opt_dfmat.txt";
%data = readmatrix(readfilepath);
%imagesc(data);

%close all
% [x,y,z]=peaks(30);
% contour3(x,y,z);
% title('山峰函数等值线图');
% xlabel('x-axis'),ylabel('y-axis '),zlabel('z-axis');
%ProName_Classic_five_uneven_peak_trap_Dim_2_SampleSize_2000002_Opt_dfcolorpoints
filename = "ProName_Classic_five_uneven_peak_trap_Dim_2";
filename = "ProName_Classic_uneven_de_maxima_Dim_2";
filename ="ProName_MMOP_CEC2013_F10_Dim_2";
filename = "ProName_BBOB_F01_Dim_2";
readfilepath = "C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_total/ProName_Classic_uneven_de_maxima_Dim_2_mesh.txt";
%readfilepath = savepath+filename+"_mesh.txt";
%savepath = "C:\Data\student/2018/diaoyiya/nbn_experiment/conProNBN_total/"
savepath = "C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_total3/";
%savepath = "Z:/student/2018/diaoyiya/nbn_experiment/conProNBN_total3/";
%readfilepath = "Z:/student/2018/diaoyiya/nbn_experiment/ProName_BBOB_F01_Dim_2_mesh.txt";
meshData=  readmatrix(savepath+ filename+"_mesh.txt");
x = meshData(:,3)';
y=meshData(:,4)';
z= meshData(:,5)';
sampleSize = size(x);
div = sqrt(sampleSize(1,2));
xx= reshape(x,div,div);
yy = reshape(y,div,div);
zz = reshape(z,div,div);
f1 = figure(1);
contour(xx,yy,zz,300);
%ProName_Classic_equal_maxima_expanded_Dim_2_SampleSize_2000025_Opt_mesh3d
%readfilepath = "Z:/student/2018/diaoyiya/nbn_experiment/conProNBN_total/ProName_BBOB_F02_Dim_2_SampleSize_2000000_mesh3d.txt";
readfilepath = savepath+filename+"_nbnFit.txt";
nbn3d_mat =readmatrix(readfilepath);

nbn3d_mat_x= nbn3d_mat(:,2)';
nbn3d_mat_y= nbn3d_mat(:,3)';
nbn3d_mat_z= nbn3d_mat(:,4)';
%ProName_Classic_equal_maxima_expanded_Dim_2_SampleSize_2000025_Opt_nbn
%readfilepath = "Z:/student/2018/diaoyiya/nbn_experiment/conProNBN_total/ProName_BBOB_F02_Dim_2_SampleSize_2000000_nbn.txt";
readfilepath = savepath+filename+"_nbn.txt";
%readfilepath="ProName_Classic_equal_maxima_expanded_Dim_2_SampleSize_2000000_nbn"

nbn_network = readmatrix(readfilepath,'NumHeaderLines',1);
from_edge = nbn_network(:,1)'+1;
to_edge = nbn_network(:,3)'+1;

nbn3d_graph = graph(from_edge,to_edge);


numColors = 1000;
ramp = linspace(-100,100, numColors);
%figure;
cform = makecform('lab2srgb');
a = repmat(ramp, [numColors 1]);           % -a on left
b = repmat(flipud(ramp'), [1 numColors]);  % -b on bottom
L = 50 * ones(numColors, numColors);  % A single L value.
Lab = cat(3, L, a, b); % A 2D image.
colormap2D = applycform(Lab, cform);
sz_nbn= size(nbn_network);
maxNBN = max(nbn_network);
minNBN = min(nbn_network);

df_x= nbn_network(:,4);
df_y = nbn_network(:,5);
df_norx= rescale(df_x,0,1);
df_nory = rescale(df_y,0,1);
df_idx = zeros(sz_nbn(1,1),1);
df_idy= zeros(sz_nbn(1,1),1);
%df_idx= uint8(df_norx.*1000+1);
%df_idy = uint8(df_nory.*1000+1);

df_color = zeros(sz_nbn(1,1),3);
for idx=1:sz_nbn(1,1)
    df_idx(idx) = int16((df_norx(idx))*1000);
    if df_idx(idx)==0
        df_idx(idx)=1;
    end
    df_idy(idx) = int16((1.0-df_nory(idx))*1000);
    if df_idy(idx)==0
        df_idy(idx)=1;
    end
    df_color(idx,:) = colormap2D(df_idy(idx),df_idx(idx), :);
end



f3 = figure(3);
scatter(df_x,df_y,1,df_color);
colormap(gca,"jet")

f4= figure(4);

c = [0 1 0; 1 0 0; 0.5 0.5 0.5; 0.6 0 1];
scatter(nbn3d_mat_x,nbn3d_mat_y,10,df_color);
colormap(gca,"jet")


%f5= figure(5);

% = [0 1 0; 1 0 0; 0.5 0.5 0.5; 0.6 0 1];
%scatter(nbn3d_mat_x,nbn3d_mat_y,[],df_color);
%colormap(gca,"jet")


%f6 = figure(6);
%surf(xx,yy,zz,df_color);
%view(2);

%scatter3()




%graph_plot = plot(nbn3d_graph,'XData',nbn3d_mat_x,'YData',nbn3d_mat_y,'ZData',nbn3d_mat_z);
%for i = 1: 
%end




%readfilepath= "Z:/student/2018/diaoyiya/nbn_experiment/ProName_MMOP_CEC2013_F12_Dim_10_SampleSize_2000000_dfcolorpoints.txt";
%data = readmatrix(readfilepath);
%colorval = data(:,3);
%x = data(:,1)';
%y=data(:,2)';
%colorval3 = ValueToColor(colorval,1e3);
%scatter(y,x,10,colorval3);