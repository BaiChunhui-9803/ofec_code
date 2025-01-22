    fileDir = 'C:/Data/student/2018/diaoyiya/nbn_experiment/corProNBN_df_cp/';
    fileDir = 'C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_df_iforest/'
dirOfDirFile= dir(fileDir);
numOfDirFiles = length(dirOfDirFile);

filenames_mat = readlines('C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBNdf_filenames.txt');
fileid = 1;
%Y = load('data.txt')
%C:\Data\student\2018\diaoyiya\nbn_experiment
outputDir= 'C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_opt_DisFit/';
mkdir(outputDir);
%filename = 'ProName_Classic_Michalewicz_Dim_2_opt_dfcolorpoints.txt';
filename = filenames_mat(fileid,1);
DirconProNBN_df_iforest='C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_df_iforest/';
Dir_conProNBN='C:/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_total3/';
%Dir_conProNBN= 'Z:\student\2018\diaoyiya\nbn_experiment\conProNBN_total3';
numColor= 1e3;
map = jet(numColor);



numcolor2d = 1000;
ramp = linspace(-100,100, numcolor2d);
%figure;
cform = makecform('lab2srgb');
a = repmat(ramp, [numcolor2d 1]);           % -a on left
b = repmat(flipud(ramp'), [1 numcolor2d]);  % -b on bottom
L = 50 * ones(numcolor2d, numcolor2d);  % A single L value.
Lab = cat(3, L, a, b); % A 2D image.
colormap2D = applycform(Lab, cform);


        disp(filename);           
        filepath=[outputDir,filename,'_opt_nbn2d.png'];
        readfilepath = DirconProNBN_df_iforest+ filename+ '_dfcolorpoints.txt';
        data = readmatrix(readfilepath);
        colorval = data(:,3);
        x = data(:,1)';
        y=data(:,2)';
        colorval3 = ValueToColor(colorval,1e3);
        scatter(y,x,10,colorval3);
        outputfilepath = outputDir + filename +'_df_density.png';
        exportgraphics(gcf,outputfilepath);
       
        readfilepath = DirconProNBN_df_iforest+ filename+ '_df_iforest.txt';
        df_iforest_mat=  readmatrix(readfilepath);
        x= df_iforest_mat(:,2);
        y= df_iforest_mat(:,3);
        colorval= df_iforest_mat(:,4);

        colorval3= zeros(sz(1,1),3);
        for id = 1:sz(1,1)
            coloridx= int32(colorval(id)*numColor)+1;
            if coloridx>numColor
                coloridx= numColor;
            end
            colorval3(id,:)= map(coloridx,:);
        end
        %colorval3 = cat(3, colorval, colorval,colorval); % A 2D image.
        %colorval3 = ValueToColor(colorval,1e3);
        scatter(y,x,10,colorval3);
        outputfilepath = outputDir + filename +'_df_iforest.png';
        exportgraphics(gcf,outputfilepath);
        
        
        histogram(colorval);
        outputfilepath = outputDir + filename +'_iforest.png';
        exportgraphics(gcf,outputfilepath);
        
        
        for id = 1:sz(1,1)
            coloridx= int32(colorval(id)*numColor)+1;
            if coloridx>numColor
                coloridx= numColor;
            end
            colorval3(id,:)= map(coloridx,:);
        end
        
        df_y= x;
        df_x = y;
        
        
        df_norx= rescale(df_x,0,1);
df_nory = rescale(df_y,0,1);
df_idx = zeros(sz(1,1),1);
df_idy= zeros(sz(1,1),1);
%df_idx= uint8(df_norx.*1000+1);
%df_idy = uint8(df_nory.*1000+1);

df_color = zeros(sz(1,1),3);
for idx=1:sz(1,1)
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
        scatter(df_x,df_y,1,df_color);
         outputfilepath = outputDir + filename +'_df_dfcolor.png';
        exportgraphics(gcf,outputfilepath);
                readfilepath = Dir_conProNBN+ filename+ '_nbnFit.txt';
        nbn_fitness= readmatrix(readfilepath);
        size_nbn = size(nbn_fitness);
        nbn_fitness(size_nbn(1,1),:)=[];
        
        x= nbn_fitness(2,:)';
        y=nbn_fitness(3,:)';
        
        scatter(nbn_fitness(:,2),nbn_fitness(:,3),10, colorval3);
                 outputfilepath = outputDir + filename +'_mesh_dfcolor.png';
        exportgraphics(gcf,outputfilepath);
        
        
        %A(2,:) = [];
%colormap(gca,"jet")
        