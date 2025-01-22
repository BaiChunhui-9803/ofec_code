fileDir = '//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_txt';
fileDir = 'Z:/student/2018/diaoyiya/nbn_experiment/conProNBN_txt/';
%dirOfDirFile= dir(fileDir);
%numOfDirFiles = length(dirOfDirFile);
picOutputDir= '//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_UAD_picture/';
mkdir(picOutputDir);
dataOutputDir=  '//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_UAD_data/';
mkdir(dataOutputDir);
coreNum =1; %设置处理器数目

%par = parpool('local', coreNum); %开启并行计算
filenames_mat = readlines('Z:/student/2018/diaoyiya/nbn_experiment/conProNBNdf_filenames.txt');
numOfDirFiles = size(filenames_mat);
numOfDirFiles = numOfDirFiles(1,1);
%numOfDirFiles=1;
DirconProNBN_df_iforest = "//172.24.24.151/Data/student/2018/diaoyiya/nbn_experiment/conProNBN_df_iforest/";

numColor= 1e3;
map = jet(numColor);

for iter_file= 1:numOfDirFiles
    filename = filenames_mat(iter_file,1);
    disp(filename);


       
        readfilepath = DirconProNBN_df_iforest+ filename+ '_df_iforest.txt';
        df_iforest_mat=  readmatrix(readfilepath);
        x= df_iforest_mat(:,2);
        y= df_iforest_mat(:,3);


        method_name='ocsvm';
      % calUADfunciton(x,y, map, numColor, picOutputDir,dataOutputDir, filename, method_name);



        sz = size(df_iforest_mat);
        data_mat = zeros(sz(1,1),2);
        data_mat(:,1)= x;
        data_mat(:,2)= y;
        %colorval= df_iforest_mat(:,4);

        [Mdl,tf,colorval] = ocsvm(data_mat);
        %histogram(scores);
        %outputfilepath = outputDir + a +'iforest.png';
        %exportgraphics(gcf,outputfilepath);

        norScore = rescale(colorval, 0,1);
        colorval3= zeros(sz(1,1),3);
        for id = 1:sz(1,1)
            coloridx= int32(norScore(id)*numColor)+1;
            if coloridx>numColor
                coloridx= numColor;
            end
            colorval3(id,:)= map(coloridx,:);
        end


        

        scatter(y,x,10,colorval3);
        outputfilepath = picOutputDir + filename +'_df_ocsvm.png';
        exportgraphics(gcf,outputfilepath);

        outputfilepath = dataOutputDir + filename +'_df_ocsvm.txt';
        writematrix(colorval,outputfilepath);


end
%delete(par) %关闭并行计算
