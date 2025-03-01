folderPath = '\\172.24.178.59\share\2018\diaoyiya\ofec-data\paper_con_realworld\AAA_algOneTrait\'; % 设置文件夹路径
outputdir = '\\172.24.6.26\f\DiaoYiya\paper_con_realworld\AlgOneTraitNBN2\'; % 设置输出目录

mkdir(outputdir);
% 获取文件夹下所有文件夹名
dirnames = dir(folderPath);
dirnames = {dirnames([dirnames.isdir]).name}; % 只保留文件夹名
dirnames(ismember(dirnames,{'.','..'})) = []; % 去除当前目录和上级目录

fid = fopen('error_info.txt','a'); % 打开或创建文件用于追加错误信息
%for idir = dirnames
for j = 1:length(dirnames)
    idir=dirnames(j);

    curfolderPath = fullfile(folderPath, idir{1})+"\";
    
    % 获取文件夹下所有文件名
    fileNames = dir(fullfile(curfolderPath, '*_network.txt'));
   curoutputdir = fullfile(outputdir, idir{1}) + "\";
   mkdir(curoutputdir);
    % 遍历文件名数组
    for i = 1:length(fileNames)
        % 获取完整文件路径
        filePath = fullfile(curfolderPath, fileNames(i).name);
      %  disp(filePath);

        curfilename = fileNames(i).name(1:end-length('_network.txt'));
    %    disp(curfilename);
        
        
        try
            drawNBN(curfolderPath, curoutputdir, curfilename);
        catch ME
            disp(idir{1});
            disp(fileNames(i).name);
            disp(ME.message);
            fprintf(fid,'%s\t',idir{1}); % 将错误信息写入文件
            fprintf(fid,'%s\t',fileNames(i).name); % 将错误信息写入文件
            fprintf(fid,'%s\n',ME.message); % 将错误信息写入文件

        end


       % drawNBN(curfolderPath, curoutputdir, curfilename);
    end
end

fclose(fid); % 关闭文件