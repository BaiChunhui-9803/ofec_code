function [tspnames] = getTspTaskNames(filenameDir)

traitFileName= filenameDir+"tspNames.txt";
formatSpec = "%s";
sizeA= [1 2]; 
%fileID = fopen(traitFileName);
%fileID = fopen('filename.txt', 'r');

% 使用 fileread 读取整个文件内容
fileContent = fileread(traitFileName);
% 关闭文件
%fclose(fileID);

% 将文件内容分割成行
lines = strsplit(fileContent, '\n');
lines(end) = [];

% 逐行处理字符串
for i = 1:length(lines)
    disp(['Line ' num2str(i) ': ' lines{i}]);
end


tspnames= lines;

    % 函数体：包含函数的代码
end