filenameDir= "//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/experiments_added/tsp_instances/";
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



