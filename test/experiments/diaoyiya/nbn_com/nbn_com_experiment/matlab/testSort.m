% 示例数据
data = [5 3; 2 8; 9 1; 4 6];

% 根据第二列的数值进行排序
sortedData = sortrows(data, 2);

% 显示排序后的数据
disp('Sorted data:');
disp(sortedData);
