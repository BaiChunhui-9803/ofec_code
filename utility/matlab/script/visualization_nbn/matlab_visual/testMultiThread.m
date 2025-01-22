numOfDirFiles= 160;
coreNum = 6;
tic
delete(par) %关闭并行计算
par = parpool('local', coreNum); %开启并行计算
parfor iter_file= 1:numOfDirFiles
    for idx=1:1e7
        a = mod(idx,3);    
        b= rand(3,3);
    end
end
delete(par) %关闭并行计算
toc