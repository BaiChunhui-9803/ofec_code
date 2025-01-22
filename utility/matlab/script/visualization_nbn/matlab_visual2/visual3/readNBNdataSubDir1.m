%fileDir= '/home/data/student/2018/diaoyiya/code/NBN_data/nbn_total_data';
fileDir = '/home/data/student/2018/diaoyiya/code/NBN_data/conData21';
%fileDir= 'E:/Diao_Yiya/code/NBN_data/testData';
%\\172.24.24.151\share\student\2018\diaoyiya\code\NBN_data\conData
fileDir = '//172.24.24.151/share/student/2018/diaoyiya/code/NBN_data/conData';
%outputDir = append('E:/Diao_Yiya/code/NBN_visual_data/test2/','NBN_visualPic_',strCurTime);
%fileDir= 'E:/Diao_Yiya/code/NBN_data/info';
curTime = datetime("now");
fmt = "yyyy-MM-dd--hh-mm-ss";
strCurTime = "continousProblem";
%outputDir = append('/home/data/student/2018/diaoyiya/code/NBN_visual_data/','NBN_visualPic_',strCurTime);
outputDir = append('E:/Diao_Yiya/code/NBN_visual_data/result/','NBN_visualPic_',strCurTime);
%mkdir(outputDir);
mkdir(outputDir);
%fileDir = '//172.24.24.151/share/student/2018/diaoyiya/code/NBN_visual_data1';
dirOfDirFile= dir(fileDir);
numOfDirFiles = length(dirOfDirFile);
format('longEng');
formatSpec = '%e';
set(0, 'DefaultFigureVisible', 'off');
%filename = dirOfFile(4).name;
%filename ='problem_name-COP_CEC2017_F02--number_of_variables-2--algorithm_name-NBN_DSA--add_optimal-1--sample_size-10002.txt';
    info = zeros(numOfDirFiles);
    ErrStr = strings(numOfDirFiles,1);

    numTask  = 3;
    taskId = 1;

for iter_file= 1:numOfDirFiles
    b = mod(iter_file,numTask);
    if(b~=taskId) 
        continue;
    end
    filename = dirOfDirFile(iter_file).name;
    if (strcmp(filename,'.'))|| (strcmp(filename,'..'))
        disp(filename);
    else
        filepath = append(fileDir,'/',filename);
        disp(filename);
        disp(filepath);       
    
    figureId= 10*(taskId+1);               
    Psize = 20;
 %   info= showNBNdata(figureId,outputDir,fileDir,filename,Psize);
 %   ErrStr(iterFile,1)  = info.ErrInfo;
 %   disp(info.ErrInfo);

    end
end




writematrix(ErrStr,'result.xls');