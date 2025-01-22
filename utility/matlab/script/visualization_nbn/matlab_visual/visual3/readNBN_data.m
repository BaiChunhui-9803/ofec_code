%fileDir= '/home/data/student/2018/diaoyiya/code/NBN_data/NBN_visual_data_2022-11-27-12-36-00';
fileDir= 'E:/Diao_Yiya/code/NBN_data/TEST_data';
curTime = datetime("now");
lastTime = datetime("now")-years(1);
fmt = "yyyy-MM-dd--hh-mm-ss";
strCurTime = string(curTime,fmt);
%outputDir = append('/home/data/student/2018/diaoyiya/code/NBN_visual_data/','NBN_visualPic_',strCurTime);
outputDir = append('E:/Diao_Yiya/code/NBN_visual_data/test2/','NBN_visualPic_',strCurTime);
%mkdir(outputDir);
mkdir(outputDir);
%mkdir(outputDir);
%fileDir = '//172.24.24.151/share/student/2018/diaoyiya/code/NBN_visual_data1';
format('longEng');
formatSpec = '%e';
%filename = dirOfFile(4).name;
%filename ='problem_name-COP_CEC2017_F02--number_of_variables-2--algorithm_name-NBN_DSA--add_optimal-1--sample_size-10002.txt';
set(0, 'DefaultFigureVisible', 'off');

n=10;
%while n > 1
    dirOfFile= dir(fileDir);
    numOfFiles = length(dirOfFile);
    info = zeros(numOfFiles);
    ErrStr = strings(numOfFiles);
    curFinishTime = datetime("now")+ minutes(10);
    finishTime = curFinishTime;
    
    
    
    
parfor iter_file= 1:numOfFiles
    filename = dirOfFile(iter_file).name;
    filetime = dirOfFile(iter_file).date;
    if filetime<=curFinishTime&&filetime>lastTime
        if (strcmp(filename,'.'))|| (strcmp(filename,'..'))
            disp(filename);
        else
            filepath = append(fileDir,'/',filename);
            disp(filename);
            disp(filepath);       
            figureId= iter_file;                  
            Psize = 20;
            info= showNBNdata(figureId,outputDir,fileDir,filename,Psize);
            ErrStr(iter_file)  = info.ErrInfo;
            disp(info.ErrInfo);
        end 
    end
 end
lastTime= finishTime;
%end 





