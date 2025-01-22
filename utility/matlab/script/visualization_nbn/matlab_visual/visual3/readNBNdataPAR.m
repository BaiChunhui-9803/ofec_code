%fileDir= '/home/data/student/2018/diaoyiya/code/NBN_data/nbn_total_data';
fileDir= 'E:/Diao_Yiya/code/NBN_data/testData';
fileDir = '//172.24.24.151/share/student/2018/diaoyiya/code/NBN_data/nbn_total_data'
%outputDir = append('E:/Diao_Yiya/code/NBN_visual_data/test2/','NBN_visualPic_',strCurTime);
%fileDir= 'E:/Diao_Yiya/code/NBN_data/info';
curTime = datetime("now");
fmt = "yyyy-MM-dd--hh-mm-ss";
strCurTime = string(curTime,fmt);
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

fileInfos = strings();
fileDirInfos = strings();
filenamesInfo = strings();


for iter_file= 1:numOfDirFiles
    filename = dirOfDirFile(iter_file).name;
    if (strcmp(filename,'.'))|| (strcmp(filename,'..'))
        disp(filename);
    else
        filepath = append(fileDir,'/',filename);
        disp(filename);
        disp(filepath);       
        dirOfFIle = dir(filepath);
        numFiles= length(dirOfFIle);
        nbnPathStr = strings(1,numFiles);
        nbnNameStr = strings(1,numFiles);
        nbnPathDirStr = strings(1,numFiles);
        for fileIter=1:numFiles
            nbnfilename = dirOfFIle(fileIter).name;
            if (strcmp(nbnfilename,'.'))|| (strcmp(nbnfilename,'..'))
            %    disp(filename);
                continue;
            end
            disp(nbnfilename);
            nbnfilepath = append(filepath,'/',nbnfilename);
            disp(nbnfilepath);
            nbnPathStr(1,fileIter)=nbnfilepath;
            nbnNameStr(1,fileIter) = nbnfilename;
            nbnPathDirStr(1,fileIter) = filepath;
        end
        
         fileInfos = [fileInfos nbnPathStr(1,:)];
         filenamesInfo = [filenamesInfo nbnNameStr(1,:)];
         fileDirInfos =[fileDirInfos nbnPathDirStr(1,:)];
    end
end


[tmp,numFile] = size(fileInfos);
    info = zeros(numFile);
    ErrStr = strings(numFile,1);
    
parfor iterFile = 2:numFile

    t = getCurrentTask;
    disp(t.ID);
    
    fileDir = fileDirInfos(iterFile);
    filename = filenamesInfo(iterFile);
    
    disp(fileDir);
    disp(filename);
    figureId= t.ID;               
    Psize = 20;
    info= showNBNdata(figureId,outputDir,fileDir,filename,Psize);
    ErrStr(iterFile,1)  = info.ErrInfo;
    disp(info.ErrInfo);



%    fh = figure(iterFile);
%    close(fh);
end


writematrix(ErrStr,'result.xls');