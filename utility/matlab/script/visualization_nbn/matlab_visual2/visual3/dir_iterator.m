

for iter_file= 1:numOfFiles
    filename = dirOfFile(iter_file).name;
    if (strcmp(filename,'.'))|| (strcmp(filename,'..'))
        disp(filename);
    else
        filepath = append(fileDir,'/',filename);
        disp(filename);
        disp(filepath);       
        fileID = fopen(filepath,'r');
        tline = fgetl(fileID);
        
        disp(tline);
        dataName = fgetl(fileID);
        disp(dataName);
    end
end 


