topDir ="D:";
readdir = topDir+ "/student/2018/DiaoYiya/tspLonGenerator/";
%savedir = dir+"";
%filename= "att532";
%filedirpah = filename+"/"+ filename;
saveDir = "E:/DiaoYiya/experiment_figure/TSP_figure/LON5/";
mkdir(saveDir);
tsp_filenames= dir(readdir);
dir_size=  size(tsp_filenames,1);
dir_size = dir_size*2;
coreNum = 40;
force_layout_iteration= 100;
lon_filename = ["_lon_filter001" "_lon_filter005"];
par = parpool('local', coreNum, 'IdleTimeout', Inf); %开启并行计算
parfor dirId = 1:dir_size
    curDirId= uint32(floor((dirId-1)/2))+1;
    curFileId = mod(dirId,2)+1;
    filename= tsp_filenames(curDirId).name;
    lonname=  lon_filename(curFileId);
    filedirpath = filename+"/"+ filename;
    filepath =  readdir+ filedirpath +lonname +".txt";
    saveFilename =filename+ lonname;
    if(tsp_filenames(curDirId).isdir&&(filename~="."&&filename~=".."))
        figureId= mod(dirId,coreNum)+1;
        display(filepath);
        try
        ShowLON(figureId, filepath, saveDir+saveFilename)
        catch ME
            display(filepath+ ME.message);
        end
       % display(saveFilename);
     %   figure(figureId+1)
    %    plot(rand(dirId*10,1));
      %  saveas(gcf,['temp' num2str(dirId) '.jpg']);
    end
end
delete(par);
