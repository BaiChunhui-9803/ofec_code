readdir=  "E:/Diao_Yiya/paper/nbn_experiment/IDEE_data/";
tspname = ["BCL380",'XQL662','PCB442','GR666'];
for idx= 1:size(tspname,2)
    readdir="E:/Diao_Yiya/paper/nbn_experiment/IDEE_data2/filterdata/";
savedir = "E:/Diao_Yiya/paper/nbn_experiment/IDEE_data2/figurei/";
drawLDEEfun(readdir,savedir,tspname(1,idx), 10);
end    



