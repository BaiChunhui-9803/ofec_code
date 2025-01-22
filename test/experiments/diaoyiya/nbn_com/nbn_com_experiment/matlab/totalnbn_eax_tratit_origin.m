readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_2/';
%readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_lon_nbndata2/';
readDriTrait= "\\172.24.242.8\share\Student\2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_4_traitdata/";
readDriTrait=readdir;
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/totalNBN_eax_figure_test/';

mkdir(outputdir);
str = ["1281" "2202" "6702" "6717" "u574","5955" ];
str = ["u574" "1281" "2202"  "5955" ];
strSize = size(str);
strSize= strSize(1,2);
formatSpec="%f";
basicNodeSize = 8;
numColor = 1000;
mapColor = jet(numColor);
solMapColor = jet(numColor);

for taskId= 1:strSize
clf;
tspname = str(1,taskId);
    
nbn_filenpath =strcat(readdir,tspname);
nbn_filenpath=strcat(nbn_filenpath,"_nbn.txt");
      %  nbn_mat  = readmatrix(nbn_filenpath);
nbn_filenpath =strcat(readdir,tspname);
nbn_filenpath=strcat(nbn_filenpath,"_network.txt");    
        
network_mat =readmatrix(nbn_filenpath, "NumHeaderLines",1);
curoutfilename= outputdir + tspname + "_eax_lkh";
curoutfilename =  outputdir + tspname + "_origin";

          drawNBN(readdir, curoutfilename, tspname);


end