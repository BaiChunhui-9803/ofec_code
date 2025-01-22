readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_nbnData_2/';
readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_100000_nbnNetwork3/';
readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final_nbnNetworks/';
readdir= '//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final_nbnNetworks/';
readdir= '//172.24.207.203/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final_nbnNetworks_remote/';


readdir= '//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/';

outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/nbn_subproblem_figures/';
outputdir= '//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/experiments_added/nbn_subproblem_figures/';
mkdir(outputdir);
str = [ "u574","5955" "2202"  "1281" "6702" "6717" ];
str = ["5955", "u574", "2202"];
str = ["2202"];
strSize = size(str);
strSize= strSize(1,2);

filenameDir= "//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/experiments_added/tsp_instances/";
tspnames = getTspTaskNames(filenameDir);
strSize = size(tspnames);
strSize= strSize(1,2);

basicNodeSize = 8;
numColor = 1000;

nbn_subfixs= ["eaxRun", "randomSamples"];
nbn_subfixs= ["randomSamples"];
nbn_subfixsSize = size(nbn_subfixs);
nbn_subfixsSize= nbn_subfixsSize(1,2);
neighborK_subfixs = ["3", "6", "12", "25", "50", "100"];

neighborK_subfixs = ["12", "25", "50"];
neighborK_subfixsSize = size(neighborK_subfixs);
neighborK_subfixsSize= neighborK_subfixsSize(1,2);

for taskId= 1:strSize
clf;
%tspname = str(1,taskId);

tspname = tspnames(taskId);
%tspname = tspname(1:end-1);  % 删除最后一个字符
tspname = regexprep(tspname, '\s+$', '');
tspname=tspname{1,1};


    for nameId= 1:nbn_subfixsSize
        filename=  tspname+"_"+ nbn_subfixs(1,nameId);
        disp(filename);
        try
         drawNBN(readdir, outputdir, filename);
       catch ME
            disp(ME.identifier);
        end
    end
    
    for nameId= 1:neighborK_subfixsSize
        filename=  tspname+"_randomSamples_neighborK_"+ neighborK_subfixs(1,nameId);
        disp(filename);
        try
          drawNBN(readdir, outputdir, filename);
        catch ME
            disp(ME.identifier);
        end
    end
    

end
