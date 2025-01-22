
readdir = "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\merge_network\";
readdir = "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers_lkh_fail\merge_network\";
readdir2= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\funnelInfo\";
readdir2= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers_lkh_fail\funnelInfo\";

%readdir2= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\nbn_data_eax_sampling_1000000_final_nbnNetworks_remote_bridgeInfo\";
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/nbn_subregion_bridgeInfo/';
outputdir ='\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/nbn_subregion_bridge_lkh_fail\';

mkdir(outputdir);
str = [ "u574","5955" "2202"  "1281" "6702" "6717" ];
str = ["5955", "u574", "2202"];
str = ["5955"];
strSize = size(str);
strSize= strSize(1,2);

basicNodeSize = 8;
numColor = 1000;

nbn_subfixs= ["eaxRun", "randomSamples"];
nbn_subfixs= ["randomSamples"];
nbn_subfixsSize = size(nbn_subfixs);
nbn_subfixsSize= nbn_subfixsSize(1,2);
neighborK_subfixs = ["3", "6", "12", "25", "50", "100"];
neighborK_subfixs = ["45"];
neighborK_subfixs=["8","10","15"];
neighborK_subfixsSize = size(neighborK_subfixs);
neighborK_subfixsSize= neighborK_subfixsSize(1,2);


numSampleSubfixs= ["10000", "100000", "1000000"];

numSampleSubfixsSize = size(numSampleSubfixs);
numSampleSubfixsSize= numSampleSubfixsSize(1,2);

for taskId= 1:strSize
clf;
tspname = str(1,taskId);
    % for nameId= 1:nbn_subfixsSize
    %     filename=  tspname+"_"+ nbn_subfixs(1,nameId);
    %     disp(filename);
    % 
    %      drawNBNbridge(readdir, readdir2,outputdir, filename);
    % 
    %     try
    %      drawNBNbridge(readdir, readdir2,outputdir, filename);
    %    catch ME
    %         disp(ME.identifier);
    %     end
    % end

        % drawNBNbridge(readdir, readdir2, outputdir, filename);
         %     drawNBN_twoBridge(readdir, readdir2, outputdir, filename);
           %   drawNBNfunnels(readdir, readdir2, outputdir, filename);
    for nameId= 1:neighborK_subfixsSize

        for sampleSizeNameId= 1:numSampleSubfixsSize
        
        filename=  tspname+"_randomSamples_neighborK_"+ neighborK_subfixs(1,nameId) + "_numSamples_" +numSampleSubfixs(1, sampleSizeNameId) ;
        disp(filename);
        try
          %  drawNBN(readdir,outputdir,filename);
          drawNBNfunnels(readdir, readdir2, outputdir, filename);
        catch ME
            disp(ME.identifier);
        end

        end
    end
    

end
