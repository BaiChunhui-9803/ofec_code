
readdir="\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers2\merge_network\";
readdir= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\nbn_data_eax_sampling_1000000_final_nbnNetworks_remote\";
readdir="\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\nbn_subregion_centers\merge_network\";
readdir ="\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers2\merge_network\";
readdir = "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\merge_network\";
readdir = "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\nbn_subregion_two_funnel2_linux\merge_network\";
%readdir= "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\nbn_subregion_two_funnel2_winServer\merge_network\";
readdir = "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\merge_network\";

readdir2= '\\172.24.242.8\share\Student\2018\YiyaDiao\code_total\data\paper_com_experiment_data\totalTsp3\nbn_subregion_twofunnels/';
readdir2= '\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers2\funnelInfo\';
readdir2="\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers2\funnelInfo\";
readdir2= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\nbn_subregion_centers\merge_network_bridgeInfo3\";
readdir2= "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers2\funnelInfo\";
readdir2= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\funnelInfo\";
readdir2= "\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\nbn_subregion_two_funnel2_linux\funnelInfo\";
readdir2= "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp3\nbn_subregion_centers\funnelInfo\";

%readdir2= "\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final\nbn_data_eax_sampling_1000000_final_nbnNetworks_remote_bridgeInfo\";
%readdir2= "\\172.24.24.151\e\DiaoYiya\code\data\ofec-data\paper_com_experiment_data\totalTsp\nbn_subregion_two_funnel2_winServer\funnelInfo\";
outputdir ='//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/total_tsp/nbn_subregion_bridgeInfo/';
outputdir ='\\172.24.24.151\f\DiaoYiya\nbn_com_exp_final_figures/nbn_center_two_funnels/';
mkdir(outputdir);
str = [ "u574","5955" "2202"  "1281" "6702" "6717" ];
str = ["5955", "u574", "2202"];
%str = ["5955"];
strSize = size(str);
strSize= strSize(1,2);

basicNodeSize = 8;
numColor = 1000;

nbn_subfixs= ["eaxRun", "randomSamples"];
nbn_subfixs= ["randomSamples"];
nbn_subfixsSize = size(nbn_subfixs);
nbn_subfixsSize= nbn_subfixsSize(1,2);
neighborK_subfixs = ["3", "6", "12", "25", "50", "100"];
neighborK_subfixs = ["35"];
neighborK_subfixsSize = size(neighborK_subfixs);
neighborK_subfixsSize= neighborK_subfixsSize(1,2);


solIds= ["0","2","3","4"];

solIds= ["0"];
solIdsSize = size(solIds);
solIdsSize= solIdsSize(1,2);

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
        for solIdIter=1:solIdsSize
          solIdName = solIds(1,solIdIter);
        filename=  tspname+"_randomSamples_neighborK_"+ neighborK_subfixs(1,nameId);
      %  filename=  tspname+"_randomSamples_solId_"+solIdName+"_neighborK_"+ neighborK_subfixs(1,nameId)+"_numSamples_1000000";
        disp(filename);

  drawNBNfunnelsAndBridge(readdir, readdir2, outputdir, filename);
        try
          %  drawNBN(readdir,outputdir,filename);
        a=3;
        catch ME
            disp(ME.identifier);
        end
        end
    end
    

end
