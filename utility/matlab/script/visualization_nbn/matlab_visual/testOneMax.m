filepath  = "//172.24.24.151/e/DiaoYiya/experiment_data/onemax/";
filename = "WModelOneMax_InstanceId_1_DummySelectRate_0_EpistasisBlockSize_0_NeutralityMu_5_RuggednessGamma_0";
display(filepath+filename+"_nbn.txt");
nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
network_mat =readmatrix(filepath+ filename + "_network.txt", "NumHeaderLines",1);