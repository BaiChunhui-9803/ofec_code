nbn_filenpath = '//172.24.207.203/share/2018/diaoyiya/paper_com_experiment_data/cwidn/cwidn_continous_figureData3/TotalResult.txt';
nbn_mat  = readmatrix(nbn_filenpath, 'NumHeaderLines',1);
x= nbn_mat(:,3);
y= nbn_mat(:,4);
boxchart(x,y,'BoxWidth',0.2);
xlabel('Distance')
ylabel('Similarity of NBN')