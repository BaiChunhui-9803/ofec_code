nbn_filenpath = '//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data/CWIDN_nbn_data/nbn_figure_dataInt4/TotalResult.txt';
nbn_mat  = readmatrix(nbn_filenpath, 'NumHeaderLines',1);
x= nbn_mat(:,2);
y= nbn_mat(:,3);
plot(x,y);
xlabel('Distance')
ylabel('Similarity of NBN')