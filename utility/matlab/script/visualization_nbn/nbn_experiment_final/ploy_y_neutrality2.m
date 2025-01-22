
filedir = '//172.24.24.151/e/DiaoYiya/a_final_NBN_result/neutrality/';
filename = 'ic_compare10';
filepath=[filedir,filename,'.txt'];
fileID = fopen(filepath,'r');

tline = fgetl(fileID);
formatSpec = '%e';
mat = fscanf(fileID,formatSpec,[9,99]);
for idx=1:8
    mat(idx,:)=rescale( mat(idx,:),0,1);
end
line_style= {'-','--',':','-.'};

x= mat(2,:);
x = rescale(x,0,1);
hmax = mat(3,:);
hmax = rescale(hmax,0,1);
eps_s = mat(4,:);
eps_s = rescale(eps_s,0,1);
eps_max = mat(5,:);
eps_max = rescale(eps_max,0,1);
eps_ratio = mat(6,:);
eps_ratio = rescale(eps_ratio,0,1);
mo = mat(7,:);
mo = rescale(mo,0,1);
nbn = mat(8,:);
nbn = rescale(nbn,0,1);
set(gca,'LineStyleOrder',{'-','--',':','-.'})
plot(x,hmax,x,eps_s,x,eps_max,x,eps_ratio,x,mo,x,nbn);

