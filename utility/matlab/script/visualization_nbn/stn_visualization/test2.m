
filedir = '//172.24.24.151/e/DiaoYiya/a_final_NBN_result/neutrality/';
filename = 'ic_compare5';
filepath=[filedir,filename,'.txt'];
fileID = fopen(filepath,'r');

tline = fgetl(fileID);
formatSpec = '%e';
mat = fscanf(fileID,formatSpec,[9,99]);
for idx=1:8
    mat(idx,:)=rescale( mat(idx,:),0,1);
end
numLS=3;
line_style= ["--",":","-."];
color_style= {	'#0072BD', '#D95319'};
clf;
hold on;
for idx=1:5
  p = plot(mat(2,:),mat(idx+2,:));
  line_id= mod(idx-1,numLS)+1;
 p.LineStyle = line_style(line_id);
 p.LineWidth = 2;
   % curplot.LineStyle = "-";
%  curplot.Color = color_stype
end
hold on;
  p = plot(mat(2,:),mat(8,:));
 p.LineStyle = "-";
 p.LineWidth = 2;
 p.Color= "black";
 
 legend_str ={	'h.max', 'eps.s', 'eps.max','mo','eps.ratio','MTNs'};
 ylim([-0 1.3])

legend(legend_str,'Location','NorthWest','NumColumns',3)
filepath=[filedir,filename,'.png'];
exportgraphics(gcf,filepath);