clf;
filedir = '//172.24.24.151/e/DiaoYiya/a_final_NBN_result/rugeness/';
filename = 'rug_compare10';
filepath=[filedir,filename,'.txt'];
fileID = fopen(filepath,'r');

[t,s] =title('D=10');
t.FontSize = 18;
s.FontAngle = 'italic';

tline = fgetl(fileID);
formatSpec = '%e';
mat = fscanf(fileID,formatSpec,[9,100]);
for idx=1:8
    mat(idx,:)=rescale( mat(idx,:),0,1);
end
numLS=3;
line_style= ["--",":","-."];
color_style= {	'#0072BD', '#D95319'};

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
 
 legend_str ={	'h.max', 'eps.s', 'eps.max','mo','eps.ratio','$I_r$'};
 ylim([-0 1.3])
set(gca,'fontsize', 15);
legend(legend_str,'Location','NorthEast','NumColumns',3, 'Interpreter','latex','fontsize', 16);
 xlal=  xlabel('$\alpha$','Interpreter','latex',"FontSize",18, "FontAngle", 'italic') ;
% xlal.FontAngle = 'italic';
ylal=ylabel('Measured value','Interpreter','latex',"FontSize",18, "FontAngle", 'italic') ;
% ylal.FontAngle = 'italic';
 
 
% legend_str ={	'h.max', 'eps.s', 'eps.max','mo','eps.ratio','MTNs'};
% ylim([-0 1.3])
%set(gca,'fontsize', 15);
%legend(legend_str,'Location','NorthWest','NumColumns',3,'fontsize', 15);

filepath=[filedir,filename,'.png'];
exportgraphics(gcf,filepath);
