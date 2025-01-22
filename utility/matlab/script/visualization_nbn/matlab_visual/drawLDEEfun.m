function[] = drawLDEEfun(readdir,savedir, tspname,fontsize)

ve_mat = readmatrix(readdir+ tspname + "_IDEE_embedded_ve.txt");
solInfo_mat = readmatrix(readdir+ tspname+"_attributes.txt");
maxDim = max(ve_mat(:,2));
[X,Y] = meshgrid(0:1:maxDim);
maxDim= maxDim+1;
sizeVe= size(ve_mat,1);
Z= zeros(maxDim, maxDim);
for idx=1:sizeVe
    Z(ve_mat(idx,2)+1,ve_mat(idx,3)+1)= solInfo_mat(ve_mat(idx,1),2);
end
colormap jet;
surf(X,Y,Z,'EdgeColor','interp');
ax=gca;
ax.FontSize = fontsize;
view(2);
filename = savedir  + tspname + "_fitness";
setExportFigureType(filename,'origin',0.15);

Z= zeros(maxDim, maxDim);
for idx=1:sizeVe
    Z(ve_mat(idx,2)+1,ve_mat(idx,3)+1)= solInfo_mat(ve_mat(idx,1),3);
end
colormap jet;
surf(X,Y,Z,'EdgeColor','interp');
view(2);
ax=gca;
ax.FontSize = fontsize;
filename = savedir  + tspname + "_generation";
setExportFigureType(filename,'origin',0.15);

Z= zeros(maxDim, maxDim);
for idx=1:sizeVe
    Z(ve_mat(idx,2)+1,ve_mat(idx,3)+1)= solInfo_mat(ve_mat(idx,1),4);
end
colormap jet;
surf(X,Y,Z,'EdgeColor','interp');
view(2);
ax=gca;
ax.FontSize = fontsize;
filename = savedir  + tspname + "_run";
setExportFigureType(filename,'origin',0.15);


Z= zeros(maxDim, maxDim);
for idx=1:sizeVe
    Z(ve_mat(idx,2)+1,ve_mat(idx,3)+1)= solInfo_mat(ve_mat(idx,1),5);
end
colormap jet;
surf(X,Y,Z,'EdgeColor','interp');
view(2);
ax=gca;
ax.FontSize = fontsize;
filename = savedir  + tspname + "_runGen";
setExportFigureType(filename,'origin',0.15);
end