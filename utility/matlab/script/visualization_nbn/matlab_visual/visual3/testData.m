

filepath = append(fileDir,'/',filename);
Psize = 20;
NBN_originData = read_NBN_DataFun(filepath);
NBN_visual_data = transfer_NBN_visual_data(NBN_originData);
fh= figure(1);
clf(fh);
nbnVisualFun(NBN_visual_data, outputDir, filename,Psize);
fh= figure(1);
clf(fh);
cut = 1000;
contour3Level=200;
contour3nbnVisualFun(NBN_visual_data, outputDir, filename,Psize,cut,contour3Level);
fh= figure(1);
clf(fh);
numColor=1000;
meshNbnVisualFun(NBN_visual_data, NBN_originData,outputDir, filename,numColor,Psize)