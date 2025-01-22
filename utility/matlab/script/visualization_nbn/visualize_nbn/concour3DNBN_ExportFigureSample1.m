% Sample script to demonstrate execution of function mesh2DconnectedEdges(n, EdgeB, EdgeA, xi, yi, zi, colorR, colorG, colorB, markerPtrs, Psize)
fontSize = 20;
numTrainingPoints = 1000; % Number of points we'll have to try to estimate our image everywhere.
imageWidth = 100;
% Create sample surface
z = peaks(imageWidth);
xi = imageWidth * rand(numTrainingPoints, 1);
yi = imageWidth * rand(numTrainingPoints, 1);
zi =imageWidth * rand(numTrainingPoints, 1);
minNum=1;
maxNum=numTrainingPoints;
EdgeA =uint32(minNum):uint32(maxNum);
EdgeB = randi([minNum maxNum],1,maxNum);
num_color = 3;
colorR =  rand(maxNum, 1, 'double');
colorG =  rand(maxNum, 1, 'double');
colorB =  rand(maxNum, 1, 'double');
n=2;
markerPtrs=[];
Psize=10;
cut = 100;
contour3Level=500;
filename= 'E:/Document/technical_report/visualize_of_fitness_landscaple/matlab_figure/show_save_img_setWindows/showMesh/imgNew202209132036';
filetype= 'pdf';
fontsize = 15;
concour3DNBN_ExportFigure(n,EdgeB,EdgeA, xi,yi,zi,colorR,colorG,colorB, cut, markerPtrs,Psize,contour3Level,filename, filetype, fontSize);
fid=fopen('concour3DNBN_ExportFigure.txt','wb');%写入'w'
fwrite(fid,a,'double');
fclose(fid);
