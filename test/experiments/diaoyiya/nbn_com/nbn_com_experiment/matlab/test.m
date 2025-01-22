filepath = "F:/code/ofec_code/out/build/x64-Debug/testMatlabInput.txt";
fileID = fopen(filepath,'r');
formatSpec = '%f';
sizeA= [1 2];
sizeM= fscanf(fileID,formatSpec,sizeA);
%sizeM = sizeM(1,2);
data= fscanf(fileID,formatSpec,sizeM);


