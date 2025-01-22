function data = inputMatrix(fileID, errInfo, rowSize)
    format('longEng');
    formatSpec = '%e';
    sizeNum=1;
    Dim = fscanf(fileID,formatSpec,sizeNum);
    data.Dim = Dim;
    data.ErrFlag = 1;
    Mat = -1;
    if Dim <= 0
    else 
       dimM= [rowSize Dim];
       Mat = fscanf(fileID,formatSpec,dimM);
       [Mat_d1,Mat_d2] = size(Mat);
       if Mat_d2~=Dim
            outputInfo = append(errInfo,' ','pos not equal');
            disp(outputInfo);
            return;
       else 
       outputInfo = append(errInfo,' ','continue');
%       disp(outputInfo);
       end
    end
    data.ErrFlag = 0;
    data.mat = Mat;
    data.fileID= fileID;
end