function NBN_originData = read_NBN_DataFun(filepath)
fileID = fopen(filepath,'r');
tline = fgetl(fileID);
disp(tline);
dataName = fgetl(fileID);
disp(dataName);
NBN_pos = inputMatrix(fileID, 'NBN_pos', 3);
NBN_originData.Dim = NBN_pos.Dim;
NBN_originData.NBN_pos = NBN_pos.mat;
NBN_originData.ErrFlag = NBN_pos.ErrFlag;

if NBN_pos.ErrFlag==1
    return
end
NBN_pos_color = inputMatrix(NBN_pos.fileID, 'NBN_pos_color', 5);
NBN_originData.NBN_pos_color = NBN_pos_color.mat;
NBN_originData.ErrFlag = NBN_pos_color.ErrFlag;
if NBN_pos_color.ErrFlag==1
    return
end
NBN_pos_zvalue = inputMatrix(NBN_pos_color.fileID, 'NBN_pos_zvalue', 1);
NBN_originData.NBN_pos_zvalue = NBN_pos_zvalue.mat;
NBN_originData.ErrFlag = NBN_pos_zvalue.ErrFlag;
if NBN_pos_zvalue.ErrFlag==1
    return
end
NBN_vitual_pos_idxs = inputMatrix(NBN_pos_zvalue.fileID, 'NBN_vitual_pos_idxs', 1);
NBN_originData.NBN_vitual_pos_idxs = NBN_vitual_pos_idxs.mat;
NBN_originData.ErrFlag = NBN_vitual_pos_idxs.ErrFlag;
if NBN_vitual_pos_idxs.ErrFlag==1
    return
end
NBN_edge = inputMatrix(NBN_vitual_pos_idxs.fileID, 'NBN_edge', 2);
NBN_originData.NBN_edge = NBN_edge.mat;
NBN_originData.ErrFlag = NBN_edge.ErrFlag;
if NBN_edge.ErrFlag==1
    return
end
NBN_edges_dis2parent  = inputMatrix(NBN_edge.fileID, 'NBN_edges_dis2parent', 1);
NBN_originData.NBN_edges_dis2parent = NBN_edges_dis2parent.mat;
NBN_originData.ErrFlag = NBN_edges_dis2parent.ErrFlag;
if NBN_edges_dis2parent.ErrFlag==1
    return
end
NBN_virtual_edges_idxs =  inputMatrix(NBN_edges_dis2parent.fileID, 'NBN_virtual_edges_idxs', 1);
NBN_originData.NBN_virtual_edges_idxs = NBN_virtual_edges_idxs.mat;
NBN_originData.ErrFlag = NBN_virtual_edges_idxs.ErrFlag;
if NBN_virtual_edges_idxs.ErrFlag==1
    return
end
NBN_opt_idxs   = inputMatrix(NBN_virtual_edges_idxs.fileID, 'NBN_opt_idxs', 1);
NBN_originData.NBN_opt_idxs = NBN_opt_idxs.mat;
NBN_originData.ErrFlag = NBN_opt_idxs.ErrFlag;
if NBN_opt_idxs.ErrFlag==1
    return
end
flag_continousFit2D   = fscanf(NBN_opt_idxs.fileID,'%d',1);
NBN_originData.flag_continousFit2D = flag_continousFit2D;

Fit2d_pos = inputMatrix(NBN_opt_idxs.fileID, 'Fit2d_pos', 3);
NBN_originData.Fit2d_pos = Fit2d_pos.mat;
NBN_originData.ErrFlag = Fit2d_pos.ErrFlag;
if Fit2d_pos.ErrFlag==1
    return
end
NBN_id_to_Fit2D_id = inputMatrix(Fit2d_pos.fileID, 'NBN_id_to_Fit2D_id', 1);
NBN_originData.NBN_id_to_Fit2D_id = NBN_id_to_Fit2D_id.mat;
NBN_originData.ErrFlag = NBN_id_to_Fit2D_id.ErrFlag;
if NBN_id_to_Fit2D_id.ErrFlag==1
    return
end
Fit2D_id_to_NBN_id = inputMatrix(NBN_id_to_Fit2D_id.fileID, 'Fit2D_id_to_NBN_id', 1);
NBN_originData.Fit2D_id_to_NBN_id = Fit2D_id_to_NBN_id.mat;
NBN_originData.ErrFlag = Fit2D_id_to_NBN_id.ErrFlag;
if Fit2D_id_to_NBN_id.ErrFlag==1
    return
end
ranges  = inputMatrix(Fit2D_id_to_NBN_id.fileID, 'ranges', 2);
NBN_originData.ranges = ranges.mat;
NBN_originData.ErrFlag = ranges.ErrFlag;
if ranges.ErrFlag==1
    return
end
mesh2D_colors = inputMatrix(ranges.fileID, 'mesh2D_colors', 3);
NBN_originData.mesh2D_colors = mesh2D_colors.mat;
NBN_originData.ErrFlag = mesh2D_colors.ErrFlag;
if mesh2D_colors.ErrFlag==1
    return
end
mesh2D_zval = inputMatrix(mesh2D_colors.fileID, 'mesh2D_zval', 1);
NBN_originData.mesh2D_zval = mesh2D_zval.mat;
NBN_originData.ErrFlag = mesh2D_zval.ErrFlag;
if mesh2D_zval.ErrFlag==1
    return
end
markers_nbn_idxs = inputMatrix(mesh2D_zval.fileID, 'markers_nbn_idxs', 3);
NBN_originData.markers_nbn_idxs = markers_nbn_idxs.mat;
NBN_originData.ErrFlag = markers_nbn_idxs.ErrFlag;
if markers_nbn_idxs.ErrFlag==1
    return
end


finish_info = append('finish',' ',dataName);
disp(finish_info);

end