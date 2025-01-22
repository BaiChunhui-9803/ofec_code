% Sample script to demonstrate execution of function NBN_Data_Info(filedir, filepath, nbn_pos, nbn_pos_color_vals, nbn_pos_zvalues, virtual_nbn_pos_idxs, m_nbn_edges, m_nbn_edges_dis2parent, m_virtual_nbn_edges_idxs, m_opt_nbn_idxs, m_flag_continousFit2D, m_fit2d_pos, m_NBN_id_to_fit2D_id, m_fit2D_id_to_NBN_id, m_ranges, m_mesh2D_colors, m_mesh2D_zval, m_markers_nbn_idxs)
filedir = 'E:/Diao_Yiya/code/OFEC/result/test'; % Initialize filedir here
filepath = 1; % Initialize filepath here
nbn_pos = [1 2 3; 1 2 3]; % Initialize nbn_pos here
nbn_pos_color_vals = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize nbn_pos_color_vals here
nbn_pos_zvalues = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize nbn_pos_zvalues here
virtual_nbn_pos_idxs = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize virtual_nbn_pos_idxs here
m_nbn_edges = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_nbn_edges here
m_nbn_edges_dis2parent = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_nbn_edges_dis2parent here
m_virtual_nbn_edges_idxs = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_virtual_nbn_edges_idxs here
m_opt_nbn_idxs = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_opt_nbn_idxs here
m_flag_continousFit2D = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_flag_continousFit2D here
m_fit2d_pos = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_fit2d_pos here
m_NBN_id_to_fit2D_id = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_NBN_id_to_fit2D_id here
m_fit2D_id_to_NBN_id = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_fit2D_id_to_NBN_id here
m_ranges = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_ranges here
m_mesh2D_colors = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_mesh2D_colors here
m_mesh2D_zval = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_mesh2D_zval here
m_markers_nbn_idxs = [0.1,0.3,0.5;0.3,0.5,0.3]; % Initialize m_markers_nbn_idxs here
NBN_Data_MatlabOutput(filedir, filepath, nbn_pos, nbn_pos_color_vals, nbn_pos_zvalues, virtual_nbn_pos_idxs, m_nbn_edges, m_nbn_edges_dis2parent, m_virtual_nbn_edges_idxs, m_opt_nbn_idxs, m_flag_continousFit2D, m_fit2d_pos, m_NBN_id_to_fit2D_id, m_fit2D_id_to_NBN_id, m_ranges, m_mesh2D_colors, m_mesh2D_zval, m_markers_nbn_idxs);
