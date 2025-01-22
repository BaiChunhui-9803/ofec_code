function []= NBN_Data_Info2(filedir, filepath,  nbn_pos,nbn_pos_color_vals,nbn_pos_zvalues, virtual_nbn_pos_idxs, m_nbn_edges,m_nbn_edges_dis2parent,m_virtual_nbn_edges_idxs, m_opt_nbn_idxs, m_flag_continousFit2D, m_fit2d_pos,m_NBN_id_to_fit2D_id,m_fit2D_id_to_NBN_id,m_ranges,m_mesh2D_colors,m_mesh2D_zval,m_markers_nbn_idxs)
filename = filedir +"/NBN_matlab_data/"+ filepath+".xls";
Fdir = filedir +"/NBN_matlab_data/"+ filepath;
mkdir(Fdir);
writematrix(nbn_pos,filename,'Sheet',1);
writematrix(nbn_pos_color_vals,filename,'Sheet',2);
writematrix(nbn_pos_zvalues,filename,'Sheet',3);
writematrix(virtual_nbn_pos_idxs,filename,'Sheet',4);

writematrix(m_nbn_edges,filename,'Sheet',5);
writematrix(m_nbn_edges_dis2parent,filename,'Sheet',6);
writematrix(m_virtual_nbn_edges_idxs,filename,'Sheet',7);

writematrix(m_flag_continousFit2D,filename,'Sheet',8);

writematrix(m_flag_continousFit2D,filename,'Sheet',9);
writematrix(m_fit2d_pos,filename,'Sheet',10);
writematrix(m_NBN_id_to_fit2D_id,filename,'Sheet',11);
writematrix(m_fit2D_id_to_NBN_id,filename,'Sheet',12);

writematrix(m_ranges,filename,'Sheet',13);
writematrix(m_mesh2D_colors,filename,'Sheet',14);
writematrix(m_mesh2D_zval,filename,'Sheet',15);
writematrix(m_markers_nbn_idxs,filename,'Sheet',16);


end