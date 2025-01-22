function[] = drawSavedf(filepath2, outputdir, filename,figureId, numColor,basicNodeSize)
mapColor = jet(numColor);
df_mat = readmatrix(filepath2+ filename + "_FitDisRadius.txt");
df_color_idx= log2(df_mat(:,4)*2)+15;
df_color_idx =df_color_idx/15*numColor;
df_color_idx = uint32(df_color_idx);
num_df_node = size(df_mat,1);
df_color =zeros(num_df_node,3);
df_node_size = zeros(num_df_node,1);
for idx=1:num_df_node
    df_node_size(idx)= basicNodeSize;
    color_idx= df_color_idx(idx);
  %  display(color_idx);
    if color_idx<=0
        color_idx=1;
    elseif color_idx>numColor
        color_idx = numColor;
    end
 %       display(color_idx);
 
    df_color(idx,:) = mapColor(color_idx,:);
%           display(df_color(color_idx));
end
f= figure(figureId);
clf(f);
scatter(df_mat(:,3),df_mat(:,2),df_node_size,df_color,'filled' );
setExportFigureType(outputdir+ filename +"_df",'origin',0.1);
end