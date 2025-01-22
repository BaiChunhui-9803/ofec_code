
%dir = "/mnt/share151/";
dirName = "D:/";
filepath =dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8/";
filepath2= dirName+"student/2018/DiaoYiya/nbnConPro_experiment_data/conProNBN_total_network8_df/";
pronamefile= dirName +"student/2018/DiaoYiya/nbnConPro_experiment_data/nbnConFilenames.txt";

filetasks = readlines(pronamefile);
%filetasks= dir(filepath);
numFile= size(filetasks,1);

nbn_name = "_nbn";
mesh_name = "_mesh";
txt_name =".txt";
picture_name ="_pic";
df_name = "_df";




outputdir = "E:/DiaoYiya/experiment_figure/conProFigure/";
mkdir(outputdir);
numThread =20; %设置处理器数目
%delete(par) %关闭并行计算
%delete(gcp('nocreate'));
par = parpool('local', numThread); %开启并行计算

%delete(par) %关闭并行计算
error_number= zeros( numFile,2);
error_info = string(error_number);
numColor = 1000;
basicNodeSize = 20;
conT = 1000;


parfor k=1:numFile
%k=309;
%k=1;

try
    figureId= mod(k,numThread)+200;
    %f = figure(figureId);
    filename=filetasks(k);
    display(filename);
    if filename~=""


        nbn_mat  = readmatrix(filepath+filename+"_nbn.txt");
        network_mat =readmatrix(filepath2+ filename + "_network.txt", "NumHeaderLines",1);
        nbn_df = readmatrix(filepath2+filename + "_FitDisRadius.txt");
        nbn_solInfo = readmatrix(filepath + filename + "_nbnSolInfo.txt", "NumHeaderLines",1);


edge_a = nbn_mat(:,1)';
edge_b = nbn_mat(:,3)';
num_edge =size(edge_a,2);
for idx= 1:num_edge
    edge_a(idx) = edge_a(idx)+1;
    edge_b(idx) = edge_b (idx) + 1;
end
nbn_graph = graph(edge_a,edge_b);

mapColor = jet(numColor);
node_fit = nbn_solInfo(:,2);
node_ColorIdx = uint32(rescale(node_fit,1,numColor));
num_node =  size(node_fit,1);
node_color=  zeros(num_node,3);
numConstrait=0;
numDynamic=0;
for idx=1:num_node
    if nbn_solInfo(idx,4)~=1
        numConstrait= numConstrait+1;
    end
    if nbn_solInfo(idx,5)~=0
        numDynamic= numDynamic +1;
    end
end
for idx=1:num_node
    node_color(idx,:) = mapColor(node_ColorIdx(idx),:);
end
node_size = zeros(num_node,1);
for idx=1:num_node
    node_size(idx)=basicNodeSize;
end
gray_color = [0.8 0.8 0.8];
black_color = [0 0 0];


if numConstrait~=0
    nodeContrait_color= node_color;
    for idx=1:num_node
    if nbn_solInfo(idx,4)~=1
    nodeContrait_color(idx,:) =black_color;
    end
    end
end

if numDynamic~=0
    node_dynamic_color = node_color;
    node_ColorIdx = uint32(rescale(nbn_solInfo(:,5),1,numColor));
    for idx=1:num_node
    node_dynamic_color(idx,:) =mapColor(node_ColorIdx(idx),:);
    end
end


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
scatter(df_mat(:,3),df_mat(:,2),df_node_size,df_color,'filled');
 ExportThreeFigures(outputdir+ filename +"_df" );
 delete(scatter);
clear("df_mat");
clear("df_node_size");
clear("df_color");


        if endsWith(filename,"_2")
            mesh_name =  replace(filename,"_Dim_2","_mesh");
            nbn_nbnfit = readmatrix(filepath + filename + "_nbnFit.txt");
            meshData = readmatrix(filepath+mesh_name+".txt");
            x = meshData(:,3)';
            y=meshData(:,4)';
            z= meshData(:,5)';
            
            sampleSize = size(x);
            div = sqrt(sampleSize(1,2));
            xx= reshape(x,div,div);
            yy = reshape(y,div,div);
            zz = reshape(z,div,div);
            
            
%f= figure(figureId);
clf(f);
            pro_contour= contour3(xx,yy,zz,conT);
            ExportThreeFigures(outputdir+ filename +"_pro_contour" );
        %    hold on;
        %    fit_nbn_plot = plot(nbn_graph,'XData',nbn_nbnfit(:,2),'YData',nbn_nbnfit(:,3),'ZData',nbn_nbnfit(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
        %    hold on;
        %    fit_node_scatter = scatter3(nbn_nbnfit(:,2),nbn_nbnfit(:,3),nbn_nbnfit(:,4),node_size,node_color,'filled');
clf(f);
delete(pro_contour);
        end



%f= figure(figureId);
clf(f);


nbn_plot = plot(nbn_graph,'XData',network_mat(:,2),'YData',network_mat(:,3),'ZData',network_mat(:,4),'NodeColor', gray_color , 'LineWidth',0.1 ,'EdgeColor',gray_color);
hold on;
nbn_points=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_color,'filled');
%ExportThreeFigures(outputdir+ filename +"_fitness_nbn" );
set(nbn_points, 'visible','off');
delete(nbn_points);


if numConstrait~=0
hold on;
nbn_points_con=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,nodeContrait_color,'filled');
ExportThreeFigures(outputdir+ filename +"_constraint_nbn" );
set(nbn_points_con, 'visible','off');
delete(nbn_points_con);
end


if numDynamic~=0
    hold on;
    nbn_points_dyn=  scatter3(network_mat(:,2),network_mat(:,3),network_mat(:,4),node_size,node_dynamic_color,'filled');
    ExportThreeFigures(outputdir+ filename +"_dynamic_nbn" );
    set(nbn_points_dyn, 'visible','off');
    delete(nbn_points_dyn);
end





%

    end

catch ME
    display(filename + ME.message);
    error_info(k,1) = filename;
    %error_info(k,2) = ME.message;
end

end
delete(par) %关闭并行计算
writematrix(error_info,'error.txt');