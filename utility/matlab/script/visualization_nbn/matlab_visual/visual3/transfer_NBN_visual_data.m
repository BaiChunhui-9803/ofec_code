function NBN_visual_data = transfer_NBN_visual_data(NBN_originData)
NBN_EdgeA  = NBN_originData.NBN_edge(1,:);
NBN_EdgeA =NBN_EdgeA';
[numEdge,k] = size(NBN_EdgeA);
for idx = 1:numEdge
     NBN_EdgeA(idx)= NBN_EdgeA(idx)+1;
end
NBN_EdgeB  = NBN_originData.NBN_edge(2,:);
NBN_EdgeB =NBN_EdgeB';
for idx = 1:numEdge
     NBN_EdgeB(idx)= NBN_EdgeB(idx)+1;
end
NBN_visual_data.NBNgraph = graph(NBN_EdgeB,NBN_EdgeA);
[k, numNBNpos] = size(NBN_originData.NBN_pos);
%disp('num of nbn pos');
%disp(numNBNpos);
NBN_visual_data.NBNpos = NBN_originData.NBN_pos';
%enum class EColorId { kNorFitness = 0, kFitnessOrder=1, 
%kDynamicFitness=2, kContraint=3, kNorId
conFitColorVal = NBN_originData.NBN_pos_color(4,:);
for idx= 1:numNBNpos
    conFitColorVal(idx)= NBN_originData.NBN_pos_color(1,idx)*NBN_originData.NBN_pos_color(4,idx);
end
NBN_visual_data.conFitColor = ValueToColor(conFitColorVal,1000);
idColorVal =  NBN_originData.NBN_pos_color(5,:);
for idx=1:numNBNpos
    idColorVal(idx) = idx;
end
idColorVal=idColorVal(randperm(length(idColorVal)));
NBN_visual_data.ids= idColorVal;
NBN_visual_data.idColor = ValueToColor(idColorVal,1000);
idConColorVal =  NBN_originData.NBN_pos_color(5,:);
for idx= 1:numNBNpos
    idConColorVal(idx)= NBN_originData.NBN_pos_color(4,idx)*idColorVal(idx);
end
NBN_visual_data.idConColor = ValueToColor(idConColorVal,1000);
dimV = size(NBN_originData.NBN_vitual_pos_idxs);
if(dimV(1)==1&&dimV(1)==1&&NBN_originData.NBN_vitual_pos_idxs==-1)
    
else
    for idx = 1:dimV
        NBN_visual_data.conFitColor(NBN_originData.NBN_vitual_pos_idxs(idx)+1,:) = [1 1 1];
        NBN_visual_data.idColor(NBN_originData.NBN_vitual_pos_idxs(idx)+1,:) = [1 1 1];
        NBN_visual_data.idConColor(NBN_originData.NBN_vitual_pos_idxs(idx)+1,:) = [1 1 1];
    end
end

NBN_visual_data.NBN_opt_idxs= NBN_originData.NBN_opt_idxs;

dimO  = size(NBN_originData.NBN_opt_idxs);
if(dimO(1)==1&&dimO(1)==1&&NBN_originData.NBN_opt_idxs==-1)
else
    for idx =1:dimO
        NBN_visual_data.NBN_opt_idxs(idx) = NBN_visual_data.NBN_opt_idxs(idx)+1;
    end
    for idx = 1:dimO
        NBN_visual_data.NBNpos(NBN_visual_data.NBN_opt_idxs(idx),3) = 1.005;
    end
end

if NBN_originData.flag_continousFit2D==0
    return;
end
Fit2d_pos  = NBN_originData.Fit2d_pos;
[k, numFitPos] = size(NBN_originData.Fit2d_pos);
%disp('num of fit pos');
%disp(numFitPos);
NBN_visual_data.numFitPos= numFitPos;
if numFitPos<=0
    return;
end
%NumPos= size(Fit2d_pos(3,:));
for idx= 1:numFitPos
   Fit2d_pos(1,idx) =  (Fit2d_pos(1,idx)+5)/10;
   Fit2d_pos(2,idx) =  (Fit2d_pos(2,idx)+5)/10;
end
for idx= 1:numFitPos
   Fit2d_pos(3,idx) =  (Fit2d_pos(3,idx)+3)/6;
end
NBN_visual_data.Fit2d_pos= Fit2d_pos';
NBN_visual_data.Fit2d_pos_z = NBN_visual_data.Fit2d_pos(:,3);
for idx= 1:numFitPos
   NBN_visual_data.Fit2d_pos_z(idx) =  NBN_visual_data.Fit2d_pos_z(idx)+0.005; 
end
    NBN_visual_data.conFitColorFit2dPos = zeros(numFitPos,3);
    NBN_visual_data.idColorFit2dPos= zeros(numFitPos,3);
    NBN_visual_data.idConColorFit2dPos= zeros(numFitPos,3);
for idx=1:numNBNpos
    NBN_originData.NBN_id_to_Fit2D_id(idx) =NBN_originData.NBN_id_to_Fit2D_id(idx)+1;
end
for idx=1:numFitPos
    NBN_originData.Fit2D_id_to_NBN_id(idx) =NBN_originData.Fit2D_id_to_NBN_id(idx)+1;
end
for idx=1:numFitPos
    NBN_visual_data.conFitColorFit2dPos(idx,:) = NBN_visual_data.conFitColor(NBN_originData.Fit2D_id_to_NBN_id(idx),:);
    NBN_visual_data.idColorFit2dPos(idx,:) = NBN_visual_data.idColor(NBN_originData.Fit2D_id_to_NBN_id(idx),:);
    NBN_visual_data.idConColorFit2dPos(idx,:) = NBN_visual_data.idConColor(NBN_originData.Fit2D_id_to_NBN_id(idx),:); 
end
numEdgeFit2D  = 0;
fit2dNBNEdges = zeros(numEdge,2);
fit2dNBNEdgesId = zeros(numEdge,1);
for idx = 1:numEdge
    nbn_from= NBN_EdgeA(idx);
    nbn_to = NBN_EdgeB(idx);
    fit2d_from =     NBN_originData.NBN_id_to_Fit2D_id(nbn_from);
    fit2d_to =   NBN_originData.NBN_id_to_Fit2D_id(nbn_to);
    fit2dNBNEdges(idx,1)= fit2d_from;
    fit2dNBNEdges(idx,2) = fit2d_to;
    if fit2d_from>=1&&fit2d_from<=numFitPos&&fit2d_to>=1&&fit2d_to<=numFitPos
        numEdgeFit2D = numEdgeFit2D+1;
        fit2dNBNEdgesId(numEdgeFit2D,1) = idx;
    end
end
fit2dEdges= zeros(numEdgeFit2D,2);
for idx = 1:numEdgeFit2D
%    disp('%d \: %d',idx,fit2dNBNEdgesId(idx,1));;
%    disp(fit2dNBNEdgesId(idx,1));
    fit2dEdges(idx,:)= fit2dNBNEdges( fit2dNBNEdgesId(idx,1),:);
end
NBN_visual_data.Fit2Dgraph = graph(fit2dEdges(:,1),fit2dEdges(:,2));
NBN_visual_data.NBN_opt_idxsFit2D = NBN_visual_data.NBN_opt_idxs;
[numOpt,tmp]= size(NBN_visual_data.NBN_opt_idxs);
for  idx =1:numOpt
    NBN_visual_data.NBN_opt_idxsFit2D =     NBN_originData.NBN_id_to_Fit2D_id(NBN_visual_data.NBN_opt_idxsFit2D(idx));
end
NBN_visual_data.Fit2d_pos_oz = zeros(numFitPos,1);
for idx =1:numFitPos
    NBN_visual_data.Fit2d_pos_oz(idx,1) = NBN_originData.NBN_pos_zvalue(1,NBN_originData.Fit2D_id_to_NBN_id(idx));
end
end


