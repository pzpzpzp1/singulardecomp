% get set of faces that are adjacent to faceseed in a regular region. will
% stop at paralell singular edges or boundary.

% fh=visualizeHmeshData(data,figure,.001); patch('faces',data.F(cutfaces,:),'vertices',data.V,'facecolor','r','facealpha',.3)

function cutfaces = extendRegularFaces(data, faceseed0, isBlockedEdge)
    if ~exist('isBlockedEdge','var'); isBlockedEdge=false(data.nE,1); end;
    cutfaces = false(data.nF,1);
    faceseed = faceseed0;
    while true
        cutfaces = cutfaces | crosscross(data, faceseed, isBlockedEdge);

        QM = getQMfromCut(data, cutfaces);
        QMbe = sort(QM.E(QM.isBoundaryEdge,:),2);
        [allf, emap] = ismember(QMbe,sort(data.E,2),'rows'); assert(all(allf));
        incompleteEdges = emap(find(~(data.isSingularEdge(emap) | data.isBoundaryEdge(emap))));
        cutfaceinds = find(cutfaces);
        continuefaceseedinds = find(any(ismember(data.F2Earray(cutfaces,:), incompleteEdges),2));
        continuefaceseeds = cutfaceinds(continuefaceseedinds);
        if numel(continuefaceseeds)~=0
            faceseed = continuefaceseeds(1);
        else
            break;
        end
    end
end

% helper propagator. creates cross at faceseed.
function cutfaces = crosscross(data, faceseed, isBlockedEdge)
    cutfaces = false(size(data.F,1),1);
    edgeseeds = data.F2Earray(faceseed,:);
    
    % for i = [1 3]
    for i = 1:4
        [logicalchord, facechord, edgechord] = extendRegularFaceChord(data, faceseed, edgeseeds(i), isBlockedEdge);
        
        for j=1:numel(facechord)
            edirs = setdiff(data.F2Earray(facechord(j),:), edgechord);
            
            logicalchord1 = extendRegularFaceChord(data, facechord(j), edirs(1), isBlockedEdge);
            logicalchord2 = extendRegularFaceChord(data, facechord(j), edirs(2), isBlockedEdge);
            
            cutfaces(logicalchord1) = true;
            cutfaces(logicalchord2) = true;
        end
        cutfaces(logicalchord) = true;
    end
    
end
    
    
%{

patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','k','facealpha',.5);
patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','y')


%}