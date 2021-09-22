% get set of faces that are adjacent to faceseed in a regular region. will
% stop at paralell singular edges or boundary.
function cutfaces = extendRegularFaces(data, faceseed)
    cutfaces = false(size(data.F,1),1);
    edgeseeds = data.F2Earray(faceseed,:);
    
    for i = [1 3]
        [logicalchord, facechord, edgechord] = extendRegularFaceChord(data, faceseed, edgeseeds(1));
        
        for j=1:numel(facechord)
            edirs = setdiff(data.F2Earray(facechord(j),:), edgechord);
            
            logicalchord1 = extendRegularFaceChord(data, facechord(j), edirs(1));
            logicalchord2 = extendRegularFaceChord(data, facechord(j), edirs(2));
            
            cutfaces(logicalchord1) = true;
            cutfaces(logicalchord2) = true;
        end
        cutfaces(logicalchord) = true;
    end
    
    
end