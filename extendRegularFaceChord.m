function [logicalchord, facechord, edgechord] = extendRegularFaceChord(data, faceseed, edgeseed)
    logicalchord = false(size(data.F,1),1);
    facechord = [];
    edgechord = [];
    
    currface = faceseed;
    curredge = edgeseed;
    while true
        logicalchord(currface) = true;
        facechord(end+1) = currface;
        edgechord(end+1) = curredge;
        
        fedges = data.F2Earray(currface,:);
        fedges2 = circshift(fedges,2);
        nextedge = fedges2(find(fedges == curredge));
        
        if data.isBoundaryEdge(nextedge) || data.isSingularEdge(nextedge) || nextedge == edgeseed
            edgechord(end+1) = nextedge;
            break;
        end
        
        nextfacecands = setdiff(find(data.E2F(nextedge,:)),currface);
        nextfacecandind = find(~any(ismember(data.F2Harray(nextfacecands,[1 3]), setdiff(data.F2Harray(currface,[1 3]),0)),2));
        nextface = nextfacecands(nextfacecandind);
        
        currface = nextface;
        curredge = nextedge;
    end
    
end