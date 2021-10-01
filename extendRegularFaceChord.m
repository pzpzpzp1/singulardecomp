%{
visualizeHmeshData(data,figure,.00001);
patch('faces',data.E(data.isSingularEdge,[1 2 1]),'vertices',data.V,'edgecolor','r','linewidth',5)
patch('faces',data.E(isBlockedEdge,[1 2 1]),'vertices',data.V,'edgecolor','r','linewidth',3)
patch('faces',data.F(faceseed,:),'vertices',data.V,'facecolor','c')
patch('faces',data.E(edgeseed,[1 2 1]),'vertices',data.V,'edgecolor','g','linewidth',5)

patch('faces',data.F(logicalchord,:),'vertices',data.V,'facecolor','c','facealpha',.5)
%}
function [logicalchord, facechord, edgechord] = extendRegularFaceChord(data, faceseed, edgeseed,isBlockedEdge)
    if ~exist('isBlockedEdge','var'); isBlockedEdge=false(data.nE,1); end;

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
        fedges2 = circshift(fedges,2); % opposite edge
        nextedge = fedges2(find(fedges == curredge));
        
        if data.isBoundaryEdge(nextedge) || data.isSingularEdge(nextedge) || nextedge == edgeseed || isBlockedEdge(nextedge)
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