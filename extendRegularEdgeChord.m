%{
edgeseed = randi(data.nE); vertseed = data.E(edgeseed,1); 
%}
%{
visualizeHmeshData(data,figure,.00001);
blockVs = data.V(data.isSingularVertex | data.isBoundaryVertex | isBlockedVert,:);
scatter3(blockVs(:,1),blockVs(:,2),blockVs(:,3),50,'r','filled')
patch('faces',data.E(edgeseed,:),'vertices',data.V,'edgecolor','c','linewidth',5)
scatter3(data.V(vertseed,1),data.V(vertseed,2),data.V(vertseed,3),50,'g','filled')

patch('faces',data.E(edgechord,:),'vertices',data.V,'edgecolor','c','linewidth',5)

%}
function [logicalchord, edgechord, vertchord] = extendRegularEdgeChord(data, edgeseed, vertseed, isBlockedVert)
    if ~exist('isBlockedVert','var'); isBlockedVert=false(data.nV,1); end;

    logicalchord = false(size(data.E,1),1);
    edgechord = [];
    vertchord = [];
    
    curredge = edgeseed;
    currvert = vertseed;
    while true
        logicalchord(curredge) = true;
        edgechord(end+1) = curredge;
        vertchord(end+1) = currvert;
        
        nextvert = setdiff(data.E(curredge,:),currvert); % next vertex
        
        if data.isBoundaryVertex(nextvert) || data.isSingularVertex(nextvert) || nextvert == vertseed || isBlockedVert(nextvert)
            vertchord(end+1) = nextvert;
            break;
        end
        
        nextedgecands = setdiff(find(data.E2V(:,nextvert)),curredge);
        nextedgeind = find(~ismember(setdiff(data.E(nextedgecands,:)',nextvert), unique(data.F(find(data.E2F(curredge,:)),:))));
        nextedge = nextedgecands(nextedgeind);
        
        curredge = nextedge;
        currvert = nextvert;
    end
    
end