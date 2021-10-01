%{
%% default args given data
hexseed = randi(data.nH); faceseed = data.H2Farray(hexseed,1); 

%% visualize starting config
visualizeHmeshData(data,figure,.00001);
patch('faces',data.F(data.H2Farray(hexseed,:),:),'vertices',data.V,'edgecolor','c','linewidth',2,'facecolor','c','facealpha',.2)
patch('faces',data.F(faceseed,:),'vertices',data.V,'edgecolor','g','linewidth',2,'facecolor','g','facealpha',1)

%% visualize current state
patch('faces',data.F(unique(data.H2Farray(hexchord,:)),:),'vertices',data.V,'edgecolor','c','linewidth',2,'facecolor','c','facealpha',.2)
%}
function [logicalchord, hexchord, facechord] = extendRegularHexChord(data, hexseed, faceseed, isBlockedFace)
    if ~exist('isBlockedFace','var'); isBlockedFace=false(data.nF,1); end;

    logicalchord = false(size(data.H,1),1);
    hexchord = [];
    facechord = [];
    
    currhex = hexseed;
    currface = faceseed;
    while true
        logicalchord(currhex) = true;
        hexchord(end+1) = currhex;
        facechord(end+1) = currface;
        
        F6 = data.H2Farray(currhex,:);
        nextface = F6(data.H2F_flip(find(F6==currface))); % opposite face
        
        if data.isBoundaryFace(nextface) || nextface == faceseed || isBlockedFace(nextface)
            facechord(end+1) = nextface;
            break;
        end
        
        nexthex = setdiff(data.F2Harray(nextface,[1 3]), currhex);
        
        currhex = nexthex;
        currface = nextface;
    end
    
end