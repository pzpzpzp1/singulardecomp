function [fh, guielems] = visualizeHmeshData(data, fh)
    if nargin < 2
        fh = figure; 
    end
    hold all; axis equal off; rotate3d on;
    
    V=data.V; F=data.F; E=data.E; 
    isSingularNode = data.isSingularNode;
    isBoundaryEdge = data.isBoundaryEdge;
    isBoundaryFace = data.isBoundaryFace;
    isSingularEdge = data.isSingularEdge;
    
    guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0);
    guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]));
    guielems{3} = patch('vertices',V,'faces',E(isSingularEdge,[1 2 1]),'linewidth',3,'edgecolor','blue');
    guielems{4} = scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled');
    
end