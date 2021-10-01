function [fh, guielems] = visualizeHmeshData(data, fh, lwfac)
    if nargin < 2
        fh = figure; 
        lwfac=1;
    end
    hold all; axis equal off; rotate3d on; axis image vis3d;
    
    V=data.V; F=data.F; E=data.E; 
    isSingularNode = data.isSingularNode;
    isBoundaryEdge = data.isBoundaryEdge;
    isBoundaryFace = data.isBoundaryFace;
    isSingularEdge = data.isSingularEdge;
    isIntSingularNode = data.isSingularNode & ~data.isBoundaryVertex;
    
    IS3 = data.efdeg==3 & data.isSingularEdge & ~data.isBoundaryEdge;
    IS5 = data.efdeg==5 & data.isSingularEdge & ~data.isBoundaryEdge;
    IS_ = (data.efdeg~=3 & data.efdeg~=5) & data.isSingularEdge & ~data.isBoundaryEdge;
    
    
    guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0);
    guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',.3,'linewidth',3*lwfac);
    guielems{3} = patch('vertices',V,'faces',E(isBoundaryEdge & isSingularEdge,[1 2 1]),'edgealpha',1,'linewidth',9*lwfac);
%     guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',0);
    guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',9*lwfac,'edgecolor','r');
    guielems{5} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',9*lwfac,'edgecolor','b');
    guielems{6} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',9*lwfac,'edgecolor','g');
    
%     guielems{4} = scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled');
    guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),500*lwfac,'k','filled');
    
    % guielems{8} = patch('vertices',V,'faces',F,'facecolor','green','facealpha',.1,'edgealpha',0);
    
end