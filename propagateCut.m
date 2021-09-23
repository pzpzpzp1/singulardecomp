function cutfaces = propagateCut(data,startnode,cutseed)
    cutfaces = false(size(data.F,1),1);
%     cutseed = cutseed(randperm(numel(cutseed))); 
%     fseed = cutseed(1);
%     hexseed = data.F2Harray(cutseed(1),1);
%     scatter3(data.cellBarycenters(hexseed,1),data.cellBarycenters(hexseed,2),data.cellBarycenters(hexseed,3),'r')
    
    %% initial cut propagation from start node
    forbiddenfaces = false(data.nF,1); % faces that cannot be cut.
    forbiddenseed = setdiff(startnode.F_adj, cutseed);
    for i=1:numel(forbiddenseed)
        if forbiddenfaces(forbiddenseed(i))==false
            forbiddenfaces(extendRegularFaces(data, forbiddenseed(i)))=true;
        end
        patch('vertices',data.V,'faces',data.F(forbiddenfaces,:),'facecolor','r','facealpha',.5);
    end
    
    cutfaces = false(size(data.F,1),1);
    patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','y')
    for i=1:numel(cutseed)
        cutfaces(extendRegularFaces(data, cutseed(i)))=true;
        
        patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','c','facealpha',.5);
    end
    %assert(~any(cutfaces(forbiddenfaces)));

    %% keep propagating cut
    resolvednodes = startnode.ind;
    while true
        % get a node that isn't fully cut
        cutnodes = intersect(unique(data.F(cutfaces,:)), find(data.isSingularNode));
        unresolvednodes = setdiff(cutnodes, resolvednodes);
        if numel(unresolvednodes)==0
            cutQM = getQMfromCut(data,cutfaces);
            if all(data.isBoundaryEdge(cutQM.HmeshCutBoundaryEdgeInds)); 
                % Cut starts and ends on boundary of hex mesh. Can try to sheet insert now.
                break;
            else
                % no singular nodes on cut boundary, but part of the cut boundary is entirely interior singular edge.
                % make up an unresolved node on the interior singular curve
                unresolvednodes = unique(data.E(data.isBoundaryEdge(cutQM.HmeshCutBoundaryEdgeInds),:));
                unresolvednodes(data.isBoundaryVertex)=[];
                
                unresolvednodes = 0;
            end
        end
        % pull off next singular node to process
        node = getNode(data, unresolvednodes(1)); 
        
        % get a cut local to node, and propagate it.
        facesToCut = selectSplit(data, node, cutfaces, forbiddenfaces);
        patch('vertices',data.V,'faces',data.F(facesToCut,:),'facecolor','y')
        for i=1:numel(facesToCut)
            cutfaces(extendRegularFaces(data, facesToCut(i)))=true;
            
            patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','c','facealpha',.5)
        end
        % mark all other faces local to node as forbidden.
        forbiddenseed = setdiff(node.F_adj, facesToCut);
        for i=1:numel(forbiddenseed)
            if forbiddenfaces(forbiddenseed(i))==false
                forbiddenfaces(extendRegularFaces(data, forbiddenseed(i)))=true;
            end
            patch('vertices',data.V,'faces',data.F(forbiddenfaces,:),'facecolor','r','facealpha',.5);
        end
        
        % this node is now resolved. it is cut and cannot be cut further.
        resolvednodes(end+1) = node.ind;
        
    end
    
end