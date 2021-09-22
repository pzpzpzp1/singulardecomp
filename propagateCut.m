function cutfaces = propagateCut(data,startnode,cutseed)
    cutfaces = false(size(data.F,1),1);
%     cutseed = cutseed(randperm(numel(cutseed))); 
%     fseed = cutseed(1);
%     hexseed = data.F2Harray(cutseed(1),1);
%     scatter3(data.cellBarycenters(hexseed,1),data.cellBarycenters(hexseed,2),data.cellBarycenters(hexseed,3),'r')
    
    %% initial cut propagation from start node
    %forbiddenfaces = false(data.nF,1); % faces that cannot be cut.
    %forbiddenfaces(setdiff(find(sum(data.F2V(:,startnode.ind),2)), cutseed))=true;
    
    cutfaces = false(size(data.F,1),1);
    for i=1:numel(cutseed)
        cutfaces(extendRegularFaces(data, cutseed(i)))=true;
    end
    %assert(~any(cutfaces(forbiddenfaces)));
%     patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','r','facealpha',.1);
%     patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','c')

    %% keep propagating cut
    resolvednodes = startnode.ind;
    while true
        % get a node that isn't fully cut
        cutnodes = intersect(unique(data.F(cutfaces,:)), find(data.isSingularNode));
        unresolvednodes = setdiff(cutnodes, resolvednodes);
        if numel(unresolvednodes)==0
            % no singular nodes are left hanging. means the cut can go through now.
            break;
        end
        node = getNode(data, unresolvednodes(1));
        
        % get a split
        facesToCut = selectSplit(data, node, cutfaces);
        for i=1:numel(facesToCut)
            cutfaces(extendRegularFaces(data, facesToCut(i)))=true;
        end
        resolvednodes(end+1) = node.ind;
        % not sure the forbidden face check here does anything or always passes.
        %forbiddenfaces(setdiff(find(sum(data.F2V(:,node.ind),2)), cutfaces))=true;
        %assert(~any(cutfaces(forbiddenfaces)));
        
%         patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','r','facealpha',.5)
%         patch('vertices',data.V,'faces',data.F(facesToCut,:),'facecolor','c')
        
        
    end
    
end