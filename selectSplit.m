% outputs cut faces. 
function cutfaces = selectSplit(data, node, alreadycutfaces)
    assert(~data.isBoundaryVertex(node.ind))
    if nargin < 3
        alreadycutfaces = []
    end
    
    if numel(alreadycutfaces)==0
        %% tetrahedral type
        if numel(node.v)==4
            % handle split on tetrahedral singularity type
            vpath = [1 2 3 4];
            [allf, einds] = ismember(sort([vpath; circshift(vpath,-1)]',2), sort(node.e,2), 'rows');
            assert(all(allf));
            cutfaces = node.e2F(einds);

            return;
        end

        signature = accumarray(data.efdeg(node.v2E),1); if numel(signature) < 6; signature(6)=0; end;
        assert(all(signature(1:2)==0)); % interior singular edges must have 3 or more faces adjacent.
        if sum(signature(6:end))==0
            % only 3,4,5 singularities. must fit into the list of 10 interior singular node types.
            if all(signature(3:5)==[2 2 2]')
                %% (2,2,2) type is a 3-5 joined perpendicularly
                error('unhandled');
            end

            if all(signature(3:5)==[0 4 4]')
                %% (0,4,4) type is a 5-5 joined perpendicularly
                error('unhandled');
            end
        end
        
        %% general case. find a cut based on the inductive simplification procedure.
        error('unhandled');
    else
        %% finish the existing cut. 
        % could be smart about this but for now lets just make a graph, and
        % fill it in arbitrarily.
        partialpath = alreadycutfaces(node.F_adj);
        vertscutcount = accumarray(reshape(node.e(partialpath,:),[],1),1);
        assert(all(ismember(vertscutcount,[0 1 2]))); % otherwise, the existing cut already invalidly cuts this node. 
        
        endpoints = find(vertscutcount == 1);
        if numel(endpoints)==0
            % alreadycutfaces partitions this singular vertex so that no
            % cut is left hanging already. can just continue. There is
            % possibility the existing cut splits triangulation into 3
            % parts though which it's unclear how to deal with. leave it as
            % an error for when trying to split the face.
            cutfaces = [];
            return;
        end
        
        % use shortest path to close up the partial path loop.
        g = graph();
        g = addedge(g, node.e(~partialpath,1), node.e(~partialpath,2));
        restofpath = shortestpath(g,endpoints(1),endpoints(2));
        
        [~,restofpath_edgeinds] = ismember(sort([restofpath(1:end-1); restofpath(2:end)]',2), node.e,'rows');
        cutfaces = node.e2F(restofpath_edgeinds);
    end
    
end



