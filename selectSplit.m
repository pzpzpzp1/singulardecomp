% 
function [cutfaces, hexesOneSide] = selectSplit(data, node)
    assert(~data.isBoundaryVertex(node.ind))
    
    %% tetrahedral type
    if numel(node.v)==4
        % handle split on tetrahedral singularity type
        vpath = [1 2 3 4];
        [allf, einds] = ismember(sort([vpath; circshift(vpath,-1)]',2), sort(node.e,2), 'rows');
        assert(all(allf));
        cutfaces = node.e2F(einds);
        
        [~, triinds] = ismember(sort([1 2 3; 1 3 4],2),sort(node.t,2),'rows');
        hexesOneSide = node.t2H(triinds);
        
        return;
    end
    
    signature = accumarray(data.efdeg(node.v2E),1); if numel(signature) < 6; signature(6)=0; end;
    assert(all(signature(1:2)==0)); % interior singular edges must have 3 or more faces adjacent.
    if sum(signature(6:end))==0
        % only 3,4,5 singularities. must fit into the list of 10 interior singular node types.
        if all(signature(3:5)==[2 2 2]')
            %% (2,2,2) type is a 3-5 joined perpendicularly
            
        end
        
        if all(signature(3:5)==[0 4 4]')
            %% (0,4,4) type is a 5-5 joined perpendicularly
            
        end
    end
    
    %% general case. find a cut based on the inductive simplification procedure.
    error('unhandled');
    
end