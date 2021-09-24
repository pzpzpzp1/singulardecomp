% outputs cut faces. 
function cutfaces = selectSplit(data, node, alreadycutfaces, forbiddenfaces)
    assert(~data.isBoundaryVertex(node.ind))
    if nargin < 3
        alreadycutfaces = false(data.nF,1);
        forbiddenfaces = false(data.nF,1);
    end
    
    if all(~alreadycutfaces)
        assert(all(~forbiddenfaces)); % no existing cuts means nothing should be off limits to cut.
        %% tetrahedral type
        if numel(node.v)==4
            %% (4,0,0) type is a 3 and 3 joined perpendicularly.
            vpath = [1 2 3 4];
            [allf, relevantEinds] = ismember(sort([vpath; circshift(vpath,-1)]',2), sort(node.e,2), 'rows');
            assert(all(allf));
            cutfaces = node.e2F(relevantEinds);
            return;
        end

        signature = accumarray(data.efdeg(node.v2E),1); if numel(signature) < 6; signature(6)=0; end;
        assert(all(signature(1:2)==0)); % interior singular edges must have 3 or more faces adjacent.
        if sum(signature(6:end))==0
            % only 3,4,5 singularities. must fit into the list of 10 interior singular node types.
            if all(signature(3:5)==[2 2 2]')
                %% (2,2,2) type is a 3 and 5 joined perpendicularly.
                % they can be pulled apart by picking all faces adjacent to one deg 3 and one deg 5 singular edge
                vdegs = data.efdeg(node.v2E);
                edegs = sort(vdegs(node.e),2);
                cutfaces = node.e2F(ismember(sort(edegs,2),[3 5],'rows'));
                return;
            end

            if all(signature(3:5)==[1 3 3]')
                %% (1,3,3) type is a 3-5 joined with a val 5.
                % the cut traverses degree 5,5,5,4 vertices in this
                % triangulation to pull a valence 5 curve off of it.
                vdegs = data.efdeg(node.v2E);
                tdegs = sort(vdegs(node.t),2);
                edegs = sort(vdegs(node.e),2);
                e55inds = find(ismember(edegs, [5 5],'rows'));
                cutfaces1 = node.e2F(e55inds(1:2)); 
                
                midp = intersect(node.e(e55inds(1),:), node.e(e55inds(2),:));
                aa = setdiff(node.e(e55inds(1),:),midp);
                bb = setdiff(node.e(e55inds(2),:),midp);
                candt = find(ismember(tdegs,[4 5 5],'rows'));
                aabb = setdiff(node.t(candt(sum(ismember(node.t(candt,:), [aa bb]),2)==2),:),[aa bb]);
                
                cutfaces2 = node.e2F(ismember(sort(node.e,2),sort([aabb aa; aabb bb],2),'rows'));
                cutfaces = [cutfaces1;cutfaces2];
                return;
            end
            
            if all(signature(3:5)==[0 4 4]')
                %% (0,4,4) type is a 5-5 joined perpendicularly
                vdegs = data.efdeg(node.v2E);
                edegs = sort(vdegs(node.e),2);
                
                cutfaces = node.e2F(ismember(sort(edegs,2),[5 5],'rows'));
                return;
            end
            
            if all(signature(3:5)==[0 3 6]')
                %% (0,3,6) type is two (0,4,4)s joined at a valence 5.
                vdegs = data.efdeg(node.v2E);
                edegs = sort(vdegs(node.e),2);
                
                % get a starting val 5 v and a val 4 v adjacent to that
                % first val 5 v.
                startv = find(vdegs==5,1);
                adjv = setdiff(node.e(sum(node.e==startv,2)~=0,:),startv);
                secondv = adjv(find(vdegs(adjv)==4,1));
                % get all triangles adj to these two.
                side1 = sum(ismember(node.t, [startv, secondv]),2)~=0;
                
                % get boundary edges of side1 tris
                side1edges = reshape(permute(reshape(node.t(side1,[1 2, 2 3, 3 1]),[],2,3),[1 3 2]),[],2);
                [side1edges_unique,ia,ic] = unique(sort(side1edges,2),'rows');
                boundarye = side1edges_unique(find(accumarray(ic,1)==1),:);
                [~,boundaryeinds] = ismember(boundarye, sort(node.e,2), 'rows')
                node.e2F(boundaryeinds);
                cutfaces = node.e2F(boundaryeinds);
                
                return;
            end
            
            if all(signature(3:5)==[2 0 6]')
                %% (2 0 6) = (1 3 3) + (1 3 3).
                vdegs = data.efdeg(node.v2E);
                edegs = sort(vdegs(node.e),2);
                
                % get a starting val 5 v and a val 3 v adjacent to that
                % first val 5 v.
                startv = find(vdegs==5,1);
                adjv = setdiff(node.e(sum(node.e==startv,2)~=0,:),startv);
                secondv = adjv(find(vdegs(adjv)==3,1));
                % get all triangles adj to these two.
                side1 = sum(ismember(node.t, [startv, secondv]),2)~=0;
                
                % get boundary edges of side1 tris
                side1edges = reshape(permute(reshape(node.t(side1,[1 2, 2 3, 3 1]),[],2,3),[1 3 2]),[],2);
                [side1edges_unique,ia,ic] = unique(sort(side1edges,2),'rows');
                boundarye = side1edges_unique(find(accumarray(ic,1)==1),:);
                [~,boundaryeinds] = ismember(boundarye, sort(node.e,2), 'rows')
                node.e2F(boundaryeinds);
                cutfaces = node.e2F(boundaryeinds);
                
                return;
            end
            
            if all(signature(3:5)==[0 2 8]')
                %% (0, 2, 8) = (0, 3, 6) + (0, 4, 4)
                vdegs = data.efdeg(node.v2E);
                edegs = sort(vdegs(node.e),2);
                
                % get a starting val 4 v and a val 5 v adjacent to that
                % first val 4 v.
                startv = find(vdegs==4,1);
                adjv = setdiff(node.e(sum(node.e==startv,2)~=0,:),startv);
                secondv = adjv(find(vdegs(adjv)==5,1));
                % get all triangles adj to these two.
                side1 = sum(ismember(node.t, [startv, secondv]),2)~=0;
                
                % get boundary edges of side1 tris
                side1edges = reshape(permute(reshape(node.t(side1,[1 2, 2 3, 3 1]),[],2,3),[1 3 2]),[],2);
                [side1edges_unique,ia,ic] = unique(sort(side1edges,2),'rows');
                boundarye = side1edges_unique(find(accumarray(ic,1)==1),:);
                [~,boundaryeinds] = ismember(boundarye, sort(node.e,2), 'rows')
                node.e2F(boundaryeinds);
                cutfaces = node.e2F(boundaryeinds);
                
                return;
            end
        end
        
        %% general case. find a cut based on the inductive simplification procedure. no conflicts to worry about.
        error('unhandled');
    else
        %% finish the existing cut. while avoiding forbidden cuts.
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
            % actually it probably works fine in that setting. but the forbidden faces should prevent this from happening.
            cutfaces = [];
            return;
        end
        
        % use shortest path to close up the partial path loop. exclude forbidden cuts.
        localforbid = forbiddenfaces(node.F_adj);
        
        g = graph();
        g = addedge(g, node.e(~partialpath & ~localforbid,1), node.e(~partialpath & ~localforbid,2));
        restofpath = shortestpath(g,endpoints(1),endpoints(2));
        
        if numel(restofpath)==0
            %{
            After removing already cut faces, and forbidden faces, there were no faces available left to bridge the endpoints. 
            Not sure if this can happen. perhaps if the singular graph wraps around itself in some way.
            %}
            error('no path found.');
        end
        
        [~,restofpath_edgeinds] = ismember(sort([restofpath(1:end-1); restofpath(2:end)]',2), node.e,'rows');
        cutfaces = node.e2F(restofpath_edgeinds);
    end
    
end



