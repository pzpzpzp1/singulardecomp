% gets vertex info from hex mesh data
function node = getNode(data, node_ind)
    %% validation. for simplificty, focus on interior singularities. easy to push boundary sings inwards with padding anyway.
    assert(~data.isBoundaryVertex(node_ind));
    
    %% get adjacencies
    E_adj = find(data.E2V(:, node_ind)); % adjacent edges
    H_adj = find(data.H2V(:, node_ind)); % adjacent hexes
    V_adj = setdiff(find(sum(data.E2V(E_adj,:))), node_ind)'; % adjacent verts
    F_adj = find(data.F2V(:, node_ind)); % adjacent faces
    
    %% build sphere triangulation
    v = setdiff(data.E(E_adj,:)',node_ind); % corresponding index in V.
    
    nt = numel(H_adj);
    t = zeros(nt,3);
    for i=1:nt
        h = H_adj(i);
        t(i,:) = intersect(data.H(h,:),V_adj);
    end
    t = bfs_orient(t);
    [allf, t] = ismember(t, v); assert(all(allf(:)));
    
    ne = size(F_adj,1);
    e = zeros(ne,2);
    for i=1:ne
        e(i,:) = intersect(setdiff(data.F(F_adj(i),:),node_ind),V_adj);
    end
    [allf, e] = ismember(e, v); assert(all(allf(:)));
    
    %% map sphere triangulation to hex mesh. lowercase is triangulation. capital is hex mesh.
    v2E = E_adj;
    t2H = H_adj;
    e2F = F_adj;
    
    % figure; hold all; rotate3d on; axis equal off; 
    % patch('vertices',data.V,'faces',T,'facecolor','red','facealpha',.1)
    
    node.v = v;
    node.t = t;
    node.e = e;
    node.v2E = v2E;
    node.t2H = t2H;
    node.e2F = e2F;
    node.E_adj = E_adj;
    node.H_adj = H_adj;
    node.V_adj = V_adj;
    node.F_adj = F_adj;
    node.ind = node_ind;
    
    %% visualize and verify
    %{
    figure; hold all; axis equal off; rotate3d on;
    scatter3(data.V(node_ind,1),data.V(node_ind,2),data.V(node_ind,3),100,'k','filled')
    scatter3(data.V(node.v,1),data.V(node.v,2),data.V(node.v,3),100,'k','filled')
    patch('vertices',data.V(node.v,:),'faces',node.t,'facecolor','y','facealpha',.3)
    patch('vertices',data.V(node.v,:),'faces',node.e(:,[1 2 1]),'facecolor','k','linewidth',2)
    patch('vertices',data.V,'faces',data.E(node.v2E,[1 2 1]),'edgecolor','r','linewidth',3)
    scatter3(data.cellBarycenters(node.t2H,1),data.cellBarycenters(node.t2H,2),data.cellBarycenters(node.t2H,3),'r','filled')
    patch('vertices',data.V,'faces',data.F(node.e2F,:),'facecolor','c')
    %}
    
end