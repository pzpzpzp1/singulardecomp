function QM = getQMfromCut(data,cut)
    cutfaceinds = find(cut);
    QM.F = data.F(cutfaceinds,:);
    QM.V = data.V;
    QM.FE = reshape(permute(reshape(QM.F(:,[1 2, 2 3, 3 4, 4 1]),[],2,4),[1 3 2]),[],4,2); % F*4, 2. faceedges
    FEflat = reshape(QM.FE, [], 2);
    [~, ia, ic] = unique(sort(reshape(QM.FE,[],2),2),'rows');
    QM.E = FEflat(ia,:);
    QM.F2Earray = reshape(ic,[],4);
    QM.nF=size(QM.F,1);
    QM.nV=size(QM.V,1);
    QM.nE=size(QM.E,1);
    QM.F2E = sparse(repmat(1:size(QM.F,1),1,4),QM.F2Earray,QM.F2Earray*0+1,QM.nF,QM.nE);
    assert(all(ismember(unique(sum(QM.F2E,1)),[1 2]))); % ensures cut surface is manifold
    QM.isBoundaryEdge = sum(QM.F2E,1)==1;
    
    [~,beind]=ismember(sort(QM.E(QM.isBoundaryEdge,:),2), sort(data.E,2), 'rows');
    % assert(all(data.isBoundaryEdge(beind))); % ensures boundary of cut surface is on the boundary of the hex mesh.
    QM.HmeshCutBoundaryEdgeInds = beind;
end