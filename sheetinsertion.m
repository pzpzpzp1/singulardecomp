% sheet insertion on hex mesh data. cut is a logical bit per face for
% whether insertion happens there or not.
function [Vnew,Hnew,hexsheetinds,VnewPreperturb]=sheetinsertion(data, cut)
    %% validation to make sure this cut can be made while preserving that the output is a hex mesh
    QM = getQMfromCut(data,cut);
    assert(all(data.isBoundaryEdge(QM.HmeshCutBoundaryEdgeInds))); % ensures boundary of cut surface is on the boundary of the hex mesh.
    
    %% start cutting
    % build new vertices. easy peasy
    oldVind = unique(data.F(cut,:));
    newV = data.V(oldVind,:);
    Vnew = [data.V; newV;];
    newVind = data.nV + (1:size(newV,1));
    
    % get list of all effected hexes in two bins on either side of the cut.
    affectedHexes = find(sum(data.H2V(:,oldVind),2));
    dualEdges = data.F2Harray(~data.isBoundaryFace & ~cut,[1 3]);
    rminds = ~all(ismember(dualEdges, affectedHexes),2);
    dualEdges(rminds,:) = [];
    g = graph; g=addedge(g,dualEdges(:,1),dualEdges(:,2));
    [ia,ib] = conncomp(g);
    seeds = data.F2Harray(find(cut,1),[1 3]);
    bin1 = ia(seeds(1));
    bin2 = ia(seeds(2));
    hexesbin1 = find(ia==bin1);
    hexesbin2 = find(ia==bin2);
    assert(numel(setdiff([hexesbin1,hexesbin2],affectedHexes))==0); % assumes sheet partitions volume into 2 pieces. ONLY 2.
    
    % orient cut surface faces by pulling from data.Fall 
    faceinds = find(cut);
    fhs = data.F2Harray(faceinds,[1 3])';
    ii = find(ismember(fhs,hexesbin1)); 
    hexinds = fhs(ii);
    
    hhff = sub2ind(size(data.H2F6), hexinds, faceinds);
    hf6 = data.H2F6(hhff);
    fallinds = sub2ind([data.nH,6],hexinds,hf6);
    orientedCutFaces = data.Fall(fallinds,:);
    QM.F = orientedCutFaces;
    
    % compute cut surface normal perturbation
    fv1 = QM.V(QM.F(:,1),:);
    fv2 = QM.V(QM.F(:,2),:);
    fv3 = QM.V(QM.F(:,3),:);
    fv4 = QM.V(QM.F(:,4),:);
    v123n = @(v1,v2,v3) cross(v2-v1,v3-v1);
    n1 = v123n(fv1,fv2,fv3);
    n2 = v123n(fv2,fv3,fv4);
    n3 = v123n(fv3,fv4,fv1);
    n4 = v123n(fv4,fv1,fv2);
    n = (n1+n2+n3+n4)/4;
    n = n./vecnorm(n,2,2);
    QM.faceNormals = n;
    xx = accumarray(repmat((1:size(QM.F,1))',4,1),QM.V(QM.F,1))/4;
    yy = accumarray(repmat((1:size(QM.F,1))',4,1),QM.V(QM.F,2))/4;
    zz = accumarray(repmat((1:size(QM.F,1))',4,1),QM.V(QM.F,3))/4;
    QM.faceBarycenters = [xx,yy,zz];
    quiver3(xx,yy,zz,n(:,1),n(:,2),n(:,3))
    
    vnx = accumarray(QM.F(:), repmat(n(:,1),4,1));
    vny = accumarray(QM.F(:), repmat(n(:,2),4,1));
    vnz = accumarray(QM.F(:), repmat(n(:,3),4,1));
    vnd = accumarray(QM.F(:), 1);
    vn = [vnx,vny,vnz]./vnd;
    vn=vn./vecnorm(vn,2,2);
    perturbdir = vn(oldVind,:);
    
    pertmag = median(data.edgelengths)/10;
    VnewPreperturb = Vnew;
    Vnew((data.nV+1):end,:) = Vnew((data.nV+1):end,:) + pertmag*perturbdir;
    
    % Create new connectivity
    vertexmap = (1:data.nV)';
    vertexmap(oldVind) = newVind;
    Hnew = data.H;
    Hnew(hexesbin2,:) = vertexmap(Hnew(hexesbin2,:));
    newH = [QM.F, vertexmap(QM.F)];
    Hnew = [Hnew; newH];
    
    hexsheetinds = (1:size(newH,1))' + data.nH;
    
    vv = data.cellBarycenters(hexesbin1,:); 
    scatter3(vv(:,1),vv(:,2),vv(:,3),'y','filled')
    vv = data.cellBarycenters(hexesbin2,:); 
    scatter3(vv(:,1),vv(:,2),vv(:,3),'r','filled')
    
    %% visualize and verify
    figure; drawnow; datanew = processhmesh(Vnew,Hnew,1);
    for i=1:size(Hnew,1)
        assert(numel(unique(Hnew(i,:)))==8)
    end
    
    
end



