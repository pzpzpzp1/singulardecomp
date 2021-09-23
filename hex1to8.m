


function [Vnew, Hnew] = hex1to8(V, H)
    if nargin==0
        file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
%         file_name = 'meshes/rockarm.vtk';
%         file_name = 'meshes/hex_sphere.vtk';
%         file_name = 'meshes/unit.vtk';
        mesh = load_vtk(file_name);
        V = mesh.points;
        H = mesh.cells;
    end
    % convenience
    data = processhmesh(V,H,0);
    % make new vertices
    newHexaVerts = data.cellBarycenters;
    newEdgeVerts = (data.V(data.E(:,1),:)+data.V(data.E(:,2),:))/2;
    newFaceVerts = data.faceBarycenters;
    Vnew = [V;newHexaVerts;newEdgeVerts;newFaceVerts];
    % vertex inds from hef index
    hv_inds = (1:data.nH)' + data.nV;
    ev_inds = (1:data.nE)' + data.nV+data.nH;
    fv_inds = (1:data.nF)' + data.nV+data.nH+data.nE;
    
    Esorted = sort(data.E,2);
    Fsorted = sort(data.F,2);
    indblock = [1 2 4 5 6 3 8 7;...
        2 6 3 1 5 7 4 8;...
        3 2 7 4 1 6 8 5 ;...
        4 1 3 8 5 2 7 6;...
        5 6 1 8 7 2 4 3;...
        6 2 5 7 3 1 8 4;...
        7 6 8 3 2 5 4 1;...
        8 5 4 7 6 1 3 2];
    Hnew = zeros(0,8);
    for i=1:8
        hexvert_vertinds = H(:,indblock(i,:));
        % get local index labelings
        v1 = sort(hexvert_vertinds(:,1),2);
        v2 = sort(hexvert_vertinds(:,2),2);
        v3 = sort(hexvert_vertinds(:,3),2);
        v4 = sort(hexvert_vertinds(:,4),2);
        v12 = sort(hexvert_vertinds(:,[1 2]),2);
        v13 = sort(hexvert_vertinds(:,[1 3]),2);
        v14 = sort(hexvert_vertinds(:,[1 4]),2);
        v1254 = sort(hexvert_vertinds(:,[1 2 5 4]),2);
        v1362 = sort(hexvert_vertinds(:,[1 3 6 2]),2);
        v1473 = sort(hexvert_vertinds(:,[1 4 7 3]),2);
        [allf, e12] = ismember(v12, Esorted,'rows'); assert(all(allf));
        [allf, e13] = ismember(v13, Esorted,'rows'); assert(all(allf));
        [allf, e14] = ismember(v14, Esorted,'rows'); assert(all(allf));
        [allf, f1] = ismember(v1254, Fsorted,'rows'); assert(all(allf));
        [allf, f2] = ismember(v1362, Fsorted,'rows'); assert(all(allf));
        [allf, f3] = ismember(v1473, Fsorted,'rows'); assert(all(allf));
        hc = (1:data.nH)';
        
        Hnew = [Hnew; v1 ev_inds(e12) fv_inds(f2) ev_inds(e13) ev_inds(e14) fv_inds(f1) hv_inds(hc) fv_inds(f3)];
        
    end
    % processhmesh(Vnew,Hnew,1);
    mesh.points = Vnew; mesh.cells = Hnew; save_vtk(mesh, 'test.vtk');
end

function verifySpiralIndex(indblock)
    uvw = [1     0     1;...
        1     1     1;...
        1     1     0;...
        1     0     0;...
        0     0     1;...
        0     1     1;...
        0     1     0;...
        0     0     0];
    figure; hold all; rotate3d on; axis equal off;
    for i=1:8
        text(uvw(i,1),uvw(i,2),uvw(i,3),['<---' num2str(i)])
    end
    
    % need to verify visually.
    for i=1:8
        try; delete(ptc); catch; end;
        spiral = uvw(indblock(i,:),:);
        ptc = plot3(spiral(:,1),spiral(:,2),spiral(:,3),'k','linewidth',2);
        pause;
    end
end