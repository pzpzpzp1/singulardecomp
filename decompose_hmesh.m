function out = decompose_hmesh(V0,H0,visualize)
    close all;  
    if nargin==0
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
        % file_name = 'meshes/rockarm.vtk';
%         file_name = 'meshes/hex_sphere.vtk';
        file_name = 'meshes/unit.vtk';
        % file_name = 'meshes/hex_tetrahedron.vtk';
%         file_name = 'meshes/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/sing1.vtk';
%         file_name = 'meshes/sing2.vtk';
%         file_name = 'meshes/sing3.vtk';
%         file_name = 'results_fmincon/sing3.vtk';
%         file_name = 'meshes/kitten.mesh';
        
        mesh = load_vtk(file_name);
%         mesh = ImportHexMesh(file_name);
        V0 = mesh.points;
        H0 = mesh.cells;
        visualize = 1;
    end
    
    %% load mesh. And preprocess with padding. extract trimesh boundary
    [dname,fname,ext]=fileparts(file_name);
    V=V0;H=H0;
    data = processhmesh(V,H,0);
    if (any(data.isSingularNode & data.isBoundaryVertex) && false) || contains(file_name,'unit.vtk')
        [V,H] = padhmesh(V,H);
    end
    data = processhmesh(V,H,visualize);
    % boundary triangle mesh for projection
    trimesh0.Vertices = data.V;
    trimesh0.Faces = [data.F(data.isBoundaryFace,[1 2 3]);  data.F(data.isBoundaryFace,[3 4 1])];
    [trimesh0.Vertices, trimesh0.Faces] = minimizeMesh(trimesh0.Vertices, trimesh0.Faces);
    % save starting mesh as index 0
    outname = sprintf('results/%s_0.vtk',fname);
    mesh.points = V; mesh.cells = H;
    save_vtk(mesh, outname)
    
    %% Begin decomposition
    Vs{1} = V; Hs{1} = H;
    iter = 2;
    while any(data.isSingularNode & ~data.isBoundaryVertex)
        %% choose random node to simplify
        singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
        selind = randi(numel(singularNodes)); 
%         selind = 1;
        node_ind = singularNodes(selind);
        
        %% build map from singular node to T(S2)
        node = getNode(data, node_ind);
        
        %% select sheet to insert in node
        cutseed = selectSplit(data,node);
        
        %% propagate sheet
        cut = propagateCut(data,node,cutseed);
        patch('vertices',data.V,'faces',data.F(cut,:),'facecolor','c')

        %% insert sheet
        [V,H]=sheetinsertion(data, cut);
        
        %% geometric simplification
        V = smoothenhmesh(V,H,trimesh0,visualize);
%{
          mesh.points = V; mesh.cells = H;
          new_mesh = repair_mesh(mesh)
          V=mesh.points; H=mesh.cells;
%}
        %% recompute data
        data = processhmesh(V,H,visualize); title(num2str(iter));
        Vs{iter} = V; Hs{iter} = H; iter = iter+1;
        outname = sprintf('results/%s_%d.vtk',fname,iter); mesh.points = V; mesh.cells = H; save_vtk(mesh, outname);
    end
    
    out.Vs=Vs;
    out.Hs=Hs;
end
%{
V=out.Vs{end}; H=out.Hs{end};
data = processhmesh(V,H,1);

figure; hold all; axis equal off; rotate3d on;
patch('vertices',data.V,'faces',data.F(:,:),'facecolor','green','facealpha',.1,'edgealpha',0)
patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]))
patch('vertices',V,'faces',E(isSingularEdge,[1 2 1]),'linewidth',3,'edgecolor','blue')
scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled')
       
%}


