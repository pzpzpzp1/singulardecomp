function out = decompose_hmesh(V,H,visualize)
    close all;  
    if nargin==0
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
        % file_name = 'meshes/rockarm.vtk';
        % file_name = 'meshes/hex_sphere.vtk';
        % file_name = 'meshes/hex_tetrahedron.vtk';
        file_name = 'meshes/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/sing1.vtk';
%         file_name = 'meshes/sing2.vtk';
%         file_name = 'meshes/sing3.vtk';
%         file_name = 'results_fmincon/sing3.vtk';
%         file_name = 'meshes/kitten.mesh';
        
        mesh = load_vtk(file_name);
%         mesh = ImportHexMesh(file_name);
        V = mesh.points;
        H = mesh.cells;
        visualize = 1;
    end

    %% load mesh. And preprocess with padding
    data = processhmesh(V,H,0);
    if any(data.isSingularNode & data.isBoundaryVertex) && false
        [V,H] = padhmesh(V,H);
    end
    data = processhmesh(V,H,visualize);
    
    %% Begin decomposition
    Vs{1} = V; Hs{1} = H;
    iter = 2;
    while any(data.isSingularNode)
        %% choose random node to simplify
        singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
        selind = randi(numel(singularNodes));
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
        
        
        %% recompute data
        data = processhmesh(V,H,visualize); title(num2str(iter));
        Vs{iter} = V; Hs{iter} = H; iter = iter+1;
    end
    
    out.Vs=Vs;
    out.Hs=Hs;
end




