function decompdata = decompose_hmesh(V0,H0,visualize)
    close all;  
    if nargin==0
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/Lpadded.vtk';
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
%         file_name = 'meshes/rockarm.vtk'; % NONMANIFOLD BOUNDARY. DONT USE.
%         file_name = 'meshes/hex_sphere.vtk';
%         file_name = 'meshes/unit.vtk';
%         file_name = 'meshes/hex_tetrahedron.vtk';
%         file_name = 'meshes/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/tetpadded.vtk';
%         file_name = 'meshes/sing1.vtk';
%         file_name = 'meshes/sing2.vtk';
%         file_name = 'meshes/sing3.vtk';
%         file_name = 'meshes/kitten.mesh';
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing222.vtk'; % same as sing3. don't need.
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing133.vtk'; % 1-3 turning point and val 5.
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing044.vtk'; % two val 5's 
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing036.vtk'; % three val 5's
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing206.vtk'; % two val 5's + one val 3
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing028.vtk'; % four val 5's
        file_name = 'extractSingularVertsFromTri/hmeshSings/sing0012.vtk'; % six val 5's
        
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
    if (any(data.isSingularNode & data.isBoundaryVertex) && false) ||...
            contains(file_name,'unit.vtk') || contains(file_name,'sing2.vtk') || contains(file_name,'bunny.vtk'); % || ...
            % contains(file_name,'Lblock.vtk')
        [V,H] = padhmesh(V,H);
        % [V,H] = hex1to8(V,H); [V,H] = hex1to8(V,H);
        % V = smoothenhmesh(V,H,[],visualize);
        % mesh.points = V; mesh.cells = H;
        % save_vtk(mesh, 'test.vtk')
    end
    data = processhmesh(V,H,visualize); title('Input mesh');
    % boundary triangle mesh for projection
    trimesh0.Vertices = data.V;
    trimesh0.Faces = [data.F(data.isBoundaryFace,[1 2 3]);  data.F(data.isBoundaryFace,[3 4 1])];
    [trimesh0.Vertices, trimesh0.Faces] = minimizeMesh(trimesh0.Vertices, trimesh0.Faces);
    
    % save starting mesh as index 0
    saveseed = randi(99);
    outdname = sprintf('results/%s_%d',fname,saveseed); mkdir(outdname);
    outname = sprintf('results/%s_%d/hmesh_1.vtk',fname,saveseed);
    mesh.points = V; mesh.cells = H;
    save_vtk(mesh, outname)
    
    Vpresmooth={};
    hexSheetInds={};
    cuts={};
    
    %% Begin decomposition
    datas{1} = data;
    iter = 2;
    selinds = [];
    while any(data.isSingularNode & ~data.isBoundaryVertex)
        %% choose random node to simplify
        singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
        selind = randi(numel(singularNodes))
        % selind = 1; 
        selinds(end+1) = selind;
        node_ind = singularNodes(selind);
        
        %% build map from singular node to T(S2)
        node = getNode(data, node_ind);
        
        %% select sheet to insert in node
        cutseed = selectSplit(data,node);
        
        %% propagate sheet
        cut = false(data.nF,1); cut(cutseed)=true;
        
        cut = propagateCut(data,node,cutseed);
        cuts{iter-1} = cut;
        ptc = patch('vertices',data.V,'faces',data.F(cut,:),'facecolor','c')
        
        %% insert sheet
        [V,H,hexSheetInds{iter-1}]=sheetinsertion(data, cut);
        Vpresmooth{iter-1} = V;
        
        %% geometric simplification
        V = smoothenhmesh(V,H,trimesh0,visualize, 1, [], 100);
%{
          mesh.points = V; mesh.cells = H;
          new_mesh = repair_mesh(mesh)
          V=mesh.points; H=mesh.cells;
%}
        %% recompute data
        data = processhmesh(V,H,visualize); title(num2str(iter));
        datas{iter} = data;
        outname = sprintf('results/%s_%d/hmesh_%d.vtk',fname,saveseed,iter); mesh.points = V; mesh.cells = H; save_vtk(mesh, outname);
        iter = iter+1;
    end
    
    decompdata.datas = datas;
    decompdata.Vpresmooth = Vpresmooth;
    decompdata.hexSheetInds = hexSheetInds;
    decompdata.cuts = cuts;
    
    outname = sprintf('results/%s_%d/decompdata.mat',fname,saveseed);
    save(outname,'decompdata');
    close all;
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


