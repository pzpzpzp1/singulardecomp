function decompdata = decompose_hmesh(V0,H0,visualize,saveres)
    close all;  
    if nargin==0
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/carter-hex.vtk'; lfac = 500; saveres = 1; % hella dense! 
%         file_name = 'meshes/blood_vessel.mesh'; lfac = 500; saveres = 1; % unhandled val 6 singular node
%         file_name = 'meshes/knob.mesh'; lfac = 500; saveres = 1;  % creates self intersections
        file_name = 'meshes/cactus.vtk'; lfac = 500; saveres = 1; % success! had to fix some nodes and add bypass though.
        
%         file_name = 'meshes/trebol.mesh'; lfac = 500; saveres = 1;  % unhandled val 6 singular node
%         file_name = 'meshes/Lpadded.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/Cblock.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/Cpadded.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
%         file_name = 'meshes/rockarm.vtk'; % NONMANIFOLD BOUNDARY. DONT USE.
%         file_name = 'meshes/hex_sphere.vtk'; lfac = 500; saveres=1;
%         file_name = 'results/tet_split_notsplit/tetnotsplit.vtk'; lfac = 500; saveres=1;
%         file_name = 'meshes/unit.vtk';
%         file_name = 'meshes/hex_tetrahedron.vtk';
%         file_name = 'meshes/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/tetpadded.vtk';
%         file_name = 'meshes/sing1.vtk'; % tet
%         file_name = 'meshes/sing2.vtk'; % tri-prism padded
%         file_name = 'meshes/sing3.vtk'; % 
%         file_name = 'meshes/kitten.mesh';
%          file_name = 'extractSingularVertsFromTri/hmeshSings/sing400.vtk'; % two val 3's
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing222.vtk'; % same as sing3. don't need.
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing133.vtk'; % 1-3 turning point and val 5.
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing044.vtk'; % two val 5's 
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing036.vtk'; % three val 5's
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing206.vtk'; % two val 5's + one val 3. % this is the one where things start to go bad. needs a 90 deg turn in the cut in a regular area.
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing028.vtk'; % four val 5's
%         file_name = 'extractSingularVertsFromTri/hmeshSings/sing0012.vtk'; % six val 5's
        
        [dname,fname,ext]=fileparts(file_name);
        if strcmp(ext,'.vtk')
            mesh = load_vtk(file_name);
        elseif strcmp(ext,'.mesh')
            mesh = ImportHexMesh(file_name);
        else
            error
        end
        V0 = mesh.points;
        H0 = mesh.cells;
        visualize = 1;
        if ~exist('saveres','var')
            saveres=0;
        end
    end
    
    %% load mesh. And preprocess with padding. extract trimesh boundary
    [dname,fname,ext]=fileparts(file_name);
    V=V0;H=H0;
    data = processhmesh(V,H,0);
    if (any(data.isSingularNode & data.isBoundaryVertex) && false) ||...
            contains(file_name,'unit.vtk') || contains(file_name,'sing2.vtk') || contains(file_name,'bunny.vtk') || ...
            contains(file_name,'tetnotsplit.vtk')
        [V,H] = padhmesh(V,H); % [V,H] = padhmesh(V,H);
        % [V,H] = hex1to8(V,H); [V,H] = hex1to8(V,H);
        % V = smoothenhmesh(V,H,[],visualize);
        % V = smoothenhmesh(V,H, [],visualize, 1, [], 100, 2, 0, 0);
        % mesh.points = V; mesh.cells = H;
        % save_vtk(mesh, 'test.vtk')
    end
    data = processhmesh(V,H,visualize); title('Input mesh');
    % boundary triangle mesh for projection
    trimesh0.Vertices = data.V;
    trimesh0.Faces = [data.F(data.isBoundaryFace,[1 2 3]);  data.F(data.isBoundaryFace,[3 4 1])];
    [trimesh0.Vertices, trimesh0.Faces] = minimizeMesh(trimesh0.Vertices, trimesh0.Faces);
    
    % save starting mesh as index 0
    saveseed = randi(99); % saveseed=64;
    outdname = sprintf('results/%s_%d',fname,saveseed); 
    outname = sprintf('results/%s_%d/hmesh_1.vtk',fname,saveseed);
    mesh.points = V; mesh.cells = H;
    if saveres
        mkdir(outdname);
        save_vtk(mesh, outname)
    end
    
    Vpresmooth={};
    hexSheetInds={};
    cuts={};
    
    %% Begin decomposition
    datas{1} = data;
    iter = 2;
%     selinds = [1 1 2 3 1 3  3   1 1 1 1 1 1 1 1 1 1 1]; % 0 0 12
    selinds = [1 6 2 5 1 1 5 5 6 5    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % cactus
%     selinds = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    selindrec = [];
    notskip = true;
    bypass = false;
    while any(data.isSingularNode & ~data.isBoundaryVertex)
        %% choose random node to simplify
        singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
        nodes={}; for i=1:numel(singularNodes)
            nodes{i} = getNode(data, singularNodes(i));
        end
        interiorsingularnodedegrees = sum(data.E2V(data.isSingularEdge,singularNodes),1);
%         selind = randi(numel(singularNodes))
        selind = selinds(iter); 
        selindrec = [selindrec selind];
        node_ind = singularNodes(selind);
        
        if iter >=11
            bypass = true;
        end
        
        %% build map from singular node to T(S2)
        node = getNode(data, node_ind);
        
        %% select sheet to insert in node
        cutseed = selectSplit(data,node,false(data.nF,1),false(data.nF,1),bypass);
        % ptc = patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','c')
        
        %% propagate sheet
        if strcmp(fname,'sing0012')
            if iter == 5
                cutseed = [65    69   128   166   168   170   172   173   174   175 169];
                notskip = false;
            elseif iter == 6
                cutseed = [21    24    54    55   109   125   160   165   166   169   170   201   202 164]
                notskip = false;
            elseif iter == 8
                cutseed = [274 275 38    39    82    83   124   143   224   225   235   236   266   269   270   272   273   276   281   282   283   284];
                notskip = false;
            end
        end
        cut = false(data.nF,1); cut(cutseed)=true;
            
        if notskip
            cut = propagateCut(data,node,cutseed);
        end
        notskip = true;
        
        cuts{iter-1} = cut;
        % data = processhmesh(V,H,visualize); title(num2str(iter));
        ptc = patch('vertices',data.V,'faces',data.F(cut,:),'facecolor','c')
        
        %% insert sheet
        [V,H,hexSheetInds{iter-1},VnewPreperturb]=sheetinsertion(data, cut);
        Vpresmooth{iter-1} = VnewPreperturb;
        
        %% geometric simplification
        lfac=1000;
        preLapSmooth=0; uniformrot = 0; if ~exist('fixb','var'); fixb = 0; end;
       V = smoothenhmesh(V,H, [],visualize, preLapSmooth, [], lfac, 2, 0, uniformrot );
       V = smoothenhmesh(V,H, [],visualize, preLapSmooth, [], lfac, 4, 0, uniformrot );
        % V = smoothenhmesh(V,H, [],visualize, 0, [], 0, 4, 0, uniformrot );

%{
          mesh.points = V; mesh.cells = H;
          new_mesh = repair_mesh(mesh)
          V=mesh.points; H=mesh.cells;
%}
        %% recompute data
        data = processhmesh(V,H,visualize); title(num2str(iter));
        datas{iter} = data;
        outname = sprintf('results/%s_%d/hmesh_%d.vtk',fname,saveseed,iter); mesh.points = V; mesh.cells = H; 
        if saveres; save_vtk(mesh, outname); end;
        iter = iter+1;
    end
    
    decompdata.datas = datas;
    decompdata.Vpresmooth = Vpresmooth;
    decompdata.hexSheetInds = hexSheetInds;
    decompdata.cuts = cuts;
    
    outname = sprintf('results/%s_%d/decompdata.mat',fname,saveseed);
    if saveres; save(outname,'decompdata'); end;
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


