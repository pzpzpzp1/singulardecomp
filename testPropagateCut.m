clear all; close all;  
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/bunny.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/Lpadded.vtk'; lfac = 500; saveres = 1;
%         file_name = 'meshes/Cblock.vtk'; lfac = 500; saveres = 1;
        file_name = 'meshes/Cpadded.vtk'; lfac = 500; saveres = 1;
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
V=V0;H=H0;
data = processhmesh(V,H,visualize); title('Input mesh');
mesh.points = V; mesh.cells = H;
Vpresmooth={};
hexSheetInds={};
cuts={};
datas{1} = data;
iter = 2;

%% test propagate singular cut
singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
nodes={}; for i=1:numel(singularNodes)
    nodes{i} = getNode(data, singularNodes(i));
end
interiorsingularnodedegrees = sum(data.E2V(data.isSingularEdge,singularNodes),1);
% selind = randi(numel(singularNodes))
%         selind = selinds(iter); 
selind=3;
node_ind = singularNodes(selind);
node = getNode(data, node_ind);
% select sheet to insert in node
% cutseed = selectSplit(data,node);  % general
cutseed = [26    30    29    25]';
fh=visualizeHmeshData(data,figure,.5); ptc = patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','c')

%% test propagate sheet chord
for i=1:numel(cutseed)
    facei = cutseed(i);
    twoEdges = intersect(data.F2Earray(facei,:), node.E_adj);
    for j=1:2
        edgej = twoEdges(j);
        isBlockedEdge = false(data.nE,1);
        [logicalchord, facechord, edgechord] = extendRegularFaceChord(data, facei, edgej,isBlockedEdge);

        patch('faces',data.F(logicalchord,:),'vertices',data.V,'facecolor','c','facealpha',.5)
    end
end

%% test propagate regular sheet
fh=visualizeHmeshData(data,figure,.5); 
for i=1:numel(cutseed)
    patch('faces',data.F(cutseed(i),:),'vertices',data.V,'facecolor','g','facealpha',.5)
    cutfaces = extendRegularFaces(data, cutseed(i));
    patch('faces',data.F(cutfaces,:),'vertices',data.V,'facecolor','y','facealpha',.5)
end

%% propagate sheet
cut = propagateCut(data,node,cutseed);
% data = processhmesh(V,H,visualize); title(num2str(iter));
ptc = patch('vertices',data.V,'faces',data.F(cut,:),'facecolor','c')
% insert sheet
[V,H,hexSheetInds,VnewPreperturb]=sheetinsertion(data, cut);
data = processhmesh(V,H,visualize); title(num2str(iter));
fh=visualizeHmeshData(data,figure,.5); 



%% mostly junk
selinds = [5 1 1 1 1 1 1 1 1 1 1 1 1 1];
while any(data.isSingularNode & ~data.isBoundaryVertex)
    %% choose random node to simplify
    singularNodes = find(data.isSingularNode & ~data.isBoundaryVertex);
    nodes={}; for i=1:numel(singularNodes)
        nodes{i} = getNode(data, singularNodes(i));
    end
    interiorsingularnodedegrees = sum(data.E2V(data.isSingularEdge,singularNodes),1);
    selind = randi(numel(singularNodes))
%         selind = selinds(iter); 
    node_ind = singularNodes(selind);

    %% build map from singular node to T(S2)
    node = getNode(data, node_ind);

    %% select sheet to insert in node
    cutseed = selectSplit(data,node);
    % ptc = patch('vertices',data.V,'faces',data.F(cutseed,:),'facecolor','c')

    %% propagate sheet
    cut = false(data.nF,1); cut(cutseed)=true;
    cut = propagateCut(data,node,cutseed);
    
    
    
    
    cuts{iter-1} = cut;
    % data = processhmesh(V,H,visualize); title(num2str(iter));
    ptc = patch('vertices',data.V,'faces',data.F(cut,:),'facecolor','c')

    %% insert sheet
    [V,H,hexSheetInds{iter-1},VnewPreperturb]=sheetinsertion(data, cut);
    Vpresmooth{iter-1} = VnewPreperturb;

    %% geometric simplification
    preLapSmooth=0; uniformrot = 0; if ~exist('fixb','var'); fixb = 0; end;
    V = smoothenhmesh(V,H, [],visualize, preLapSmooth, [], lfac, 2, 0, uniformrot );
    V = smoothenhmesh(V,H, [],visualize, preLapSmooth, [], lfac, 4, 0, uniformrot );
  
%{
      mesh.points = V; mesh.cells = H;
      new_mesh = repair_mesh(mesh)
      V=mesh.points; H=mesh.cells;
%}
    %% recompute data
    data = processhmesh(V,H,visualize); title(num2str(iter));
    datas{iter} = data;
    iter = iter+1;
end

decompdata.datas = datas;
decompdata.Vpresmooth = Vpresmooth;
decompdata.hexSheetInds = hexSheetInds;
decompdata.cuts = cuts;
close all;

%{
V=out.Vs{end}; H=out.Hs{end};
data = processhmesh(V,H,1);

figure; hold all; axis equal off; rotate3d on;
patch('vertices',data.V,'faces',data.F(:,:),'facecolor','green','facealpha',.1,'edgealpha',0)
patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]))
patch('vertices',V,'faces',E(isSingularEdge,[1 2 1]),'linewidth',3,'edgecolor','blue')
scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled')
       
%}


