% loads hex mesh, and tries to maximize its min scaled jacobian
clear all; close all;

%% canon singularities
% file_name = 'results/sing1_59/hmesh_1.vtk'; 
% file_name = 'results/sing1_59/hmesh_2.vtk'; 
% file_name = 'results/sing3_92/hmesh_1.vtk'; 
% file_name = 'results/sing3_92/hmesh_2.vtk'; 
% file_name = 'results/sing133_42/hmesh_1.vtk'; 
% file_name = 'results/sing133_42/hmesh_3.vtk'; 
% file_name = 'results/sing044_73/hmesh_1.vtk'; 
% file_name = 'results/sing044_73/hmesh_2.vtk'; 
% file_name = 'results/sing036_32/hmesh_1.vtk'; 
% file_name = 'results/sing036_32/hmesh_4.vtk'; 
% file_name = 'results/sing206_69/hmesh_1.vtk'; 
% file_name = 'results/sing206_69/hmesh_8.vtk'; 
% file_name = 'results/sing028_39/hmesh_1.vtk'; 
% file_name = 'results/sing028_39/hmesh_6.vtk'; 
% file_name = 'results/sing0012_65/hmesh_1.vtk'; 
% file_name = 'results/sing0012_65/hmesh_8.vtk'; 

%% not canon singularities
% file_name = 'results/sing2_88/hmesh_1.vtk'; % triprism padded
% file_name = 'results/sing2_88/hmesh_4.vtk'; 
% file_name = 'results/hex_ellipsoid_coarse_78/hmesh_1.vtk'; % ellipsoid
% file_name = 'results/hex_ellipsoid_coarse_78/hmesh_5.vtk';
% file_name = 'results/tet_split_notsplit/tetnotsplit.vtk'; % tet. 1to8 twice
file_name = 'results/tetnotsplit_64/hmesh_1.vtk'; % tet. padded
% file_name = 'results/tet_split_notsplit/tetsplit.vtk'; 
% file_name = 'results/tetpadded_16/hmesh_1.vtk'; skipL = 1; % tet padded
% file_name = 'results/tetpadded_16/hmesh_4.vtk'; 
% file_name = 'results/unit_15/hmesh_1.vtk'; skipL = 1;
% file_name = 'results/unit_15/hmesh_5.vtk'; 
% file_name = 'results/unit_70/hmesh_1.vtk';  skipL = 1;
% file_name = 'results/unit_70/hmesh_5.vtk'; 
% file_name = 'meshes/hex_sphere.vtk'; fixb = 1; skipL=1;
% file_name = 'results/hex_sphere_64/hmesh_5.vtk'; fixb = 0; skipL=1;

mesh = load_vtk(file_name);
V0 = mesh.points;
H = mesh.cells;

if ~exist('skipL','var')
    skipL = 0;
end
if ~exist('fixb','var')
    fixb = 0;
end
if ~exist('unifrot','var')
    unifrot = 0;
end

if ~skipL
    [V, out] = smoothenhmesh(V0, H, [], 1, 1, [], 0, 2,fixb,unifrot);
else
    V=V0;
end

[V, out] = smoothenhmesh(V, H, [], 1, 0, [], 0, 2, fixb,unifrot);
[V, out] = smoothenhmesh(V, H, [], 1, 0, [], 0, 4, fixb,unifrot);
[V, out] = smoothenhmesh(V, H, [], 1, 0, [], 0, 8, fixb,unifrot);
mesh.points = V;

[dname,fname,ext] = fileparts(file_name);
dname(ismember(dname,'/'))='_';

outname = sprintf('bestembedding/%s_%s.vtk',dname,fname);
save_vtk(mesh,outname)    
    