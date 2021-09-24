% clear all; close all; 
% files = dir('outputm/mesh*');
% nfiles = numel(files);
%  for i = 1:nfiles
%      file = files(i);
%      read the face list
%      F = readmatrix([file.folder '\' file.name])+1;
% nverts = max(F(:));
% nfaces = size(F,1);
%F = readmatrix('outputm/meshV06_I00.txt')+1;
%generate the adj matrice7745 2217 8411s
%S = sparse(F);
function V = embedtrimesh(F)
% if nargin == 0
%     filename = 'C:\Users\judy8\Documents\GitHub\special-hex\matlab_decompose_singularities\outputm\meshV05_I00.txt';
% end   
% F = readmatrix(filename)+1;
A = adjacency_matrix(F);
%find the shortest path distance
D = graphallshortestpaths(A);

%find the embedding. center and scale it.
Dp = D+randn(size(D))*.2; % add noise. too much symmetry results in colocated points in embedding. not that it matters too much because i'll have to post process the trimesh anyway :'(
Dp=(Dp+Dp')/2; Dp(1:(size(Dp,1)+1):end)=0; % symmetrize and zero diag.
V = mdscale(Dp,3); 
V=V-mean(V);
V=V./vecnorm(V,2,2);

% V = smoothSphereTriEmbedding(V,F);

%{
figure; axis equal; hold all; rotate3d on;
patch('Faces', F, 'Vertices', V, 'facecolor',  'blue', 'facealpha', 0.1);
scatter3(V(:,1), V(:,2), V(:,3), 'k', 'filled');
scatter3(0, 0, 0, 'r', 'filled');
%}
end