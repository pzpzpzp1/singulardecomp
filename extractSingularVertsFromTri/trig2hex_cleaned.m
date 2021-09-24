function [Vh,H,centerind] = trig2hex_cleaned(V,F)

v0 = [0,0,0]; %reference point
E = edges(F);
E = sort(E,2);
E = unique(E, 'rows'); %edge list
%number of vertices, edges, trig faces
n_v = size(V,1);
n_e = size(E,1);
n_f = size(F,1);
%First, we build the vertex list
Ea = V(E(:,1),:)+ V(E(:,2),:);
Ta = V(F(:,1),:)+ V(F(:,2),:)+V(F(:,3),:);
% for i = 1:n_e
%     Ea = [Ea; V(E(i,1),:)+ V(E(i,2),:)];%average of points on each edge
% end
% for i = 1:n_f
%     Ta = [Ta; V(F(i,1),:)+ V(F(i,2),:)+V(F(i,3),:)];%average of points on each triangular face
% end

Vh = [V; v0; Ea; Ta]; %vertex list for hex

%Then, we build the hex list
%index of v0
centerind = (n_v+1);
%index of v1, v2, v3 are just the index of the face matrix
%index of v12, v23, v13
[~,e_12] = ismember(sort(F(:,[1 2]),2), E, 'rows'); 
[~,e_13] = ismember(sort(F(:,[1 3]),2), E, 'rows'); 
[~,e_23] = ismember(sort(F(:,[2 3]),2), E, 'rows'); 

v_edge = [e_12 e_13 e_23]+n_v+1;

%index of v123
v_123 = (n_v+1+n_e+(1:n_f))';
%hex list
H = [repmat(centerind,n_f,1), F, v_edge, v_123];
%reorder hex list to make it correct
H = H(:,[1 2 5 3 4 6 8 7]);
H = H(:,[5 6 7 8 1 2 3 4]);

% figure; title('original'); axis equal; hold all; rotate3d on;
% patch('Faces', Fh, 'Vertices', Vh, 'facecolor', 'blue', 'facealpha', 0.1);
% scatter3(Vh(:,1), Vh(:,2), Vh(:,3), 'k', 'filled');
% figure; title('MDS');axis equal; hold all; rotate3d on;
% patch('Faces', Fh, 'Vertices', Vh2, 'facecolor', 'blue', 'facealpha', 0.1);
% scatter3(Vh2(:,1), Vh2(:,2), Vh2(:,3), 'k', 'filled');
%tsurf(Fh,Vh); axis equal;
%shading interp;
%colorbar;
end