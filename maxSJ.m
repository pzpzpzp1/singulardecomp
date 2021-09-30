%% this doc keeps track of maximum SJ bounds for some singular node/vertex types.

%% (4,0,0) -> 4/sqrt(3)/3
V = [1 1 1; -1 -1 1; -1 1 -1; 1 -1 -1];
V=V - mean(V);
T = [1 2 3 4];
F = [1 2 3; 2 3 4; 3 4 1; 1 2 4];

v0=V(F(:,1),:);
v1=V(F(:,2),:);
v2=V(F(:,3),:);
normals = cross(v1-v0, v2-v0);
normals = normals./vecnorm(normals,2,2);
SJ = abs(dot(cross(normals(1,:),normals(2,:)),normals(3,:)));

figure; hold all; axis equal; rotate3d on; scatter3(V(:,1),V(:,2),V(:,3))
patch('faces',F,'vertices',V,'facecolor','green','facealpha',.5)
 
%% (0,5,2) -> sin(2*pi/3) = .86603
a = 0; v1 = [cos(a) sin(a) 0];
a = 2*pi/3; v2 = [cos(a) sin(a) 0];
SJ = dot(cross(v1,v2),[0 0 1]);

%% (2,3,0) -> sin(2*pi/5) = .95106
a = 0; v1 = [cos(a) sin(a) 0];
a = 2*pi/5; v2 = [cos(a) sin(a) 0];
SJ = dot(cross(v1,v2),[0 0 1]);

%% (0,0,12) -> sqrt(2*(5+sqrt(5)))/5
% syms t real;
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));
% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20
v1 = v(f(1,1),:);
v2 = v(f(1,2),:);
v3 = v(f(1,3),:);
SJ = dot(cross(v1,v2),v3);

 figure; hold all; axis equal; rotate3d on; scatter3(v(:,1),v(:,2),v(:,3))
patch('faces',f,'vertices',v,'facecolor','green','facealpha',.5)


