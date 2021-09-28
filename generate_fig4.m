clear all; close all;
% file_name = 'extractSingularVertsFromTri/hmeshSings/sing206.vtk'; 
% file_name = 'extractSingularVertsFromTri/hmeshSings/sing133.vtk'; 
file_name = 'extractSingularVertsFromTri/hmeshSings/sing028.vtk'; 
% file_name = 'extractSingularVertsFromTri/hmeshSings/sing400.vtk'; 
mesh = load_vtk(file_name);
V0 = mesh.points;
H0 = mesh.cells;
[dname,fname,ext]=fileparts(file_name);
V=V0;H=H0;
[V, out] = smoothenhmesh(V, H, [], 0, 0, [], 0, 2, 0);
[V, out] = smoothenhmesh(V, H, [], 0, 0, [], 0, 4, 0);
[V, out] = smoothenhmesh(V, H, [], 0, 0, [], 0, 8, 0);

data = processhmesh(V,H,0);

fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
lwfac=1;
hold all; axis equal off; rotate3d on; axis image vis3d;
V=data.V; F=data.F; E=data.E; 
isSingularNode = data.isSingularNode;
isBoundaryEdge = data.isBoundaryEdge;
isBoundaryFace = data.isBoundaryFace;
isSingularEdge = data.isSingularEdge;
isIntSingularNode = data.isSingularNode & ~data.isBoundaryVertex;
IS3 = data.efdeg==3 & data.isSingularEdge & ~data.isBoundaryEdge;
IS5 = data.efdeg==5 & data.isSingularEdge & ~data.isBoundaryEdge;
IS_ = (data.efdeg~=3 & data.efdeg~=5) & data.isSingularEdge & ~data.isBoundaryEdge;
purple = [62.7, 12.5, 94.1]/100;
guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','blue','facealpha',.1,'edgealpha',.2);
% guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',3*lwfac,'edgecolor','r');
% guielems{4} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',3*lwfac,'edgecolor','b');
% guielems{4} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',3*lwfac,'edgecolor','g');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),200*lwfac,'k','filled');
c = V(isIntSingularNode,:);
r=.4;
[xx,yy,zz]=sphere(50);
surf(r*xx+c(1),r*yy+c(2),r*zz+c(3),'facecolor','yellow','edgealpha',0,'facealpha',1)

node = getNode(data, find(isIntSingularNode));
cutfaces = false(data.nF,1);
cutfaces(selectSplit(data, node, [], []))=true; 
patch('vertices',V,'faces',F(cutfaces,:),'facecolor','r','facealpha',.5,'edgealpha',0);
% guielems{1} = patch('vertices',V,'faces',F(~isBoundaryFace & ~cutfaces,:),'facecolor',purple,'facealpha',.2,'edgealpha',0,'linewidth',1);
v0ind = find(isIntSingularNode);
v0 = data.V(v0ind,:);
for i=1:size(node.e,1)
    Einds = node.v2E(node.e(i,:))
    Vinds = setdiff(data.E(Einds,:),v0ind);
    vs = data.V(Vinds,:);
    
    xi = linspace(vs(1,1),vs(2,1),10);
    yi = linspace(vs(1,2),vs(2,2),10);
    zi = linspace(vs(1,3),vs(2,3),10);
    
    xi = (xi-v0(1));
    yi = (yi-v0(2));
    zi = (zi-v0(3));
    xyzn = vecnorm([xi;yi;zi],2,1);
    xi=xi./xyzn;
    yi=yi./xyzn;
    zi=zi./xyzn;
    
    
    xi = (xi)*r + v0(1);
    yi = (yi)*r + v0(2);
    zi = (zi)*r + v0(3);
    plot3(xi,yi,zi,'k','linewidth',2)
    
    if ismember(node.e2F(i), find(cutfaces))
        plot3(xi,yi,zi,'r','linewidth',5)
    end
    
    % plot3(vs(:,1),vs(:,2),vs(:,3),'k','linewidth',2)
end
% campos([  -39.821      -14.611     -0.22972]);
campos([-31.441       3.4051       4.9599])
campos([-28.989      -7.9055       11.013])
pic=getframe(gcf);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig4/fig4a.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);


%% do split
visualize=0;
[V,H,inds]=sheetinsertion(data, cutfaces);
V = smoothenhmesh(V,H,[],1, 1, [], 100, 2, 0);
data = processhmesh(V,H,visualize); 


% visualize split nodes
fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
lwfac=1;
hold all; axis equal off; rotate3d on; axis image vis3d; 
V=data.V; F=data.F; E=data.E; 
isSingularNode = data.isSingularNode;
isBoundaryEdge = data.isBoundaryEdge;
isBoundaryFace = data.isBoundaryFace;
isSingularEdge = data.isSingularEdge;
isIntSingularNode = data.isSingularNode & ~data.isBoundaryVertex;
IS3 = data.efdeg==3 & data.isSingularEdge & ~data.isBoundaryEdge;
IS5 = data.efdeg==5 & data.isSingularEdge & ~data.isBoundaryEdge;
IS_ = (data.efdeg~=3 & data.efdeg~=5) & data.isSingularEdge & ~data.isBoundaryEdge;
purple = [62.7, 12.5, 94.1]/100;
guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','blue','facealpha',.1,'edgealpha',.2);
% guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',1*lwfac,'edgecolor','k');
% guielems{4} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',1*lwfac,'edgecolor','k');
% guielems{4} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',1*lwfac,'edgecolor','k');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),200*lwfac,'k','filled');


sinds = find(isIntSingularNode);
sind1 = sinds(1); sind2 = sinds(2);
c = V(sind1,:);[xx,yy,zz]=sphere(50);surf(r*xx+c(1),r*yy+c(2),r*zz+c(3),'facecolor','yellow','edgealpha',0,'facealpha',1)
c = V(sind2,:);[xx,yy,zz]=sphere(50);surf(r*xx+c(1),r*yy+c(2),r*zz+c(3),'facecolor','yellow','edgealpha',0,'facealpha',1)
newhexes = inds;
newfaces = unique(data.H2Farray(inds,:));
oldfaces = unique(data.H2Farray(setdiff(1:data.nH,inds),:));
newwithoutold = setdiff(newfaces,[oldfaces; find(data.isBoundaryFace)]);
% newandold = intersect(newfaces,oldfaces);
patch('vertices',V,'faces',data.F(newwithoutold,:),'facecolor','blue','facealpha',.1,'edgealpha',1,'linewidth',2);
% patch('vertices',V,'faces',data.F(newandold,:),'facecolor','red','facealpha',.1,'edgealpha',1);

node = getNode(data, sind1);
v0ind = sind1;
v0 = data.V(v0ind,:);
for i=1:size(node.e,1)
    Einds = node.v2E(node.e(i,:))
    Vinds = setdiff(data.E(Einds,:),v0ind);
    vs = data.V(Vinds,:);    
    xi = linspace(vs(1,1),vs(2,1),10);
    yi = linspace(vs(1,2),vs(2,2),10);
    zi = linspace(vs(1,3),vs(2,3),10);
    xi = (xi-v0(1));
    yi = (yi-v0(2));
    zi = (zi-v0(3));
    xyzn = vecnorm([xi;yi;zi],2,1);
    xi=xi./xyzn;
    yi=yi./xyzn;
    zi=zi./xyzn;
    xi = (xi)*r + v0(1);
    yi = (yi)*r + v0(2);
    zi = (zi)*r + v0(3);
    plot3(xi,yi,zi,'k','linewidth',2)
    if ismember(node.e2F(i), newfaces)
        plot3(xi,yi,zi,'b','linewidth',5)
    end
end
node = getNode(data, sind2);
v0ind = sind2;
v0 = data.V(v0ind,:);
for i=1:size(node.e,1)
    Einds = node.v2E(node.e(i,:))
    Vinds = setdiff(data.E(Einds,:),v0ind);
    vs = data.V(Vinds,:);    
    xi = linspace(vs(1,1),vs(2,1),10);
    yi = linspace(vs(1,2),vs(2,2),10);
    zi = linspace(vs(1,3),vs(2,3),10);
    xi = (xi-v0(1));
    yi = (yi-v0(2));
    zi = (zi-v0(3));
    xyzn = vecnorm([xi;yi;zi],2,1);
    xi=xi./xyzn;
    yi=yi./xyzn;
    zi=zi./xyzn;
    xi = (xi)*r + v0(1);
    yi = (yi)*r + v0(2);
    zi = (zi)*r + v0(3);
    plot3(xi,yi,zi,'k','linewidth',2)
    if ismember(node.e2F(i), newfaces)
        plot3(xi,yi,zi,'b','linewidth',5)
    end
end
campos([-31.441       3.4051       4.9599])
campos([-28.989      -7.9055       11.013])

% campos([  -39.821      -14.611     -0.22972]);
pic=getframe(gcf);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig4/fig4b.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);
