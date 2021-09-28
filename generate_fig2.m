clear all; close all;
% file_name = 'extractSingularVertsFromTri/hmeshSings/sing206.vtk'; 
% file_name = 'extractSingularVertsFromTri/hmeshSings/sing133.vtk'; 
file_name = 'extractSingularVertsFromTri/hmeshSings/sing028.vtk'; 
mesh = load_vtk(file_name);
V0 = mesh.points;
H0 = mesh.cells;
[dname,fname,ext]=fileparts(file_name);
V=V0;H=H0;


[V, out] = smoothenhmesh(V, H, [], 1, 1, [], 0, 2, 0);
[V, out] = smoothenhmesh(V, H, [], 1, 0, [], 0, 4, 0);
[V, out] = smoothenhmesh(V, H, [], 1, 0, [], 0, 8, 0);

data = processhmesh(V,H,0);
fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
lwfac=1;

hold all; axis equal off; rotate3d on;
    
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
guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','blue','facealpha',.1,'edgealpha',0);
guielems{1} = patch('vertices',V,'faces',F(~isBoundaryFace,:),'facecolor',purple,'facealpha',.2,'edgealpha',.5,'linewidth',2);
guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',3*lwfac,'edgecolor','r');
guielems{4} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',3*lwfac,'edgecolor','b');
guielems{4} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',3*lwfac,'edgecolor','g');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),200*lwfac,'k','filled');

c = V(isIntSingularNode,:);
r=.8;
[xx,yy,zz]=sphere(50);
surf(r*xx+c(1),r*yy+c(2),r*zz+c(3),'facecolor','yellow','edgealpha',0)

colors = {'','','red','black','green','blue'};
singedgeinds = find(~isBoundaryEdge);
for i=1:numel(singedgeinds)
    deg = data.efdeg(singedgeinds(i));
    c = colors{deg};
    vi = data.E(singedgeinds(i),:);
    
    
    v0 = find(isIntSingularNode);
    v1 = setdiff(vi,v0);
    e12 = data.V(v1,:)-data.V(v0,:);
    e12=e12/norm(e12);
    vp = data.V(v0,:) + e12*r;
    
    
    scatter3(vp(1),vp(2),vp(3),100,'k','filled')
    scatter3(vp(1),vp(2),vp(3),50,c,'filled')
    
    
    %vp = data.V(v1,:);
    %scatter3(vp(1),vp(2),vp(3),200,c,'filled')
    
end

sind1 = find(data.isSingularNode & ~data.isBoundaryVertex);
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
end


% campos([7.4758      -41.097      -1.8722])
campos([12.748      -39.352      -6.0979])

pic=getframe(gcf);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig2/fig2_2.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);
