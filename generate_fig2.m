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

guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','blue','facealpha',.1,'edgealpha',0);
guielems{1} = patch('vertices',V,'faces',F(~isBoundaryFace,:),'facecolor','green','facealpha',.2,'edgealpha',.5,'linewidth',2);
guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',3*lwfac,'edgecolor','r');
guielems{4} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',3*lwfac,'edgecolor','b');
guielems{4} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',3*lwfac,'edgecolor','g');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),200*lwfac,'k','filled');

c = V(isIntSingularNode,:);
r=.37;
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
campos([7.4758      -41.097      -1.8722])


pic=getframe(gcf);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig2/fig2.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);
