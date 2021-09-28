clear all; close all;
file_name = 'meshes/hex_sphere.vtk';
mesh = load_vtk(file_name);
V0 = mesh.points;
H0 = mesh.cells;
[dname,fname,ext]=fileparts(file_name);
V=V0;H=H0;
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

guielems{1} = patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0);
% guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',.3,'linewidth',3*lwfac);
% guielems{3} = patch('vertices',V,'faces',E(isBoundaryEdge & isSingularEdge,[1 2 1]),'edgealpha',1,'linewidth',9*lwfac);
%     guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',0);
guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',3*lwfac,'edgecolor','r');


%     guielems{4} = scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),200*lwfac,'k','filled');
guielems{7} = scatter3(V(data.isSingularVertex,1),V(data.isSingularVertex,2),V(data.isSingularVertex,3),50*lwfac,'b','filled');

pic=getframe(gcf);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig1/fig1.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);
