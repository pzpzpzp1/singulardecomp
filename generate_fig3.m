
file_name = 'meshes/hex_ellipsoid_coarse.vtk';
mesh = load_vtk(file_name);
V0 = mesh.points;
H0 = mesh.cells;
V=V0;H=H0;
visualize=1;

V = smoothenhmesh(V,H,[],visualize, 1, [], 100, 2, 0);

data = processhmesh(V,H,visualize); 

cutseed = 263;
cutfaces = false(size(data.F,1),1); cutfaces(cutseed)=true;
cutfaces(extendRegularFaces(data, cutseed))=true;


f1=figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); campos([ -6.1338      -3.6235     -0.38119]);
cla; hold all; axis equal off; rotate3d on;
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
%guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',.3,'linewidth',3*lwfac);
%guielems{3} = patch('vertices',V,'faces',E(isBoundaryEdge & isSingularEdge,[1 2 1]),'edgealpha',1,'linewidth',9*lwfac);
%     guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',0);
guielems{4} = patch('vertices',V,'faces',E(IS3,[1 2 1]),'linewidth',4,'edgecolor','r');
%guielems{5} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',9*lwfac,'edgecolor','b');
%guielems{6} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',9*lwfac,'edgecolor','g');
%     guielems{4} = scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled');
guielems{7} = scatter3(V(isIntSingularNode,1),V(isIntSingularNode,2),V(isIntSingularNode,3),100,'k','filled');
ptc = patch('vertices',data.V,'faces',data.F(cutfaces,:),'facecolor','r','facealpha',.5)

%% insert sheet
[V2,H2,inds]=sheetinsertion(data, cutfaces);

V2 = smoothenhmesh(V2,H2,[],visualize, 1, [], 100, 2, 0);

data2 = processhmesh(V2,H2,visualize); 

f2=figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); campos([ -6.1338      -3.6235     -0.38119]);
cla; hold all; axis equal off; rotate3d on;
V2=data2.V; F2=data2.F; E2=data2.E; 
isSingularNode = data2.isSingularNode;
isBoundaryEdge = data2.isBoundaryEdge;
isBoundaryFace = data2.isBoundaryFace;
isSingularEdge = data2.isSingularEdge;
isIntSingularNode = data2.isSingularNode & ~data2.isBoundaryVertex;
IS3 = data2.efdeg==3 & data2.isSingularEdge & ~data2.isBoundaryEdge;
IS5 = data2.efdeg==5 & data2.isSingularEdge & ~data2.isBoundaryEdge;
IS_ = (data2.efdeg~=3 & data2.efdeg~=5) & data2.isSingularEdge & ~data2.isBoundaryEdge;
guielems{1} = patch('vertices',V2,'faces',F2(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0);
guielems{1} = patch('vertices',V2,'faces',F2(data2.H2Farray(inds,:),:),'facecolor','blue','facealpha',.2,'edgealpha',1);
%guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',.3,'linewidth',3*lwfac);
%guielems{3} = patch('vertices',V,'faces',E(isBoundaryEdge & isSingularEdge,[1 2 1]),'edgealpha',1,'linewidth',9*lwfac);
%     guielems{2} = patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]),'edgealpha',0);
guielems{4} = patch('vertices',V2,'faces',E2(IS3,[1 2 1]),'linewidth',4,'edgecolor','r');
%guielems{5} = patch('vertices',V,'faces',E(IS_,[1 2 1]),'linewidth',9*lwfac,'edgecolor','b');
%guielems{6} = patch('vertices',V,'faces',E(IS5,[1 2 1]),'linewidth',9*lwfac,'edgecolor','g');
%     guielems{4} = scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled');
guielems{7} = scatter3(V2(isIntSingularNode,1),V2(isIntSingularNode,2),V2(isIntSingularNode,3),100,'k','filled');


pic=getframe(f1);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig3/fig3a.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);

pic=getframe(f2);
pic2=RemoveWhiteSpace(pic.cdata);
oname = 'figures/fig3/fig3b.png';
dname = fileparts(oname);
mkdir(dname);
imwrite(pic2, oname);



        
        
        
        
        
        
        
        
        
        
        