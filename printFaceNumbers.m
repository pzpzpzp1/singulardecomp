
% label all faces. works surprisingly well.
for i=1:data.nF
    if ~data.isBoundaryFace(i)
        bc = data.faceBarycenters(i,:);
        text(bc(1),bc(2),bc(3),['<-' num2str(i)]);
    end
end


% extract trimesh from hexmesh
iv = find(~data.isBoundaryVertex);
av = setdiff(unique(data.E(any(data.E==iv,2),:)),iv);
tf = [];
for i=1:data.nH
    tf = [tf; intersect(data.H(i,:),av)'];
end
figure; hold all; axis equal; rotate3d on;
patch('vertices',data.V,'faces',tf,'facecolor','green','facealpha',.5)



