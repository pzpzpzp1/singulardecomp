function [Vo,Ho] = padhmesh(V,H)
    data = processhmesh(V,H,0);
    nV = size(V,1);

    %% build new vertices and topology
    movedist = min(data.edgelengths)/2;
    Vo=[V;V]; % double verts. lower half are unreferenced currently.
    Vo(data.isBoundaryVertex,:) = V(data.isBoundaryVertex,:) - movedist*data.bvertNormals; % move old boundary inwards.
    
    oldboundaryF = data.F(data.isBoundaryFace,:);
    newboundaryF = oldboundaryF + nV;
    newH = [oldboundaryF newboundaryF]; % pad layer of hexes.
    Ho = [H; newH];
    
    %% remove unreferenced vertices
    [Vo,Ho] = minimizeMesh(Vo,Ho);
    
    %% vis and verify
%     data = processhmesh(Vo,Ho,1);

%     mesh.points = Vo; mesh.cells = Ho;
%     save_vtk(mesh, 'test.vtk')
    
end