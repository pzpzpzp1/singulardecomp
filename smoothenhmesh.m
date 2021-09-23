% minimizes dirichlet energy for interior vertices to get smooth equi-edge length mesh
% also minimizes scaled jacobian to ensure no hex element gets horribly skewed compared to any other
% DOES NOT PRESERVE BOUNDARY AT ALL. used for topological analysis. boundary isn't much of a concern.
function [V, out] = smoothenhmesh(V0, H, trimesh, visualize)
    if nargin==0
        file_name = 'results/sing1_59/hmesh_2.vtk';
%         file_name = 'results/hex_ellipsoid_coarse_78/hmesh_5.vtk';
%         file_name = 'results/sing3_59/hmesh_2.vtk';
        mesh = load_vtk(file_name);
        V0 = mesh.points;
        H = mesh.cells;
        visualize = 1;
%         [V0,H] = hex1to8(V0,H); [V0,H] = hex1to8(V0,H); 
    end
    
    if visualize
        figure; hold all; axis equal; rotate3d on; % axis off;
        F = hex2face(H);
    end
    data = processhmesh(V0,H,0); L = data.graphlaplacian;
    centroid = (min(V0)+max(V0))/2;
    BBTR = max(V0) - centroid;
    
    maxiters=1000;
    V=V0; p=1; 
    E = hex2edge(H);
    elens = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
    dt = min(elens)/1;
    %% just dir E to start. unwinds edges that are too short.
    for i=1:10
        Vs{i}=V;
        if visualize
            try; delete(ptc); catch; end;
            ptc = patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');
            drawnow;
        end
        [Elap, grad_lap] = dirE(L, V); 
        
        % store/accumulate energies
        grad = grad_lap;
        
        % line search. max 50 iterations to get ls to work.
        for j=1:50
            candV = V-dt*grad;
            elap = dirE(L, candV);
            etot = elap;
            if etot > Elap
                dt = dt/2;
            else
                break
            end
        end
        dt = dt*(2^(1/4));
        dts(i) = dt;
        
        % update V. and recenter/rescale
        V = V - dt*grad;
        newC = (min(V)+max(V))/2;
        newBBTR = max(V)-newC;
        V=V-newC;
        V=V.*BBTR./newBBTR;
        V=V+newC;
    end
    
    %% joint minimization
    E = hex2edge(H);
    elens = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
    dt = min(elens)/300;    
    energybreakdown = [];
    lfac = 5000; lfac = 10;
    for i=1:maxiters
        Vs{i}=V;
        if visualize
            try; delete(ptc); catch; end;
            ptc = patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');
            drawnow;
        end
        
        [Esj, grad_sj, out] = scaledJacobian_hmesh(V,H,p);
        [Elap, grad_lap] = dirE(L, V); grad_lap(data.isBoundaryVertex,:) = 0;
        
        % store/accumulate energies
        if i==1
            norm1 = Esj; norm2 = Elap;
        end
        energybreakdown(i,1) = Esj/norm1;
        energybreakdown(i,2) = Elap/norm2;
        energy(i) = Esj/norm1 + lfac*Elap/norm2;
        grad = grad_sj/norm2 + lfac*grad_lap/norm2;
        
        % line search. max 50 iterations to get ls to work.
        for j=1:50
            candV = V-dt*grad;
            esj = scaledJacobian_hmesh(candV,H,p);
            elap = dirE(L, candV);
            etot = esj/norm1 + lfac*elap/norm2;
            if etot > energy(i)
                dt = dt/2;
            else
                break
            end
        end
        dt = dt*(2^(1/4));
        dts(i) = dt;
        
        % update V
        V = V - dt*grad;
    end
    out.Vs = Vs;
    out.energy = energy;
    out.dts = dts;
    out.energybreakdown=energybreakdown;
    
    
    
    
    
    
    
    
    
    %{
    c0 = (max(V0)+min(V0))/2; 
    V0 = V0 - c0;
    BBd0 = norm([max(V0)-min(V0);]);
    trimesh.Vertices = trimesh.Vertices - c0;
    
    V=V0; 
    data = processhmesh(V,H,0);
    L = data.graphlaplacian;
    
    figure; hold all; rotate3d on; axis equal off;
    ptc1 = patch('faces',trimesh.Faces,'vertices',trimesh.Vertices,'facealpha',.1,'facecolor','red','edgealpha',0);
    for i=1:10
        V = (speye(data.nV) + .1*L)\V;
        energy(end+1) = norm((L*V).*V,'fro')^2/2
    
        % recenter rescale
        V=V- [max(V)+min(V)]/2;
        BBd = norm([max(V)-min(V)]);
        V=V/BBd*BBd0;
        
%         project boundary back to boundary
%         [dist, projpts, faces2, vertices2, corresponding_vertices_ID, new_faces_ID] = point2trimesh('Faces',trimesh.Faces,'Vertices',trimesh.Vertices,'QueryPoints',V(data.isBoundaryVertex,:));
%         V(data.isBoundaryVertex,:) = projpts;
        V(data.isBoundaryVertex,:) = V0(data.isBoundaryVertex,:);
        
        try; delete(ptc); catch; end;
        ptc = patch('faces',data.F,'vertices',V,'facealpha',.1,'facecolor','green');
        drawnow;
    end
    %}
    
    %{
     mesh.points = V; mesh.cells = H;
    save_vtk(mesh, 'test2.vtk')     
    save_vtk(mesh, 'tetnotsplit.vtk')
     save_vtk(mesh, 'tetsplit.vtk')
     save_vtk(mesh, 'tetpadded.vtk')
    
    figure; hold all; 
    plot(energybreakdown(:,1),'r'); 
    plot(energybreakdown(:,2),'g'); 
    legend('sj','laplacian');
    plot(energy,'b')
    %}

end


function [energy,grad] = dirE(L, V)
    grad = (L*V);
    energy = sum(grad.*V,[1 2])/2;
    
end