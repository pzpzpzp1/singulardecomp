% minimizes dirichlet energy for interior vertices to get smooth equi-edge length mesh
% also minimizes scaled jacobian to ensure no hex element gets horribly skewed compared to any other
% DOES NOT PRESERVE BOUNDARY AT ALL. used for topological analysis. boundary isn't much of a concern.
function [V, out] = smoothenhmesh(V0, H, trimesh, visualize, preLapSmooth,fixednodes,lfac,p,fixb,uniformrot)
    if nargin==0
        file_name = 'results/sing1_59/hmesh_2.vtk';
%         file_name = 'results/hex_ellipsoid_coarse_78/hmesh_5.vtk';
%         file_name = 'results/sing3_59/hmesh_2.vtk';
        mesh = load_vtk(file_name);
        V0 = mesh.points;
        H = mesh.cells;
        visualize = 1;
        preLapSmooth = 1; % this essentially jostles the mesh.
%         [V0,H] = hex1to8(V0,H); [V0,H] = hex1to8(V0,H); 
        lfac = 500; % scaling factor on laplacian
        p = 2; % p norm
        fixb = false; % fixed boundary;
        uniformrot = false;
    end
    
    if visualize
        figure; hold all; axis equal; rotate3d on; % axis off;
        F = hex2face(H);
    end
    data = processhmesh(V0,H,0); L = data.graphlaplacian;
    centroid = (min(V0)+max(V0))/2;
    BBTR = max(V0) - centroid;
    
    V=V0;
    %% just dir E to start. unwinds edges that are too short.
    if preLapSmooth
        E = hex2edge(H);
        elens = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
        dt = min(elens)/10;
        for i=1:10
            Vs{i}=V;
            if visualize
                try; delete(ptc); catch; end;
                ptc = patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');
                drawnow;
            end
            [Elap, grad_lap] = dirE(L, V); 
            Elaps(i) = Elap;
            % store/accumulate energies
            grad = grad_lap;
            grad(fixednodes,:)=0;
            if fixb
                grad(data.isBoundaryVertex,:)=0;
            end

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
            randR = axang2rotm([randn(1,3) rand*2*pi]);
            V=V-newC;
            if uniformrot; V=V*randR; end; % rotate so bounding box applies uniformly instead of just on diagonal direcs
            newBBTR = max(V)-newC; % get BBTR on this rotation. this wont work well for very long skinny meshes. It'll try to force it into not long/skinny.
            V=V.*BBTR./newBBTR;
            if uniformrot; V=V*randR'; end; % unrotate
            V=V+newC;
            
        end
    end
    
    %% joint minimization
    
    maxiters=100;
    E = hex2edge(H);
    elens = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
    dt = min(elens)/300;    
    energybreakdown = []; minSJ=[]; maxSJ=[];
%     lfac = 5000; 
%     lfac = 500;
%     lfac = 0;
    for i=1:maxiters
        Vs{i}=V;
        if visualize
            try; delete(ptc); catch; end;
            ptc = patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');
            drawnow;
        end
        
        [Esj, grad_sj, out] = scaledJacobian_hmesh(V,H,p);
        [Elap, grad_lap] = dirE(L, V); grad_lap(data.isBoundaryVertex,:) = 0;
        minSJ(i) = min(out.SJ);
        maxSJ(i) = max(out.SJ);
        
        % mesh.points = V; mesh.cells = H; save_vtk(mesh, 'test.vtk');
        
        % store/accumulate energies
        if i==1
            norm1 = Esj; norm2 = Elap;
        end
        energybreakdown(i,1) = Esj/norm1;
        energybreakdown(i,2) = Elap/norm2;
        energy(i) = Esj/norm1 + lfac*Elap/norm2;
        grad = grad_sj/norm2 + lfac*grad_lap/norm2;
        grad(fixednodes,:)=0;
        if fixb
            grad(data.isBoundaryVertex,:)=0;
        end
        
        %{
        % remove gradients with norms that are too large.
        gradnorms = vecnorm(grad,2,2);
        unstablegrad = find(gradnorms > (mean(gradnorms)+std(gradnorms)*3));
        if numel(unstablegrad)~=0
            stc = scatter3(V(unstablegrad,1),V(unstablegrad,2),V(unstablegrad,3),'r','filled');
            qvr = quiver3(V(unstablegrad,1),V(unstablegrad,2),V(unstablegrad,3),grad(unstablegrad,1),grad(unstablegrad,2),grad(unstablegrad,3));
            delete(stc);
            delete(qvr);
        end
        grad(unstablegrad,:)=0;
        %}
        
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
        % hindsight adaptive timestep.
%         if i~=1
%             if energy(i) > energy(i-1)
%                 dt=dt/2;
%             else
%                 dt=dt*sqrt(2);
%             end
%         end
        dts(i) = dt;
        if dt < 1e-7
            display('timestep too small.');
            break;
        end
    
        % update V
        V = V - dt*grad;
        
        if i>30
            if energy(i-30)-energy(i) < .002
                display('converged');
                break;
            end
        end
        
    end
    out.minSJ = minSJ;
    out.maxSJ = maxSJ;
    out.Vs = Vs;
    out.energy = energy;
    out.dts = dts;
    out.energybreakdown=energybreakdown;
    % mesh.points = V; mesh.cells = H; save_vtk(mesh, 'test.vtk');
    
    
    
    
    
    
    
    
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