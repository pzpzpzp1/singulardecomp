%{ 
yeesh this is harder than i expected. 
minimizing squared volume degenerates to a flat config wtih 0 volume.
maximizing volume results in n-cover wrapovers.
%}
function V = smoothSphereTriEmbedding(V,F,vis)
if nargin < 3
    vis=1;
end

V=V-mean(V);
V=V./vecnorm(V,2,2);
nF = size(F,1);
nV = size(V,1);

if vis; figure; axis equal; hold all; rotate3d on; scatter3(0, 0, 0, 'r', 'filled'); end;
dt = .1;
maxiters=  100;
spherevol = pi;
% spherevol = 0;
p=2;
for i=1:maxiters
    
    n34 = zeros(nF,3,4);
    n34(:,:,2:end) = permute(reshape(V(F',:),3,nF,3),[2 3 1]);
    [svols, svolgrad] = face_planarity(n34); % Want all to be as sas possible.
    
    % add spherevol to make all values positive. use p norm to minimize
    % worst case tet.
    energy(i) = sum((svols+spherevol).^p);
    pscaledgrad = (p*(svols+spherevol).^(p-1)).*svolgrad; 
    
    % hindsight adaptive timestep
    if i~=1
        if energy(i) > energy(i-1)
            dt=dt/2;
        else
            dt=dt*sqrt(2);
        end
    end
    
    svolgrad_reduced = pscaledgrad(:,:,2:end);
    xx = accumarray(F(:), reshape(svolgrad_reduced(:,1,:),[],1), [nV, 1]);
    yy = accumarray(F(:), reshape(svolgrad_reduced(:,2,:),[],1), [nV, 1]);
    zz = accumarray(F(:), reshape(svolgrad_reduced(:,3,:),[],1), [nV, 1]);
    Vgrad = [xx yy zz];
    
    
    % advance and project onto centered sphere
    V=V-Vgrad*dt;
    V=V-mean(V);
    V=V./vecnorm(V,2,2);
    
    % vis
    if vis
        try; delete(ptc); catch; end;
        ptc = patch('Faces', F, 'Vertices', V, 'facecolor',  'blue', 'facealpha', 0.1);
        drawnow;
    end
    
    
end



end