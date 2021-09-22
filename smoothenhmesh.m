function [V,H] = smoothenhmesh(V0,H)
    
    V=V0;
    data = processhmesh(V,H,1);
    L = data.graphlaplacian;
    
    figure; hold all; rotate3d on; axis equal off;
    onorm = norm(V,'fro');
    for i=1:100
        V=V-.01*L*V;
        V=V/norm(V,'fro');
        V=V*onorm;
        
        try; delete(ptc); catch; end;
        ptc = patch('faces',data.F,'vertices',V,'facealpha',.1,'facecolor','green')
        drawnow;
    end
    
    

end