function replayDecomp(decompdata)
    if nargin==0
        close all; clc;
%         decompdata = decompose_hmesh; % save('decompdata2.mat','decompdata');
        load results/unit_45.mat;
        load results/unit_36.mat;
    end
    
    datas = decompdata.datas;
    Vpresmooth = decompdata.Vpresmooth;
    cuts = decompdata.cuts;
    
    fh = figure; guielems={};
    for i=1:numel(cuts)
        guielems = deleteElems(guielems); % clean slate
        
        data = datas{i};
        %% plot current state and selected cut sheet
        [~, guielems] = visualizeHmeshData(data, fh);
        guielems{end+1} = patch('vertices',data.V,'faces',data.F(cuts{i},:),'facecolor','c','facealpha',.9);
        title(sprintf('cut %d',i))
        pause;
        
        %% plot transition
        tempdata = datas{i+1};
        Vold = Vpresmooth{i};
        Vnew = tempdata.V;
        ts = linspace(0,1,30);
        for ti=1:numel(ts)
            guielems = deleteElems(guielems);
            Vi = Vold*(1-ts(ti)) + Vnew*ts(ti);
            tempdata.V = Vi;
            [~, guielems] = visualizeHmeshData(tempdata,fh);
            drawnow;
        end
        title(sprintf('post sheet insertion %d',i))
%         pause;
    end
    
    
    
end

function guielems = deleteElems(guielems)
    try;
    for i=1:numel(guielems)
        delete(guielems{i});
    end
    catch;
    end;
    guielems = {};
end