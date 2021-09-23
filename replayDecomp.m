function replayDecomp(decompdata)
    if nargin==0
        close all; clc;
%         decompdata = decompose_hmesh; % save('decompdata2.mat','decompdata');
%         load results/unit_70/decompdata.mat;
%         load results/unit_15/decompdata.mat;
%         load results/sing1_59/decompdata.mat;
%         load results/tetpadded_16/decompdata.mat;
%         load results/hex_ellipsoid_coarse_78/decompdata.mat;
%         load results/sing2_88/decompdata.mat;
        load results/sing3_92/decompdata.mat;
    end
    
    datas = decompdata.datas;
    Vpresmooth = decompdata.Vpresmooth;
    cuts = decompdata.cuts;
    hexSheetInds = decompdata.hexSheetInds;
    
    fh = figure; guielems={};
    for i=1:numel(cuts)
        guielems = deleteElems(guielems); % clean slate
        
        data = datas{i};
        %% plot current state and selected cut sheet
        [~, guielems] = visualizeHmeshData(data, fh);
        guielems{end+1} = patch('vertices',data.V,'faces',data.F(cuts{i},:),'facecolor','r','facealpha',.9);
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
            
            insertionlayerfaces = unique(tempdata.H2Farray(hexSheetInds{i},:));
            guielems{end+1} = patch('vertices',tempdata.V,'faces',tempdata.F(insertionlayerfaces,:),'facecolor','c','facealpha',.4);
            
            drawnow;
        end
        title(sprintf('post sheet insertion %d',i))
        pause;
    end
    
    guielems = deleteElems(guielems);
    [~, guielems] = visualizeHmeshData(tempdata,fh);
    
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