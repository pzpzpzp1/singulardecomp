function replayDecomp(decompdata)
    if nargin==0
        close all; clc;
%         decompdata = decompose_hmesh; % save('decompdata2.mat','decompdata');
%         filename='results/unit_70/decompdata.mat'; lwfac = .5; cpos=[-6.7853      -10.909       5.0505]; cpos = [ -0.59946     -0.79364       10.618];
%         filename='results/unit_15/decompdata.mat'; lwfac = .5; cpos=[-6.7853      -10.909       5.0505]; cpos=[-1.8606      -11.936      0.53508];
%         filename='results/sing1_59/decompdata.mat';
%         filename='results/tetpadded_16/decompdata.mat';  lwfac = .5; 
%         filename='results/hex_ellipsoid_coarse_78/decompdata.mat'; lwfac = .5; cpos=[-8.5277       5.9426       2.1422];
%         filename='results/Lpadded_77/decompdata.mat'; lwfac = .5; 
        filename='results/Cpadded_64/decompdata.mat'; lwfac = .5; cpos=[ 0.58797      -6.9005      -7.4385];
%         filename='results/hex_sphere_64/decompdata.mat'; lwfac = .5; 
%          filename='results/sing400_68/decompdata.mat'; cpos = [0.27478 -8.7551       17.758]; 
%          filename='results/sing400_68/decompdata.mat'; cpos=[17.376          -0.39907      -5.3769]; % view 2. not a great angle.
%         filename='results/sing2_88/decompdata.mat'; lwfac=.5; cpos=[15.132       -22.09       1.9178];
%         filename='results/sing3_92/decompdata.mat'; lwfac=.8; cpos=[ 17.899      -8.3873       19.894];
%         filename='results/sing133_42/decompdata.mat'; 
%         filename='results/sing044_73/decompdata.mat'; cpos=[0.17802      -6.5468       30.996];
%         filename='results/sing036_32/decompdata.mat'; cpos=[-0.059017    -0.086719        41.25];
%         filename='results/sing206_69/decompdata.mat'; cpos=[-11.68      -34.268       2.9365]; lwfac=.7;
%         filename='results/sing028_39/decompdata.mat'; lwfac=.5; cpos=[ 2.1639      -29.849       30.217];
%         filename='results/sing0012_65/decompdata.mat'; lwfac=.5; 

        load(filename);
    end
    if ~exist('lwfac','var'); lwfac=1; end;
    
    imc = 0;
    [dname,fname,ext]=fileparts(filename);
    dname(ismember(dname,'/'))='_';
    dname=[dname '_2']; % used for a second view.
%     dname = 'test'; % used for something we can overwrite.
    outname = ['figures/',dname,'/', dname,'_%d','.png'];
    
    
    datas = decompdata.datas;
    % get singularities and signatures
    %{
    for i=1:numel(datas)
        d = datas{i};
        snodes = find(d.isSingularNode & ~d.isBoundaryVertex);
        for j=1:numel(snodes)
            sigs{i,j} = getNode(d,snodes(j)).signature
        end
    end
    %}
    Vpresmooth = decompdata.Vpresmooth;
    cuts = decompdata.cuts;
    hexSheetInds = decompdata.hexSheetInds;
    
    fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); axis image vis3d;
    if exist('cpos','var'); campos(cpos); end;
    [~, guielems] = visualizeHmeshData(datas{1}, fh, lwfac);
    pause; campos
    imc=localsh(outname,imc);
    for i=1:numel(cuts)
        guielems = deleteElems(guielems); % clean slate
        
        data = datas{i};
        %% plot current state and selected cut sheet
        [~, guielems] = visualizeHmeshData(data, fh, lwfac);
        guielems{end+1} = patch('vertices',data.V,'faces',data.F(cuts{i},:),'facecolor','r','facealpha',.5);
%         title(sprintf('cut %d',i))
        pause;
        imc=localsh(outname,imc);
        
        
        %% plot transition
        tempdata = datas{i+1};
        Vold = Vpresmooth{i};
        Vnew = tempdata.V;
        ts = linspace(0,1,30);
        for ti=1:numel(ts)
            guielems = deleteElems(guielems);
            Vi = Vold*(1-ts(ti)) + Vnew*ts(ti);
            tempdata.V = Vi;
            [~, guielems] = visualizeHmeshData(tempdata,fh, lwfac);
            
            insertionlayerfaces = unique(tempdata.H2Farray(hexSheetInds{i},:));
            guielems{end+1} = patch('vertices',tempdata.V,'faces',tempdata.F(insertionlayerfaces,:),'facecolor','b','facealpha',.2);
            
            drawnow;
        end
%         title(sprintf('post sheet insertion %d',i))
        pause;        
        imc=localsh(outname,imc);
        
    end
    
    guielems = deleteElems(guielems);
    [~, guielems] = visualizeHmeshData(tempdata,fh, lwfac);
    imc=localsh(outname,imc);
    
end

function imc=localsh(outname,imc)
    pic=getframe(gcf);
    pic2=RemoveWhiteSpace(pic.cdata);
    oname = sprintf(outname,imc);
    dname = fileparts(oname);
    mkdir(dname);
    imwrite(pic2, oname);
    imc=imc+1;
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