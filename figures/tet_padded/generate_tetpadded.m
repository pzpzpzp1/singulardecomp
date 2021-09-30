

% save('temp.mat','data','cut','data2')
clear all; close all;
load temp.mat;

fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); axis image vis3d;
[~, guielems] = visualizeHmeshData(data, fh, .5);
patch('vertices',data.V,'faces',data.F,'facecolor','g','edgealpha',0,'facealpha',.02);
campos([-15.148      -7.7379       3.0292]);

hexsheetinds = [];
[Vnew,Hnew,hexsheetindsi,VnewPreperturb]=sheetinsertion(data, cut);
hexsheetinds=[hexsheetinds; hexsheetindsi];
cutf1234 = sort(data.F(cut,:),2);
data = processhmesh(Vnew,Hnew,1); 
newcut = ismember(sort(data.F,2), cutf1234,'rows')
patch('vertices',data.V,'faces',data.F(newcut,:),'facecolor','b','facealpha',.2);
[V] = smoothenhmesh(Vnew, Hnew, [], 1, 1,[],500,2,0,0);

for i=1:4
    [Vnew,Hnew,hexsheetindsi,VnewPreperturb]=sheetinsertion(data, newcut);
    hexsheetinds=[hexsheetinds; hexsheetindsi];
    cutf1234 = sort(data.F(newcut,:),2);
    data = processhmesh(Vnew,Hnew,1);
    newcut = ismember(sort(data.F,2), cutf1234,'rows')
    patch('vertices',data.V,'faces',data.F(newcut,:),'facecolor','b','facealpha',.2);
    [Vnew] = smoothenhmesh(Vnew, Hnew, [], 1, 1,[],500,2,0,0);
end


[V] = smoothenhmesh(Vnew, Hnew, [], 1, 1,[],500,2,0,0);
% [V] = smoothenhmesh(V, Hnew, [], 1, 1,[],500,2,0,0);
data = processhmesh(V,Hnew,1);

fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); axis image vis3d;
[~, guielems] = visualizeHmeshData(data2,fh, .3);
patch('vertices',data2.V,'faces',data2.F,'facecolor','g','edgealpha',0,'facealpha',.02);
finds = unique(data.H2Farray(unique(hexsheetinds),:));
patch('vertices',data.V,'faces',data.F(finds,:),'facecolor','b','edgealpha',1,'facealpha',.1);
campos([-15.148      -7.7379       3.0292]);

%% tet padded side
file_name = 'results/tet_split_notsplit/tetnotsplit.vtk'; lfac = 500; saveres=1;
mesh = load_vtk(file_name);
V0 = mesh.points;
H0 = mesh.cells;
visualize = 1;

[dname,fname,ext]=fileparts(file_name);
V=V0;H=H0;
data = processhmesh(V,H,0);
[V,H] = padhmesh(V,H); % [V,H] = padhmesh(V,H);
V = smoothenhmesh(V,H, [],visualize, 1, [], 100, 2, 0, 0);
data = processhmesh(V,H,visualize);

fh = figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w'); axis image vis3d;
[~, guielems] = visualizeHmeshData(data, fh, .5);
patch('vertices',data.V,'faces',data.F,'facecolor','g','edgealpha',0,'facealpha',.06);
campos([9.1342       12.352       1.9076]);

%% crop and re-save
files=dir('*.png');
for i=1:numel(files)
    pic = imread(files(i).name);
    pic2=RemoveWhiteSpace(pic);
    imwrite(pic2, ['redo_' files(i).name]);
end

