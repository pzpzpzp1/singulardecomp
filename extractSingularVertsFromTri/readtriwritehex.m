clear all; close all;

indir = 'trimeshes';
outdir = 'hmeshSings';

files = dir([indir '/*.txt']);

for i=1:numel(files)
    [~,fname,ext] = fileparts(files(i).name);
    outname = [outdir '/' fname '.vtk'];
    if ~exist(outname,'file')
        T = readmatrix([files(i).folder '/' files(i).name])+1;
        V = embedtrimesh(T);

        [Vh,H,centerind] = trig2hex_cleaned(V,T);
        [Vh, out] = smoothenhmesh(Vh, H, [], 1, 0,[],0);
        if out.minSJ(end) > .4
            mesh.points = Vh; mesh.cells = H;
            save_vtk(mesh, outname);
        end
    end
end
