 function [] = eg_point_cloud_curve_skeleton(filename, outfilename, extension, delta, bsave)
% extract curve skeleton from a point cloud or triangular mesh
% update: 2010-8-19
% update: 2010-7-12
% create: 2009-4-26
% by: JJCAO, deepfish @ DUT
%
%% setting
% clear;clc;close all;
path('toolbox',path);
options.USING_POINT_RING = GS.USING_POINT_RING;
options.WC = GS.POSITION_CONSTRAINT_WEIGHT;
% extension='.off';
% extension='.ply';
%% Step 0: read file (point cloud & local feature size if possible), and
% normalize the modle.
% filename = '../data/9HandleTorus';% which file we should run on
% filename = '../../data/airplane/test/airplane_0667.2048';% which file we should run on
% middlename = 'airplane/';
tic
P.filename = [filename extension];% point set
name = filename(max(strfind(filename, '/'))+1 : end);
if(contains(name, '.'))
    name = name(1:max(strfind(name, '.'))-1);
end
P.name = name;

if extension == '.ply'
    ptCloud = pcread(P.filename);
    Location = sortrows(ptCloud.Location, 1, 'descend');
    Location = sortrows(ptCloud.Location, 2);
    
    P.pts = double(Location*delta);
else
    [P.pts,P.faces] = read_mesh(P.filename);
end

% if exist([filename '_fe.txt'],'file') % result of Tamal K Dey's NormFet
%     P.radis = load([filename '_fe.txt']);
% else
%     P.radis = ones(P.npts,1);
% end

if size(P.pts,1) > 4000
    ptCloud= pointCloud(P.pts);
    gridStep = 0.035;
    ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
    P.pts = ptCloudA.Location;
    figure;
    pcshow(ptCloudA,'MarkerSize', 36);
end
P.npts = size(P.pts,1);

P.pts = GS.normalize(P.pts);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);
disp(sprintf('read point set:'));
toc
%% Step 0: build local 1-ring
% build neighborhood, knn?
tic
P.k_knn = GS.compute_k_knn(P.npts);
if options.USING_POINT_RING
    P.rings = compute_point_point_ring(P.pts, P.k_knn, []);
else    
    P.frings = compute_vertex_face_ring(P.faces);
%     P.rings = compute_vertex_ring(P.faces, P.frings);
    P.rings = compute_vertex_ring(P.faces);
    options.USING_POINT_RING = 1;
end
disp(sprintf('compute local 1-ring:'));
toc

%% Step 1: Contract point cloud by Laplacian
%point cloud contraction_by_mesh_laplacian
P.delta = delta
tic
[P.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options);
disp(sprintf('Contraction:'));
toc
%% step 2: Point to curve – by cluster ROSA2.0
%RECOVER CONNECTIVITY after each edge collapse
tic
P.sample_radius = P.diameter*0.04;  

% P.sample_radius = 0.01;  
P.sample_ratio = 0.01;%暂时没有用
P = rosa_lineextract(P, 1);
disp(sprintf('to curve:'));


%% show results
figure('Name','Original point cloud and its contraction');movegui('northeast');set(gcf,'color','white')
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),30,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
scatter3(P.cpts(:,1),P.cpts(:,2),P.cpts(:,3),30,'.r'); axis off;axis equal;set(gcf,'Renderer','OpenGL');
camorbit(0,0,'camera'); axis vis3d; view(0,90);view3d rot;

figure('Name','Original point cloud and its skeleton'); movegui('center');set(gcf,'color','white');
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),20,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
showoptions.sizep=400;showoptions.sizee=2;
plot_skeleton(P.spls, P.spls_adj, showoptions);
axis off;axis equal;set(gcf,'Renderer','OpenGL');view(0,90);view3d rot;
%% save results
if bsave == 1
    % 删掉突出的点
    spls = P.spls./delta;
    spls_adj = P.spls_adj;
    idx = spls(:,2)<P.bbox(5);
    P.spls = spls(idx, :);
    spls_adj = spls_adj(idx,:);
    spls_adj = spls_adj(:,idx);
%     spls_adj(33,39) = 0;
%     spls_adj(39,33) = 0;
    P.spls_adj = spls_adj;
    
    default_filename = sprintf('%s%s_skeleton.mat',  outfilename, P.name);
    save(default_filename,'P');
    % filename_Ske = [outfilename, middlename, P.name, '.ply'];
    filename_Ske = [outfilename, P.name, '.ply'];
    ptCon = pointCloud(P.cpts); 


    toc
    ptSke = pointCloud(P.spls);
    % pcwrite(ptCon, filename_Con);
    pcwrite(ptSke, filename_Ske);
end

 end