function [cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options)
% point cloud contraction_by_mesh_laplacian
% refer to Skeleton Extraction by Mesh Contraction 08
% mesh could be simplified by MeshLab/Filter/Clustering Decimation, if it
% is too larged.
% 
% inputs:
%   P.pts
%   P.faces
%   P.npts
%   P.k_knn: k of knn
%   P.rings:
%   options.WL = 3;% initial contraction weight
%   options.WC = 1;% initial attraction weight
%   options.sl: scalar for increasing WL in each iteration
%   options.tc: Termination Conditions for total area ratio
%   options.iterate_time = 10; max iteration steps
%
% outputs:
%   cpts: contracted vertices
%
% notes:
%   compute_point_laplacian() is slow. Not changing weights in iteration
%   and using combinatorial weights both lead to failed to contraction
%   enough.
% @author: deepfish, jjcao
% @create-data:     2009-4-10
% @modify-data:     2009-4-23
% @modify-data:     2009-5-19
% @modify-data:     2009-9-03
% @modify-data:     2010-8-19
% profile on;

if nargin < 1
    clear;clc;close all;
    P.filename = '../data/simplejoint_v4770.off'; 
    options.USING_POINT_RING = GS.USING_POINT_RING;
    options.iterate_time = 10;
[P.pts,P.faces] = read_mesh(P.filename);
P.npts = size(P.pts,1);
P.pts = GS.normalize(P.pts);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);

P.k_knn = GS.compute_k_knn(P.npts);  
atria = nn_prepare(P.pts); 
[P.knn_idx, P.knn_dist] = nn_search(P.pts, atria, P.pts, P.k_knn); 
P.rings = compute_point_point_ring(P.pts,P.k_knn, P.knn_idx);
% P.rings = compute_vertex_face_ring(P.faces);
end

%##########################################################################
%% setting
%##########################################################################
% visual debug conditions
RING_SIZE_TYPE = 2;%1:min, 2:mean, 3:max
SHOW_CONTRACTION_PROGRESS = true;
Laplace_type = 'conformal';%conformal%combinatorial%spring%mvc

% setting
tc = getoptions(options, 'tc', GS.CONTRACT_TERMINATION_CONDITION); 
iterate_time = getoptions(options, 'iterate_time', GS.MAX_CONTRACT_NUM); 

initWL = getoptions(options, 'WL', GS.compute_init_laplacian_constraint_weight(P,Laplace_type)); 
% set initial attraction weights according to different type of discrete
% Laplacian
if strcmp(Laplace_type,'mvc')
    WC = getoptions(options, 'WC', GS.POSITION_CONSTRAINT_WEIGHT)*10;
elseif strcmp(Laplace_type,'conformal')
    WC = getoptions(options, 'WC', 1);
else
    WC = getoptions(options, 'WC', 1);
end
WH = ones(P.npts, 1)*WC; % 初始约束权
sl = getoptions(options, 'sl', GS.LAPLACIAN_CONSTRAINT_SCALE); % scale factor for WL in each iteration! in original paper is 2;
s2 = getoptions(options, 's2', GS.POSITION_CONSTRAINT_SCALE); 
minWL = GS.MIN_POSITION_CONSTRAINT_WEIGHT;
WL = initWL;%*sl;

sprintf(['1) k of knn: %d\n 2) termination condition: %f \n 3)' ...
    'Init Contract weight: %f\n 4) Init handle weight: %f\n 5) Contract scalar: %f\n' ...
    '6) Max iter steps: %d'], P.k_knn, tc, initWL, WC, sl,iterate_time)

%% init iteration
t = 1; % current iteration step
%% left side of the equation
tic
if options.USING_POINT_RING
    L = -compute_point_laplacian(P.pts,Laplace_type,P.rings, options);%conformal;%spring
else
%     L = -compute_point_laplacian(P.pts,Laplace_type,P.rings, options);%conformal;%spring
    L = -compute_mesh_laplacian(P.pts,P.faces,Laplace_type,options);
end
A = [L.*WL;sparse(1:P.npts,1:P.npts, WH)];
% right side of the equation
b = [zeros(P.npts,3);sparse(1:P.npts,1:P.npts, WH)*P.pts];
cpts = (A'*A)\(A'*b); 
% newVertices = A\b; % this is slow than above line
disp(sprintf('solve equation:'));
toc

if SHOW_CONTRACTION_PROGRESS
    tic
    figure();
    movegui('northeast');axis off;axis equal;set(gcf,'Renderer','OpenGL'); view3d rot;hold on;set(gcf,'color','white');
    camorbit(0,0,'camera'); axis vis3d; view(-90,0);    
    h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,'b','filled');
    h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,'r','filled');
    %legend('orignal points','contracted points');
    title(['iterate ',num2str(t),' time(s)'])    
    disp(sprintf('draw mesh:'));
    toc
end
%%
tic
if options.USING_POINT_RING
    sizes = GS.one_ring_size(P.pts, P.rings, RING_SIZE_TYPE);   % min radius of 1-ring
    size_new = GS.one_ring_size(cpts, P.rings, RING_SIZE_TYPE);
    a(t) = sum(size_new)/sum(sizes);
else
%     sizes = GS.one_ring_area(P.pts,P.rings);
%     size_new = GS.one_ring_area(cpts,P.rings);

%     sizes = area_1_face_ring(P.pts, P.faces, P.frings);
%     size_new =area_1_face_ring(cpts, P.faces, P.frings);    
    ratio_new = area_ratio_1_face_ring(P.pts, cpts, P.faces, P.frings);
    ratio = ones(size(ratio_new));
    a(t) = mean(ratio_new);%sum(ratio_new)/sum(ratio);
end

disp(sprintf('compute area:'));
toc

% mwrite(['A' num2str(t) '.txt'], A);
% mwrite(['b' num2str(t) '.txt'], b);
% mwrite(['sizes' num2str(t) '.txt'], sizes);
% mwrite(['size_new' num2str(t) '.txt'], size_new);
% mwrite(['newVertices' num2str(t) '.txt'], cpts);
npts = size(P.pts, 1);
pts = P.pts;
pStable = ones(npts,1);
%pdists = pdist2(pts, pts);
[knn_idx, pdists] = knnsearch(cpts, cpts, 'K', P.k_knn*2);
while t<iterate_time
    if options.USING_POINT_RING
        % Jerry： 更新ring
        [rings, flag] = compute_point_point_ring(cpts, P.k_knn, []);
        if flag == 1
            emp = cellfun(@isempty, rings);
            idx = find(emp == 0);
            P.rings(idx) = rings(idx);
        end
        L = -compute_point_laplacian(cpts,Laplace_type, P.rings, options);%conformal
        
        %计算局部稳定的状态
        pStable = ones(npts,1);
        if t > 2
            pCoefs = double(zeros(npts,1));
            for i = 1:npts
                h = P.diameter*0.08;
                pdist = pdists(i,:);
                idx_nei = knn_idx(i,pdist<h);
                pts_i = cpts(idx_nei,:);
                %pts_i = pts( knn_idx(i,1:5), :);
                [coefs,~,latent] = pca(pts_i);  
                pCoefs(i) = latent(1)/sum(latent);
            end
            for i = 1:npts
                sigma = mean(pCoefs(knn_idx(i,:)));
                if sigma > 0.93
                    pStable(i) = 0;
                end
            end
            if sum(pStable) < 0.1*npts
                sprintf('most points stabled, break\n')
                break;
            end
            sprintf('number of points stabled: %d\n', npts-sum(pStable))
            L = L.* pStable';
            L = L.* pStable;
        end
    else
%         L = -compute_point_laplacian(cpts,Laplace_type,P.rings, options);%conformal;%spring
        L = -compute_mesh_laplacian(cpts,P.faces,Laplace_type,options);
    end
    
    WL = sl*WL;    
    WC = WC/s2;
    WC = max(WC, minWL);
    if WL>GS.MAX_LAPLACIAN_CONSTRAINT_WEIGHT,WL=GS.MAX_LAPLACIAN_CONSTRAINT_WEIGHT;end; % from Oscar08's implementation, 2048
    if options.USING_POINT_RING
        if strcmp(Laplace_type,'mvc')
            WH = WC.*(sizes./size_new)*10;% 初始约束权
        else
            WH = WC.*(sizes./size_new);
        end
    else    
%         WH = WC.*(sizes./size_new);
%         WH = WC.*((sizes./size_new).^0.5); % a little diff with Oscar08
        WH = WC*(ratio_new.^(-0.5));% from Oscar08, if the size is area
    end
    
    WH(WH>GS.MAX_POSITION_CONSTRAINT_WEIGHT) = GS.MAX_POSITION_CONSTRAINT_WEIGHT;% from Oscar08's implementation, 10000

    A = real([WL*L;sparse(1:P.npts,1:P.npts, WH)]);
    
    % update right side of the equation
    b(P.npts+1:end, :) = sparse(1:P.npts,1:P.npts, WH)*cpts;
    tmp = (A'*A)\(A'*b);
    
    if options.USING_POINT_RING
        size_new = GS.one_ring_size(tmp, P.rings, RING_SIZE_TYPE);  
        a(end+1) = sum(size_new)/sum(sizes);
    else
%         size_new = GS.one_ring_area(cpts,P.rings);
%         size_new =area_1_face_ring(tmp, P.faces, P.frings);
        ratio = ratio_new;
        ratio_new = area_ratio_1_face_ring(P.pts, tmp, P.faces, P.frings);
        a(end+1) = mean(ratio_new);%sum(ratio_new)/sum(ratio);
    end    

    vec = tmp-cpts;
    vec = sum(vec.*vec,2);
    r = P.diameter*0.08;
    idx_wrong = vec > r*r; 
    tmp(idx_wrong,:) = cpts(idx_wrong,:)
    %终止的条件
%     if a(t)-a(end)<tc || isnan(a(end))
%         break;
%     else 
%         cpts = tmp;
%     end
    cpts = tmp;
    t = t+1
    
    if SHOW_CONTRACTION_PROGRESS
        % 显示前后点云     
        delete(h1);delete(h2);
        h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,WH,'filled');
        h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,ones(P.npts,1)*WL,'filled');
        %legend('orignal points','contracted points');
        title(['iterate ',num2str(t),' time(s)']); drawnow; 
    end
end
clear tmp;

if SHOW_CONTRACTION_PROGRESS
    figure;
    plot(a);xlabel('Iteration times');ylabel('Ratio of original and current volume');
    if t<iterate_time
        if a(end) < tc 
            sprintf('exist iteration for termination condition meets!')
        else
             warning('exist iteration unnormally!');
        end
    end
end
% profile off;
% profile viewer;

%%
% default_filename = sprintf('result\\%s_contract_t(%d)_nn(%d)_WL(%f)_WH(%f)_sl(%f).off',...
%     P.filename(1:end-4), t, P.k_knn, initWL, WC, sl);
% [FileName,PathName,FilterIndex] = uiputfile( {'*.off';'*.obj';'*.wrl';'*.*'}, 'Save as',default_filename);
% if FilterIndex~=0
%     write_mesh([PathName FileName],cpts, P.faces);
% end
