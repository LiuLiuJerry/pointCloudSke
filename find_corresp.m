%  @author:Jerry
%   2018-8-22
function [spls, corresp] = find_corresp(pts, sample_ratio)
    ptCloudOut = pcdownsample(pointCloud(pts), 'random', sample_ratio);
%     pcshow(ptCloudOut)
    spls = ptCloudOut.Location;
    [Idx, Dist] = knnsearch(spls, pts);
    A = 1 : length(spls);
    corresp = Idx(:, 1);

    x = corresp;
    x = unique(x);
    out = spls(setdiff(A, x));
    
    if size(out, 2) == 0 || size(out, 1)==0
        return;
    end
    out_near = knnsearch(pts, out, 'K', 1);
    for i = 1:length(out)
        idx = out_near(i, 1);
        corresp(idx) = i;
    end
end