function [] = plot_skeleton(pts, A, options)

if nargin<3
    options.sizep=400;
end
sizep = getoptions(options, 'sizep', 6); 
sizee = getoptions(options, 'sizee', 2); 
colorp = getoptions(options, 'colorp', [.0, .8, .8]); 
colore = getoptions(options, 'colore', [.8, .0, .0]); 

scatter3(pts(:,1),pts(:,2),pts(:,3),sizep,'.','MarkerEdgeColor', colorp);  hold on;
plot_connectivity(pts, A, sizee, colore);