clear all;
close all;
idx = 0;
%·Ö¸îºóµÄ
formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%.4d';
% input_file = 'Octopus\octopus_seg';
input_file = 'Octopus\octopus2048';
output_file = 'Octopus\skeletons';
%input_file = 'Chen\samba2048_ply'
%output_file = 'Chen\skeletons_samba'
extension1 = '.off';
extension2 = '.ply';
for ii = 1
     filename = sprintf(formatSpec, input_file, ii);
     outfilename = sprintf(formatSpec, output_file, ii);
%     filename = 'legs'
%     pc_legs = pcread(filename);
%     xyz = pc_legs.Location;
%     idx = kmeans(xyz, K, 'Distance', 'cosine' , 'MaxIter', 5000);
%     figure();
%     c = jet(K);
%     for i = 1:8
%         leg = xyz(idx==i,:);
%         pcshow(leg, c(i,:)); hold on;
%     end
    eg_point_cloud_curve_skeleton(filename, outfilename, extension2, 1, 1);
%     clear;clc;close all;
end