clear all
close all
clc

NpixelsWidth = 2560; 
NpixelsHeight = 2160;
clip_height_top = 350;
clip_height_bottom = 100;
clipImageHeight = NpixelsHeight-clip_height_top-clip_height_bottom+1;
clip_lateral_edges = 250;
threshold_range = [100 200];
boundary_width = 20;
N_steelballs = 8;

L_noise = 2;
r_noise = 1;
A_noise = 1; 

r_bearing = 20;
L_bearing = 18;
%% plot filters
noise_filter = gaussian2D(L_noise, r_noise,1);
bearing_filter = gaussian2D(L_bearing, r_bearing, 1);

figure(90);
subplot(1,2,1); imagesc(noise_filter); colormap gray; axis equal;
title('noise filter');
subplot(1,2,2); imagesc(bearing_filter); colormap gray; axis equal;
title('bearing filter');
%% load image

params = struct('NpixelsWidth', NpixelsWidth, ...
                'NpixelsHeight', NpixelsHeight, ...
                'clip_height_top', clip_height_top, ...
                'clip_height_bottom', clip_height_bottom, ...
                'clip_lateral_edges', clip_lateral_edges, ...
                'boundary_width', boundary_width, ...
                'threshold_range', threshold_range, ...
                'N_steelballs', N_steelballs, ...
                'L_noise', L_noise, ...
                'r_noise', r_noise, ...
                'L_bearing', L_bearing, ...
                'r_bearing', r_bearing);

folder_root = 'H:\Calibration\112613\';

u = zeros(N_steelballs, 180);
v = zeros(N_steelballs, 180);

for n = 1:1 %180
    filename = strcat(folder_root, '100kV_2sec_200uA_P1_', num2str(n-1), '.ct');
    [u(:,n), v(:,n), img] = find_centroidfcn(params, filename);
end

% save(strcat(folder_root, 'uvdata.mat'), 'u', 'v', 'params');
% save(strcat(folder_root, 'uvdata_new_extraction.mat'), 'u', 'v', 'params');