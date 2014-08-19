% clear all
% close all
% clc
% 
% NpixelsWidth = 2560; 
% NpixelsHeight = 2160;
% clip_height_top = 380;
% clip_height_bottom = 100;
% clipImageHeight = NpixelsHeight-clip_height_top-clip_height_bottom+1;
% threshold_range = [100 200];
% boundary_width = 20;
% N_steelballs = 8;
% 
% L_noise = 2;
% r_noise = 1;
% A_noise = 1; 
% 
% r_bearing = 20;
% L_bearing = 15;
% % plot filters
% noisefilter = gaussian2D(L_noise, r_noise,1);
% gaussianfilter2d = gaussian2D(L_bearing, r_bearing, 1);
% 
% figure;
% subplot(1,2,1); imagesc(noisefilter); colormap gray; axis equal;
% title('noise filter');
% subplot(1,2,2); imagesc(gaussianfilter2d); colormap gray; axis equal;
% title('bearing filter');
% % load image
% 
% params = struct('NpixelsWidth', NpixelsWidth, ...
%                 'NpixelsHeight', NpixelsHeight, ...
%                 'clip_height_top', clip_height_top, ...
%                 'clip_height_bottom', clip_height_bottom, ...
%                 'boundary_width', boundary_width, ...
%                 'threshold_range', threshold_range, ...
%                 'N_steelballs', N_steelballs, ...
%                 'L_noise', L_noise, ...
%                 'r_noise', r_noise, ...
%                 'L_bearing', L_bearing, ...
%                 'r_bearing', r_bearing);
% 
% folder_root = 'C:\Users\Ryan-Helen\Documents\Calibration\112613\';
% 
% filename = strcat(folder_root, '100kV_2sec_200uA_P1_', num2str(35), '.ct');

function [u, v, img] = find_centroid2(params, filename, noisefilter, gaussianfilter2d)

%% load image
dat = readBinary(filename, params.NpixelsHeight*params.NpixelsWidth+1, 'uint16');
dat = reshape(dat(2:end), params.NpixelsWidth, params.NpixelsHeight);
    figure(1); imagesc(dat, [0 500]); colormap gray; axis equal; 

%% clip the original image
clip_dat = dat(params.clip_height_top:end-params.clip_height_bottom,params.clip_width_left:end-params.clip_width_right);
figure; imagesc(clip_dat, [0 500]); colormap gray; colorbar; axis equal;

clear dat;
dat_noedge = zeros(size(clip_dat));
dat_noedge(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width) = ...
           clip_dat(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width); 
figure; imagesc(dat_noedge, [0 500]); colormap gray; axis equal;       
%% simple threshold 
n = find(dat_noedge(:)<= params.threshold_range(2) & dat_noedge(:)>= params.threshold_range(1));
img_simple_threshold = zeros(size(clip_dat));
img_simple_threshold(n) = clip_dat(n);
figure; imagesc(img_simple_threshold); colormap gray; colorbar; axis equal;
  
%% find x and y coordinates of filter1
clipImageHeight = params.NpixelsHeight-params.clip_height_top-params.clip_height_bottom+1;
xcoord_noedge = floor(n/clipImageHeight)+1;
ycoord_noedge = mod(n, clipImageHeight);

%% create a small gaussian filter
img_filter_gnoise = zeros(size(img_simple_threshold));
ycoor_filtered_gnoise = zeros(length(xcoord_noedge),1);
xcoor_filtered_gnoise = zeros(length(ycoord_noedge),1);
count = 0;
for index = 1:length(xcoord_noedge)
    ycoord2 = ycoord_noedge(index)-params.L_noise:ycoord_noedge(index)+params.L_noise;
    xcoord2 = xcoord_noedge(index)-params.L_noise:xcoord_noedge(index)+params.L_noise;
    temp = sum(sum(img_simple_threshold(ycoord2, xcoord2).*noisefilter));
    
    if (temp > params.threshold_range(2)*2)
        count = count + 1;
        img_filter_gnoise(ycoord_noedge(index), xcoord_noedge(index)) = img_simple_threshold(ycoord_noedge(index), xcoord_noedge(index));    
        ycoor_filtered_gnoise(count) = ycoord_noedge(index);
        xcoor_filtered_gnoise(count) = xcoord_noedge(index);
    end       
    
end
figure; imagesc(img_filter_gnoise); colormap gray; colorbar;
%% centroid x histogram
% [Nxelements,ycenters] = hist(xcoor_filtered_gnoise(1:count), 3);
% figure; hist(xcoor_filtered_gnoise(1:count), 3);

% dy = ycenters(2)-ycenters(1);
% xhistogram_edges = [ycenters-dy/2; ycenters+dy/2]';
% 
% Nxelements_index = find(Nxelements > 150);
% xpoints = zeros(sum(Nxelements(Nxelements_index)),1);
% xpoints_index = zeros(length(Nxelements_index),2);
% xpoints_index(1,1) = 1;
% xpoints_index(1,2) = Nxelements(Nxelements_index(1));
% for i = 1:length(Nxelements_index)-1
%     xpoints_index(i+1,1) = xpoints_index(i,2)+1;
%     xpoints_index(i+1,2) = xpoints_index(i+1,1) + Nxelements(Nxelements_index(i+1))-1;       
% end
% 
% for i = 1:length(Nxelements_index)
%    xpoints(xpoints_index(i,1):xpoints_index(i,2)) = find(xcoor_filtered_gnoise(1:count) > xhistogram_edges(Nxelements_index(i),1) & xcoor_filtered_gnoise(1:count) < xhistogram_edges(Nxelements_index(i),2)); 
%    
% end
% xhist_image = zeros(size(img_filter_gnoise));
% xhist_image(ycoor_filtered_gnoise(xpoints), xcoor_filtered_gnoise(xpoints)) = img_filter_gnoise(ycoor_filtered_gnoise(xpoints), xcoor_filtered_gnoise(xpoints));
% figure; imagesc(xhist_image); colormap gray; colorbar;
 
%% centroid y histogram
% xcoord = xcoor_filtered_gnoise(1:count);
% ycoord = ycoor_filtered_gnoise(1:count);
% 
% [~,ycenters] = hist(ycoord, params.N_steelballs);
% dy = ycenters(2)-ycenters(1);
% 
% z = zeros(size(img_filter_gnoise));
% 
% index = zeros(params.N_steelballs,1);
% R2 = 1;
% for n = 1:params.N_steelballs
%    bin_index = find(ycoord > (ycenters(n)-dy/2) & ycoord < (ycenters(n)+dy/2));
%    filtered_point = zeros(length(bin_index),1);
%    
%    for m = 1:length(bin_index)
%        x = xcoord(bin_index(m));
%        y = ycoord(bin_index(m));       
%        filtered_point(m) = sum(sum(img_filter_gnoise(y-params.L_bearing:y+params.L_bearing, x-params.L_bearing:x+params.L_bearing).*gaussianfilter2d));
%    end
%    
%    M = find(filtered_point == max(filtered_point));
%    index(n) = bin_index(round(median(M)));
%    z(ycoord(index(n))-R2:ycoord(index(n))+R2, xcoord(index(n))-R2:xcoord(index(n))+R2) = 300;
% end
% % figure; imagesc(z+img_filter_gnoise); colormap gray; colorbar; axis equal;
% figure(80); imagesc(z + clip_dat, [0 500]); colormap gray; colorbar;
% 
% img = z + clip_dat;
% u = xcoord(index);
% v = ycoord(index);
u = 0;
v = 0;
img = 1;