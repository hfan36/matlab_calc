clear all
close all
clc

%% set up pathes
currentpath = cd;
newpath = 'C:\Users\Ryan-Helen\Documents\CT data\042414';

% change to new path
cd(newpath);
%%

filename = strcat('callibration_50kV_5000ms_400uA_', num2str(0), '.ct');

dat = readBinary(filename, 2160*2560+1, 'uint16');
dat = dat(2:end);
dat = reshape(dat, 2560, 2160);

figure; imagesc(dat, [0 400]); colormap gray; axis equal; colorbar;
%%

NpixelsWidth = 2160; 
NpixelsHeight = 2560;
clip_height_top = 650;
clip_height_bottom = 80;
clip_width_left = 120;
clip_width_right = 80;
clipImageHeight = NpixelsHeight-clip_height_top-clip_height_bottom+1;
clipImageWidth = NpixelsWidth-clip_width_left-clip_width_right+1;
threshold_range = [90 130];
boundary_width = 20;
N_steelballs = 8;

L_noise = 2;
r_noise = 1;
A_noise = 1; 

r_bearing = 20;
L_bearing = 19;
%% plot filters
noise_filter = gaussian2D(L_noise, r_noise,1);
bearing_filter = gaussian2D(L_bearing, r_bearing, 1);

% figure;
% subplot(1,2,1); imagesc(noise_filter); colormap gray; axis equal;
% title('noise filter');
% subplot(1,2,2); imagesc(bearing_filter); colormap gray; axis equal;
% title('bearing filter');
%% load image

params = struct('NpixelsWidth', NpixelsWidth, ...
                'NpixelsHeight', NpixelsHeight, ...
                'clip_height_top', clip_height_top, ...
                'clip_height_bottom', clip_height_bottom, ...
                'clip_width_left', clip_width_left, ...
                'clip_width_right', clip_width_right, ...
                'boundary_width', boundary_width, ...
                'threshold_range', threshold_range, ...
                'N_steelballs', N_steelballs, ...
                'L_noise', L_noise, ...
                'r_noise', r_noise, ...
                'L_bearing', L_bearing, ...
                'r_bearing', r_bearing);

folder_root = 'C:\Users\Ryan-Helen\Documents\CT data\042414\';

u = zeros(N_steelballs, 180);
v = zeros(N_steelballs, 180);

n = 1;
filename = strcat(folder_root, 'callibration_50kV_5000ms_400uA_', num2str(n-1), '.ct');
% [u(:,n), v(:,n), img] = find_centroid2(params, filename, noise_filter, bearing_filter);


% load image
dat = readBinary(filename, params.NpixelsHeight*params.NpixelsWidth+1, 'uint16');
dat = reshape(dat(2:end), params.NpixelsHeight, params.NpixelsWidth);
figure; imagesc(dat, [0 500]); colormap gray; axis equal; colorbar;
xlim([1 NpixelsWidth]);
ylim([1 NpixelsHeight]);

binary_mask = zeros(size(dat));
binary_mask(params.clip_height_top:end-params.clip_height_bottom,params.clip_width_left:end-params.clip_width_right) = 1;

masked_image = binary_mask.*dat;
figure; imagesc(masked_image, [0 500]); colormap gray; colorbar; axis equal;
xlim([1 NpixelsWidth]);
ylim([1 NpixelsHeight]);

clear dat;
% simple threshold 
n = find(masked_image(:)<= params.threshold_range(2) & masked_image(:)>= params.threshold_range(1));
img_simple_threshold = zeros(size(masked_image));
img_simple_threshold(n) = masked_image(n);
figure; imagesc(img_simple_threshold); colormap gray; colorbar; axis equal;
  
img_simple_threshold = medfilt2(img_simple_threshold);
figure; imagesc(img_simple_threshold, [90 130]); colormap gray; colorbar; axis equal;
 
% n = find(img_simple_threshold(:)<= params.threshold_range(2) & img_simple_threshold(:)>= params.threshold_range(1));
% 
% % find x and y coordinates of filter1
% clipImageHeight = params.NpixelsHeight-params.clip_height_top-params.clip_height_bottom+1;
% xcoord_noedge = floor(n/clipImageHeight)+1;
% ycoord_noedge = mod(n, clipImageHeight);
% 
% % % create a small gaussian filter
% % img_filter_gnoise = zeros(size(img_simple_threshold));
% % ycoor_filtered_gnoise = zeros(length(xcoord_noedge),1);
% % xcoor_filtered_gnoise = zeros(length(ycoord_noedge),1);
% % count = 0;
% % 
% A = zeros(size(img_simple_threshold));
% for index = 1:length(xcoord_noedge)
%     A(ycoord_noedge(index), xcoord_noedge(index)) = 1;
% end
% figure; imagesc(A); colormap gray; axis equal; colorbar;
% title('binary');
% 
% M = [115 1806 801 806 802 802 804 806 806 802; 1765 1776 1687 1455 1241 1013 788 571 350 120]';
% figure; plot(xcoord_noedge, ycoord_noedge, 'o');
% [idx, ctrs] = kmeans([xcoord_noedge ycoord_noedge], 10, 'start', M);
% hold on;
% plot(ctrs(:,1), ctrs(:,2), 'g.', 'MarkerSize', 18); hold off;
% 



% for index = 1:length(xcoord_noedge)
%     ycoord2 = ycoord_noedge(index)-params.L_noise:ycoord_noedge(index)+params.L_noise;
%     xcoord2 = xcoord_noedge(index)-params.L_noise:xcoord_noedge(index)+params.L_noise;
%     temp = sum(sum(img_simple_threshold(ycoord2, xcoord2).*noise_filter));
%     
%     if (temp > params.threshold_range(2)*2)
%         count = count + 1;
%         img_filter_gnoise(ycoord_noedge(index), xcoord_noedge(index)) = img_simple_threshold(ycoord_noedge(index), xcoord_noedge(index));    
%         ycoor_filtered_gnoise(count) = ycoord_noedge(index);
%         xcoor_filtered_gnoise(count) = xcoord_noedge(index);
%     else
%     end       
%     
% end
% figure; imagesc(img_filter_gnoise); colormap gray; colorbar;


% % centroid x histogram
% [Nxelements,ycenters] = hist(xcoor_filtered_gnoise(1:count), 3);
% figure; hist(xcoor_filtered_gnoise(1:count), 3);
% 
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
%  
% % centroid y histogram
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
% u = 0;
% v = 0;
% img = 1;

%% change back to the matlab file path;
cd(currentpath);