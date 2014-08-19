%older code, don't use, it' probably don't work or not well
function [u, v] = find_centroid(params, filename, noisefilter, gaussianfilter2d)

%% load image
dat = readBinary(filename, params.NpixelsHeight*params.NpixelsWidth, 'uint16');
dat = reshape(dat, params.NpixelsWidth, params.NpixelsHeight)';
%     figure(1); imagesc(dat, [0 500]); colormap gray; axis equal; 

%% clip the original image
clip_dat = dat(params.clip_height_top:end-params.clip_height_bottom,:);
% figure; imagesc(clip_dat, [0 500]); colormap gray; colorbar; axis equal;
clear dat;
dat_noedge = zeros(size(clip_dat));
dat_noedge(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width) = ...
           clip_dat(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width); 
       
%% simple threshold 
n = find(dat_noedge(:)<= params.threshold_range(2) & dat_noedge(:)>= params.threshold_range(1));
img_simple_threshold = zeros(size(clip_dat));
img_simple_threshold(n) = clip_dat(n);
figure; imagesc(img_simple_threshold); colormap gray; colorbar;
  
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
% figure; imagesc(img_filter_gnoise); colormap gray; colorbar;

%% centroid x histogram
[Nxelements,ycenters] = hist(xcoor_filtered_gnoise(1:count), 8);
% figure; bar(ycenters,Nxelements);
dy = ycenters(2)-ycenters(1);
xhistogram_edges = [ycenters-dy/2; ycenters+dy/2]';

Nxelements_index = find(Nxelements > 150);
xpoints = zeros(sum(Nxelements(Nxelements_index)),1);
xpoints_index = zeros(length(Nxelements_index),2);
xpoints_index(1,1) = 1;
xpoints_index(1,2) = Nxelements(Nxelements_index(1));
for i = 1:length(Nxelements_index)-1
    xpoints_index(i+1,1) = xpoints_index(i,2)+1;
    xpoints_index(i+1,2) = xpoints_index(i+1,1) + Nxelements(Nxelements_index(i+1))-1;       
end

for i = 1:length(Nxelements_index)
   xpoints(xpoints_index(i,1):xpoints_index(i,2)) = find(xcoor_filtered_gnoise > xhistogram_edges(Nxelements_index(i),1) & xcoor_filtered_gnoise <= xhistogram_edges(Nxelements_index(i),2)); 
end
xhist_image = zeros(size(img_filter_gnoise));
xhist_image(ycoor_filtered_gnoise(xpoints), xcoor_filtered_gnoise(xpoints)) = img_filter_gnoise(ycoor_filtered_gnoise(xpoints), xcoor_filtered_gnoise(xpoints));
% figure; imagesc(xhist_image); colormap gray; colorbar;
 
%% centroid y histogram
xcoord = xcoor_filtered_gnoise(xpoints);
ycoord = ycoor_filtered_gnoise(xpoints);

[~,ycenters] = hist(ycoord, params.N_steelballs);
dy = ycenters(2)-ycenters(1);

z = zeros(size(xhist_image));

index = zeros(params.N_steelballs,1);
R2 = 1;
for n = 1:params.N_steelballs
   bin_index = find(ycoord > (ycenters(n)-dy/2) & ycoord < (ycenters(n)+dy/2));
   filtered_point = zeros(length(bin_index),1);
   
   for m = 1:length(bin_index)
       x = xcoord(bin_index(m));
       y = ycoord(bin_index(m));       
       filtered_point(m) = sum(sum(xhist_image(y-params.L_bearing:y+params.L_bearing, x-params.L_bearing:x+params.L_bearing).*gaussianfilter2d));
   end
   
   M = find(filtered_point == max(filtered_point));
   index(n) = bin_index(round(median(M)));
   z(ycoord(index(n))-R2:ycoord(index(n))+R2, xcoord(index(n))-R2:xcoord(index(n))+R2) = 300;
end
% figure; imagesc(z+xhist_image); colormap gray; colorbar; axis equal;
figure(80); imagesc(z + clip_dat, [0 500]); colormap gray; colorbar; axis equal;

u = xcoord(index);
v = ycoord(index);
