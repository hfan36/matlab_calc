clear all
close all
clc

currentpath = cd;
newpath = 'H:\CT data\042414';
% change to new path
cd(newpath);
%%
%%
NpixelsWidth = 2160; 
NpixelsHeight = 2560;
clip_height_top = 732;
clip_height_bottom = 2370;
clip_width_left = 870;
clip_width_right = 990;
threshold_range = [100 130];
N_steelballs = 8;

L_noise = 2;
r_noise = 1;
A_noise = 1; 

r_bearing = 15;
L_bearing = 20;

xc = zeros(8,180);
yc = zeros(8,180);

M = [923 926 923 920 923 923 923 920; 765 989 1216 1439 1660 1881 2108 2322]';
%% plot filters
noise_filter = gaussian2D(L_noise, r_noise,1);
bearing_filter = gaussian2D(L_bearing, r_bearing, 1);

binary_mask = zeros(NpixelsHeight, NpixelsWidth);
binary_mask(clip_height_top:clip_height_bottom,clip_width_left:clip_width_right) = 1;

mask_w = clip_width_right-clip_width_left;
mask_h = clip_height_bottom-clip_height_top;

cmass_x = zeros(3,1);
cmass_y = zeros(3,1);
%%
for image_index = 1:4
    filename = strcat('callibration_50kV_5000ms_400uA_', num2str(image_index-1), '.ct');
    dat = readBinary(filename, 2160*2560+1, 'uint16');
    dat = dat(2:end);
    dat = reshape(dat, NpixelsHeight, NpixelsWidth);

    %%
    % load image
    dat = readBinary(filename, NpixelsHeight*NpixelsWidth+1, 'uint16');
    dat = reshape(dat(2:end), NpixelsHeight, NpixelsWidth);
    figure; imagesc(dat, [0 500]); colormap gray; axis equal; colorbar;
    xlim([1 NpixelsWidth]);
    ylim([1 NpixelsHeight]);

    masked_image = medfilt2(dat.*binary_mask, [5 5]);
    figure(1); 
    imagesc(masked_image, [0 500]); colormap gray; axis equal; colorbar;
    xlim([1 NpixelsWidth]);
    ylim([1 NpixelsHeight]);

    n = find(masked_image(:)>= threshold_range(1) & masked_image(:) <= threshold_range(2));
    threshold = zeros(size(masked_image));
    threshold(n) = 1./(masked_image(n)/threshold_range(2));

    figure(2); imagesc(threshold);  colormap gray; axis equal; colorbar;
    xlim([1 NpixelsWidth]);
    ylim([1 NpixelsHeight]);

    xcoord_noedge = floor(n/NpixelsHeight)+1;
    ycoord_noedge = mod(n, NpixelsHeight);
    
    [idx center] = kmeans([xcoord_noedge, ycoord_noedge], 8, 'start', M);

    figure(3); plot(xcoord_noedge, ycoord_noedge, '.b'); hold on;
    plot(center(:,1), center(:,2), '.r', 'MarkerSize', 14);
    xlim([1 NpixelsWidth]);
    ylim([1 NpixelsHeight]);

    cmap = hsv(8);
    B = sum(sum(bearing_filter));

    for m = 1:8

        x = xcoord_noedge(idx == m);
        y = ycoord_noedge(idx == m);

        B = 0;   
        for a = 1:length(x)
            A = sum(sum(threshold(y(a)-L_bearing:y(a)+L_bearing, x(a)-L_bearing:x(a)+L_bearing).*bearing_filter));

            if (A >= B)
                xc(m, image_index) = x(a);
                yc(m, image_index) = y(a);
                B = A;
            end
        end    

    end

    figure(4); 
    plot(xcoord_noedge, ycoord_noedge, '.b'); hold on;
    plot(center(:,1), center(:,2), '.r', 'MarkerSize', 14);
    plot(xc(:,1:image_index), yc(:,1:image_index), '.g', 'MarkerSize', 14); hold off;
    xlim([1 NpixelsWidth]);
    ylim([1 NpixelsHeight]);

    M = [xc(:, image_index) yc(:, image_index)];
    
    cmass_x(image_index) = mean(xc(:, image_index));
    cmass_y(image_index) = mean(yc(:, image_index));
    

    
    clip_width_left = floor(cmass_x(image_index))-floor(mask_w/2);
    clip_height_top = floor(cmass_y(image_index))-floor(mask_h/2);
    
    binary_mask = zeros(NpixelsHeight, NpixelsWidth);
    binary_mask(clip_height_top:clip_height_top+mask_h,clip_width_left:clip_width_left+mask_w) = 1;
   
end


%% change back to the matlab file path;
cd(currentpath);