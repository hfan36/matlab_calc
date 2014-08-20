clear all
close all
clc

%% read data file

%switch directory
% data_directory = 'H:\CT data\051314\cleaned_data';
% cd(data_directory);

image_data = readBinary('H:\CT data\051314\cleaned_data\clean_phantom_60kV_400uA_10sec_0.ct', 2160*2560, 'float');
image_data = reshape(image_data, 2560, 2160);
figure; imagesc(image_data); colormap gray; axis equal; colorbar;





small_image_data = readBinary('H:\CT data\051314\cleaned_data\C1024_clean_phantom_60kV_400uA_10sec_0.ct', 1024*1024, 'float');
small_image_data = reshape(small_image_data, 1024, 1024);
figure; imagesc(small_image_data'); colormap gray; axis equal; colorbar; title('center before');


image_data2 = image_data(769:769+1023, 569:569+1023);
figure; imagesc(10 - log(image_data2)); colormap gray; axis equal; colorbar; title('centered today');


load('H:\CT data\051514\cleaned_data\clean_mean_variance.mat');

%     t = (log(flood_mean) - log(dat));
%     t(t < 0.08) = 0;
t = log(flood_mean) - log(image_data(:));
t(t < 0.08) = 0;




sim_image_data = readBinary('H:\CT data\051314\simulation\fp_phantom_voxelvolume_10.bin', 1024*1024, 'float');
sim_image_data = reshape(sim_image_data, 1024, 1024);
figure; imagesc(sim_image_data'); colormap gray; axis equal; colorbar; title('simulation with calibration parameters');

%% read phantom voxelvolume
phantom = readBinary('H:\CT data\051314\simulation\phantom_voxelvolume.bin', 64*64*64, 'float');
phantom = reshape(phantom, 64, 64, 64);

% figure; imagesc(phantom(:,:,32)); colormap gray; axis equal; colorbar;

% cd ..