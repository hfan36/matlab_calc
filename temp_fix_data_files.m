clear all
close all
clc

%% read data file

%%%% load average flood data 
% load('H:\CT data\051514\cleaned_data\clean_mean_variance.mat');

%%%% This section is used to clip the clean phantom images to a certain size,
%%%% subtract flood correction and use them to test the reconstruction
%%%% algorithm, or actually test to see how the reconstruction looks like
%%%% because I think the actual reconstruction algorithm works just fine....
% for n = 1:180
%    image_data_folder = 'H:\CT data\051314\cleaned_data\';
%    image_data_filename = strcat('clean_phantom_60kV_400uA_10sec_', num2str(n-1), '.ct');
%    image_data = readBinary(strcat(image_data_folder, image_data_filename), 2160*2560, 'float');
%    flood_corrected = log(flood_mean) - log(image_data);
%    flood_corrected (flood_corrected < 0.08) = 0;
%    flood_corrected = reshape(flood_corrected, 2560, 2160);
%    dat = flood_corrected(769:769+1023, 569:569+1023)';
%    
%    save_file_name = strcat(image_data_folder, 'C1024_2_', image_data_filename);
%    writeBinary(dat(:), save_file_name, 'float');
% end

% figure; imagesc(reshape(dat, 1024, 1024)); colormap gray; axis equal;

%% reading the reconstructed images
xdim = 512;
ydim = 512;
zdim = 10;
slice = 4;
recon_data = readBinary('H:\CT data\051314\cleaned_data\recon_phantom_VOX015_512_s10_9.bin', xdim*ydim*zdim, 'float');
recon_data = reshape(recon_data, xdim, ydim, zdim);
figure; imagesc(recon_data(:,:,slice)'); colormap gray; axis equal; colorbar;

%% 
% test = readBinary('H:\Visual Studio 2010\CTSolution\Siddon\data\recon_trimmed_data1024_9.bin', 128*128*128, 'float');
% test = reshape(test, 128, 128, 128);
% figure; imagesc(test(:,:,3)); colormap gray; axis equal; colorbar;

%% just looking at the simulation images
% sim_image_data = readBinary('H:\CT data\051314\simulation\fp_phantom_voxelvolume_10.bin', 1024*1024, 'float');
% sim_image_data = reshape(sim_image_data, 1024, 1024);
% figure; imagesc(sim_image_data'); colormap gray; axis equal; colorbar; title('simulation with calibration parameters');

%% read phantom voxelvolume
% phantom = readBinary('H:\CT data\051314\simulation\phantom_voxelvolume.bin', 64*64*64, 'float');
% phantom = reshape(phantom, 64, 64, 64);
% figure; imagesc(phantom(:,:,32)); colormap gray; axis equal; colorbar;



