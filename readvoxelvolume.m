clear all
close all
clc

% phantom1_1.raw
currentpath = cd;
cd('H:\CT data\051314\simulation\');


% S = 128;
% dat = readBinary('phantom1_4.raw', S^3, 'int8');
% dat = reshape(dat, S, S, S);
% % 
% dat_planar = shiftdim(dat(:,2,:),2);

% figure; imagesc(dat_planar); axis equal;
% dat_small = dat_planar(1:38, 1:38);


% figure; imagesc(dat_small); axis equal;

% out = padarray(dat_small, [13 13]);
% figure; imagesc(-out); axis equal;


% v = zeros(64,64,64);
% for n = 1:64
%     v(:,:,n) = -out;
% end
% figure; imagesc(v(:,:,1)); axis equal;

% writeBinary(v, 'phantom_voxelvolume.bin', 'float');

% 
fp = readBinary('fp_sim_phantom_0.bin', 512*512, 'float');
fp = reshape(fp, 512, 512);
figure; imagesc(fp'); axis equal; colormap gray;

% d = readBinary('d_detector_value.bin', 512*512, 'float');
% d = reshape(d, 512,512);
% figure; imagesc(d); axis equal; colormap gray; colorbar;

cd('H:\CT data\051314\simulation');

bp = readBinary('recon_raw_5.bin', 64*64*64, 'float');
bp = reshape(bp, 64,64,64);

% figure;
% imagesc(shiftdim(bp(32,:,:), 1)); axis equal; colormap gray; colorbar 
% 
% figure;
% imagesc(shiftdim(bp(:,32,:), 2)); axis equal; colormap gray; colorbar 
% 
% figure;
% imagesc(shiftdim(bp(:,:,32), 3)); axis equal; colormap gray; colorbar 


% for n = 1:1
% fp = readBinary(strcat('fp_sim_phantom_', num2str(n-1), '.bin'), 512*512, 'float');
% fp = reshape(fp, 512,512);
% 
% figure(12); imagesc(fp'); axis equal; colormap gray;
% end

cd('H:\CT data\051314\cleaned_data');
raw_fp = readBinary('C1024_clean_phantom_60kV_400uA_10sec_0.ct', 1024^2, 'float');
raw_fp = reshape(raw_fp, 1024, 1024);
figure; imagesc(raw_fp'); axis equal; colormap gray; 
% 
%%
cd('H:\CT data\051514\cleaned_data');
% flood_sum = zeros(2560*2160,1);
% 
% for n = 1:100
%     filename = strcat('clean_Flood_60kV_400uA_10sec_', num2str(n-1), '.seq');
%     flood = readBinary(filename, 2560*2160+1, 'float');
%     flood = flood(1:end);
%     flood_sum(:,n) = flood;
% end
% 
% flood_mean = mean(flood_sum,2);
% flood_var = var(flood_sum, 0 , 2);
% save('clean_mean_variance.mat', 'flood_mean', 'flood_var');

load('clean_mean_variance.mat');
% 
flood_mean = reshape(flood_mean, 2560, 2160);
% flood_var = reshape(flood_var, 2560, 2160);
% 
% figure;
% imagesc(flood_mean, [0 600]); colormap gray; axis equal;
% 
% % figure;
% % imagesc(flood_var, [0 600]); colormap gray; axis equal;
% 
% 
dat = readBinary('H:\CT data\051314\cleaned_data\clean_phantom_60kV_400uA_10sec_0.ct', 2560*2160, 'float');
dat = reshape(dat, 2560, 2160);

% flood1 = readBinary('H:\CT data\051514\Flood_60kV_400uA_10sec_0.seq', 2560*2160+1, 'uint16');
% flood1 = flood1(2:end);
% flood1 = reshape(flood1, 2560, 2160);

% figure;
% imagesc(flood1, [0 600]); colormap gray; axis equal;


% figure;
% imagesc(dat, [0 600]); colormap gray; axis equal;
% 


%%

for n = 1:1
    dat = readBinary(strcat('H:\CT data\051314\cleaned_data\clean_phantom_60kV_400uA_10sec_', num2str(n-1), '.ct'), 2560*2160, 'float');
    dat = reshape(dat, 2560, 2160);
    t = (log(flood_mean) - log(dat));
    t(t < 0.08) = 0;
% 
    figure;
    imagesc(t); colormap gray; axis equal; colorbar; title('t');
% 
%     figure;
%     t_small = t(1280-512:1280+511, 1080-512:1080+511);
%     imagesc(t_small); colormap gray; axis square; colorbar;
%     writeBinary(t_small(:), strcat('H:\Visual Studio 2010\CTSolution\Siddon\data\trimmed_data1024_', num2str(n-1), '.bin'), 'float');
end

% figure;
% imagesc(dat-flood1,  [-300 200]); colormap gray; axis equal;
%% 
% dat = readBinary('H:\CT data\051314\phantom_60kV_400uA_10sec_0.ct', 2560*2160, 'uint16');
% dat = reshape(dat, 2560, 2160);
% 
% figure; imagesc(dat, [0 800]); colormap gray; axis equal;
% 
% dat = readBinary('H:\CT data\051514\Flood_60kV_400uA_10sec_0.seq', 2560*2160, 'uint16');
% dat = reshape(dat, 2560, 2160);
% figure; imagesc(dat, [0 800]); colormap gray; axis equal;


cd(currentpath);
