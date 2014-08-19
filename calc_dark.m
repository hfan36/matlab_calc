clear all
close all
clc
%% set up pathes
currentpath = cd;
newpath = 'H:\CT data\042414';
%making new changes on this part, it doesn't exist yet

% change to new path
cd(newpath);
%%
% pixels_width = 2160;
% pixels_height = 2560;
time_ms = [100 500 1000 1500 2000 3000 5000];
% N = 100;
% 
% data = zeros(pixels_width*pixels_height, N);
% var_data = zeros(pixels_width*pixels_height, length(time_ms));
% mean_data = zeros(pixels_width*pixels_height, length(time_ms));
% 
% for time_it = 1:length(time_ms)
%     for img_it = 1:N
%         
%         filename = strcat('dark_', num2str(time_ms(time_it)), 'ms_', num2str(img_it-1), '.spool');
%         dat = readBinary(filename, 2160*2560+1, 'uint16');
%         data(:,img_it) = dat(2:end);
%     end
%     var_data(:, time_it) = var(data, 0, 2);
%     mean_data(:, time_it) = mean(data, 2);
% end
% save('dark_variance_mean.mat', 'var_data', 'mean_data', 'pixels_width', 'pixels_height', 'time_ms', 'N');
%%
% load('dark_variance_mean.mat');
% window_x = 1000:1100;
% 
% m_data_window = zeros(length(window_x), length(window_x), length(time_ms));
% v_data_window = zeros(length(window_x), length(window_x), length(time_ms));
% for time_iteration = 1:length(time_ms)
%     data = reshape(mean_data(:,time_iteration), pixels_width, pixels_height);
%     m_data_window(:,:,time_iteration) = data(window_x, window_x);
%     
%     data = reshape(var_data(:, time_iteration), pixels_width, pixels_height);
%     v_data_window(:,:,time_iteration) = data(window_x, window_x);
% end
% save('window_dark_variance_mean.mat', 'm_data_window', 'v_data_window', 'window_x');

load('window_dark_variance_mean.mat');
m_data = zeros(length(window_x).^2, length(time_ms));
v_data = zeros(length(window_x).^2, length(time_ms));

for n = 1:length(time_ms)
   m_data(:,n) = reshape(m_data_window(:,:,n), length(window_x).^2,1); 
   v_data(:,n) = reshape(v_data_window(:,:,n), length(window_x).^2,1);
end


figure(1);
for n = 1:length(time_ms)
    subplot(3,3,n);
    imagesc(m_data_window(:,:,n)); colormap gray; axis equal; colorbar;
end
% 
% figure(2);
% for n = 1:length(time_ms)
%     subplot(3,3,n);
%     imagesc(v_data_window(:,:,n)); colormap gray; axis equal; colorbar;
% end
% 
% figure(3);
% imagesc(mean(v_data_window,3));  colormap gray; axis equal; colorbar;


figure;
plot(time_ms, m_data(509,:));


% figure; plot(mean(m_data,1), mean(v_data,1), '-x')
% 
% p = fit(mean(m_data,1)', mean(v_data,1)', 'poly1');
% 
% figure; plot(p, mean(m_data,1)', mean(v_data,1)')


%% change back to the matlab file path;
cd(currentpath);