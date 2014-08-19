clear all
close all
clc
%% set up pathes
currentpath = cd;
newpath = 'H:\CT data\042414';

% change to new path
cd(newpath);
%% load original image and getting the part that I need
dat = readBinary('ruler.single', 2160*2560+1, 'uint16');
dat = reshape(dat(2:end), 2160, 2560);
dat_small = dat(1091:1684, 983:1550);

% display of this section 
% figure; imagesc(dat); colormap gray; axis equal; colorbar;
% title('Original image of the ruler');
% 
% figure; imagesc(dat_small); colormap gray; axis equal; colorbar;
% title('Clipped image of the ruler');
% 
% figure; plot(dat_small(100:550,:)');
% title('cross section of the important part');

%% fit, getting the angular deviation of the ruler
%this section first extract the two dips created by the 2 vertical lines on
%the ruler, then I fit a 6 order polynomial to the 'dip'/parabola and find
%the minimum of the 'dip' as well as it's location.  The location is
%recorded in xmin_left and xmin_right as locations, then I fitted another
%linear line to the deviation of the minimum as it is tracked over a bunch
%of vertical pixels, this deviation is then presented as an angle in
%degrees at the end.  It tells how much of an angle the ruler is placed
%with respect to the detector. As a result I can estimate the optical
%magnification (next section) a little bit more accurately.  Turns out if
%the deviation is too small, then it doesn't matter too much...

% fit_dat_left = dat_small(100:550, 237:256)';
% fit_dat_right = dat_small(100:550, 386:405)';
% fit_dat_x = 0:1:size(fit_dat_left,1)-1;
% 
% xmin_left = zeros(size(fit_dat_left,2),1);
% xmin_right = zeros(size(fit_dat_right,2),1);
%  
% for n = 1:size(fit_dat_left,2)
%    [xmin_left(n), xExa] = getPolyMin(fit_dat_x', fit_dat_left(:,n));
%    [xmin_right(n), xExb] = getPolyMin(fit_dat_x', fit_dat_right(:,n));
% end
% 
% pixel_index = (0:1:length(fit_dat_left)-1)';
% fit_line_left = fit(pixel_index, xmin_left, 'poly1');
% fit_line_right = fit(pixel_index, xmin_right, 'poly1');
% 
% angle_deg =( atan(fit_line_left.p1) + atan(fit_line_right.p1) )./2 * 180/pi;
% fprintf('deviation = %.3f degrees \n', angle_deg);
% save('angular_deviation.mat', 'angle_deg', 'xmin_left', 'xmin_right', 'fit_line_left', 'fit_line_right', 'fit_dat_left', 'fit_dat_right', 'pixel_index', 'fit_dat_x');

load('angular_deviation.mat');
% display part of this section
% figure;
% subplot(1,2,1);
% plot(fit_dat_left);
% subplot(1,2,2);
% plot(fit_dat_right);
% hold off;
% title('cross section of the 2 black lines');

figure;
subplot(1,2,1);
plot(fit_line_left, pixel_index, xmin_left);
subplot(1,2,2);
plot(fit_line_right, pixel_index, xmin_right);
title('fitted data');
%% getting the period of the ruler
% load('angular_deviation.mat');
rulerimage_mm = dat_small(:, 460:487);
figure; plot(rulerimage_mm);

N = size(rulerimage_mm,2);

ruler = zeros(size(rulerimage_mm,1), N);
for n = 1:N
   ruler(:,n) = rulerimage_mm(:, n);
end
ruler_ave = mean(ruler,2);
%Fourier Transform, doesn't work so hot...
% FT_ruler = fftshift(fft(ruler_ave));
% [FT_peaks, FT_locs] = findpeaks(abs(FT_ruler));
% FT = [FT_peaks FT_locs];
% [B, IX] = sort(FT, 1, 'descend');
% p = (IX(3,1) - IX(2,1))/2;
% figure;
% plot(FT);
% index = (1:length(ruler_ave));
% figure;
% plot(index, abs(FT_ruler));


x = linspace(0, 359, length(ruler_ave))*pi/180;
[pks, locs] = findpeaks(ruler_ave);

figure;
plot(x, ruler_ave); hold on;
plot(x(locs), pks, 'k^', 'markerfacecolor', [1 0 0]);
hold off;
%need to extract double peaks?

period_locs = zeros(length(locs)-1, 1);
for n = 1:length(locs)-1
    period_locs(n) = x(locs(n+1))-x(locs(n));
end

m = 2*pi./mean(period_locs);
approx2 = 500*cos(m*x-57)+6000;
figure;
plot(x(locs), 6500, 'k^', 'markerfacecolor', [1 0 0]); hold on;
plot(x, approx2, 'b'); hold off;
ylim([5000 7000]);
xlim([0 3]);

demag1 = 1./((1/m*cos(angle_deg*pi/180))*length(ruler_ave)*6.5e-3);
demag2 = 1./((1/m)*length(ruler_ave)*6.5e-3);
fprintf('demagnification = %.5f \n', demag1);
%% test FFT
% x_rad = (-720:1:720)*pi/180; %radians
% x_deg = x_rad*180/pi;
% 
% y = cos(3*x_rad);
% Y = fftshift(fft(y));
% 
% figure; plot(y);
% X = (1:1441)-721;
% figure; plot(X, abs(Y));
% 
% index_deg = linspace(0, 359, length(y));
% index_rad = index_deg*pi/180;
% 
% test = cos(12.025*index_rad);
% figure; plot(index_rad, test);
% 
% figure;
% plot(index_rad, y); hold on;
% plot(index_rad, test, 'r'); hold off;
% 
% disp(sum( (test-y).^2));
%% change back to the matlab file path;
cd(currentpath);