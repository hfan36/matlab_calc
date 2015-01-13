clear all
close all
clc
%% set up pathes
currentpath = cd;
% newpath = 'H:\CT data\042414'; %original

newpath = 'H:\Images\010915';

% change to new path
cd(newpath);
%% load original image and getting the part that I need
% dat = readBinary('ruler.single', 2160*2560+1, 'uint16'); %original
% dat = reshape(dat(2:end), 2160, 2560); %original
% dat_small = dat(1091:1684, 983:1550); %original

dat = readBinary('farthest_away_ruler.seq', 2160*2560, 'uint16');
dat = reshape(dat, 2160, 2560);


% display of this section 
figure(88); imagesc(dat, [0 2000]); colormap gray; axis equal; colorbar;
title('Original image of the ruler');

dat_small = dat(900:1150, 800:1170)';
figure; imagesc(dat_small); colormap gray; axis equal; colorbar;
title('Clipped image of the ruler');

% figure; plot(dat_small(100:550,:)'); %original
% title('cross section of the important part'); %original

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

% fit_dat_left = dat_small(100:550, 237:256)'; %original
% fit_dat_right = dat_small(100:550, 386:405)'; %original
% fit_dat_x = 0:1:size(fit_dat_left,1)-1; %original

fit_dat_left = dat_small(50:250, 72:95)';
fit_dat_right = dat_small(50:250, 162:185)';
fit_dat_x = 0:1:size(fit_dat_left,1)-1;
 
xmin_left = zeros(size(fit_dat_left,2),1);
xmin_right = zeros(size(fit_dat_right,2),1);
 
for n = 1:size(fit_dat_left,2)
   [xmin_left(n), xExa] = getPolyMin_edit(fit_dat_x', fit_dat_left(:,n), 4);
   [xmin_right(n), xExb] = getPolyMin_edit(fit_dat_x', fit_dat_right(:,n), 4);
end

pixel_index = (0:1:length(fit_dat_left)-1)';
fit_line_left = fit(pixel_index, xmin_left, 'poly1');
fit_line_right = fit(pixel_index, xmin_right, 'poly1');

angle_deg =( atan(fit_line_left.p1) + atan(fit_line_right.p1) )./2 * 180/pi;
fprintf('deviation = %.3f degrees \n', angle_deg);

% save('angular_deviation.mat', 'angle_deg', 'xmin_left', 'xmin_right', 'fit_line_left', 'fit_line_right', 'fit_dat_left', 'fit_dat_right', 'pixel_index', 'fit_dat_x');

% load('angular_deviation.mat');
% display part of this section
figure(1);
subplot(1,2,1);
plot(fit_dat_left);
subplot(1,2,2);
plot(fit_dat_right);
hold off;
title('cross section of the 2 black lines');

figure(2);
subplot(1,2,1);
plot(fit_line_left, pixel_index, xmin_left);
subplot(1,2,2);
plot(fit_line_right, pixel_index, xmin_right);
title('fitted data');

% figure(2);
% set(gcf, 'Position', [100 100 340 250]);
% plot(fit_line_left, pixel_index, xmin_left);
% set(gca, 'FontSize', 8, 'xlim', [0 460], 'ylim', [6.5 10.5]);
% xlabel('shift in x direction (mm)');
% ylabel('shift in y direction (mm)');
% legend('extracted peak locations', 'Fitted line');
% set(gcf, 'Color', 'white');
% export_fig(gcf, ...      % figure handle
%     'calc_opticalmag',... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', '-eps', ...           % file format
%     '-r100' );             % resolution in dpi

% figure(1);
% set(gcf, 'Position', [100 100 450 350]);
% plot(Dwb/D_50, P_wb(1,:), '-k'); hold on;
% plot(Dwb/D_50, P_wb(2,:), '--k');
% plot(Dwb/D_50, P_wb(3,:), ':k'); hold off;
% set(gca, 'FontSize', 12, 'ylim', [0 1], 'xlim', [0 2], 'XTick', 0:0.4:2, 'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.5 1.2 1]);
% ylabel('Probability of a tissue-reaction effect, P'); 
% xlabel('Dose, D/D_{50}');
% legend('V = 2', 'V = 5', 'V = 10', 'Location', 'NorthWest');
% text(0.04, 0.02, 'T');
% text(0.2, 0.02, 'T');
% text(0.3, 0.02, 'T');
% name1 = 'Weibull-P_vs_D';
% set(gcf, 'Color', 'white'); % white bckgr
% export_fig(gcf, ...      % figure handle
%     name1,... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', '-eps', ...           % file format
%     '-r100' );             % resolution in dpi
% % saveas(gcf, strcat(name1, '.fig'));

%% getting the period of the ruler; mm side
close all; clc;
% load('angular_deviation.mat');
% rulerimage_mm = dat_small(50:end, 215:220);

rulerimage_mm = dat(1113:1114, 836:1242)';
% rulerimage_mm = dat_small(:, 460:487); %original 

figure; plot(rulerimage_mm); xlim([0 350]);

N = size(rulerimage_mm,2);

ruler = zeros(size(rulerimage_mm,1), N);
ft_ruler = zeros(size(ruler));
x_ruler = zeros(size(ruler));
for n = 1:N
   ruler(:,n) = rulerimage_mm(:, n);
   x_ruler(:,n) = (0:1:size(rulerimage_mm,1)-1) + (n-1)*tan(angle_deg*pi/180);
   ft_ruler(:,n) = fftshift(fft(ruler(:,n)));
end

ft_ruler = abs(ft_ruler);
figure; plot(ft_ruler);

delta = zeros(N,1);
for n = 1:N
    [FT_peaks, FT_locs] = findpeaks(ft_ruler(:,n));
    A = sortrows([FT_peaks FT_locs]);
    delta(n) = A(end-1, 2) - A(end, 2);
end

delta = mean(delta);
k = size(ruler,1)/(delta);

x_ft = linspace(0, 359, size(x_ruler,1))*pi/180;
figure; plot(x_ft, ruler(:,2));
hold on;
y_ft = 100*cos(2*pi*k*x_ft) + 1200;
plot(x_ft, y_ft, 'k', 'LineWidth', 2);
plot(x_ft, y_ft, 'k');
hold off;

period = k*6.5e-3;
fprintf('period = %.3f mm \n', period);
fprintf('demagnification = %.3f \n', 1/period);



% % %Fourier Transform, doesn't work so hot...
% % FT_ruler = fftshift(fft(ruler_ave));
% % [FT_peaks, FT_locs] = findpeaks(abs(FT_ruler));
% % FT = [FT_peaks FT_locs];
% % [B, IX] = sort(FT, 1, 'descend');
% % p = (IX(3,1) - IX(2,1))/2;
% % figure;
% % plot(FT);
% % index = (1:length(ruler_ave));
% % figure;
% % plot(index, abs(FT_ruler));
% % 
% % 
% x = linspace(0, 359, length(ruler_ave))*pi/180;
% % [pks, locs] = findpeaks(ruler_ave);
% % 
% % figure;
% % plot(x, ruler_ave); hold on;
% % plot(x(locs), pks, 'k^', 'markerfacecolor', [1 0 0]);
% % hold off;
% % %need to extract double peaks?
% % 
% % period_locs = zeros(length(locs)-1, 1);
% % for n = 1:length(locs)-1
% %     period_locs(n) = x(locs(n+1))-x(locs(n));
% % end
% % 
% % m = 2*pi./mean(period_locs);
% m = 2*pi*322/46;
% approx2 = 500*cos(m*x-57)+1200;
% figure;
% % plot(x(locs), 6500, 'k^', 'markerfacecolor', [1 0 0]); hold on;
% plot(linspace(0, max(x), length(x)), ruler_ave, 'r'); hold on;
% plot(x, approx2, 'b'); hold off;
% % ylim([5000 7000]);
% % xlim([0 3]);
% % 
% demag1 = 1./((1/m*cos(angle_deg*pi/180))*length(ruler_ave)*6.5e-3);
% demag2 = 1./((1/m)*length(ruler_ave)*6.5e-3);
% fprintf('demagnification = %.5f \n', demag1);
%% test FFT
% close all
% clc
% 
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