% assume everything is in millimeters!!!!
% assume magnification is about 10 right now
% 2560 x 2160 at 6.5 microns pixel size

clear all
close all
clc
%%  folder, and plot raw data
%this is the raw data from experiment after centroid has been extracted
% folder_root = 'H:\Calibration\112613\';
% load(strcat(folder_root, 'uvdata.mat'));
load('H:\CT data\042414\calibration_042414.mat');
% load('C:\Users\Ryan-Helen\Documents\CT data\matlab_calc\calibration_042414.mat');
% u = u(:,:,2);
% v = -v(:,:,2)+2560;

u = u(:,1:2:end,2);
v = v(:,1:2:end,2);

figure;
plot(u', v', '.');

figure; 
v = -v+2560;
plot(u', v', '.');

% title('Raw data from experiment after centroid extraction');
% xlim([1 2160]);
% ylim([1 2560]);
%% Enter Detector values
DetectorWidth = 2560;
DetectorHeight = 2160;
PixelSize = 6.5e-3; % millimeters
Mag = 13.7; % magnification

%% contracting grid to find uc and vc for centering the image
uc0 = 1208; %1096; %detector shift in pixels (width)
vc0 = 1080; %996; %detector shift in pixels (height)
delta_uc = 5; %tolerance in width
delta_vc = 5; %tolderance in height
N_grid = 6; %number of grid points
Rate = 1.15; %contraction rate
N_loop = 50; %number of contractions
alpha = 2*pi-(0:2:358)*pi/180; %angles
alpha = alpha(1:2:end);
cg_var = make_CG_variable(N_grid, Rate, N_loop);
det_var = make_det_variable(DetectorWidth, DetectorHeight, PixelSize, Mag);
[u0, v0, result, ucenter, vcenter] = find_imagecenter(u, v, delta_uc, delta_vc, uc0, vc0, alpha, cg_var, det_var);

%% contracting grid to find Rf
%% first lets shorten u0 and v0 a bit since it's so long
theta0 = 0;
eta0 = 0; 
phi0 = 0;
% 
Nball = 8;
R_phantom = 60;
N_grid = 4;
dz0 = 0;
% 
xi = R_phantom*ones(Nball,1);
yi = sqrt(R_phantom^2-xi.^2);
Rf0 = round(15/abs(result.z0_Rf(5)-result.z0_Rf(4)));
z = 0:15:15*(Nball-1);
zi0 = (z+(result.z0_Rf(5)*Rf0 - z(5)))';
% 
dz0_Rf = zeros(Nball-1, 1);
for n = 1:Nball-1
   dz0_Rf(n) = result.z0_Rf(n+1)-result.z0_Rf(n);
end
Rf0 = 15/mean(abs(dz0_Rf));
% Rf0 = 898.5;
zi = round(result.z0_Rf*Rf0);
zi0 = (min(zi):15:15*(Nball-1)+min(zi))';

ri_Rf = [xi yi zi0]'/Rf0;
% 
% p = make_phantom_variable(0, 0, 0, 60, 0, 15, 8, result.z0_Rf);
 
%% setting up the tolerance for each variables
delta_dz = 10;
delta_Rf = 10;
delta_eta = 10 * pi/180;
delta_theta = 10 *pi/180;
delta_phi = 10 *pi/180;
%% find initial phi (rotation like CT)
% m = zeros(360,1);
% for phi0 = 1:360
%     [u_contract, v_contract] = CalibrationPlot_smekal(result, ri_Rf, theta0, eta0, phi0*pi/180, alpha);
%     figure(1);
%     title('find initial phi');
%     subplot(1,2,1);
%     for j = 1:Nball
%         plot(u0(j, :), v0(j, :), '.b'); hold on;
%         plot(u_contract(j,:), v_contract(j,:), '.r');
%     end
%     hold off
%     
%     subplot(1,2,2);
%     m(phi0) = CalibrationMSE(u0, v0, u_contract, v_contract);
%     plot(phi0, m(phi0), '.b'); hold on;
% end
% phi0 = find(min(m) == m)*pi/180;

%% find initial eta values
% eta_guessArray = getgridarray(0, 0.2*pi/180, 50);
% m_eta = zeros(length(eta_guessArray), 1);
% for eta_index = 1:length(eta_guessArray)
%     [u_contract, v_contract] = CalibrationPlot_smekal(result, ri_Rf, theta0, eta_guessArray(eta_index), phi0, alpha);
%     
%     figure(2);
%     title('find initial eta');
%     subplot(1,2,1);
%     plot(u0', v0', '.b'); hold on;
%     plot(u_contract', v_contract', '.r'); hold off;
%     
%     subplot(1,2,2);
%     m_eta(eta_index) = CalibrationMSE(u0, v0, u_contract, v_contract);
%     plot(eta_index, m_eta(eta_index), '.b'); hold on;
% end
% eta0 = eta_guessArray(min(m_eta) == m_eta);

%% find initial theta values
% theta_guessArray = getgridarray(0, pi/180, 360);
% m_theta = zeros(length(theta_guessArray),1);
% for theta_index = 1:length(theta_guessArray)
%    [u_contract, v_contract] = CalibrationPlot_smekal(result, ri_Rf, theta_guessArray(theta_index), eta0, phi0, alpha); 
%     
%    figure(3);
%    title('find initial theta');
%    subplot(1,2,1);
%    plot(u0', v0', '.b'); hold on;
%    plot(u_contract', v_contract', '.r'); hold off;
%    
%    subplot(1,2,2);
%    m_theta(theta_index) = CalibrationMSE(u0, v0, u_contract, v_contract);
%    plot(theta_index, m_theta(theta_index), '.b'); hold on; 
% end
% theta0 = theta_guessArray(min(m_theta) == m_theta);

%%

theta0 = 180.02 *pi/180;
eta0 = -1.3 * pi/180;
phi0 = 260 *pi/180;


%% actual contracting grid to find Rf
N_loop = 100;
rate = 1.05;
dz_guessArray = getgridarray(dz0, delta_dz, N_grid);
phi_guessArray = getgridarray(phi0, delta_phi, N_grid);
theta_guessArray = getgridarray(theta0, delta_theta, N_grid);
eta_guessArray = getgridarray(eta0, delta_eta, N_grid);
Rf_guessArray = getgridarray(Rf0, delta_Rf, N_grid);
mse_array = zeros(N_loop,1);

% fps= 0.5;
% outfile = sprintf('%s','contracting_grid_movie.avi');
% mov = VideoWriter(outfile, 'Motion JPEG AVI', 'FrameRate',fps,'Quality',10);

% nframes = 100;
% mov(1:nframes) = struct('cdata', [], 'colormap', []);

for loop_index = 1:N_loop
    zi = zi0 + dz_guessArray(1);
    ri_Rf0 = [xi yi zi]'./Rf_guessArray(1);
    [u_contract, v_contract] = CalibrationPlot_smekal(result, ri_Rf0, ...
    theta_guessArray(1), eta_guessArray(1), phi_guessArray(1), alpha);
    min_mse = CalibrationMSE(u0, v0, u_contract, v_contract);
    
    for dz_index = 1:N_grid
        for phi_index = 1:N_grid
            for theta_index = 1:N_grid
                for eta_index = 1:N_grid
                    for Rf_index = 1:N_grid
                        
                        zi = zi0 + dz_guessArray(dz_index);
                        ri_Rf = [xi yi zi]'/Rf_guessArray(Rf_index);
                        [u_contract, v_contract] = CalibrationPlot_smekal(result, ...
                         ri_Rf, theta_guessArray(theta_index), ...
                         eta_guessArray(eta_index), ...
                         phi_guessArray(phi_index), alpha); 
                     
                         arg = CalibrationMSE(u0, v0, u_contract, v_contract);
                         
                         
                             if arg <= min_mse
                                min_mse = arg;
                                min_Rf_index = Rf_index;
                                min_eta_index = phi_index;
                                min_theta_index = theta_index;
                                min_phi_index = phi_index;
                                min_dz_index = dz_index;
                                mse_array(loop_index) = min_mse;
                             end %end if statement
                             
                    end % end Rf_index
                end % end eta_index
            end %end theta_index
        end %end phi_index
    end %end dz_index
    zi = zi0 + dz_guessArray(min_dz_index);
    ri_Rf = [xi yi zi]'/Rf_guessArray(min_Rf_index);
    [u_min, v_min] = CalibrationPlot_smekal(result, ri_Rf, theta_guessArray(min_theta_index), ...
        eta_guessArray(min_eta_index), phi_guessArray(min_phi_index), alpha);
    
    figure(88)
    title('contracting grid result');
    subplot(1,2,1);
    plot(u0', v0', '.b'); hold on
    plot(u_min', v_min', '.r'); hold off;
    subplot(1,2,2);
    plot(loop_index, mse_array(loop_index), '.'); hold on;

%------------------------------------------------------------------------------------------------
% only used for making GIF video
%     figure(88)
%     set(gcf, 'Position', [100 100 1200 600]);
%     subplot(1,2,1);
%     plot(u0', v0', '.b'); hold on;
%     plot(u_min', v_min', '.r'); hold off;
%     xlim([-90 90]); ylim([-90 90]);
%     xlabel('mm'); ylabel('mm');
%     
%     subplot(1,2,2);
%     plot(loop_index, mse_array(loop_index), '.'); hold on;
%     xlim([0 N_loop]);
%     xlabel('number of iterations');
%     ylabel('sum of squared errors');
    
    
%     f = getframe(88);
%     im = frame2im(f);
%     [A, map] = rgb2ind(im, 256);
%     if loop_index == 1;
%         imwrite(A, map, 'contracting_grid.gif', 'gif', 'LoopCount', 1, 'DelayTime', 2);
%     else
%         imwrite(A, map, 'contracting_grid.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
%     end
%     
%     mov(loop_index) = f;
%------------------------------------------------------------------------------------------------    
    
    min_dz = dz_guessArray(min_dz_index);
    min_phi = phi_guessArray(min_phi_index);
    min_theta = theta_guessArray(min_theta_index);
    min_eta = eta_guessArray(min_eta_index);
    min_Rf = Rf_guessArray(min_Rf_index);
    
    dz_guessArray = getgridarray(min_dz, delta_dz/(rate^loop_index), N_grid);
    phi_guessArray = getgridarray(min_phi, delta_phi/(rate^loop_index) , N_grid);
    theta_guessArray = getgridarray(min_theta, delta_theta/(rate^loop_index) , N_grid);
    eta_guessArray = getgridarray(min_eta, delta_eta/(rate^loop_index), N_grid);
    Rf_guessArray = getgridarray(min_Rf, delta_Rf/(rate^loop_index) , N_grid); 
    fprintf('dz = %.2f (mm), phi = %.2f (deg), theta = %.2f (deg), eta = %.2f (deg) Rf = %.2f (mm), mse is %.5f \n', ...
        min_dz, min_phi*180/pi, min_theta*180/pi, min_eta*180/pi, min_Rf, mse_array(loop_index));
    
end
% mov = close(mov);
% movie2avi(mov, 'test.avi', 'compression', 'Cinepak', 'fps', 3, 'quality', 10);

figure;
plot(u0', v0', '.b'); hold on;
plot(u_min', v_min', '.r'); hold off;

save('calibration_values_dissertation.mat', 'min_Rf', 'min_eta', 'min_theta', 'min_phi', 'min_dz', 'u0', 'v0', 'u_min', 'v_min');
% % % save(strcat(folder_root, 'calibration_values.mat'), 'min_Rf', 'min_eta', 'min_theta', 'min_phi', 'min_dz', 'u0', 'v0');

%% display pretty figure for dissertation
clear all
close all
clc

load('calibration_values_dissertation.mat');

figure;
plot(u0(:), v0(:), '.b'); 
hold on;
plot(u_min(:), v_min(:), '.r'); hold off; axis square;
set(gca, 'FontSize', 12, 'ylim', [-80 80], 'xlim', [-90 90], 'XTick', -90:20:90, 'XMinorTick', 'off', 'PlotBoxAspectRatio', [1.5 1.5 1]);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 5, 5]); % white bckgr
xlabel('phosphor screen lateral dimension (mm)'); ylabel('phosphor screen vertical dimension (mm)');
legend('points from experiment', 'points using calibration parameters', ...
       'Location', 'northoutside');
name = 'calibration_plot';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi
% saveas(gcf, strcat(name, '.fig'));