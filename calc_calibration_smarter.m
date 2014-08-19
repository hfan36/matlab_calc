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
Bearing_window = 120;
N_bearing = 8;

%%
L_noise = 2;
r_noise = 1;
A_noise = 1; 

r_bearing = 10;
L_bearing = 20;

M = [923 926 923 920 923 923 923 920; 765 989 1216 1439 1660 1881 2108 2322]';
bearing_filter = gaussian2D(L_bearing, r_bearing, 1);

Points = zeros(Bearing_window+1, Bearing_window+1, N_bearing);


Q = 180;
xc = zeros(N_bearing,Q);
yc = zeros(N_bearing,Q);

u = zeros(8, Q, 2);
v = zeros(8, Q, 2);
%%

% for iteration = 1:2
%     
%     for image_index = 1:Q
%         
%     filename = strcat('callibration_50kV_5000ms_400uA_', num2str(image_index-1), '.ct');
%     dat = readBinary(filename, 2160*2560+1, 'uint16');
%     dat = dat(2:end);
%     dat = reshape(dat, NpixelsHeight, NpixelsWidth);
% 
%     %%
%     % load image
%     dat = readBinary(filename, NpixelsHeight*NpixelsWidth+1, 'uint16');
%     dat = reshape(dat(2:end), NpixelsHeight, NpixelsWidth);
%     figure(1); imagesc(dat, [0 500]); colormap gray; axis equal; colorbar;
%     xlim([1 NpixelsWidth]);
%     ylim([1 NpixelsHeight]);
% 
%     % pick the window
%     
%     if iteration == 1
%     else
%         M = [u(:, image_index) v(:, image_index)];
%     end
% %         figure(2);
%         for n = 1:N_bearing
%             Points(:,:,n) = medfilt2(dat( M(n,2)-floor(Bearing_window/2): M(n,2)+floor(Bearing_window/2), M(n,1)-floor(Bearing_window/2):M(n,1)+floor(Bearing_window/2) ));
% %             subplot(3,3,n);
% %             imagesc(Points(:,:,n)); axis equal;
% %             xlim([1 Bearing_window]); ylim([1 Bearing_window]);
%         end
%         hold off;
% 
%         pad_p = zeros(Bearing_window+2*L_bearing+1, Bearing_window+2*L_bearing+1);
%         pad_p(L_bearing+1:L_bearing+1+Bearing_window, L_bearing+1:L_bearing+1+Bearing_window) = Points(:,:,1);
% 
%         test = zeros(size(Points));
%         for n = 1:N_bearing
%             pad_p = zeros(Bearing_window+2*L_bearing+1, Bearing_window+2*L_bearing+1);
%             pad_p(L_bearing+1:L_bearing+1+Bearing_window, L_bearing+1:L_bearing+1+Bearing_window) = Points(:,:,n)-mean(mean(Points(:,:,n)));
% 
%             T = 0;
%             for mx = 1:Bearing_window+1
%                 for my = 1:Bearing_window+1
% 
%                     x = L_bearing + mx;
%                     y = L_bearing + my;
%                     test(my, mx, n) = sum(sum(pad_p(y-L_bearing:y+L_bearing, x-L_bearing:x+L_bearing).*-bearing_filter));
% 
%                     if (test(my,mx,n) >= T)
%                         yc(n, image_index) = my;
%                         xc(n, image_index) = mx;
%                         T = test(my,mx,n);
%                     end
% 
%                 end
%             end
% 
%             dx = xc(:,2)-xc(:,1);
%             dy = yc(:,2)-yc(:,1);   
%         end  %end of 1:N_bearing
% 
%         M = M-floor(Bearing_window/2) + [xc(:,image_index) yc(:,image_index)]-1;
%         
%         u(:, image_index, iteration) = M(:,1);
%         v(:, image_index, iteration) = M(:,2);
% 
% 
%         A = zeros(size(Points));
%         A_large = zeros(size(dat));
% 
% %         figure(3);
%         for n = 1:N_bearing
%             A(yc(n,image_index), xc(n,image_index), n) = 1000;
%             A_large(M(n,2), M(n,1)) = 200;
% %             subplot(3,3,n);
% %             imagesc(test(:,:,n)+A(:,:,n)); axis equal;
% %             xlim([1 Bearing_window]); ylim([1 Bearing_window]);
%         end
%         hold off;
% 
%         figure(4); 
%         for n = 1:N_bearing
%            subplot(3,3,n);
%            imagesc(Points(:,:,n)+A(:,:,n)/20); axis equal;
%            xlim([1 Bearing_window]); ylim([1 Bearing_window]);
%         end
% 
%         fprintf('image index is = %d \n', image_index);
%         
% 
%     end
%     
%     
%     
% end %end of iteration

% save('calibration_042414.mat', 'u', 'v');
%%
%% change back to the matlab file path;
filename = strcat('callibration_50kV_5000ms_400uA_0.ct');
img = readBinary(filename, 2160*2560+1, 'uint16');
img = img(2:end);
img = reshape(img, NpixelsHeight, NpixelsWidth);

figure; imagesc(img, [0 500]); colormap gray; axis equal;

load('calibration_042414.mat');


u = u(:,:,2);
v = v(:,:,2);
%% Enter Detector values
DetectorWidth = 2560;
DetectorHeight = 2160;
PixelSize = 6.5e-3; % millimeters
Mag = 13.10; % magnification
Nball = 8;
% 
%% Smekal... this method sucks!
uc0 = 997.61; %detector shift in pixels (width)
vc0 = 1497.97; %detector shift in pixels (height)

u = (u - uc0)*PixelSize*Mag;
v = (v - vc0)*PixelSize*Mag;

figure; plot(u', v', '.'); title('scaled experimental raw image');



%%
rotation_angle_alpha = 2*pi-(0:2:358)*pi/180;
z = (0:15:8*15-15)-50;
Nangles = length(rotation_angle_alpha);
%% Fourier Fit using the paper
U0 = 2/(length(rotation_angle_alpha))*sum(u,2);
U1 = 2/(length(rotation_angle_alpha))*u*cos(rotation_angle_alpha');
U2 = 2/(length(rotation_angle_alpha))*u*cos(2*rotation_angle_alpha');
U3 = 2/(length(rotation_angle_alpha))*u*cos(3*rotation_angle_alpha');
U1_t = 2/(length(rotation_angle_alpha))*u*sin(rotation_angle_alpha');
U2_t = 2/(length(rotation_angle_alpha))*u*sin(2*rotation_angle_alpha');
U3_t = 2/(length(rotation_angle_alpha))*u*sin(3*rotation_angle_alpha');

V0 = 2/(length(rotation_angle_alpha))*sum(v,2);
V1 = 2/(length(rotation_angle_alpha))*v*cos(rotation_angle_alpha');
V2 = 2/(length(rotation_angle_alpha))*v*cos(2*rotation_angle_alpha');
V3 = 2/(length(rotation_angle_alpha))*v*cos(3*rotation_angle_alpha');
V1_t = 2/(length(rotation_angle_alpha))*v*sin(rotation_angle_alpha');
V2_t = 2/(length(rotation_angle_alpha))*v*sin(2*rotation_angle_alpha');
V3_t = 2/(length(rotation_angle_alpha))*v*sin(3*rotation_angle_alpha');

%% see if the U and V's are a good fit
u_FT_fit = zeros(length(z), length(rotation_angle_alpha));
v_FT_fit = u_FT_fit;
 
for n = 1:Nball
    u_FT_fit(n,:) = U0(n)/2 +  U1(n)*cos(rotation_angle_alpha) + U1_t(n)*sin(rotation_angle_alpha) + ...
                           U2(n)*cos(2*rotation_angle_alpha) + U2_t(n)*sin(2*rotation_angle_alpha) + ...
                           U3(n)*cos(3*rotation_angle_alpha) + U3_t(n)*sin(3*rotation_angle_alpha);
    v_FT_fit(n,:) = V0(n)/2 +  V1(n)*cos(rotation_angle_alpha) + V1_t(n)*sin(rotation_angle_alpha) + ...
                           V2(n)*cos(2*rotation_angle_alpha) + V2_t(n)*sin(2*rotation_angle_alpha) + ...
                           V3(n)*cos(3*rotation_angle_alpha) + V3_t(n)*sin(3*rotation_angle_alpha);
end
% 


figure;
plot(u', v', 'b');
hold on;
plot(u_FT_fit', v_FT_fit', 'r');
hold off;
title('data fit using Fourier series method (paper)');
fprintf('fit = %.5f \n', (sum(sum((u-u_FT_fit).^2 + (v-v_FT_fit).^2))));

%% the d's
d22_rho = (U3.^2 - U1.^2 + U3_t.^2 - U1_t.^2)/2;
d20_p = ((U1-U3).*U2 - (U3_t-U1_t).*U2_t)./d22_rho;
d21_p = ((U1+U3).*U2_t - (U3_t+U1_t).*U2)./d22_rho;

d00_p = ((U0+U2).*d20_p + U2_t.*   d21_p + 2*U1)./2;
d01_p = (U2_t.*   d20_p + (U0-U2).*d21_p + 2*U1_t)./2;
d02_p = (U1.*     d20_p + U1_t   .*d21_p + U0)./2;

d10_p = ((V0+V2).*d20_p + V2_t.*   d21_p + 2*V1)./2;
d11_p = (V2_t.*   d20_p + (V0-V2).*d21_p + 2*V1_t)./2;
d12_p = (V1.*     d20_p + V1_t.*   d21_p + V0)./2; 
 
%% the c's
c00_p = zeros(length(z),1);
c01_p = zeros(length(z),1);
c10_p = zeros(length(z),1);
c11_p = zeros(length(z),1);

for n = 1:length(z)
    denominator = 1./(d20_p(n).^2 + d21_p(n).^2);
    matrix = [d20_p(n) d21_p(n); d21_p(n) -d20_p(n)];
    c1 = 1/2*[U2(n) U2_t(n); U2_t(n) -U2(n)]*[d20_p(n); d21_p(n)] + [U1(n); U1_t(n)];
    c2 = 1/2*[V2(n) V2_t(n); V2_t(n) -V2(n)]*[d20_p(n); d21_p(n)] + [V1(n); V1_t(n)];
    c00_c01 = denominator*matrix*c1 + [U0(n)/2; 0];
    c10_c11 = denominator*matrix*c2 + [V0(n)/2; 0];
    
    c00_p(n) = c00_c01(1);
    c01_p(n) = c00_c01(2);
    c10_p(n) = c10_c11(1);
    c11_p(n) = c10_c11(2);
end

%% eta
eta = atan2(-c11_p./c01_p,1);

%% A, B, C, E, and F
A = sin(eta).*c00_p + cos(eta).*c10_p;
B = cos(eta).*c01_p - sin(eta).*c11_p;
C = cos(eta).*c00_p - sin(eta).*c10_p;
E = sin(eta).*d02_p + cos(eta).*d12_p;
F = cos(eta).*d02_p - sin(eta).*d12_p;

%% theta equation 48
sin_theta_kj = zeros(Nball, Nball);
sin_theta = zeros(Nball,1);
for K = 1:8
    for k = 1:8
        for j = 1:8
            sin_theta_kj(k, j) = B(K)*(F(k)-F(j))/(E(k)-E(j)*(C(K)-F(j)) - (F(k)-F(j))*(A(K)-E(j)));
        end
        sin_theta_kj(k,k) = 0;
    end
    
    sin_theta(K) = sum(sum(sin_theta_kj))./(Nball.^2-Nball);
end    
sin_theta = sum(sum(sin_theta_kj,2),1)/(Nball-1)^2;
theta = shiftdim(real(asin(sin_theta)),2);

%% phi
arctan_phi = zeros(length(z), length(z));
 for j = 1:length(z)
     arctan_phi(:,j) = (C - F(j))./(sin(theta).*(A - E(j))+B);
 end
phi = atan(mean(real(arctan_phi),2));

%% dx, dz, R
dx_all = zeros(length(z), length(z));
for j = 1:length(z)
    dx_all(:,j) = sin(phi).*sin(theta)*E(j) - cos(phi).*F(j);
end
dx = mean(dx_all,2);
dz = -cos(theta).*A;
R = sqrt(A.^2 + B.^2 + C.^2 + 2*sin(theta).*A.*B);
Ry_p = sqrt(R.^2 - dx.^2 - dz.^2);
 
%% Rf scales
eta_mean = mean(eta);
phi_mean = mean(phi);
theta_mean = mean(theta);
O = RM(theta_mean, eta_mean, phi_mean);
 
u_plus = u(:,1:Nangles/2) + u(:,Nangles/2+1:end);
u_minus = u(:,1:Nangles/2) - u(:,Nangles/2+1:end);
 
v_plus = v(:,1:Nangles/2) + v(:,Nangles/2+1:end);
v_minus = v(:,1:Nangles/2) - v(:,Nangles/2+1:end);

a =    O(3,1)*u_plus + O(3,3)*v_plus + 2*dz*ones(1,Nangles/2);
d = -( O(1,1)*u_plus + O(1,3)*v_plus + 2*dx*ones(1,Nangles/2) );
e = -( O(2,1)*u_plus + O(2,3)*v_plus + 2*Ry_p*ones(1,Nangles/2) );

b = (  O(3,1)*u_minus + O(3,3)*v_minus).*(ones(length(z),1)*cos(rotation_angle_alpha(1:Nangles/2)));

c = (  O(3,1)*u_minus + O(3,3)*v_minus).*(ones(length(z),1)*sin(rotation_angle_alpha(1:Nangles/2)));
 
f = -(O(1,1)*u_minus + O(1,3)*v_minus).*(ones(length(z),1)*cos(rotation_angle_alpha(1:Nangles/2))) ...
    +(O(2,1)*u_minus + O(2,3)*v_minus).*(ones(length(z),1)*sin(rotation_angle_alpha(1:Nangles/2)));
 
g = -(O(2,1)*u_minus + O(2,3)*v_minus).*(ones(length(z),1)*cos(rotation_angle_alpha(1:Nangles/2))) ...
    -(O(1,1)*u_minus + O(1,3)*v_minus).*(ones(length(z),1)*sin(rotation_angle_alpha(1:Nangles/2)));
 
Rf_p = a.*e + c.*f - b.*g;
x0_Rf = mean((c.*e + a.*f - b.*d)./Rf_p,2);
y0_Rf = mean(-(b.*e + c.*d - a.*g)./Rf_p,2);
z0_Rf = mean((-a.^2 + b.^2 + c.^2)./Rf_p,2);

%% plot things back
R_mean = mean(R);
Ry_p_mean = mean(Ry_p);
dz_mean = mean(dz);
dx_mean = mean(dx);

O = RM(theta_mean, eta_mean, phi_mean);
x_alpha_Rf = x0_Rf*cos(rotation_angle_alpha) + y0_Rf*sin(rotation_angle_alpha);
y_alpha_Rf = -x0_Rf*sin(rotation_angle_alpha) + y0_Rf*cos(rotation_angle_alpha);

u_id_Rf = (R_mean*x_alpha_Rf)./(y_alpha_Rf + 1);
v_id_Rf = (R_mean*z0_Rf*ones(1, length(rotation_angle_alpha)))./(y_alpha_Rf + 1);


u_idtest = u_id_Rf;
v_idtest = v_id_Rf; 

%% faster
det_Q = (O(1,1)-O(2,1)*u_idtest/R_mean).*(O(3,3)-O(2,3)*v_idtest/R_mean) - ...
        (O(1,3)-O(2,3)*u_idtest/R_mean).*(O(3,1)-O(2,1)*v_idtest/R_mean);

u_idptest = Ry_p_mean/R_mean *u_idtest;
v_idptest = Ry_p_mean/R_mean *v_idtest;

Ma = O(3,3)-O(2,3)*v_idtest/R_mean;
Mb = -(O(1,3)-O(2,3)*u_idtest/R_mean);
Mc = -(O(3,1)-O(2,1)*v_idtest/R_mean);
Md = O(1,1)-O(2,1)*u_idtest/R_mean;
Me = (u_idptest) - dx_mean;
Mf = (v_idptest) - dz_mean;

u_test = (1./det_Q).*(Ma.*Me + Mb.*Mf);
v_test = (1./det_Q).*(Mc.*Me + Md.*Mf);

figure; 
plot(u', v', 'b');
hold on;
plot(u_test', v_test', 'r');
hold off;

fit_2 = CalibrationMSE(u', v', u_test', v_test');
%%
cd(currentpath);