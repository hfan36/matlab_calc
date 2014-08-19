% folder_root = 'F:\Calibration\112613\';
% load(strcat(folder_root, 'smekal_results.mat'));
% theta_obj = 0;
% eta_obj = 0; 
% phi_obj = 0;
% alpha = (0:2:358)*pi/180;
% Nball = 8;
% R_phantom = 60;
% xi = 30*ones(Nball,1);
% yi = sqrt(R_phantom^2-xi.^2);
% Rf_guess0 = round(15/abs(result.z0_Rf(5)-result.z0_Rf(4)));
% z = 0:15:15*(Nball-1);
% zi = (z+(result.z0_Rf(4)*Rf_guess0 - z(4)))';
% ri_Rf = [xi yi zi]'/Rf_guess0;
% smekal = result;
% close all;

function[u, v] = CalibrationPlot_smekal(smekal, ri_Rf, theta_obj, eta_obj, phi_obj, alpha)

r0_Rf = RM(theta_obj, eta_obj, phi_obj)*ri_Rf;
O = RM(smekal.theta_mean, smekal.eta_mean, smekal.phi_mean);
x_alpha_Rf = r0_Rf(1,:)'*cos(alpha) + r0_Rf(2,:)'*sin(alpha);
y_alpha_Rf = -r0_Rf(1,:)'*sin(alpha) + r0_Rf(2,:)'*cos(alpha);

u_id = (smekal.R_mean*x_alpha_Rf)./(y_alpha_Rf + 1);
v_id = (smekal.R_mean*r0_Rf(3,:)'*ones(1,length(alpha)))./(y_alpha_Rf + 1);

%% this the the original form in loops, I've speed it up a bit by multiplying the matrix
%terms separately


% alphaL = length(alpha);
% u = zeros(Nball, alphaL);
% v = zeros(Nball, alphaL);
% for n = 1:Nball
%     for m = 1:alphaL
%         Q_test = [O(1,1)-O(2,1)*u_id(n,m)/smekal.R_mean O(1,3)-O(2,3)*u_id(n,m)/smekal.R_mean; ...
%                   O(3,1)-O(2,1)*v_id(n,m)/smekal.R_mean O(3,3)-O(2,3)*v_id(n,m)/smekal.R_mean];
% 
%         u_id_p = smekal.Ry_p_mean/smekal.R_mean *u_id(n,m);
%         v_id_p = smekal.Ry_p_mean/smekal.R_mean *v_id(n,m);
% 
%         det_real_test = 1/det(Q_test) * [O(3,3)-O(2,3)*v_id(n,m)/smekal.R_mean -(O(1,3)-O(2,3)*u_id(n,m)/smekal.R_mean);...
%                                        -(O(3,1)-O(2,1)*v_id(n,m)/smekal.R_mean)  O(1,1)-O(2,1)*u_id(n,m)/smekal.R_mean]* ...
%                                         [u_id_p - smekal.dx_mean; v_id_p - smekal.dz_mean];
%         u(n,m) = det_real_test(1);
%         v(n,m) = det_real_test(2);      
%     end
% end
% toc

det_Q = (O(1,1)-O(2,1)*u_id/smekal.R_mean).*(O(3,3)-O(2,3)*v_id/smekal.R_mean) - ...
        (O(1,3)-O(2,3)*u_id/smekal.R_mean).*(O(3,1)-O(2,1)*v_id/smekal.R_mean);
u_id_p = smekal.Ry_p_mean/smekal.R_mean *u_id;
v_id_p = smekal.Ry_p_mean/smekal.R_mean *v_id;

% form (1/det_Q).*[a b; c d]*[e; f];
a = O(3,3)-O(2,3)*v_id/smekal.R_mean;
b = -(O(1,3)-O(2,3)*u_id/smekal.R_mean);
c = -(O(3,1)-O(2,1)*v_id/smekal.R_mean);
d = O(1,1)-O(2,1)*u_id/smekal.R_mean;
e = (u_id_p) - smekal.dx_mean;
f = (v_id_p) - smekal.dz_mean;

u = (1./det_Q).*(a.*e + b.*f);
v = (1./det_Q).*(c.*e + d.*f);

% figure;
% plot(u', v', '.b');

 


