function [u_test, v_test, out, u_FT_fit, v_FT_fit] = smekal_method(u, v, Nball, z, rotation_angle_alpha)
% this is the actual fitting part!
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
v_FT_fit =  zeros(length(z), length(rotation_angle_alpha));
for n = 1:length(z)
    u_FT_fit(n,:) = U0(n)/2 +  U1(n)*cos(rotation_angle_alpha) + U1_t(n)*sin(rotation_angle_alpha) + ...
                           U2(n)*cos(2*rotation_angle_alpha) + U2_t(n)*sin(2*rotation_angle_alpha) + ...
                           U3(n)*cos(3*rotation_angle_alpha) + U3_t(n)*sin(3*rotation_angle_alpha);
    v_FT_fit(n,:) = V0(n)/2 +  V1(n)*cos(rotation_angle_alpha) + V1_t(n)*sin(rotation_angle_alpha) + ...
                           V2(n)*cos(2*rotation_angle_alpha) + V2_t(n)*sin(2*rotation_angle_alpha) + ...
                           V3(n)*cos(3*rotation_angle_alpha) + V3_t(n)*sin(3*rotation_angle_alpha);
end

% figure;
% plot(u', v', '-x');
% hold on;
% plot(u_FT_fit', v_FT_fit', 'o');
% hold off;
% title('data fit using Fourier series method (paper)');
% fprintf('fit = %.5f \n', (sum(sum((u-u_FT_fit).^2 + (v-v_FT_fit).^2))));

% figure;
% plot(u(3,:), 'or'); hold on;
% plot(u_FT_fit(3,:)); hold off;
% 
% figure;
% plot(v(3,:), 'or'); hold on;
% plot(v_FT_fit(3,:)); hold off;

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
for k = 1:8
    for j = 1:8
        sin_theta_kj(k, j) = B(k)*(F(k)-F(j))/(E(k)-E(j)*(C(k)-F(j)) - (F(k)-F(j))*(A(k)-E(j)));
    end
    sin_theta_kj(k,k) = 0;
    
end

sin_theta = sum(sin_theta_kj,2)/(Nball-1);
theta = real(asin(sin_theta));

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

%% slower
% u_test = zeros(length(z), length(rotation_angle_alpha));
% v_test = zeros(length(z), length(rotation_angle_alpha));
% for n = 1:Nball
%     for m = 1:length(rotation_angle_alpha)
%         Q_test = [O(1,1)-O(2,1)*u_idtest(n,m)/R_mean O(1,3)-O(2,3)*u_idtest(n,m)/R_mean; ...
%                   O(3,1)-O(2,1)*v_idtest(n,m)/R_mean O(3,3)-O(2,3)*v_idtest(n,m)/R_mean];        
%         u_idptest = Ry_p_mean/R_mean*u_idtest(n,m);
%         v_idptest = Ry_p_mean/R_mean*v_idtest(n,m);
%         
%         det_real_test = 1/det(Q_test) * [O(3,3)-O(2,3)*v_idtest(n,m)/R_mean -(O(1,3)-O(2,3)*u_idtest(n,m)/R_mean);...
%                                        -(O(3,1)-O(2,1)*v_idtest(n,m)/R_mean)  O(1,1)-O(2,1)*u_idtest(n,m)/R_mean]* ...
%                                         [u_idptest - dx(n); v_idptest - dz_mean];
%                  
%         u_test(n,m) = det_real_test(1);
%         v_test(n,m) = det_real_test(2);      
%     end
% end

%% faster
det_Q = (O(1,1)-O(2,1)*u_idtest/R_mean).*(O(3,3)-O(2,3)*v_idtest/R_mean) - ...
        (O(1,3)-O(2,3)*u_idtest/R_mean).*(O(3,1)-O(2,1)*v_idtest/R_mean);

u_idptest = Ry_p_mean/R_mean *u_idtest;
v_idptest = Ry_p_mean/R_mean *v_idtest;

% form (1/det_Q).*[Ma Mb; Mc Md]*[Me; Mf];

Ma = O(3,3)-O(2,3)*v_idtest/R_mean;
Mb = -(O(1,3)-O(2,3)*u_idtest/R_mean);
Mc = -(O(3,1)-O(2,1)*v_idtest/R_mean);
Md = O(1,1)-O(2,1)*u_idtest/R_mean;
Me = (u_idptest) - dx_mean;
Mf = (v_idptest) - dz_mean;

u_test = (1./det_Q).*(Ma.*Me + Mb.*Mf);
v_test = (1./det_Q).*(Mc.*Me + Md.*Mf);


out = make_smekal_variables(eta, theta, phi, dx, dz, R, Ry_p, Rf_p, x0_Rf, y0_Rf, z0_Rf);











