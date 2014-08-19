function [out] = make_smekal_variables(eta, theta, phi, dx, dz, R, Ry_p, Rf_p, x0_Rf, y0_Rf, z0_Rf)

out.eta = eta;
out.theta = theta;
out.phi = phi;
out.dx = dx;
out.dz = dz;
out.R = R;
out.Ry_p = Ry_p;
out.Rf_p = Rf_p;
out.x0_Rf = x0_Rf;
out.y0_Rf = y0_Rf;
out.z0_Rf = z0_Rf;
out.eta_mean = mean(eta);
out.theta_mean = mean(theta);
out.phi_mean = mean(phi);
out.dx_mean = mean(dx);
out.dz_mean = mean(dz);
out.R_mean = mean(R);
out.Ry_p_mean = mean(Ry_p);
out.Rf_p_mean = mean(Rf_p);


