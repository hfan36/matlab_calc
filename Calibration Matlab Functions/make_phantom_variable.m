function [out] = make_phantom_variable(theta, eta, phi, R_phantom, dz, delta_z_ball, Nballs, z0_Rf)

out.theta0 = theta;
out.eta0 = eta;
out.phi0 = phi;
out.dz0 = dz;
out.xi0 = R_phantom*ones(Nballs,1);
out.yi0 = sqrt(R_phantom.^2 - ones(Nballs,1));
z = 0:delta_z_ball:(delta_z_ball*Nballs-1);
out.Rf0 = round(delta_z_ball./abs(z0_Rf(round(Nballs/2)) - z0_Rf(round(Nballs/2)+1)));
out.zi0 = z + z0_Rf(round(Nballs/2))*out.Rf0 - z(round(Nballs/2));

