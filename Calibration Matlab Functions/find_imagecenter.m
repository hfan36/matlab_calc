%function, it's actually a contracting grid algorithm to center the image
%(U, V values from pixels to mm size).  As long as I get close, dx and dz
%will take care of the rest. 

function[u0, v0, result, uc, vc] = find_imagecenter(u, v, delta_uc, delta_vc, uc0, vc0, alpha, cgvar, detvar)

Nball = size(u,1);
z = (0:15:15*(Nball-1)); % in millimeters, maybe, probably doesn't matter

mse = zeros(cgvar.N_loops,1);

u_center = getgridarray(uc0, delta_uc, cgvar.N_grid);
v_center = getgridarray(vc0, delta_vc, cgvar.N_grid);
mse_all = zeros(cgvar.N_grid, cgvar.N_grid, cgvar.N_loops);

for loopindex = 1:cgvar.N_loops
    
%     u0 = (u-u_center(1))/(detvar.DetectorWidth*detvar.PixelSize*detvar.Mag);
%     v0 = (v-v_center(1))/(detvar.DetectorHeight*detvar.PixelSize*detvar.Mag);
    u0 = (u-u_center(1))*(detvar.PixelSize*detvar.Mag);
    v0 = (v-v_center(1))*(detvar.PixelSize*detvar.Mag);

    [u_test, v_test, ~] = smekal_method(u0, v0, Nball, z, alpha);
    min_arg = CalibrationMSE(u0, v0, u_test, v_test);
%     fprintf('min_arg = %f \n', min_arg);
    
    for u_index = 1:cgvar.N_grid
        for v_index = 1:cgvar.N_grid
%             u_grid = (u-u_center(u_index))/(detvar.DetectorWidth*detvar.PixelSize*detvar.Mag);
%             v_grid = (v-v_center(v_index))/(detvar.DetectorHeight*detvar.PixelSize*detvar.Mag);

            u_grid = (u-u_center(u_index))*detvar.PixelSize*detvar.Mag;
            v_grid = (v-v_center(v_index))*detvar.PixelSize*detvar.Mag;

            [u_test, v_test, ~] = smekal_method(u_grid, v_grid, Nball, z, alpha);
            
            arg = CalibrationMSE(u_grid, v_grid, u_test, v_test);
%             fprintf('arg = %f \n', arg);
            mse_all(u_index, v_index, loopindex) = arg;
            
            if arg <= min_arg
                uc = u_center(u_index);
                vc = v_center(v_index);
                min_arg = arg;
                mse(loopindex) = min_arg; 
            else
                uc = u_center(1);
                vc = v_center(1);
            end
            
            
        end
    end
    
%     u0 = (u-uc)/(detvar.DetectorWidth*detvar.PixelSize*detvar.Mag);
%     v0 = (v-vc)/(detvar.DetectorWidth*detvar.PixelSize*detvar.Mag);

     u0 = (u-uc)*(detvar.PixelSize*detvar.Mag);
     v0 = (v-vc)*(detvar.PixelSize*detvar.Mag);
    
    [u_fit, v_fit, result, ~] = smekal_method( u0, v0, Nball, z, alpha);
    
    u_center = getgridarray(uc, delta_uc/(cgvar.Rate^loopindex), cgvar.N_grid);
    v_center = getgridarray(vc, delta_vc/(cgvar.Rate^loopindex), cgvar.N_grid);   
    
    
    figure(1000);
    subplot(1,3,1);
    plot(u0', v0', 'b'); hold on;
    plot(u_fit', v_fit', 'r'); hold off;
    
    subplot(1,3,2);
    plot(1:loopindex, mse(1:loopindex), '.b');
    
    subplot(1,3,3);
    imagesc(mse_all(:,:,loopindex)); colormap gray;
    
    fprintf('uc = %f, vc = %f, arg = %.10f mm \n', uc, vc, min_arg);
end
        

figure(2000);
plot(u0', v0', 'b');   hold on;
plot(u_fit', v_fit', 'r');
