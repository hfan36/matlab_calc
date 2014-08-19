clear all
close all
clc
%% set up pathes
currentpath = cd;
newpath = 'H:\CT data\051314';
% change to new path
cd(newpath);
% phantom_60kV_400uA_10sec_0.ct
%% clean up image and save it
for image_index = 1:2
    filename = strcat('phantom_60kV_400uA_10sec_', num2str(image_index-1), '.ct');
    dat = readBinary(filename, 2560*2160+1, 'uint16');
    dat = dat(2:end);
    dat = reshape(dat, 2560, 2160);
    
    figure(1);
    imagesc(dat, [0 600]); colormap gray; 
    axis equal; xlim([1 2160]);
    colorbar;
    
    dat = padarray(dat, [10 10], 'both');


    [y_spot, x_spot] = find( (dat >= 600) );

    dr = zeros(3, length(x_spot));
    for n = 1:3
         dr(n,:) = sqrt( (x_spot-circshift(x_spot, -2 + n)).^2 + (y_spot-circshift(y_spot, -2 + n)).^2 );
    end

    index = zeros(length(x_spot),1);
    count = 1;

    threshold = 5;
    for m = 1:size(dr,2)

        if (dr(1,m) > threshold) || (dr(3,m) > threshold)
            index(count) = m;
            count = count + 1;
        end

    end

    index = index(index~=0);

    tic
    pts2 = zeros(size(dr, 2),1);
    count = 1;

    threshold2 = 5;
    for m = 1:size(dr,2)

        if ((dr(1,m) < threshold2) && (dr(3,m) > threshold2))
            pts2(count) = m;
            count = count + 1;
        elseif (dr(1,m) > threshold2) && (dr(3,m) > threshold2)
            pts2(count) = m;
            count = count + 1;
        end  
    end 
    pts2 = pts2(pts2~=0);
    toc

%     figure(2);
%     plot(x_spot, y_spot, '.-');
%     hold on;
%     plot(x_spot(index), y_spot(index), '+r');
%     plot(x_spot(pts2), y_spot(pts2), 'sk');
%     hold off;

    row = zeros(length(pts2),1);
    col = zeros(length(pts2),1);
    dat_2 = dat;
    A = zeros(11,11,length(pts2));

    for j = 1:length(pts2)

        A(:,:,j) = dat(y_spot(pts2(j))-5:y_spot(pts2(j))+5, x_spot(pts2(j))-5:x_spot(pts2(j))+5);
        [row_a, col_a] = find(A(:,:,j) == max(max(A(:,:,j))));
    % 
        row = y_spot(pts2(j))-5+row_a(1) - 1;
        col = x_spot(pts2(j))-5+col_a(1) - 1;

        ave = sum(sum(dat_2(row-7:row+7, col-7:col+7) - padarray(dat_2(row-3:row+3, col-3:col+3), [4 4]),1),2)/176;
        dat_2(row-3:row+3, col-3:col+3) = ave + 2*(rand(7,7)-0.5)*ave*0.08;

    end
    
    dat_2 = dat_2(11:end-10, 11:end-10);
    clean_data = medfilt2(dat_2, [3 3]);
    
    %replace surrounding pixels
     [rowz, colz] = find(clean_data == 0);
     

      temp = (circshift(clean_data, [1, 0])  + ...
              circshift(clean_data, [-1, 0]) + ...
              circshift(clean_data, [0, 1])  + ...
              circshift(clean_data, [0, -1]) + ...
              circshift(clean_data, [1, 1])  + ...
              circshift(clean_data, [-1, 1]) + ...
              circshift(clean_data, [-1,-1]) + ...
              circshift(clean_data, [1, 1])  + ...
              circshift(clean_data, [1, -1]))./8;       
     
    clean_data(rowz, colz) = temp(rowz, colz);
          
%     figure(3);
%     imagesc(dat_2, [0 600]); colormap gray; 
%     axis equal; xlim([1 2160]);
%     colorbar;
    
    figure(4);
    imagesc(clean_data, [0 600]); colormap gray; 
    axis equal; xlim([1 2160]);
    colorbar;
    
%     mkdir(newpath, 'cleaned_data');
    newdata_filename = strcat(newpath, '\cleaned_data\', 'clean_', filename);
%     writeBinary(clean_data, newdata_filename, 'float');
    disp(image_index);
end
%%
% newpath = strcat(newpath, '\cleaned_data\');
% cd(newpath);
% for n = 1:180
%     filename = strcat('clean_phantom_60kV_400uA_10sec_', num2str(n-1), '.ct');
% 
%     dat = readBinary(filename, 2560*2160, 'float');
%     dat = reshape(dat, 2560, 2160);
% 
% %     figure; imagesc(dat, [0 500]); colormap gray; axis equal; colorbar;
% 
%     dat2 = dat(996-511:996+512, 1096-511:1096+512);
%     dat2 = 10-log(dat2');
%     figure(2); imagesc(dat2); colormap gray; axis equal; colorbar;
%     filename2 = strcat('C1024_', filename);
%     writeBinary(dat2(:), filename2, 'float');
% end
%%
% cd('H:\CT data\051314\simulation\');
% dat = readBinary('bp2.bin', 64^3, 'float');
% dat = reshape(dat, 64,64,64);
% 
% figure;
% imagesc((shiftdim(dat(16,:,:),1))); axis equal; colormap gray; colorbar; 
% figure;
% imagesc((shiftdim(dat(32,:,:),1))); axis equal; colormap gray; colorbar; 
% figure; 
% imagesc((shiftdim(dat(40,:,:),1))); axis equal; colormap gray; colorbar; 

%% change back to the matlab file path;
cd(currentpath);