clear all
close all
clc
%% set up pathes
currentpath = cd;
newpath = 'H:\CT data\042414';

% change to new path
cd(newpath);
%%



filename = 'dat_70kV_5000ms_400uA_50.ct';
dat = readBinary(filename, 2560*2160+1, 'uint16');
dat = dat(2:end);
dat = reshape(dat, 2560, 2160);
dat = padarray(dat, [10 10], 'both');

figure;
imagesc(dat, [0 600]); colormap gray; 
axis equal; xlim([1 2160]);
colorbar;
% set(gca, 'YDir', 'normal');

dat2 = medfilt2(dat, [3 3]);

% figure;
% imagesc(dat2, [0 600]); colormap gray;
% axis equal; colorbar;

%use the median filter
% hotspot_index = find((dat2(:) >= 500));
[y_spot, x_spot] = find( (dat >= 510) );

% figure; plot(x_spot, y_spot, '.');
% mask = (dat2(:) >= 500);
% mask = reshape(mask, size(dat2,1), size(dat2,2));
% 
% % figure; imagesc(mask); colormap gray; colorbar;
% 
% % x = floor(hotspot_index/size(dat2,1))+1;
% % y = rem(hotspot_index, size(dat2,1));
% 
% % Area = zeros(5,5,length(y));
% 
% % for n = 1:length(hotspot_index)
% %     Area(:,:,n) = dat(y(n)-2:y(n)+2, x(n)-2:x(n)+2); 
% %     replacement = dat(y(n)-4:y(n)+4, x(n)-4:x(n)+4); 
% %     area = padarray(Area(:,:,n), [2,2], 'both');
% %     m = sum(sum(replacement-area))/(size(replacement,1)^2-size(Area(:,:,1),1).^2);
% %     dat(y(n)-2:y(n)+2, x(n)-2:x(n)+2) = m;
% % end
% % 
% % figure
% % imagesc(dat, [0 600]); colormap gray; axis equal; colorbar;
% % 
% % fin = medfilt2(dat);
% % figure;
% % imagesc(fin, [0 600]); colormap gray; axis equal; colorbar;
% 
% %%
% a = [];
% rd = zeros(length(y_spot)-1,1);
% for n = 1:length(y_spot)-1
%     
%     rd(n) = sqrt((x_spot(1)-x_spot(n+1)).^2 + (y_spot(1) - y_spot(n+1)).^2);  
%     
%     if (n ~= 1)
%         if ( abs(rd(n)-rd(n-1)) > 5 )
%             a = [n a];
%         end
%     else
%         
%     end
% end
% 
% A = zeros(11,11, length(a));
% r_spot_single = zeros(length(a), 2);
% for i = 1:length(a)
%     r_spot_single(i,:) = [x_spot(a(i)); y_spot(a(i))];
%     
%     
% %     A(:,:,i) = dat(y_spot(a(i))-5:y_spot(a(i))+5, x_spot(a(i))-5:x_spot(a(i))+5);    
%     A(:,:,i) = dat(y_spot(a(i))-5:y_spot(a(i))+5, x_spot(a(i))-5:x_spot(a(i))+5);    
% end
% 
% figure; plot(x_spot, y_spot, '.-');
% hold on;
% plot(r_spot_single(:,1), r_spot_single(:,2), 'or');
% hold off;
% 
% %display
% % for f = 1: ceil(length(a)/9)   
% %         figure(f+1);
% %         for i = 1:9
% %             subplot(3,3,i);
% %             if ((f-1)*(9)+i) <=length(a)
% %                 imagesc(A(:,:,(f-1)*(9)+i)); colormap gray;
% %             end
% %         end
% % end
% 
% row = zeros(length(a),1);
% col = zeros(length(a),1);
% for j = 1:3
%     [row_a, col_a] = find (A(:,:,j) == max(max(A(:,:,j))));
% %     row = y_spot(a(j)-5+row_a - 1);
% %     col = x_spot(a(j)-5+col_a) - 1;
% 
%     row(j) = y_spot(a(j))-5+row_a - 1;
%     col(j) = x_spot(a(j))-5+col_a - 1;
% 
%     dat_2 = dat;
%     %four corners
%     dat_2(row(j)-1, col(j)-1) = round( sum(sum(dat(row(j)-1-2: row(j)-1, col(j)-1-2:col(j)-1),1),2)/9 );
%     dat_2(row(j)-1, col(j)+1) = round( sum(sum(dat(row(j)-1-2: row(j)-1, col(j)+1:col(j)+1+2),1),2)/9 );
%     dat_2(row(j)+1, col(j)-1) = round( sum(sum(dat(row(j)+1:row(j)+1+2, col(j)-1-2:col(j)-1), 1),2)/9 );
%     dat_2(row(j)+1, col(j)+1) = round( sum(sum(dat(row(j)+1:row(j)+1+2, col(j)+1:col(j)+1+2), 1),2)/9 );
%     %four crosses
%     dat_2(row(j)-1, col(j)) = round( sum(sum(dat(row(j)-1-2:row(j)-1, col(j)-1:col(j)+1),1),2)/9);%up
%     dat_2(row(j)+1, col(j)) = round( sum(sum(dat(row(j)+1:row(j)+1+2, col(j)-1:col(j)+1),1),2)/9);%down
%     dat_2(row(j), col(j)-1) = round( sum(sum(dat(row(j)-1:row(j)+1, col(j)-1-2:col(j)-1),1),2)/9);%left
%     dat_2(row(j), col(j)+1) = round( sum(sum(dat(row(j)-1:row(j)+1, col(j)+1:col(j)+1+2),1),2)/9);%right
%     %center
%     dat_2(row(j), col(j)) = round( (sum(sum(dat_2(row(j)-1:row(j)+1, col(j)-1:col(j)+1),1),2) - dat_2(row(j), col(j)))/8 ) ;
%     
% end
% 
% % figure;
% 
% % imagesc(dat(row-5:row+5, col-5:col+5));
% % 
% % figure;
% % imagesc(dat_2(row-5:row+5, col-5:col+5));
% 
% figure;
% imagesc((dat_2), [0 600]); colormap gray; axis equal; colorbar;
% xlim([1 2160]);
% % figure; plot(x_spot, y_spot, '.-');
% % hold on;
% % plot(r_spot_single(:,1), r_spot_single(:,2), 'or');
% % % hold off;
% % plot(col, row, 'dg');
% % hold off;

%%
% figure;
% plot(x_spot, y_spot, '.-');
% xlim([1 2160]);
% ylim([1 2560]);
% r_spot = sqrt(x_spot.^2 + y_spot.^2);

% figure;
% plot(r_spot, '.-');

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

% tic
% m = 1;
% pts = [];
% while m < size(dr,2)
%     if ((dr(1,m) < threshold) && (dr(3,m) > threshold))
%         pts = [m pts];
%         m = m+1;
%     elseif (dr(1,m) > threshold) && (dr(3,m) > threshold)
%         pts = [m pts];    
%     end
%     m = m+1;
% end
% toc

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

figure;
plot(x_spot, y_spot, '.-');
hold on;
plot(x_spot(index), y_spot(index), '+r');
% plot(x_spot(pts), y_spot(pts), 'dm');
plot(x_spot(pts2), y_spot(pts2), 'sk');
hold off;

% figure;
% plot(x_spot(index), y_spot(index), 'or');
% hold on;
% plot(x_spot(pts), y_spot(pts), '.b');
% hold off;



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

% method 1: filling in corners first, then crosses, then center of a 3x3
% area centered at the maximum point

% %     four corners
%     dat_2(row-1, col-1) = round( sum(sum(dat(row-1-2: row-1, col-1-2:col-1),1),2)/9 );
%     dat_2(row-1, col+1) = round( sum(sum(dat(row-1-2: row-1, col+1:col+1+2),1),2)/9 );
%     dat_2(row+1, col-1) = round( sum(sum(dat(row+1:row+1+2, col-1-2:col-1), 1),2)/9 );
%     dat_2(row+1, col+1) = round( sum(sum(dat(row+1:row+1+2, col+1:col+1+2), 1),2)/9 );
% %     four crosses
%     dat_2(row-1, col) = round( sum(sum(dat(row-1-2:row-1, col-1:col+1),1),2)/9);%up
%     dat_2(row+1, col) = round( sum(sum(dat(row+1:row+1+2, col-1:col+1),1),2)/9);%down
%     dat_2(row, col-1) = round( sum(sum(dat(row-1:row+1, col-1-2:col-1),1),2)/9);%left
%     dat_2(row, col+1) = round( sum(sum(dat(row-1:row+1, col+1:col+1+2),1),2)/9);%right
% %     center
%     dat_2(row, col) = round( (sum(sum(dat_2(row-1:row+1, col-1:col+1),1),2) - dat_2(row, col))/8 ) ;
    

% method 2: interpolate and just remove the center
    
    ave = sum(sum(dat_2(row-7:row+7, col-7:col+7) - padarray(dat_2(row-3:row+3, col-3:col+3), [4 4]),1),2)/176;
    dat_2(row-3:row+3, col-3:col+3) = ave + 2*(rand(7,7)-0.5)*ave*0.08;



end


figure;
imagesc(dat_2, [0 600]); colormap gray; 
axis equal; xlim([1 2160]);
colorbar;

figure;
imagesc(medfilt2(dat_2, [3 3]), [0 600]); colormap gray; 
axis equal; xlim([1 2160]);
colorbar;


%% change back to the matlab file path;
cd(currentpath);