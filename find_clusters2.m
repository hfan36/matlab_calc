clear all
close all
clc


% for t = 1:1%180
for t = 1;
    close all
    filename = strcat('temp_', num2str(t-1), '.mat');
    load(filename);
    
    params.threshold_range(1) = 100;
    params.threshold_range(2) = 200;
    
    clip_dat_medfilt = medfilt2(clip_dat);
    
	figure(8); imagesc(clip_dat_medfilt', [0 500]); colormap gray; axis square;
    
    [ycoord_noedge, xcoord_noedge] = find(clip_dat_medfilt <= params.threshold_range(2) & clip_dat_medfilt >= params.threshold_range(1));
    test = find(clip_dat_medfilt(:) <= params.threshold_range(2) & clip_dat_medfilt(:) >= params.threshold_range(1));
    figure(1); plot(xcoor_filtered_gnoise, ycoor_filtered_gnoise, '.');
    figure(2); plot(xcoord_noedge, ycoord_noedge, '.');
    
%     x = xcoor_filtered_gnoise;
%     y = ycoor_filtered_gnoise;
    x = xcoord_noedge;
    y = ycoord_noedge;
    %%
    index_array = 1:length(x);
    radius = 25;

    x_sum = 0;
    y_sum = 0;
    x_com = x(1);
    y_com = y(1);

    loop_iteration = 1;
    n_points = 1;

    key = zeros(size(x));

    n = 1;
    clusterN = 1;
    l = length(index_array);

    while isempty(index_array) ~= 1 % <= length(x)
    % while l >= 1300;
    
        %finding the distance between 1 and all other points
        d = sqrt( (x(index_array(1)) - x(index_array)).^2 + (y(index_array(1)) - y(index_array)).^2 );
        
        %find the index of distances, d that is less than 2*radius and it's not the point
        %itself
        valid_index = find(d <= 2*radius & d ~= 0);

        %if the validpoint clusters are less than 100, and valid_index is
        %not empty -> a point hanging out by itself
        %->delete those index_arrays, because the cluster is too small
        if length(valid_index) < 100 && ~isempty(valid_index) 
            index_array(valid_index) = [];
        %if the point is hanging out by itself, then delete it, because
        %it's not one of the interested clusters
        elseif isempty(valid_index) %is empty
            index_array(1) = [];

        %if the clusters are greater than 100 points, then we'll see
        elseif length(valid_index) >= 100 && ~isempty(valid_index)

           xsum = sum(x(index_array(valid_index)));
           ysum = sum(y(index_array(valid_index)));
           %finding the center-of-mass points in this cluster
           xcom = xsum/length(valid_index);
           ycom = ysum/length(valid_index);
           fprintf('xcom = %.f, ycom = %.f \n', xcom, ycom);

           %given the COM point (xcom, ycom), lets look at all points
           %again, and find those points that are within radius of this
           %cluster
           d2 = sqrt( (x(index_array) - xcom).^2 + (y(index_array) - ycom).^2 );
           index2 = find(d2 <= radius);

           %look at the range, and std of the x, and y values
           %if the cluster is roundish (ystd approximately equals to xstd,
           %xrange approximately equal to yrange)
           %then it's probably one of the clusters we are looking for.
           xr = range(x(index_array(index2)));
           yr = range(y(index_array(index2)));
           xstd = std(x(index_array(index2)));
           ystd = std(y(index_array(index2)));

           if ( abs(xr-yr) <= 5) && (abs(xstd - ystd) <= 5 && length(index2) >= 500 ) %checking to make sure this cluster is somewhat roundish
                key(index_array(index2)) = clusterN;
                figure(5); subplot(3,3,clusterN); plot(x(index_array(index2)), y(index_array(index2)), '.'); axis equal;
                xlim([xcom-30 xcom+30]); ylim([ycom-30 ycom+30]);
                index_array(index2) = [];
                clusterN = clusterN + 1;

           elseif ( abs(xr-yr) > 10 ) || (abs(xstd - ystd) > 5) %if it's not really round, then just delete them from the array...
               index_array(index2) = [];
           else
               index_array(index2) = [];
           end
        end

        l = length(index_array);

    %     disp(l);
%         n = n+1;
%         figure; plot(x(index_array), y(index_array), '.');
    %     pause();

    end
% disp(clusterN);
fprintf('t = %d, clusters = %d \n', t, clusterN-1);
% pause();
end

%%
% for q = 1:8
%     n1 = find(key == q);
%     xn = x(n1);
%     yn = y(n1);
%     dat = zeros(size(xn));
% 
%     for a = 1:length(xn)
%        dat(a) = clip_dat(yn(a), xn(a));  
%     end
%     fprintf('mean = %.3f, std = %.3f \n', mean(dat), std(dat));
% end
%%

%%
% t = 1;
% filename = strcat('temp_', num2str(t-1), '.mat');
% load(filename);
% clip_dat_medfilt = medfilt2(clip_dat);  
% % figure(8); imagesc(clip_dat_medfilt', [0 500]); colormap gray; axis square;
% [ycoord, xcoord] = find(clip_dat_medfilt <= params.threshold_range(2) & clip_dat_medfilt >= params.threshold_range(1));
% index = ycoord + (xcoord-1)*2061;
% 
% img = zeros(size(clip_dat));
% img(index) = clip_dat(index);
% figure; imagesc(img'); colormap gray; axis square;
% data = [ycoord, xcoord];
% 
% % randa = datasample(data, 1);
% randa = [366 1631];
% disp(randa);
% disp(img(randa(1), randa(2)));
% 
% 
% d = sqrt( (ycoord - randa(1)).^2 + (xcoord - randa(2)).^2 );
% plot(d);
% 
% x_com = randa(2);
% y_com = randa(1);
% 
% index_array = 1:length(ycoord);
% d = sqrt( (ycoord(index_array) -y_com).^2 + (xcoord(index_array) - x_com).^2 );
% B_d = d;
% s = 2;
% x_sum = 0;
% y_sum = 0;
% len = 0;
% IX = 1:length(ycoord);
% key = zeros(size(ycoord));
% 
% while s <= 25
%     
%     B_d = sqrt( (ycoord(IX) -y_com).^2 + (xcoord(IX) - x_com).^2 );
%     [B_d, IX] = sort(B_d, 'ascend');
%     small_index = IX(B_d <= s);
%     y_sum = y_sum + sum(ycoord(small_index));
%     x_sum = x_sum + sum(xcoord(small_index));
%     
%     key(small_index) = 1;
%     len = len + length(small_index);
%     
%     y_com = y_sum/len;
%     x_com = x_sum/len;
%     s = s + 1;
%     
%     IX(B_d <= s) = [];
%     B_d(B_d <= s) = [];
% end
% 
% figure; plot(key);
% 
% % b = [1 4 2 5 6 3 9 7];
% % ix = 1:length(b);
% % 
% % [b, ix] = sort(b, 'ascend');
% % ix(b <= 3) = [];
% % b(b <= 3) = [];




