clear all
close all
clc

load('temp_34.mat');
figure(1); plot(xcoor_filtered_gnoise, ycoor_filtered_gnoise, '.');
x = xcoor_filtered_gnoise;
y = ycoor_filtered_gnoise;
%%
% index_array = 1:length(x);
% radius = 30; %if the points are within this radius, centered at COM, then sure
% x0 = x(1);
% y0 = y(1);
% 
% x_sum = 0; %initial, x0
% y_sum = 0; %initial, y0
% 
% x_com = x0; %initialize x_com
% y_com = y0; %initialize y_com
% 
% n = 1;
% n_points = 1;
% key = zeros(size(x));
% 
% clusterN = 1;
% index_limit = 0;
% % figure;
% 
% while clusterN < 10
%     l = length(index_array); %preset so it does not loop over and over again
%     while n <= l
%         d = sqrt( (x(index_array(1)) - x_com)^2 + (y(index_array(1))-y_com)^2 );
%         if (d <= radius)
%             %this point can be used in COM calculation
%             x0 = x(index_array(1));
%             y0 = y(index_array(1));
% 
%             index_array(1) = [];
%                 
%             key(index_array(1)) = clusterN;
%             %calculate COM
%             x_sum = x_sum + x0;
%             y_sum = y_sum + y0;     
%             x_com = x_sum / n_points; %center of mass in x
%             y_com = y_sum / n_points; %center of mass in y
%                 
%             n_points = n_points + 1;
%             
% %                 %display to figure  
% %                 figure(1);
% %                 subplot(1,2,1); plot(x0, y0, 'x'); hold on;
% %                 h = get(gca);
% %                 subplot(1,2,2); plot(x_com, y_com, 'r.'); hold on; 
% %                 set(gca, 'xlim', h.XLim, 'ylim', h.YLim);
% 
%         end
%         n = n+1;  
%     end
%     
%     x_sum = 0;
%     y_sum = 0;
%     x_com = x(index_array(1));
%     y_com = y(index_array(1));
%     n_points = 1;
%     clusterN = clusterN + 1;    
%     n = 1;
%     figure;
%     plot(x(index_array), y(index_array), '.');
% end

%%
index_array = 1:length(x);
radius = 25; %if the points are within this radius, centered at COM, then sure
x0 = x(1);
y0 = y(1);

x_sum = 0; %initial, x0
y_sum = 0; %initial, y0

x_com = x0; %initialize x_com
y_com = y0; %initialize y_com

n = 1;
n_points = 1;
key = zeros(size(x));

clusterN = 1;
index_limit = 0;
% figure;
dd = 40;
cluster_index_array = 1:10;
cluster_temp = 10;
emptyTF = 0;

% while emptyTF ~= 1 && cluster_index_array(1) <= 8
    
    l = length(index_array); %preset so it does not loop over and over again
    
    while n <= l && emptyTF ~= 1
        d = sqrt( (x(index_array(1)) - x_com)^2 + (y(index_array(1))-y_com)^2 );

        if (d <= radius)
            %this point can be used in COM calculation
            x0 = x(index_array(1));
            y0 = y(index_array(1));
            
%             fprintf('x0 = %.f, y0 = %.0f \n', x0, y0);               
            key(index_array(1)) = cluster_temp;
            %calculate COM
            x_sum = x_sum + x0;
            y_sum = y_sum + y0;     
            x_com = x_sum / n_points; %center of mass in x
            y_com = y_sum / n_points; %center of mass in y
                
            n_points = n_points + 1;
            
                %display to figure  
%                 figure(2);
%                 subplot(1,2,1); plot(x0, y0, 'x'); hold on;
%                 h = get(gca);
%                 subplot(1,2,2); plot(x_com, y_com, 'r.'); hold on; 
%                 set(gca, 'xlim', h.XLim, 'ylim', h.YLim);

            index_array(1) = [];
        end
        n = n+1;  
        index_array(1) = [];
        emptyTF = isempty(index_array);
    end
    
    %this part is to double check everything in this cluster is valid
    index_cluster = find(key == cluster_temp);
    x_cluster = x(index_cluster);
    y_cluster = y(index_cluster); 
    rx = range(x_cluster);
    ry = range(y_cluster);
    x_std = std(x_cluster);
    y_std = std(y_cluster);
    if abs(rx - ry) < 10 && abs( x_std - y_std ) < 5 && length(index_cluster) > 500
       key(index_cluster) = cluster_index_array(1); 
       cluster_index_array(1) = [];
    end
    
    figure; plot(x_cluster, y_cluster, 'x');
    hold on;
    plot(x_com, y_com, 'ro');    
    
    x_sum = 0;
    y_sum = 0;
    x_com = x(index_array(1));
    y_com = y(index_array(1));
    n_points = 1;

    n = 1;
%     figure;
%     plot(x(index_array), y(index_array), '.');
    disp(emptyTF);
% % end
% figure;
% plot(x(index_array), y(index_array), '.');