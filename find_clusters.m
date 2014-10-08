clear all
close all
clc

load('temp.mat');
figure(1); plot(xcoor_filtered_gnoise, ycoor_filtered_gnoise, '.');
x = xcoor_filtered_gnoise;
y = ycoor_filtered_gnoise;
index_array = 1:length(x);
%%

radius = 30; %if the points are within this radius, centered at COM, then sure
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

while clusterN < 10
    l = length(index_array); %preset so it does not loop over and over again
    while n <= l
        d = sqrt( (x(index_array(1)) - x_com)^2 + (y(index_array(1))-y_com)^2 );
        if (d <= radius)
            %this point can be used in COM calculation
            x0 = x(index_array(1));
            y0 = y(index_array(1));

            index_array(1) = [];
                
            key(index_array(1)) = clusterN;
            %calculate COM
            x_sum = x_sum + x0;
            y_sum = y_sum + y0;     
            x_com = x_sum / n_points; %center of mass in x
            y_com = y_sum / n_points; %center of mass in y
                
            n_points = n_points + 1;
            
%                 %display to figure  
%                 figure(1);
%                 subplot(1,2,1); plot(x0, y0, 'x'); hold on;
%                 h = get(gca);
%                 subplot(1,2,2); plot(x_com, y_com, 'r.'); hold on; 
%                 set(gca, 'xlim', h.XLim, 'ylim', h.YLim);

        end
        n = n+1;  
    end
    
    x_sum = 0;
    y_sum = 0;
    x_com = x(index_array(1));
    y_com = y(index_array(1));
    n_points = 1;
    clusterN = clusterN + 1;    
    n = 1;
    figure;
    plot(x(index_array), y(index_array), '.');
end
%%
% 
% while clusterN < 3
%     l = length(index_array); %preset so it does not loop over and over again
%     while n <= 3
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
% %     x_sum = 0;
% %     y_sum = 0;
% %     x_com = x(index_array(1));
% %     y_com = y(index_array(1));
% %     n_points = 1;
% %     clusterN = clusterN + 1;    
% %     
% %     figure;
% %     plot(x(index_array), y(index_array), '.');
% % end


%%
% x2 = x;
% y2 = y;
% figure; plot(x2, y2, '.');
% 
% x0 = x2(1);
% y0 = y2(1);
% x_sum = x0;
% y_sum = y0;
% x_com = x0;
% y_com = y0;
% n = 1;
% n_points = 1;
% clusterN = 2;
% 
% 
%     while n <= length(x)+1
%         d = sqrt( (x0 - x_com)^2 + (y0-y_com)^2 );
%         
%         if (n < (length(x) + 1))
%             if (d <= radius)
%                 %next point can be used in COM calculation
%                 x0 = x(n+1);
%                 y0 = y(n+1);
% 
%                 key(index_limit + n+1) = clusterN;
%                 n_points = n_points + 1;
% 
%                 %calculate COM
%                 x_sum = x_sum + x0;
%                 y_sum = y_sum + y0;     
%                 x_com = x_sum / n_points; %center of mass in x
%                 y_com = y_sum / n_points; %center of mass in y
% 
%                 %display to figure  
% %                 figure(2);
% %                 subplot(1,2,1); plot(x0, y0, 'x'); hold on;
% %                 h = get(gca);
% %                 subplot(1,2,2); plot(x_com, y_com, 'r.'); hold on; 
% %                 set(gca, 'xlim', h.XLim, 'ylim', h.YLim);
%             end 
%         elseif n == length(x)+1
%                 d = sqrt( (x(1) - x_com)^2 + (y(1)-y_com)^2 );
%                 if (d <= radius)
%                     key(index_limit + 1) = clusterN;
%                 end
%         end
%         n = n+1;  
% %         disp(n);
% %         pause(0.01);
%     end
% %     index = find(key == clusterN);
% %     index_limit = min(index)-1;
% %     
% %     %reset everything again...
% %     x = xcoor_filtered_gnoise(index);
% %     y = ycoor_filtered_gnoise(index);
%%
% x3 = x;
