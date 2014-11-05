clear all
close all
clc

for t = 1:180; %only selecting the first image, '0'
%     disp(t);
    close all;
    filename = strcat('temp_', num2str(t-1), '.mat');
    load(filename);
    
    params.threshold_range(1) = 100;
    params.threshold_range(2) = 200;
    
    clip_dat_medfilt = medfilt2(clip_dat);
    
% 	figure(8); imagesc(clip_dat_medfilt', [0 500]); colormap gray; axis square;
    
    
    [ycoord_noedge, xcoord_noedge] = find(clip_dat_medfilt <= params.threshold_range(2) & clip_dat_medfilt >= params.threshold_range(1));
%     test = find(clip_dat_medfilt(:) <= params.threshold_range(2) & clip_dat_medfilt(:) >= params.threshold_range(1));
%     figure(1); plot(xcoor_filtered_gnoise, ycoor_filtered_gnoise, '.');
    figure(2); plot(xcoord_noedge, ycoord_noedge, '.'); xlim([0 1800]); ylim([0 2000]);  
    
%     x = xcoor_filtered_gnoise;
%     y = ycoor_filtered_gnoise;
    x = xcoord_noedge;
    y = ycoord_noedge;

 %% try it again similar to what was described in the dissertation
 key = zeros(size(x));
 radius = 35;
 counter = 0;
 
 Ncluster = [];
 
 index_array = 1:length(x);
 keyvalue = 1;

%%
g = cell(30,7);
%the final clusters are stored in a cell format, I've allocated 30 maximum
%clusters, but it can change, though data has never exceeded 30 so far.
%They are organized as follows:
%{xclusterpoints, yclusterpoints, xcom, ycom, xrange, yrange, Ncluster}
%you can plot xclusterpoints and yclusterpoints straight, these are the
%coordinates from original image so if these points were used on the
%original image, then you can retrieve the image values at these
%coordinates.
%xcom and ycom are the center-of-mass points calculated using the cluster
%points
%xrange and yrange are basically range calculated using xclusterpoints,
%yclusterpoints
%Ncluster = number of points in this cluster, it's also the length of
%xclusterpoints and yclusterpoints
while size(index_array)~= 0
    
     %setting up stuff
     xcom = x(index_array(1));
     ycom = y(index_array(1));
     COM_xsum = 0;
     COM_ysum = 0;
     counter = 0;
     index = [];
     
     for n = 1:length(index_array)
        %calculate the distance of the next point away from xcom and ycom
        distance = sqrt( ( xcom-x(index_array(n)) ).^2 + ( ycom-y(index_array(n)) ).^2 );
        
        
        if (distance <= radius)
            %key basically iterates through all the identified clusters
            key(index_array(n)) = keyvalue;
            index = [index n];
            counter = counter + 1;
            COM_xsum = COM_xsum + x(index_array(n));
            COM_ysum = COM_ysum + y(index_array(n));
            xcom = COM_xsum/counter;
            ycom = COM_ysum/counter;
            
%             figure(5);
%             plot(xcom, ycom, 'or'); hold on;
%             plot(x(index_array(index)), y(index_array(index)), '.'); hold off;
%             xlim([xcom-30 xcom+30]); ylim([ycom-30 ycom+30]);
        end
     end
     
     figure(5);
     plot(x(index_array(index)), y(index_array(index)), '.', 'Color', [rand rand rand]); hold on;
     plot(xcom, ycom, 'o', 'MarkerSize', 20);
     axis equal;
     
     %calculate x and y range values
     xrange = range(x(index_array(index)));
     yrange = range(y(index_array(index)));
     xcell = x(index_array(index));
     ycell = y(index_array(index));
    
     %enter values into cell array
     g(keyvalue,:) = {xcell, ycell, xcom, ycom, xrange, yrange, counter};  
     
     %record and destroy
     Ncluster = [Ncluster counter];
     index_array(index) = [];
     keyvalue = keyvalue + 1;
end
    
    %delete those cells that are empty then reshape back into the original
    %cell array format, only shorter
    g(cellfun(@isempty, g))=[];
    l = length(g)/7;
    g = reshape(g, l, 7);

    LL = size(g,1);   
    %finding the clusters that are roundish
    for mm = 1:LL
        if ( (cell2mat(g(mm,5)) >= 32) && (cell2mat(g(mm,5)) <= 60) && ...
             (cell2mat(g(mm,6)) >= 32) && (cell2mat(g(mm,6)) <= 60) && ...
           (abs(cell2mat(g(mm,5))-cell2mat(g(mm,6))) <= 12) && cell2mat(g(mm,7)) >= 800)
       
%             g(mm,:) = g(mm,:);
        else
            g(mm,:) = {[]}; %if they don't fit the criteria, then it's emptied out
        end  
        
    end
    
    %delete those that are empty
    g(cellfun(@isempty, g))=[];
    g = reshape(g, length(g)/7, 7);
    
    %plot subplots of clusters
    figure(1000);
    for j = 1:size(g,1)
       subplot(2,4,j); 
       plot( cell2mat(g(j,1)), cell2mat(g(j,2)),  '.'); 
       xlim([cell2mat(g(j,3))-30 cell2mat(g(j,3))+30]); ylim([cell2mat(g(j,4))-30 cell2mat(g(j,4))+30]);
       axis square;
       hold on;
       plot( cell2mat(g(j,3)), cell2mat(g(j,4)), 'or');
       hold off;
    end
        
%     pause();

end
