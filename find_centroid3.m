clear all
close all
clc

NpixelsWidth = 2560; 
NpixelsHeight = 2160;
clip_height_top = 350;
clip_height_bottom = 100;
clipImageHeight = NpixelsHeight-clip_height_top-clip_height_bottom+1;
clip_lateral_edges = 250;
threshold_range = [100 200];
boundary_width = 20;
N_steelballs = 8;

L_noise = 2;
r_noise = 1;
A_noise = 1; 

r_bearing = 20;
L_bearing = 18;


% % plot filters
noise_filter = gaussian2D(L_noise, r_noise,1);
bearing_filter = gaussian2D(L_bearing, r_bearing, 1);

% figure(90);
% subplot(1,2,1); imagesc(noise_filter); colormap gray; axis equal;
% title('noise filter');
% subplot(1,2,2); imagesc(bearing_filter); colormap gray; axis equal;
% title('bearing filter');
% load image
params = struct('NpixelsWidth', NpixelsWidth, ...
                'NpixelsHeight', NpixelsHeight, ...
                'clip_height_top', clip_height_top, ...
                'clip_height_bottom', clip_height_bottom, ...
                'clip_lateral_edges', clip_lateral_edges, ...
                'boundary_width', boundary_width, ...
                'threshold_range', threshold_range, ...
                'N_steelballs', N_steelballs, ...
                'L_noise', L_noise, ...
                'r_noise', r_noise, ...
                'L_bearing', L_bearing, ...
                'r_bearing', r_bearing);

folder_root = 'H:\Calibration\112613\';


for t = 9
filename = strcat(folder_root, '100kV_2sec_200uA_P1_', num2str(t-1), '.ct');
close all;
% function [u, v, img] = find_centroid3(params, filename)

%% load image
dat = readBinary(filename, params.NpixelsHeight*params.NpixelsWidth+1, 'uint16');
dat = reshape(dat(2:end), params.NpixelsWidth, params.NpixelsHeight);
figure(1); imagesc(dat', [0 500]); colormap gray; axis equal; 

%% clip the original image
clip_dat = zeros(size(dat));
clip_dat(params.clip_lateral_edges:end-params.clip_lateral_edges,params.clip_height_top:end-params.clip_height_bottom)...
   = dat(params.clip_lateral_edges:end-params.clip_lateral_edges,params.clip_height_top:end-params.clip_height_bottom);
% figure(2); imagesc(clip_dat', [0 500]); colormap gray; colorbar; axis equal;
% 
% clear dat;
% 
% dat_noedge = zeros(size(clip_dat));
% dat_noedge(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width) = ...
%            clip_dat(params.boundary_width:end-params.boundary_width, params.boundary_width:end-params.boundary_width); 
% figure(3); imagesc(dat_noedge', [0 500]); colormap gray; axis equal;       
%% simple threshold 
n = find(clip_dat(:)<= params.threshold_range(2) & clip_dat(:)>= params.threshold_range(1));
img_simple_threshold = zeros(size(clip_dat));
img_simple_threshold(n) = clip_dat(n);
% figure(4); imagesc(img_simple_threshold'); colormap gray; colorbar; axis equal;
% clear n;

img_filtered = medfilt2(img_simple_threshold);
% figure(5); imagesc(img_filtered'); colormap gray; colorbar; axis equal;

[ycoord, xcoord] = find(img_filtered <= params.threshold_range(2) & img_filtered >= params.threshold_range(1));
%% find clusters
key = zeros(size(xcoord));
radius = 35;
counter = 0;
Ncluster = [];
index_array = 1:length(xcoord);
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
while size(index_array,2)~= 0

     %setting up stuff
     xcom = xcoord(index_array(1));
     ycom = ycoord(index_array(1));
     COM_xsum = 0;
     COM_ysum = 0;
     counter = 0;
     index = [];

     for n = 1:length(index_array)
        %calculate the distance of the next point away from xcom and ycom
        distance = sqrt( ( xcom-xcoord(index_array(n)) ).^2 + ( ycom-ycoord(index_array(n)) ).^2 );


        if (distance <= radius)
            %key basically iterates through all the identified clusters
            key(index_array(n)) = keyvalue;
            index = [index n];
            counter = counter + 1;
            COM_xsum = COM_xsum + xcoord(index_array(n));
            COM_ysum = COM_ysum + ycoord(index_array(n));
            xcom = COM_xsum/counter;
            ycom = COM_ysum/counter;

%             figure(5);
%             plot(xcom, ycom, 'or'); hold on;
%             plot(xcoord(index_array(index)), ycoord(index_array(index)), '.'); hold off;
%             xlim([xcom-30 xcom+30]); ylim([ycom-30 ycom+30]);
        end
     end

%          figure(50);
%          plot(xcoord(index_array(index)), ycoord(index_array(index)), '.', 'Color', [rand rand rand]); hold on;
%          plot(xcom, ycom, 'o', 'MarkerSize', 20);
%          axis equal;

     %calculate x and y range values
     xrange = range(xcoord(index_array(index)));
     yrange = range(ycoord(index_array(index)));
     xcell = xcoord(index_array(index));
     ycell = ycoord(index_array(index));

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



img = zeros(size(dat));




% plot subplots of clusters
figure(1000);

for j = 1:size(g,1)
    
    x = cell2mat(g(j,1));
    y = cell2mat(g(j,2));
    
    for test_index = 1:length(x)
    img(y(test_index), x(test_index)) = dat(y(test_index), x(test_index));
    end
     
    subplot(2,4,j); 
    plot( x, y,  '.'); 
    xlim([cell2mat(g(j,3))-30 cell2mat(g(j,3))+30]); ylim([cell2mat(g(j,4))-30 cell2mat(g(j,4))+30]);
    axis square;
    hold on;
    plot( cell2mat(g(j,3)), cell2mat(g(j,4)), 'or');
    hold off;
    
end

v = cell2mat(g(:,3));
u = cell2mat(g(:,4));

% pause();     


end

