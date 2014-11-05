clear all
close all
clc

cd('H:\Calibration\112613');

dat = readBinary('100kV_2sec_200uA_P1_0.ct', 2560*2160+1, 'uint16');
dat = dat(2:end);
dat = reshape(dat, 2560, 2160);
dat = dat';
figure; imagesc(dat, [0 550]); colormap gray; axis equal;
% set(gca, 'FontSize', 12, 'ylim', [1 2160], 'xlim', [1 2560], 'XMinorTick', 'off', 'PlotBoxAspectRatio', [1.5 1.5 1.5]);
% set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 5, 5]); % white bckgr
% xlabel('pixels'); ylabel('pixels');
% name = 'calibration_projection';
% export_fig(gcf, ...      % figure handle
%     name,... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', '-eps', ...           % file format
%     '-r100' );             % resolution in dpi
%%
x = [851 1931 1641 245 146];
y = [1372 184 1100 509 633];
s = 10;

% xlabel('pixels'); ylabel('pixels');
% name = 'calibration_projection';
% export_fig(gcf, ...      % figure handle
%     name,... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', '-eps', ...           % file format
%     '-r100' );             % resolution in dpi

a1 = dat(y(1)-s:y(1)+s, x(1)-s:x(1)+s);
a2 = dat(y(2)-s:y(2)+s, x(2)-s:x(2)+s);
a3 = dat(y(3)-s:y(3)+s, x(3)-s:x(3)+s);
a4 = dat(y(4)-s:y(4)+s, x(4)-s:x(4)+s);
a5 = dat(y(5)-s:y(5)+s, x(5)-s:x(5)+s);
figure; imagesc(a1, [0 550]); colormap gray; axis equal; axis off
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 4.4, 4.4]); % white bckgr
name = 'a1';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi

figure; imagesc(a2, [0 550]); colormap gray; axis equal; axis off
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 4.4, 4.4]); % white bckgr
name = 'a2';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi


figure; imagesc(a3, [0 550]); colormap gray; axis equal; axis off
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 4.4, 4.4]); % white bckgr
name = 'a3';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi


figure; imagesc(a4, [0 550]); colormap gray; axis equal; axis off
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 4.4, 4.4]); % white bckgr
name = 'a4';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi


figure; imagesc(a5, [0 550]); colormap gray; axis equal; axis off
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white', 'Units', 'inches', 'Position', [4 4 4.4, 4.4]); % white bckgr
name = 'a5';
export_fig(gcf, ...      % figure handle
    name,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', '-eps', ...           % file format
    '-r100' );             % resolution in dpi



