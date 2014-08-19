%% written by Helen Fan
% outputs a 2D image of gaussian function
% L = used to create the x and y indices of the gaussian function,
% i.e. x = linspace(-L, L, 2*L+1), 
% r = radius
% A = max amplitude, usually 1... 
% good luck!

function out = gaussian2D(L, r, A)

[X,Y] = meshgrid(linspace(-L, L, 2*L+1),linspace(-L, L, 2*L+1));
out = A*exp( -( X.^2+Y.^2)/(2*r^2) );