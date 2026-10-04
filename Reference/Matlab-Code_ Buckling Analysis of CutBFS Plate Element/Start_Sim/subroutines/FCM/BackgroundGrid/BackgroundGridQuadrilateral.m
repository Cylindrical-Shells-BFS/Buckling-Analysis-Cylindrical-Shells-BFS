%% Nodal distribution in a quadrilateral
%-------------------------------------------------------------------------%
%
% (c) 01.09.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Compute the global coordinates of a fine element background grid.
%
%-------------------------------------------------------------------------%
%% Input
%
% n:        Number of points in a local direction
%
%-------------------------------------------------------------------------%
%% Output
%
% [xi,eta]:	Nodal coordinates
%
%-------------------------------------------------------------------------%
function [xi,eta] = BackgroundGridQuadrilateral(n)

% Number of points in the background grid (interior points)
NumberOfPointsBG = (n-2)*(n-2);

% Local coordinates of the background grid
coor = linspace(-1,1,n);

% Divide the background grid into two parts: (1) points on the boundary of
% the element and (2) points in the interior of the element
% Boundary vertices
xi.Boundary  = [ coor,      ones(1,n-1), -coor(2:n), -ones(1,n-2)];
eta.Boundary = [-ones(1,n), coor(2:n),    ones(1,n-1), -coor(2:n-1)];

% Generate a two-dimensional grid
[Xi,Eta] = meshgrid(coor(2:n-1),coor(2:n-1));

% Interior vertices
xi.Interior   = reshape(Xi,1,NumberOfPointsBG);
eta.Interior  = reshape(Eta,1,NumberOfPointsBG);

end